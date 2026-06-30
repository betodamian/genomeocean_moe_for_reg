#!/usr/bin/env python3
"""
Build the CONFOUND-FREE promoter-detection test (within-class dinucleotide-shuffle).

Question: does GO-MoE encode promoter-specific structure that no composition, region,
or naturalness confound can explain? Compare each class to its OWN dinucleotide-
preserving shuffle, on UPSTREAM-ONLY windows, and subtract the intergenic baseline:

    signal = MCC(P vs P') - MCC(I vs I')

  P  = real promoter window, upstream-only (-80..+1; window_seq[70:151], TSS at idx 150)
  P' = dinucleotide-preserving shuffle of P   (same length/GC/dinuc; motif destroyed)
  I  = same upstream slice of the intergenic decoys
  I' = dinucleotide-preserving shuffle of I

Confounds killed: downstream-genic (upstream-only), GC + dinucleotide (shuffle preserves
both -> invisible in X-vs-X'), naturalness (cancels in the P−I subtraction), region-type
(never compare P to I). See experiments/phase0/PROMOTER_DETECTION_TEST.md.

Output: data/phase0/promoter_dinuc_ALL.tsv  (id, organism, split, klass, seq)
        klass in {P, P_shuf, I, I_shuf}
Run: .venv/bin/python pipeline/build_promoter_dinuc_test.py
"""
import csv, os
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ALL   = os.path.join(ROOT, "data/datasets/promoters/ALL.tsv")
SPLIT = os.path.join(ROOT, "splits/promoters/split_80_20.tsv")
OUT   = os.path.join(ROOT, "data/phase0/promoter_dinuc_ALL.tsv")
LO, HI = 70, 151          # -80..+1 relative to TSS (index 150); 81 bp ~ 16.5 tokens


def dinuc_shuffle(seq, rng):
    """Altschul-Erikson dinucleotide-preserving shuffle (verified: exact dinuc counts)."""
    seq = seq.upper()
    if len(seq) < 3:
        return seq
    last = seq[-1]
    edges = {}
    for i in range(len(seq) - 1):
        edges.setdefault(seq[i], []).append(seq[i + 1])
    while True:
        shuf = {v: [edges[v][i] for i in rng.permutation(len(edges[v]))] for v in edges}
        last_edge = {v: shuf[v][-1] for v in shuf if v != last}
        ok = True
        for v in last_edge:
            cur, seen = v, set()
            while cur != last:
                if cur in seen or cur not in last_edge:
                    ok = False; break
                seen.add(cur); cur = last_edge[cur]
            if not ok:
                break
        if ok:
            break
    out = [seq[0]]; idx = {v: 0 for v in shuf}; cur = seq[0]
    for _ in range(len(seq) - 1):
        nxt = shuf[cur][idx[cur]]; idx[cur] += 1
        out.append(nxt); cur = nxt
    return "".join(out)


def dinuc_counts(s):
    from collections import Counter
    return Counter(s[i:i+2] for i in range(len(s) - 1))


def main():
    split_of = {}
    with open(SPLIT) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            split_of[r["window_id"]] = r["split"]

    rows = []
    with open(ALL) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            seq = r["window_seq"].upper()
            if len(seq) < HI:
                continue
            sl = seq[LO:HI]
            if set(sl) - set("ACGT"):            # skip windows with N etc.
                continue
            klass = "P" if r["label"] == "positive" else "I"
            rows.append((r["id"], r["organism"], split_of.get(r["id"], "none"), klass, sl))

    n_bad = 0
    out_rows = []
    for i, (wid, org, sp, klass, sl) in enumerate(rows):
        rng = np.random.default_rng((i * 2654435761) & 0xFFFFFFFF)
        sh = dinuc_shuffle(sl, rng)
        if dinuc_counts(sh) != dinuc_counts(sl) or len(sh) != len(sl):
            n_bad += 1
            continue
        out_rows.append((f"{wid}__real", org, sp, klass, sl))
        out_rows.append((f"{wid}__shuf", org, sp, f"{klass}_shuf", sh))

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["id", "organism", "split", "klass", "seq"])
        w.writerows(out_rows)

    from collections import Counter
    kc = Counter(k for _, _, _, k, _ in out_rows)
    print(f"slice -80..+1 ({HI-LO} bp); dinuc-shuffle mismatches: {n_bad}")
    print(f"wrote {len(out_rows):,} rows -> {os.path.relpath(OUT, ROOT)}")
    print("  class counts:", dict(kc))
    # sanity: GC identical real vs shuf, motif destroyed at population level
    import statistics
    def gc(s): return (s.count("G")+s.count("C"))/len(s)
    P  = [r[4] for r in out_rows if r[3] == "P"]
    Ps = [r[4] for r in out_rows if r[3] == "P_shuf"]
    print(f"  GC(P)={statistics.mean(map(gc,P)):.4f}  GC(P_shuf)={statistics.mean(map(gc,Ps)):.4f} (must match)")
    print(f"  TATAAT in P: {sum('TATAAT' in s for s in P)}  in P_shuf: {sum('TATAAT' in s for s in Ps)}")


if __name__ == "__main__":
    main()
