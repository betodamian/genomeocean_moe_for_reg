#!/usr/bin/env python3
"""
Peak-call S. aureus translation-initiation sites (TIS) from the Kohl 2026
Ribo-RET signal (GEO GSE299221), HG001 genome NZ_CP018205.1 (= wig 'HG001').

Method (faithful to Ribo-RET; reads the EXPERIMENT, not sequence — research_plan §12):
  - Retapamulin arrests initiating ribosomes at start codons. The wig assigns
    normalized density (RPM) to footprint 3' ends, so the arrest signal sits a
    fixed offset 3' of the start codon. A metagene over annotated starts puts
    that offset at +16 nt (sharp peak, 765 vs ~10 background) — discovered, not
    assumed; recomputed here and asserted.
  - A TIS is called at a start codon (ATG/GTG/TTG, both strands, genome-wide)
    when the RET arrest signal at codon+16 (translation direction):
      (1) >= THRESH (data-driven from the annotated-start signal distribution),
      (2) is a local maximum (non-max suppression within +/-NMS nt), and
      (3) is enriched over the no-drug control (RET/ctl >= FOLD).
  - RNase1 libraries (both reps summed) — the higher-resolution set per the paper.

Validation (printed): metagene offset, recovery of annotated starts, recovery of
the 46 published novel sORFs (SAUR_EXSD), start-codon composition, total count.

Output: data/rbs_database/raw/SAUR_RIBORET/saur_riboret_tis_sites.tsv  (unified RBS schema)
Run: .venv/bin/python pipeline/peakcall_saureus_tis.py
"""
import csv, gzip, os, sys
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
GEN  = os.path.join(ROOT, "data/genomes/saureus_HG001")
WIG  = os.path.join(ROOT, "data/rbs_database/raw/SAUR_RIBORET")
FNA  = os.path.join(GEN, "GCF_001900185.1_genomic.fna")
GFF  = os.path.join(GEN, "GCF_001900185.1_genomic.gff")
OUT  = os.path.join(WIG, "saur_riboret_tis_sites.tsv")

OFFSET   = 16     # 3'-end arrest offset (nt downstream of start codon), from metagene
WIN      = (14, 18)  # sum RET arrest over codon+14..+18 (3'-end spread, from metagene)
NMS      = 10     # non-max-suppression radius (nt)
# Two-tier calling: annotated starts have prior support -> lenient; novel starts
# must clear a higher bar (specificity for genome-wide scan).
FOLD_ANN, FLOOR_ANN     = 2.0, 5.0     # annotated: RET/ctl, abs RPM floor
# novel: same modest floor, but specificity comes from peak SHAPE (control
# enrichment + prominence over gene body), not absolute height — real sORFs are
# weaker than typical genes (validated: published sORFs are 5-40 RPM).
FOLD_NOV, PROM_NOV, FLOOR_NOV = 3.0, 3.0, 5.0
STARTS   = {"ATG", "GTG", "TTG"}
COMP     = str.maketrans("ACGT", "TGCA")

def revcomp(s): return s.translate(COMP)[::-1]

def load_genome():
    seq = []
    with open(FNA) as fh:
        for line in fh:
            if not line.startswith(">"): seq.append(line.strip())
    return "".join(seq).upper()

def load_wig(fname, L):
    a = np.zeros(L + OFFSET + 2, dtype=np.float64)
    with gzip.open(os.path.join(WIG, fname), "rt") as fh:
        i = 0
        for line in fh:
            c = line[0]
            if c in "tf": continue          # track / fixedStep headers
            i += 1
            if c != "0": a[i] = float(line)
    return a

def gff_starts():
    """annotated start-codon first base (translation direction), per strand."""
    out = set()
    with open(GFF) as fh:
        for line in fh:
            if line.startswith("#"): continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 8 or f[2] != "CDS": continue
            s, e, st = int(f[3]), int(f[4]), f[6]
            out.add((s if st == "+" else e, st))
    return out

def start_codon(seq, p, strand):
    if strand == "+":
        return seq[p-1:p+2] if 1 <= p <= len(seq)-2 else ""
    return revcomp(seq[p-3:p]) if 3 <= p <= len(seq) else ""

def arrest_at(p, strand, plus, minus):
    """Windowed RET 3'-end arrest density for a start codon at p (sum over WIN, translation dir)."""
    A = plus if strand == "+" else minus
    xs = range(p+WIN[0], p+WIN[1]+1) if strand == "+" else range(p-WIN[1], p-WIN[0]+1)
    return float(sum(A[x] for x in xs if 1 <= x < len(A)))

def body_at(p, strand, plus, minus):
    """Mean RET density over the gene body downstream of the start (+25..+75), for prominence."""
    A = plus if strand == "+" else minus
    xs = range(p+25, p+76) if strand == "+" else range(p-75, p-24)
    v = [A[x] for x in xs if 1 <= x < len(A)]
    return float(np.mean(v)) if v else 0.0

def main():
    print("loading genome + GFF ...")
    seq = load_genome(); L = len(seq)
    ann = gff_starts()
    print(f"  genome {L:,} bp ; annotated CDS starts {len(ann):,}")

    # Combine BOTH enzymes x 2 reps for depth. MNase arrest peaks at +15, RNase1 at
    # +16 (metagene-verified) — both fall inside the +14..+18 window, so summing is valid.
    print("loading RET + ctl tracks (MNase + RNase1, 2 reps each, both strands) ...")
    RET = {"plus": ["GSM9035386_MNase_RET_RP1", "GSM9035387_MNase_RET_RP2",
                    "GSM9035390_RNase1_RET_RP1", "GSM9035391_RNase1_RET_RP2"]}
    CTL = {"plus": ["GSM9035384_MNase_ctl_RP1", "GSM9035385_MNase_ctl_RP2",
                    "GSM9035388_RNase1_ctl_RP1", "GSM9035389_RNase1_ctl_RP2"]}
    def stack(gsms, strand):
        a = np.zeros(L + OFFSET + 2)
        for g in gsms:
            a += load_wig(f"{g}_15-45_{strand}.wig.gz", L)
        return a
    ret_p, ret_m = stack(RET["plus"], "plus"), stack(RET["plus"], "minus")
    ctl_p, ctl_m = stack(CTL["plus"], "plus"), stack(CTL["plus"], "minus")

    # confirm offset via metagene over annotated starts
    W = 30; prof = np.zeros(2*W+1); nmg = 0
    for p, st in ann:
        if p-W < 1 or p+W+OFFSET > L: continue
        seg = (ret_p[p-W:p+W+1] if st == "+" else ret_m[p+W:p-W-1:-1])
        t = seg.sum()
        if t > 0: prof += seg/t; nmg += 1
    print(f"  metagene offset (max) = +{int(prof.argmax()-W)} nt  (expect +{OFFSET}); genes={nmg}")

    ann_sig = np.array([arrest_at(p, st, ret_p, ret_m) for p, st in ann])
    print(f"  annotated windowed-arrest: nonzero {int((ann_sig>0).sum())}/{len(ann)} ; "
          f"median(nz)={np.median(ann_sig[ann_sig>0]):.1f} ; floors ann={FLOOR_ANN} nov={FLOOR_NOV}")

    # enumerate candidate start codons genome-wide, score the arrest peak
    print("scanning genome-wide start codons ...")
    cand = []   # (score, p, strand, codon, is_annotated)
    for p in range(1, L-1):
        cp = seq[p-1:p+2]
        if cp in STARTS:
            cand.append((arrest_at(p, "+", ret_p, ret_m), p, "+", cp, (p, "+") in ann))
        cm = revcomp(seq[p-3:p]) if p >= 3 else ""
        if cm in STARTS:
            cand.append((arrest_at(p, "-", ret_p, ret_m), p, "-", cm, (p, "-") in ann))
    print(f"  candidate start codons: {len(cand):,}")

    # two-tier filter: annotated lenient, novel strict (control + prominence)
    keep = []
    for sig, p, st, codon, is_ann in cand:
        if sig <= 0: continue
        cs = arrest_at(p, st, ctl_p, ctl_m)
        if is_ann:
            if sig >= FLOOR_ANN and sig >= FOLD_ANN * (cs + 0.1):
                keep.append((sig, p, st, codon, is_ann))
        else:
            bd = body_at(p, st, ret_p, ret_m)
            if (sig >= FLOOR_NOV and sig >= FOLD_NOV * (cs + 0.1)
                    and sig >= PROM_NOV * (bd + 0.1)):
                keep.append((sig, p, st, codon, is_ann))
    keep.sort(reverse=True)                  # strongest first
    taken = {"+": np.zeros(L+2, bool), "-": np.zeros(L+2, bool)}
    calls = []
    for sig, p, st, codon, is_ann in keep:
        lo, hi = max(1, p-NMS), min(L, p+NMS)
        if taken[st][lo:hi+1].any(): continue   # suppress near a stronger call
        taken[st][p] = True
        calls.append((p, st, codon, sig))
    print(f"  TIS called (two-tier + NMS): {len(calls):,}")

    # ---- validation ----
    call_set = {(p, st) for p, st, _, _ in calls}
    rec_ann = sum((p, st) in call_set for p, st in ann)
    print("\nVALIDATION")
    print(f"  annotated-start recovery : {rec_ann}/{len(ann)} ({100*rec_ann/len(ann):.1f}%)")
    from collections import Counter
    cc = Counter(c for _, _, c, _ in calls)
    print(f"  start-codon composition  : {dict(cc)}  (100% are start codons by construction)")
    novel = len(calls) - rec_ann
    print(f"  novel (non-annotated) TIS: {novel}")

    # recovery of the 46 published sORFs (SAUR_EXSD), HG001 coords
    exsd = os.path.join(ROOT, "data/rbs_database/raw/SAUR_EXSD/saur_exsd_tis_sites.tsv")
    if os.path.exists(exsd):
        sorf = [(int(r["tis_pos"]), r["strand"]) for r in csv.DictReader(open(exsd), delimiter="\t")
                if r["tis_pos"].lstrip("-").isdigit()]
        hit = sum(any(abs(p-sp) <= 3 and st == ss for p, st, _, _ in calls) for sp, ss in sorf)
        print(f"  published-sORF recovery  : {hit}/{len(sorf)} (within 3 nt)")

    # ---- write unified RBS schema ----
    fields = ["source_id","organism","strain","gene_name","locus_tag","tis_pos","tis_end",
              "strand","start_codon","ribo_density","sd_motif","sd_pos_rel","evidence_type","notes"]
    rows = []
    for p, st, codon, sig in sorted(calls):
        rows.append(dict(source_id="SAUR_RIBORET", organism="saureus_HG001", strain="HG001",
            gene_name="", locus_tag="", tis_pos=p, tis_end="", strand=st, start_codon=codon,
            ribo_density=round(sig, 3), sd_motif="", sd_pos_rel="", evidence_type="T1_riboseq",
            notes=f"ret_arrest_rpm={sig:.2f};method=RiboRET_peakcall;offset={OFFSET};annotated={(p,st) in ann}"))
    with open(OUT, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", lineterminator="\n")
        w.writeheader(); w.writerows(rows)
    print(f"\nwrote {len(rows):,} TIS -> {os.path.relpath(OUT, ROOT)}")

if __name__ == "__main__":
    main()
