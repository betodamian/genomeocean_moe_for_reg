#!/usr/bin/env python3
"""
Compute biophysical Shine-Dalgarno strength for every RBS TIS window as the
hybridization free energy (ΔG) between the upstream region and the organism's own
16S rRNA 3' tail (the anti-SD). research_plan §5d:

  "SD class is assigned with organism-specific anti-SD sequences (the 16S rRNA 3'
   tail of each genome) and organism-specific ΔG thresholds, not one global cutoff."

This REPLACES the circular regex SD labels (build_rbs_windows.py sd_label = "does the
upstream contain a GGAGG-like substring") that made Phase-0's SD-vs-UNSD task a
substring-rediscovery test (a k-mer counter beat the model; PHASE0_FINDINGS.md).
ΔG is non-circular (thermodynamic base-pairing to anti-SD), GC-aware (G:C vs A:U
weighted), and organism-specific — so the SD task becomes a fair §9 routing probe.

Method:
  - anti-SD = last ANTISD_LEN nt of each genome's 16S rRNA (3' end, transcription dir),
    as RNA. (E. coli -> GAUCACCUCCUUA, the canonical CCUCC core; validated.)
  - For each TIS window, take the upstream region window[U0:U1] (positions ~-20..-3
    relative to the start codon at index 150) and compute the minimum-energy duplex
    with the anti-SD via ViennaRNA RNA.duplexfold. More negative ΔG = stronger SD.
  - Canonical binary call (interpretability only): SD if ΔG <= SD_DG_THRESH.

Output: data/datasets/rbs/sd_deltaG.tsv  (window_id, organism, dG, sd_call)
Run:    .venv/bin/python pipeline/compute_sd_deltaG.py
"""
import csv, os, sys
import RNA

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
GEN  = os.path.join(ROOT, "data/genomes")
ALL  = os.path.join(ROOT, "data/datasets/rbs/ALL.tsv")
OUT  = os.path.join(ROOT, "data/datasets/rbs/sd_deltaG.tsv")

ANTISD_LEN   = 13          # nt of the 16S 3' tail to use as anti-SD
WIN_HALF     = 150         # TIS is at window index 150 (0-based) — start codon
U0, U1       = 130, 148    # upstream region scanned for SD: ~ -20 .. -3 nt
SD_DG_THRESH = -3.4        # kcal/mol; canonical E. coli/Salis cut (interpretability only)

# organism -> genome dir (same panel as build_rbs_windows.py)
ORG_GENOME = {
    "ecoli_K12_MG1655":    "ecoli_K12_MG1655",
    "ecoli_BL21_DE3":      "ecoli_BL21_DE3",
    "bsubtilis_168":       "bsubtilis_168",
    "mtuberculosis_H37Rv": "mtuberculosis_H37Rv",
    "ccrescentus_NA1000":  "ccrescentus_NA1000",
    "saureus_HG001":       "saureus_HG001",
    "hvolcanii_DS2":       "hvolcanii_DS2",
}
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
def rc(s): return s.translate(COMP)[::-1]


def load_fasta(path):
    seqs, cid, buf = {}, None, []
    with open(path) as fh:
        for ln in fh:
            if ln.startswith(">"):
                if cid: seqs[cid] = "".join(buf).upper()
                cid = ln[1:].split()[0]; buf = []
            else:
                buf.append(ln.strip())
    if cid: seqs[cid] = "".join(buf).upper()
    return seqs


def is_16s(attrs):
    a = attrs.lower()
    return ("16s" in a) and ("ribosomal rna" in a or "rrna" in a or "ribosomal-rna" in a
                             or "rna-16s" in a)


def extract_anti_sd(gdir):
    """Return the anti-SD RNA (3' tail of the 16S rRNA) for a genome dir, or None."""
    d = os.path.join(GEN, gdir)
    fna = next(os.path.join(d, f) for f in os.listdir(d) if f.endswith(".fna"))
    gff = next(os.path.join(d, f) for f in os.listdir(d) if f.endswith(".gff"))
    seqs = load_fasta(fna)
    with open(gff) as fh:
        for ln in fh:
            if ln.startswith("#") or "\t" not in ln:
                continue
            c = ln.rstrip("\n").split("\t")
            if len(c) < 9 or c[2] != "rRNA":
                continue
            if not is_16s(c[8]):
                continue
            s, e, st = int(c[3]), int(c[4]), c[6]
            gene = seqs[c[0]][s-1:e] if st == "+" else rc(seqs[c[0]][s-1:e])
            return gene[-ANTISD_LEN:].replace("T", "U")     # 3' end of mature 16S
    return None


def main():
    print("extracting per-organism anti-SD (16S 3' tail) ...")
    anti = {}
    for org, gdir in ORG_GENOME.items():
        a = extract_anti_sd(gdir)
        if a is None:
            print(f"  !! {org}: no 16S rRNA feature found — skipping", file=sys.stderr)
            continue
        anti[org] = a
        print(f"  {org:22s} anti-SD = {a}")

    print("\ncomputing upstream:anti-SD ΔG per TIS window ...")
    rows = []
    n_sd = 0
    with open(ALL) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if not r["label"].startswith("positive_"):
                continue                      # only TIS windows have an SD region
            org = r["organism"]
            if org not in anti:
                continue
            win = r["window_seq"]
            if len(win) < U1:
                continue
            upstream = win[U0:U1].replace("T", "U")
            if "N" in upstream:
                continue
            dup = RNA.duplexfold(anti[org], upstream)
            dg = round(dup.energy, 3)
            call = "SD" if dg <= SD_DG_THRESH else "UNSD"
            if call == "SD":
                n_sd += 1
            rows.append((r["id"], org, dg, call))

    with open(OUT, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["window_id", "organism", "dG", "sd_call"])
        w.writerows(rows)

    print(f"\nwrote {len(rows):,} ΔG labels ({n_sd:,} SD / {len(rows)-n_sd:,} UNSD "
          f"at ΔG<={SD_DG_THRESH}) -> {os.path.relpath(OUT, ROOT)}")
    # per-organism ΔG summary
    import statistics as st
    by = {}
    for _id, org, dg, _c in rows:
        by.setdefault(org, []).append(dg)
    print("\nper-organism ΔG (more negative = stronger SD):")
    for org in sorted(by):
        v = by[org]
        sd_frac = sum(1 for x in v if x <= SD_DG_THRESH) / len(v)
        print(f"  {org:22s} n={len(v):5d}  median={st.median(v):+.2f}  "
              f"mean={st.mean(v):+.2f}  SD_frac={sd_frac:.2f}")


if __name__ == "__main__":
    main()
