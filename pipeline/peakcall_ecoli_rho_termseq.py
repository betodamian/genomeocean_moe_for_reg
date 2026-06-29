#!/usr/bin/env python3
"""
Call high-confidence bp-resolution Rho-dependent termination sites for E. coli K-12
by intersecting two orthogonal in vivo datasets (research_plan §5c, §5d):

  1. BCM BSTs (Peters 2012)   — regions of BCM-induced RNA readthrough  → Rho-specificity
  2. Term-seq 3'-ends (NAR 2018, GEO GSE109766) — exact RNA 3'-end positions → bp resolution

A site is called when the dominant Term-seq 3'-end WITHIN a BCM BST clears the
reproducibility gate (sum across 3 reps >= MIN_SUM, ≥ MIN_REPS reps with signal).
Sites outside any BCM BST are NOT called — BCM concordance is the Rho-specificity filter.

This avoids centering 300-bp windows on large BST regions (median 557 bp, max 6.5 kb)
where the actual termination point could be anywhere inside the region.

Thresholds (data-driven):
  MIN_SUM  = 10  → sum of 3 reps at the best position within the BST
  MIN_REPS = 2   → at least 2 of 3 replicates must have reads

Validation: print how many BSTs yield a site, signal distribution, strand balance.

Output: data/rho_database/raw/ECOLI_BCM_RHO/ecoli_rho_termseq_sites.tsv
Run:    .venv/bin/python pipeline/peakcall_ecoli_rho_termseq.py
"""
import csv, gzip, os
import numpy as np
from collections import defaultdict

ROOT     = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TERMSEQ  = os.path.join(ROOT, "data/rho_database/raw",
                        "GSE109766_ecoli_termseq_rep1-3.counts_per_position.txt.gz")
BCM_TSV  = os.path.join(ROOT, "data/rho_database/raw/ECOLI_BCM_RHO",
                        "ecoli_bcm_rho_sites.tsv")
OUT      = os.path.join(ROOT, "data/rho_database/raw/ECOLI_BCM_RHO",
                        "ecoli_rho_termseq_sites.tsv")

MIN_SUM  = 10   # sum of 3 Term-seq reps at peak position
MIN_REPS = 2    # min replicates with at least 1 read

RHO_FIELDS = ["source_id", "organism", "strain", "site_id",
              "site_start", "site_end", "strand", "site_class",
              "terminates_or_silences", "sequence", "evidence_type", "notes"]


def load_termseq(path):
    """Load per-position Term-seq counts: (pos, strand) -> [rep1, rep2, rep3]."""
    data = defaultdict(lambda: [0, 0, 0])
    with gzip.open(path, "rt") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            pos = int(r["Genomic_position"].strip())
            st  = r["strand"].strip()
            data[(pos, st)][0] += int(r["reads_rep1"])
            data[(pos, st)][1] += int(r["reads_rep2"])
            data[(pos, st)][2] += int(r["reads_rep3"])
    return data


def load_bsts(path):
    """Load BCM BSTs: list of (site_id, start, end, strand)."""
    bsts = []
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            bsts.append((r["site_id"], int(r["site_start"]), int(r["site_end"]), r["strand"]))
    return bsts


def main():
    print("loading Term-seq 3'-end counts ...")
    ts_data = load_termseq(TERMSEQ)
    print(f"  {len(ts_data):,} positions with signal")

    print("loading BCM BSTs ...")
    bsts = load_bsts(BCM_TSV)
    print(f"  {len(bsts):,} BST regions")

    # Build per-BST index of Term-seq positions
    # Efficiency: sort positions, then scan BSTs
    pos_by_strand = {"+": {}, "-": {}}
    for (pos, st), reps in ts_data.items():
        pos_by_strand[st][pos] = reps

    calls = []
    n_no_signal = n_below_thresh = 0

    for bst_id, s, e, st in bsts:
        # find all Term-seq positions inside this BST on matching strand
        strand_data = pos_by_strand.get(st, {})
        candidates = {p: reps for p, reps in strand_data.items() if s <= p <= e}

        if not candidates:
            n_no_signal += 1
            continue

        # Pick the position with the highest total signal
        best_pos = max(candidates, key=lambda p: sum(candidates[p]))
        reps = candidates[best_pos]
        total = sum(reps)
        n_nz  = sum(1 for r in reps if r > 0)

        if total < MIN_SUM or n_nz < MIN_REPS:
            n_below_thresh += 1
            continue

        calls.append(dict(
            source_id="ECOLI_BCM_RHO",
            organism="ecoli_K12_MG1655",
            strain="MG1655",
            site_id=f"RhoTS_{bst_id}",
            site_start=best_pos,
            site_end=best_pos,
            strand=st,
            site_class="rho_termseq_concordance",
            terminates_or_silences="T",
            sequence=None,
            evidence_type="T1_invivo",
            notes=(f"bcm_bst={bst_id};bst_range={s}-{e};"
                   f"termseq_sum={total};rep_counts={reps[0]},{reps[1]},{reps[2]};"
                   f"method=BCM_concordance_termseq;MIN_SUM={MIN_SUM};MIN_REPS={MIN_REPS}"),
        ))

    # ── validation ────────────────────────────────────────────────────────────
    print(f"\nVALIDATION")
    print(f"  BSTs with no Term-seq signal       : {n_no_signal}/{len(bsts)}")
    print(f"  BSTs with signal but below threshold: {n_below_thresh}")
    print(f"  High-confidence bp-resolution sites : {len(calls)}")

    if calls:
        totals  = [int(c["notes"].split("termseq_sum=")[1].split(";")[0]) for c in calls]
        arr = np.array(totals)
        print(f"  Signal (sum 3 reps) — median={np.median(arr):.0f}  "
              f"p10={np.percentile(arr,10):.0f}  p90={np.percentile(arr,90):.0f}  "
              f"max={arr.max()}")
        n_plus  = sum(1 for c in calls if c["strand"] == "+")
        n_minus = len(calls) - n_plus
        print(f"  Strand: +={n_plus}  -={n_minus}")

    # ── write ─────────────────────────────────────────────────────────────────
    with open(OUT, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=RHO_FIELDS, delimiter="\t",
                           lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        w.writerows(calls)
    print(f"\nwrote {len(calls):,} sites -> {os.path.relpath(OUT, ROOT)}")


if __name__ == "__main__":
    main()
