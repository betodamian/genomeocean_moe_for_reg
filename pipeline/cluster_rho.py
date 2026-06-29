#!/usr/bin/env python3
"""
§5f data-census gate + §5a/§5b sequence-dissimilar split, for the Rho windows.

Only T1 in-vivo positives (positive_rho_t1, positive_rho_t1_rsr) enter clustering
and the census gate — T2 in-vitro sites (positive_rho_t2) are not co-mingled with
in-vivo evidence per research_plan §5d.

Census gate (pre-registered, Rho-specific):
  - cross-genome unseen claim: ≥ 30 clusters spanning ≥ 2 in-vivo phyla
                                (E. coli Gammaproteobacteria ↔ MTB Actinobacteria)
  - intra-genome claim:        ≥ 30 clusters per organism
  The 2-phylum scope is stated honestly: this is a "2-organism transfer", not a
  broad multi-phylum sweep (research_plan §5d honest-scope note).

T2 positives and all decoy classes receive a simple per-organism 80/20 random split.

Outputs (committed): splits/rho/{clusters.tsv, split_80_20.tsv, folds_loo.tsv, census.md}
Run: .venv/bin/python pipeline/cluster_rho.py
"""
import csv, os, random, subprocess
from collections import defaultdict, Counter

random.seed(20260625)
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DS   = os.path.join(ROOT, "data/datasets/rho")
MM   = os.path.join(DS, "mmseqs")
OUT  = os.path.join(ROOT, "splits/rho")
ALL  = os.path.join(DS, "ALL.tsv")

MIN_ID        = 0.60
VAL_FRAC      = 0.20
GATE_CLUSTERS = 30
GATE_PHYLA    = 2    # Rho: 2-phylum in-vivo panel (honest scope)

T1_LABELS = {"positive_rho_t1", "positive_rho_t1_rsr"}

def load_rows():
    with open(ALL) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))

def write_fasta(rows, path):
    """Write T1 positives only for clustering."""
    with open(path, "w") as fh:
        for r in rows:
            if r["label"] in T1_LABELS:
                fh.write(f">{r['id']}\n{r['window_seq']}\n")

def run_mmseqs(fasta, prefix):
    os.makedirs(MM, exist_ok=True)
    tmp = os.path.join(MM, "tmp")
    cmd = ["mmseqs", "easy-cluster", fasta, prefix, tmp,
           "--min-seq-id", str(MIN_ID), "-c", "0.8", "--cov-mode", "0",
           "-v", "1", "--threads", str(os.cpu_count() or 4)]
    print("  $", " ".join(cmd))
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL)
    return prefix + "_cluster.tsv"

def parse_clusters(tsv):
    m = {}
    with open(tsv) as fh:
        for ln in fh:
            rep, mem = ln.rstrip("\n").split("\t")
            m[mem] = rep
    return m

def main():
    os.makedirs(OUT, exist_ok=True)
    rows   = load_rows()
    t1_pos = [r for r in rows if r["label"] in T1_LABELS]
    org_of = {r["id"]: r["organism"] for r in rows}
    phy_of = {r["id"]: r["phylum"]   for r in rows}

    fasta = os.path.join(DS, "t1_positives.fasta")
    write_fasta(rows, fasta)
    print(f"Clustering {len(t1_pos):,} T1 Rho positives at ≤{int(MIN_ID*100)}% identity ...")
    ctsv    = run_mmseqs(fasta, os.path.join(MM, "rho"))
    mem2rep = parse_clusters(ctsv)

    clus   = defaultdict(list)
    for mid, rep in mem2rep.items():
        clus[rep].append(mid)
    cl_org = {c: set(org_of[m] for m in members) for c, members in clus.items()}
    cl_phy = {c: set(phy_of[m] for m in members) for c, members in clus.items()}

    n_clusters     = len(clus)
    singletons     = sum(1 for c in clus if len(clus[c]) == 1)
    cross_org_clus = sum(1 for c in cl_org if len(cl_org[c]) > 1)
    all_phyla      = sorted(set(phy_of[m] for m in mem2rep))
    n_phyla        = len(all_phyla)

    per_org_clusters = Counter()
    for c, orgs in cl_org.items():
        for o in orgs: per_org_clusters[o] += 1

    cross_ok = n_clusters >= GATE_CLUSTERS and n_phyla >= GATE_PHYLA
    intra_ok = {o: per_org_clusters[o] >= GATE_CLUSTERS for o in per_org_clusters}

    # ── split T1 positives at cluster level ───────────────────────────────────
    cl_primary_org = {c: Counter(org_of[m] for m in clus[c]).most_common(1)[0][0]
                      for c in clus}
    by_org = defaultdict(list)
    for c in clus: by_org[cl_primary_org[c]].append(c)

    split_of_cluster = {}
    for o, cs in by_org.items():
        random.shuffle(cs)
        nval = max(1, round(len(cs) * VAL_FRAC))
        for i, c in enumerate(cs):
            split_of_cluster[c] = "val" if i < nval else "train"

    # ── write output files ────────────────────────────────────────────────────
    with open(os.path.join(OUT, "clusters.tsv"), "w") as fh:
        fh.write("window_id\tcluster_rep\torganism\tphylum\tlabel\n")
        for mid, rep in sorted(mem2rep.items()):
            fh.write(f"{mid}\t{rep}\t{org_of[mid]}\t{phy_of[mid]}\t"
                     f"{next(r['label'] for r in rows if r['id']==mid)}\n")

    with open(os.path.join(OUT, "split_80_20.tsv"), "w") as fh:
        fh.write("window_id\tlabel\torganism\tsplit\n")
        # T1 positives: cluster-level split
        for r in t1_pos:
            fh.write(f"{r['id']}\t{r['label']}\t{r['organism']}\t"
                     f"{split_of_cluster[mem2rep[r['id']]]}\n")
        # T2 positives + all decoys: per-organism random 80/20
        other = [r for r in rows if r["label"] not in T1_LABELS]
        by_org_other = defaultdict(list)
        for r in other: by_org_other[r["organism"]].append(r)
        for o, rs in by_org_other.items():
            random.shuffle(rs)
            nval = round(len(rs) * VAL_FRAC)
            for i, r in enumerate(rs):
                fh.write(f"{r['id']}\t{r['label']}\t{o}\t"
                         f"{'val' if i < nval else 'train'}\n")

    # LOO: E. coli and MTB are the T1 in-vivo phyla; B. subtilis is T2 only
    with open(os.path.join(OUT, "folds_loo.tsv"), "w") as fh:
        fh.write("# leave-one-organism-out for 2-phylum Rho in-vivo cross-genome test (P5, P8)\n")
        fh.write("# Only T1 in-vivo organisms are in the LOO; B. subtilis is T2 in-vitro only\n")
        fh.write("held_out_organism\tphylum\tlabel_types\tn_t1_clusters\n")
        for o in sorted(by_org):
            phy = next(phy_of[m] for c in by_org[o] for m in clus[c])
            lbls = ";".join(sorted({next(r["label"] for r in rows if r["id"]==m)
                                    for c in by_org[o] for m in clus[c]}))
            fh.write(f"{o}\t{phy}\t{lbls}\t{len(by_org[o])}\n")

    # ── census report ─────────────────────────────────────────────────────────
    n_t1_ecoli = sum(1 for r in t1_pos if r["organism"] == "ecoli_K12_MG1655")
    n_t1_mtb   = sum(1 for r in t1_pos if r["organism"] == "mtuberculosis_H37Rv")
    n_t2_total = sum(1 for r in rows if r["label"] == "positive_rho_t2")
    n_dec_i    = sum(1 for r in rows if r["label"] == "decoy_intrinsic")
    n_dec_g    = sum(1 for r in rows if r["label"] == "decoy_intergenic")

    lines = []
    lines.append("# Rho data-census gate (§5f)\n\n")
    lines.append(f"- clustering: MMseqs2 easy-cluster, min-seq-id={MIN_ID}, -c=0.8\n")
    lines.append(f"- T1 in-vivo positives: **{len(t1_pos):,}** "
                 f"(E. coli={n_t1_ecoli}, MTB={n_t1_mtb})\n")
    lines.append(f"  -> **{n_clusters:,}** ≤{int(MIN_ID*100)}%-identity clusters "
                 f"({singletons:,} singletons, {cross_org_clus:,} cross-organism)\n")
    lines.append(f"- in-vivo phyla: **{n_phyla}** ({', '.join(all_phyla)})\n")
    lines.append(f"- T2 in-vitro positives: {n_t2_total:,} (not in census, flagged separately)\n")
    lines.append(f"- Tier-2 decoy_intrinsic: {n_dec_i:,}   "
                 f"Tier-1 decoy_intergenic: {n_dec_g:,}\n\n")
    lines.append("| organism | phylum | T1 windows | clusters | intra-genome gate (≥30) |\n")
    lines.append("|---|---|---|---|---|\n")
    for o in sorted(by_org):
        npos = sum(1 for r in t1_pos if r["organism"] == o)
        phy  = next(phy_of[m] for c in by_org[o] for m in clus[c])
        ok   = "PASS" if per_org_clusters[o] >= GATE_CLUSTERS else "**DEMOTE -> case study**"
        lines.append(f"| {o} | {phy} | {npos:,} | {per_org_clusters[o]:,} | {ok} |\n")
    lines.append(f"\n**Cross-genome unseen claim "
                 f"(≥{GATE_CLUSTERS} clusters, ≥{GATE_PHYLA} in-vivo phyla): "
                 f"{'PASS' if cross_ok else 'FAIL -> intra-genome only'}**\n")
    lines.append(f"({n_clusters:,} clusters; {n_phyla} in-vivo phyla)\n\n")
    lines.append("## Honest scope\n\n")
    lines.append("Rho cross-genome claim = 2-organism transfer "
                 "(E. coli BCM/Term-seq ↔ MTB RhoDUC), not a broad multi-phylum sweep.\n")
    lines.append("B. subtilis H-SELEX (T2 in-vitro) is supporting evidence only.\n")
    lines.append("Primary Tier-2 decoy: Rho vs intrinsic terminators (TERMITe, E. coli + B. subtilis).\n")

    report = "".join(lines)
    with open(os.path.join(OUT, "census.md"), "w") as fh:
        fh.write(report)
    print("\n" + report)
    print(f"Committed artifacts -> {OUT}/  "
          f"(clusters.tsv, split_80_20.tsv, folds_loo.tsv, census.md)")

if __name__ == "__main__":
    main()
