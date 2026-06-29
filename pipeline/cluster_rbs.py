#!/usr/bin/env python3
"""
§5f data-census gate + §5a/§5b sequence-dissimilar split, for the RBS windows.

1. Cluster all TIS windows (positive_SD + positive_UNSD) at ≤60% identity
   (MMseqs2 easy-cluster) so near-duplicate / homologous starts fall in one cluster.
2. CENSUS GATE (§5f): count independent ≤60%-identity clusters per organism / phylum.
   Pre-registered thresholds (RBS-specific):
     - cross-genome unseen claim: ≥ 30 held-out clusters spanning ≥ 4 bacterial phyla
                                   + an archaeal domain-transfer fold (H. volcanii)
     - intra-genome claim:        ≥ 30 clusters per organism
   Below threshold → demoted to "case study, no generalization claim".
3. SPLIT (committed, §14): 80/20 at the CLUSTER level (so a validation window never
   has a ≤60% near-duplicate in train), stratified by organism.
4. LOO FOLDS for the cross-genome unseen test (P5, P8).
   NOTE: E. coli K-12 (ecoli_K12_MG1655) and BL21 (ecoli_BL21_DE3) share
   near-identical core genes and will fall in the same MMseqs2 clusters.
   They MUST be co-held-out in the Gammaproteobacteria LOO fold to prevent leakage.

Tier-1 decoys (decoy_intergenic) are split by a simple per-organism 80/20 random
split (same as the promoter pipeline — negatives are less leakage-critical).

Outputs (committed): splits/rbs/{clusters.tsv, split_80_20.tsv, folds_loo.tsv, census.md}
Run: .venv/bin/python pipeline/cluster_rbs.py
"""
import csv, os, random, subprocess
from collections import defaultdict, Counter

random.seed(20260625)
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DS   = os.path.join(ROOT, "data/datasets/rbs")
MM   = os.path.join(DS, "mmseqs")
OUT  = os.path.join(ROOT, "splits/rbs")
ALL  = os.path.join(DS, "ALL.tsv")

MIN_ID       = 0.60
VAL_FRAC     = 0.20
GATE_CLUSTERS = 30   # pre-registered minimum per §5f
# RBS gate: ≥4 bacterial phyla + archaea (H. volcanii) = 5 domain groups.
# The census checks n_phyla in the cluster label, where phylum includes "Archaea".
GATE_PHYLA   = 4     # minimum distinct bacterial phyla (archaea counted separately)

# Organisms that must be co-held-out (near-identical core genes → leakage risk)
ECOLI_GROUP  = {"ecoli_K12_MG1655", "ecoli_BL21_DE3"}

def load_rows():
    with open(ALL) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))

def write_fasta(rows, path):
    """Write all TIS positives (SD + UNSD) as FASTA for clustering."""
    with open(path, "w") as fh:
        for r in rows:
            if r["label"].startswith("positive_"):
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
    """member_id -> cluster_rep_id"""
    m = {}
    with open(tsv) as fh:
        for ln in fh:
            rep, mem = ln.rstrip("\n").split("\t")
            m[mem] = rep
    return m

def main():
    os.makedirs(OUT, exist_ok=True)
    rows   = load_rows()
    pos    = [r for r in rows if r["label"].startswith("positive_")]
    org_of = {r["id"]: r["organism"] for r in rows}
    phy_of = {r["id"]: r["phylum"]   for r in rows}

    fasta = os.path.join(DS, "positives.fasta")
    write_fasta(rows, fasta)
    print(f"Clustering {len(pos):,} RBS positives at ≤{int(MIN_ID*100)}% identity ...")
    ctsv    = run_mmseqs(fasta, os.path.join(MM, "rbs"))
    mem2rep = parse_clusters(ctsv)

    clus   = defaultdict(list)
    for mid, rep in mem2rep.items():
        clus[rep].append(mid)
    cl_org = {c: set(org_of[m] for m in members) for c, members in clus.items()}
    cl_phy = {c: set(phy_of[m] for m in members) for c, members in clus.items()}

    n_clusters       = len(clus)
    singletons       = sum(1 for c in clus if len(clus[c]) == 1)
    multi            = n_clusters - singletons
    cross_org_clus   = sum(1 for c in cl_org if len(cl_org[c]) > 1)

    all_phyla = sorted(set(phy_of.values()))
    bact_phyla = [p for p in all_phyla if p != "Archaea"]
    n_phyla   = len(all_phyla)
    has_archaea = "Archaea" in all_phyla

    per_org_clusters = Counter()
    for c, orgs in cl_org.items():
        for o in orgs:
            per_org_clusters[o] += 1

    # ── census gate ──────────────────────────────────────────────────────────
    cross_ok = (n_clusters >= GATE_CLUSTERS and
                len(bact_phyla) >= GATE_PHYLA and has_archaea)
    intra_ok = {o: per_org_clusters[o] >= GATE_CLUSTERS for o in per_org_clusters}

    # ── splits ───────────────────────────────────────────────────────────────
    # E. coli K-12 and BL21 are grouped together (same Gammaproteobacteria fold).
    # For stratified splitting, treat them as a single group ("ecoli_group").
    def split_org_key(org):
        return "ecoli_group" if org in ECOLI_GROUP else org

    cl_primary_org = {c: Counter(org_of[m] for m in clus[c]).most_common(1)[0][0]
                      for c in clus}
    by_group = defaultdict(list)
    for c in clus:
        by_group[split_org_key(cl_primary_org[c])].append(c)

    split_of_cluster = {}
    for grp, cs in by_group.items():
        random.shuffle(cs)
        nval = max(1, round(len(cs) * VAL_FRAC))
        for i, c in enumerate(cs):
            split_of_cluster[c] = "val" if i < nval else "train"

    # Write clusters.tsv
    with open(os.path.join(OUT, "clusters.tsv"), "w") as fh:
        fh.write("window_id\tcluster_rep\torganism\tphylum\tlabel\n")
        for mid, rep in sorted(mem2rep.items()):
            row = next(r for r in rows if r["id"] == mid)
            fh.write(f"{mid}\t{rep}\t{org_of[mid]}\t{phy_of[mid]}\t{row['label']}\n")

    # Write split_80_20.tsv
    with open(os.path.join(OUT, "split_80_20.tsv"), "w") as fh:
        fh.write("window_id\tlabel\torganism\tsplit\n")
        for r in pos:
            fh.write(f"{r['id']}\t{r['label']}\t{r['organism']}\t"
                     f"{split_of_cluster[mem2rep[r['id']]]}\n")
        # Tier-1 decoys: simple per-organism 80/20
        dec_by_org = defaultdict(list)
        for r in rows:
            if r["label"] == "decoy_intergenic":
                dec_by_org[r["organism"]].append(r["id"])
        for o, ids in dec_by_org.items():
            random.shuffle(ids)
            nval = round(len(ids) * VAL_FRAC)
            for i, did in enumerate(ids):
                fh.write(f"{did}\tdecoy_intergenic\t{o}\t{'val' if i < nval else 'train'}\n")

    # Write folds_loo.tsv
    # LOO grouped: E. coli K-12 + BL21 held out together as "ecoli_group"
    loo_groups = {}
    for org in per_org_clusters:
        grp = split_org_key(org)
        loo_groups.setdefault(grp, {"orgs": set(), "phyla": set(), "n_clusters": 0})
        loo_groups[grp]["orgs"].add(org)
        loo_groups[grp]["phyla"].update(
            phy_of[m] for c in by_group[grp] for m in clus[c]
        )
        loo_groups[grp]["n_clusters"] = len(by_group[grp])

    with open(os.path.join(OUT, "folds_loo.tsv"), "w") as fh:
        fh.write("# leave-one-organism-out folds for the cross-genome unseen test (P5, P8)\n")
        fh.write("# E. coli K-12 + BL21 must be co-held-out (near-identical core genes → leakage)\n")
        fh.write("loo_group\theld_out_organisms\tphyla\tn_clusters\n")
        for grp, info in sorted(loo_groups.items()):
            orgs  = ";".join(sorted(info["orgs"]))
            phyla = ";".join(sorted(info["phyla"]))
            fh.write(f"{grp}\t{orgs}\t{phyla}\t{info['n_clusters']}\n")

    # ── census report ────────────────────────────────────────────────────────
    lines = []
    lines.append("# RBS data-census gate (§5f)\n\n")
    lines.append(f"- clustering: MMseqs2 easy-cluster, min-seq-id={MIN_ID}, -c=0.8\n")
    lines.append(f"- positives: **{len(pos):,}** TIS windows "
                 f"-> **{n_clusters:,}** independent ≤{int(MIN_ID*100)}%-identity clusters "
                 f"({singletons:,} singletons, {multi:,} multi-member, "
                 f"{cross_org_clus:,} cross-organism)\n")
    lines.append(f"- phyla / domains: **{n_phyla}** ({', '.join(all_phyla)})\n")
    lines.append(f"  - bacterial phyla: {len(bact_phyla)} ({', '.join(bact_phyla)})\n")
    lines.append(f"  - archaea domain: {'present' if has_archaea else 'absent'} "
                 f"(H. volcanii DS2)\n\n")
    lines.append("| organism | phylum | TIS windows | clusters | intra-genome gate (≥30) |\n")
    lines.append("|---|---|---|---|---|\n")
    for o in sorted(by_group):
        if o == "ecoli_group":
            for real_org in sorted(ECOLI_GROUP):
                npos = sum(1 for r in pos if r["organism"] == real_org)
                phy  = next(phy_of[m] for c in by_group[o] for m in clus[c]
                            if org_of[m] == real_org)
                ok   = "PASS (grouped with ecoli_group)" if per_org_clusters[real_org] >= GATE_CLUSTERS \
                       else "**DEMOTE**"
                lines.append(f"| {real_org} | {phy} | {npos:,} | "
                              f"{per_org_clusters[real_org]:,} | {ok} |\n")
        else:
            npos = sum(1 for r in pos if r["organism"] == o)
            phy  = next(phy_of[m] for c in by_group[o] for m in clus[c])
            ok   = "PASS" if per_org_clusters[o] >= GATE_CLUSTERS else "**DEMOTE -> case study**"
            lines.append(f"| {o} | {phy} | {npos:,} | {per_org_clusters[o]:,} | {ok} |\n")
    lines.append(f"\n**Cross-genome unseen claim "
                 f"(≥{GATE_CLUSTERS} clusters, ≥{GATE_PHYLA} bacterial phyla + archaea): "
                 f"{'PASS' if cross_ok else 'FAIL -> intra-genome only'}**\n")
    lines.append(f"({n_clusters:,} clusters; {len(bact_phyla)} bacterial phyla; "
                 f"archaea: {'yes' if has_archaea else 'no'})\n\n")
    lines.append("## SD / UNSD breakdown\n\n")
    n_sd   = sum(1 for r in pos if r["label"] == "positive_SD")
    n_unsd = sum(1 for r in pos if r["label"] == "positive_UNSD")
    lines.append(f"- positive_SD   : {n_sd:,} ({100*n_sd/len(pos):.1f}%)\n")
    lines.append(f"- positive_UNSD : {n_unsd:,} ({100*n_unsd/len(pos):.1f}%)\n")
    lines.append("- Tier-2 SD-detection task: positive_SD (pos) vs positive_UNSD (neg)\n")
    lines.append("- TIS-detection task: both SD+UNSD (pos) vs decoy_intergenic (neg)\n")

    report = "".join(lines)
    with open(os.path.join(OUT, "census.md"), "w") as fh:
        fh.write(report)
    print("\n" + report)
    print(f"Committed artifacts -> {OUT}/  "
          f"(clusters.tsv, split_80_20.tsv, folds_loo.tsv, census.md)")

if __name__ == "__main__":
    main()
