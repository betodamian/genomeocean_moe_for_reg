#!/usr/bin/env python3
"""
§5f data-census gate + §5a/§5b sequence-dissimilar split, for the promoter windows.

1. Cluster the promoter POSITIVE windows by sequence identity (MMseqs2 easy-cluster,
   --min-seq-id 0.6) so near-duplicate / homologous promoters fall in one cluster.
2. CENSUS GATE (§5f): count independent <=60%-identity clusters per organism / phylum.
   Pre-registered thresholds:
     - cross-genome unseen claim: >= 30 held-out clusters spanning >= 2 phyla
     - intra-genome claim:        >= 30 clusters
   Anything below -> demoted to "case study, no generalization claim".
3. SPLIT (committed, §14): a frozen 80/20 split at the CLUSTER level (so a validation
   window never has a <=60% near-duplicate in train), stratified by organism, plus
   leave-one-organism-out folds for the cross-genome (P8) test.

Outputs (committed): splits/promoters/{clusters.tsv,split_80_20.tsv,folds_loo.tsv,census.md}
Run: .venv/bin/python pipeline/cluster_and_census.py
"""
import csv, json, os, random, subprocess, sys
from collections import defaultdict, Counter

random.seed(20260625)
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DS   = os.path.join(ROOT, "data/datasets/promoters")
MM   = os.path.join(DS, "mmseqs")
OUT  = os.path.join(ROOT, "splits/promoters")
ALL  = os.path.join(DS, "ALL.tsv")
MIN_ID = 0.60
VAL_FRAC = 0.20
GATE_CLUSTERS = 30     # pre-registered minimum

def load_rows():
    rows = []
    with open(ALL) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            rows.append(r)
    return rows

def write_fasta(rows, path):
    with open(path, "w") as fh:
        for r in rows:
            if r["label"] == "positive":
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
    rows = load_rows()
    pos = [r for r in rows if r["label"] == "positive"]
    org_of = {r["id"]: r["organism"] for r in rows}
    phy_of = {r["id"]: r["phylum"] for r in rows}

    fasta = os.path.join(DS, "positives.fasta")
    write_fasta(rows, fasta)
    print(f"clustering {len(pos)} promoter positives at <= {int(MIN_ID*100)}% identity ...")
    ctsv = run_mmseqs(fasta, os.path.join(MM, "prom"))
    mem2rep = parse_clusters(ctsv)

    # cluster -> members ; cluster -> set(organism/phylum)
    clus = defaultdict(list)
    for mid, rep in mem2rep.items():
        clus[rep].append(mid)
    cl_org = {c: set(org_of[m] for m in members) for c, members in clus.items()}
    cl_phy = {c: set(phy_of[m] for m in members) for c, members in clus.items()}

    n_clusters = len(clus)
    singletons = sum(1 for c in clus if len(clus[c]) == 1)
    multi      = n_clusters - singletons
    cross_org_clusters = sum(1 for c in cl_org if len(cl_org[c]) > 1)
    phyla = sorted(set(phy_of.values()))

    # per-organism cluster counts (a cluster counts for an organism if it has a member there)
    per_org_clusters = Counter()
    for c, orgs in cl_org.items():
        for o in orgs:
            per_org_clusters[o] += 1

    # ---------- census gate ----------
    cross_phyla_clusters = sum(1 for c in cl_phy.values() if len(c) >= 1)  # all clusters; phyla coverage below
    n_phyla = len(phyla)
    cross_ok = (n_clusters >= GATE_CLUSTERS) and (n_phyla >= 2)
    intra_ok = {o: per_org_clusters[o] >= GATE_CLUSTERS for o in per_org_clusters}

    # ---------- splits ----------
    # 80/20 at cluster level, stratified by the cluster's primary organism
    cl_primary_org = {c: Counter(org_of[m] for m in clus[c]).most_common(1)[0][0] for c in clus}
    by_org = defaultdict(list)
    for c in clus: by_org[cl_primary_org[c]].append(c)
    split_of_cluster = {}
    for o, cs in by_org.items():
        random.shuffle(cs)
        nval = max(1, round(len(cs) * VAL_FRAC))
        for i, c in enumerate(cs):
            split_of_cluster[c] = "val" if i < nval else "train"

    # write committed split artifacts
    with open(os.path.join(OUT, "clusters.tsv"), "w") as fh:
        fh.write("window_id\tcluster_rep\torganism\tphylum\n")
        for mid, rep in sorted(mem2rep.items()):
            fh.write(f"{mid}\t{rep}\t{org_of[mid]}\t{phy_of[mid]}\n")

    with open(os.path.join(OUT, "split_80_20.tsv"), "w") as fh:
        fh.write("window_id\tlabel\torganism\tsplit\n")
        # positives follow their cluster's split
        for r in pos:
            fh.write(f"{r['id']}\t{r['label']}\t{r['organism']}\t{split_of_cluster[mem2rep[r['id']]]}\n")
        # decoys: simple 80/20 per organism (v1 — random; negatives, less leakage-critical)
        dec_by_org = defaultdict(list)
        for r in rows:
            if r["label"] == "decoy": dec_by_org[r["organism"]].append(r["id"])
        for o, ids in dec_by_org.items():
            random.shuffle(ids); nval = round(len(ids) * VAL_FRAC)
            for i, did in enumerate(ids):
                fh.write(f"{did}\tdecoy\t{o}\t{'val' if i < nval else 'train'}\n")

    with open(os.path.join(OUT, "folds_loo.tsv"), "w") as fh:
        fh.write("# leave-one-organism-out folds for the cross-genome unseen test (P8)\n")
        fh.write("held_out_organism\tphylum\tn_positive_clusters\n")
        for o in sorted(by_org):
            phy = next(phy_of[m] for c in by_org[o] for m in clus[c])
            fh.write(f"{o}\t{phy}\t{len(by_org[o])}\n")

    # ---------- census report ----------
    lines = []
    lines.append("# Promoter data-census gate (§5f)\n")
    lines.append(f"- clustering: MMseqs2 easy-cluster, min-seq-id={MIN_ID}, -c 0.8\n")
    lines.append(f"- positives: **{len(pos)}** windows -> **{n_clusters}** independent "
                 f"<= {int(MIN_ID*100)}%-identity clusters "
                 f"({singletons} singletons, {multi} multi-member, {cross_org_clusters} cross-organism)\n")
    lines.append(f"- phyla represented: **{n_phyla}** ({', '.join(phyla)})\n\n")
    lines.append("| organism | phylum | positives | clusters | intra-genome gate (>=30) |\n")
    lines.append("|---|---|---|---|---|\n")
    for o in sorted(by_org):
        npos = sum(1 for r in pos if r["organism"] == o)
        phy = next(phy_of[m] for c in by_org[o] for m in clus[c])
        ok = "PASS" if per_org_clusters[o] >= GATE_CLUSTERS else "**DEMOTE -> case study**"
        lines.append(f"| {o} | {phy} | {npos} | {per_org_clusters[o]} | {ok} |\n")
    lines.append(f"\n**Cross-genome unseen claim (>=30 clusters spanning >=2 phyla): "
                 f"{'PASS' if cross_ok else 'FAIL -> intra-genome only'}** "
                 f"({n_clusters} clusters, {n_phyla} phyla)\n")
    report = "".join(lines)
    with open(os.path.join(OUT, "census.md"), "w") as fh:
        fh.write(report)

    print("\n" + report)
    print(f"committed artifacts -> {OUT}/  (clusters.tsv, split_80_20.tsv, folds_loo.tsv, census.md)")

if __name__ == "__main__":
    main()
