# Promoter data-census gate (§5f)
- clustering: MMseqs2 easy-cluster, min-seq-id=0.6, -c 0.8
- positives: **14049** windows -> **9454** independent <= 60%-identity clusters (7039 singletons, 2415 multi-member, 0 cross-organism)
- phyla represented: **3** (Archaea, Firmicutes, Gammaproteobacteria)

| organism | phylum | positives | clusters | intra-genome gate (>=30) |
|---|---|---|---|---|
| bsubtilis_168 | Firmicutes | 690 | 629 | PASS |
| ecoli_K12_MG1655 | Gammaproteobacteria | 8614 | 4710 | PASS |
| hvolcanii_DS2 | Archaea | 4745 | 4115 | PASS |

**Cross-genome unseen claim (>=30 clusters spanning >=2 phyla): PASS** (9454 clusters, 3 phyla)
