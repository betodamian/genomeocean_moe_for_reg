# Rho data-census gate (§5f)

- clustering: MMseqs2 easy-cluster, min-seq-id=0.6, -c=0.8
- T1 in-vivo positives: **893** (E. coli=151, MTB=742)
  -> **883** ≤60%-identity clusters (873 singletons, 0 cross-organism)
- in-vivo phyla: **2** (Actinobacteria, Gammaproteobacteria)
- T2 in-vitro positives: 4,789 (not in census, flagged separately)
- Tier-2 decoy_intrinsic: 3,456   Tier-1 decoy_intergenic: 5,682

| organism | phylum | T1 windows | clusters | intra-genome gate (≥30) |
|---|---|---|---|---|
| ecoli_K12_MG1655 | Gammaproteobacteria | 151 | 151 | PASS |
| mtuberculosis_H37Rv | Actinobacteria | 742 | 732 | PASS |

**Cross-genome unseen claim (≥30 clusters, ≥2 in-vivo phyla): PASS**
(883 clusters; 2 in-vivo phyla)

## Honest scope

Rho cross-genome claim = 2-organism transfer (E. coli BCM/Term-seq ↔ MTB RhoDUC), not a broad multi-phylum sweep.
B. subtilis H-SELEX (T2 in-vitro) is supporting evidence only.
Primary Tier-2 decoy: Rho vs intrinsic terminators (TERMITe, E. coli + B. subtilis).
