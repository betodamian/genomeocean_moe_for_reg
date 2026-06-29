# RBS data-census gate (§5f)

- clustering: MMseqs2 easy-cluster, min-seq-id=0.6, -c=0.8
- positives: **20,696** TIS windows -> **17,695** independent ≤60%-identity clusters (15,413 singletons, 2,282 multi-member, 1,648 cross-organism)
- phyla / domains: **5** (Actinobacteria, Alphaproteobacteria, Archaea, Firmicutes, Gammaproteobacteria)
  - bacterial phyla: 4 (Actinobacteria, Alphaproteobacteria, Firmicutes, Gammaproteobacteria)
  - archaea domain: present (H. volcanii DS2)

| organism | phylum | TIS windows | clusters | intra-genome gate (≥30) |
|---|---|---|---|---|
| bsubtilis_168 | Firmicutes | 4,344 | 4,279 | PASS |
| ccrescentus_NA1000 | Alphaproteobacteria | 3,848 | 3,724 | PASS |
| ecoli_BL21_DE3 | Gammaproteobacteria | 3,152 | 2,669 | PASS (grouped with ecoli_group) |
| ecoli_K12_MG1655 | Gammaproteobacteria | 2,164 | 1,928 | PASS (grouped with ecoli_group) |
| hvolcanii_DS2 | Archaea | 1,542 | 1,485 | PASS |
| mtuberculosis_H37Rv | Actinobacteria | 3,508 | 3,387 | PASS |
| saureus_HG001 | Firmicutes | 2,138 | 1,871 | PASS |

**Cross-genome unseen claim (≥30 clusters, ≥4 bacterial phyla + archaea): PASS**
(17,695 clusters; 4 bacterial phyla; archaea: yes)

## SD / UNSD breakdown

- positive_SD   : 12,866 (62.2%)
- positive_UNSD : 7,830 (37.8%)
- Tier-2 SD-detection task: positive_SD (pos) vs positive_UNSD (neg)
- TIS-detection task: both SD+UNSD (pos) vs decoy_intergenic (neg)
