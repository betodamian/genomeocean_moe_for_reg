# Data Acquisition Manifest

**Last updated:** 2026-06-24
**Scope:** every dataset is sourced strictly from `context/research_plan.md` (§5d, §5f) and the
per-element catalogs (`data/rbs_database/sources.tsv`, `data/rho_database/sources.tsv`).
**Reproducibility (research_plan §14):** data is fetched by the committed scripts in
`data/download/` and pinned by `data/download/checksums.txt`. Bulk data is `.gitignore`d
(kept locally); the scripts + checksums + this manifest are committed.

Status legend: **PULLED** = downloaded locally & checksummed · **PIPELINE** = raw signal pulled
or available, needs the Week-1 TIS/peak-calling pipeline to become labels · **MANUAL** = behind a
SPA/journal SI, fetch by hand (URL given).

---

## 1. Genome sequences + structural annotation (NCBI RefSeq)
Source: **NCBI Datasets API v2** (https://www.ncbi.nlm.nih.gov/datasets/). Role per §5f:
sequence substrate + structural scaffold (NOT a regulatory-label source). Script:
`data/download/download_genomes.sh`.

| Accession | Organism | Role (research_plan) | Status |
|---|---|---|---|
| GCF_000005845.2 | *E. coli* K-12 MG1655 | core; promoters/RBS/Rho | **PULLED** |
| GCF_000009045.1 | *B. subtilis* 168 | core; promoters/RBS; Rho (in vitro) | **PULLED** |
| GCF_000195955.2 | *M. tuberculosis* H37Rv | core; RBS leaderless; **Rho 2nd phylum** | **PULLED** |
| GCF_000006765.1 | *P. aeruginosa* PAO1 | core; RBS (PGAP-derived T2) | **PULLED** |
| GCF_000013425.1 | *S. aureus* NCTC 8325 | core; RBS extended-SD | **PULLED** |
| GCF_000022005.1 | *C. crescentus* NA1000 | RBS 4th phylum (Alphaproteobacteria) | **PULLED** |
| GCF_000025685.1 | *H. volcanii* DS2 | RBS archaeal domain-transfer | **PULLED** |

Each: `<accession>_*_genomic.fna` + `<accession>_genomic.gff` in `data/genomes/<label>/`.

## 2. σ-Promoters / TSS (§5d)
Script: `data/download/download_promoters.sh`.

| Source | Citation / accession | Provides | Status |
|---|---|---|---|
| **PPD** (Prokaryotic Promoter Database) | Liu et al. 2021, *J Mol Biol* 433:166860; http://lin-group.cn/database/ppd/ | **129,148** experimental promoters, 29 species CSVs pulled (incl. E. coli, B. subtilis, S. aureus MW2, H. volcanii, Campylobacter, Synechocystis…); cols: TSSPosition, Strand, 81-bp PromoterSeq | **PULLED** |
| RegulonDB v12 | Tierrafría et al. 2024, *NAR* 52:D255 (PMC10767902); https://regulondb.ccg.unam.mx | E. coli σ-factor sub-typing + strong-evidence cross-check | **MANUAL** (SPA; Downloads tab) |
| dRNA-seq TSS catalogs | MTB Cortes 2013 (PMC3898074); Salmonella Kröger 2013 (PMC3963941); Campylobacter Dugar 2013 (PMC3656092); H. pylori Sharma 2010; Synechocystis; B. subtilis BSGatlas | supplement PPD per organism | **MANUAL** (per-paper SRA/SI) |
| *excluded* CDBProm | Pérez-Rueda 2024 (predicted) | — software-predicted, barred by §12 | n/a |

> **§12 anti-circularity:** before use, filter PPD to experimental evidence and drop any entry
> traceable to a benchmarked tool (BPROM, bTSSfinder, iPromoter-2L, G4PromFinder).

## 3. RBS / Shine-Dalgarno (§5d; full catalog `data/rbs_database/sources.tsv`)
Script: `data/download/download_labels.sh`. Per-gene TIS *class* tables (SD/UNSD/leaderless)
live in journal SI → **PIPELINE/MANUAL**; GEO supplementary signal pulled where small.

| Organism | Source / accession | Status |
|---|---|---|
| *E. coli* Ribo-RET | Meydan et al. 2019 *Mol Cell*, PMID 30904393, PMC7115971 | Tables S1/S2 xlsx downloaded manually 2026-06-25; re-parsed 2026-06-26 into unified schema → `ecoli_riboret_tis_sites.tsv` (**2,296** BW25113 K-12 TIS). NOTE: dropped 3,608 BL21 rows (E. coli B strain — coords do not map to NC_000913.3); BW25113 coords need gene-name liftover to MG1655 (flagged per row). **PULLED** |
| *E. coli* TetRP | Nakahigashi 2016, DDBJ **PRJDB2960** | **PIPELINE** (raw reads) |
| *E. coli* ΔaSD | Saito 2020, GEO **GSE135906** | wig signal **PULLED** (`GSE135906_62wigfiles.tar.gz`) → PIPELINE |
| *M. tuberculosis* Ribo-seq | Sawyer et al. 2021 *Cell Reports*, PMID 33535039, PMC7856553. DOI 10.1016/j.celrep.2021.108695 | mmc2/mmc3/mmc4.xlsx auto-fetched via EuropePMC. **PULLED** (note: prior entry listed Zhu 2021 / PMID 33513356 — wrong paper, corrected 2026-06-25) |
| *M. tuberculosis* dRNA-seq leaderless | Cortes 2013, SRA **SRP028740** | **PIPELINE** |
| *B. subtilis* Ribo-seq | Lalanne et al. 2018 *Cell*, PMID 29606352, PMC5978003. DOI 10.1016/j.cell.2018.03.007 | TableS1 xlsx downloaded manually 2026-06-25; re-parsed 2026-06-26 into unified schema → `bsub_riboseq_tis_sites.tsv` (4,176 TIS positions, strand-aware tis_pos, low-reads flag). **PULLED** (note: prior entry listed Lalanne 2017 / PMID 29144454 — wrong paper, corrected 2026-06-25; file renamed from bsub_riboret→bsub_riboseq) |
| *B. subtilis* sporulation | Iwańska et al. 2024 *Nat Commun*, PMCID PMC11339384. DOI 10.1038/s41467-024-51654-6 | MOESM xlsx auto-fetched via EuropePMC. **PULLED** (note: prior entry listed Bhatt 2024 / PMID 39179838 — wrong paper, corrected 2026-06-25) |
| *S. aureus* extended-SD (sORF subset) | Kohl 2026, *Nat Commun* PMID 41680142, PMC13009471. DOI 10.1038/s41467-026-69079-8 | MOESM3.xlsx auto-fetched 2026-06-25 → `saur_exsd_tis_sites.tsv` (46 novel sORF TIS, HG001). **PULLED** (subset only) |
| *S. aureus* Ribo-RET (full TIS signal) | Kohl 2026 [Ribo-Seq], GEO **GSE299221** (+RNA-seq GSE299222) | `GSE299221_RAW.tar` (8.8 MB, 16 wig tracks: RNase1/MNase × Ret/ctl × 2 reps) pulled 2026-06-26 to `raw/SAUR_RIBORET/`. fixedStep wig, **chrom=HG001** (2,819,767 bp ≈ NCTC 8325 2,821,361, ~1.6 kb diff → near-trivial liftover). Translatome-wide TIS (~2,700 expected) requires Week-1 peak-calling on Ret-vs-ctl at start codons. **PIPELINE** |
| *C. crescentus* | Schrader 2014, GEO **GSE54883** | **PIPELINE** |
| *H. volcanii* (archaea) | Gelsinger 2020, *NAR* 48:5201, PMID 32382758, PMC7261190. DOI 10.1093/nar/gkaa210 | Ribo_MS_TableS1_final.xlsx auto-fetched via EuropePMC 2026-06-25; parsed → `hvolc_riboseq_tis_sites.tsv` (1,555 rows: 1,413 annotated TIS + 142 novel/extension). **PULLED** |
| *P. aeruginosa* | PGAP-derived from GCF_000006765.1 GFF (T2) | derivable from pulled genome |

## 4. Rho-dependent terminators (§5d; full catalog `data/rho_database/sources.tsv`)

| Organism | Source / accession | Status |
|---|---|---|
| *E. coli* Term-seq 3′-ends | NAR 2018 46:6797, GEO **GSE109766** | per-position signal **PULLED** (`GSE109766_…counts_per_position.txt.gz`) → PIPELINE |
| *E. coli* BCM Rho | Peters et al. 2012, *Genes Dev* 26:2621, PMID 23207917, PMC3521622. DOI 10.1101/gad.196741.112; GEO **GSE41936** (ChIP-chip) + **GSE41939** (RNA-seq BCM) | `Supplemental_tables.xls` downloaded manually 2026-06-26; parsed → `ecoli_bcm_rho_sites.tsv` (1,264 BST sites: BCM-significant transcripts = Rho-dependent regions, incl. sense + antisense). NOTE: research plan incorrectly cited as "PNAS 109:15584" — correct journal is Genes & Development. **PULLED** |
| *M. tuberculosis* RhoDUC | Botella et al. 2022/2023 *iScience*, PMID 37096044, PMC10122055. DOI 10.1016/j.isci.2023.106465; ArrayExpress **E-MTAB-11753** (raw BAMs only) | mmc5.xlsx downloaded manually 2026-06-26; parsed → `mtb_rhoduc2_sites.tsv` (439 high-confidence RD-TTS: 299 True + 125 Cond + 15 SecCond). NOTE: research plan said "1,385" — mmc5 has 439 high-confidence; mmc4 has 802 total non-intrinsic TTS. Also have Botella 2017 (303 RSRs) as separate source. **PULLED** |
| *B. subtilis* H-SELEX (in vitro) | PMC12350095 (600 BsRho rut) | Auto-fetched via EuropePMC 2026-06-24; parsed → `bsub_hselex_sites.tsv` (4,789 rows). **PULLED** |
| Intrinsic decoy (TERMITe) | PMC12207403, Zenodo + GitHub (13 species) | Auto-fetched via EuropePMC 2026-06-24; parsed → `intrinsic_termite_sites.tsv` (5,646 rows). **PULLED** |
| RhoTermPredict (cross-check only) | PMC6407284 (algorithm-derived; never headline, §12) | Auto-fetched via EuropePMC 2026-06-24; parsed → `rhotermpredict_sites.tsv` (23,976 rows). **PULLED** |
| *M. tuberculosis* Rho (Botella 2017) | Botella et al. 2017 *Nat Commun*, PMID 28348398, PMC5379054. DOI 10.1038/ncomms14731 | Auto-fetched via EuropePMC 2026-06-24; parsed → `mtb_rhoduc_sites.tsv` (303 RSRs). **PULLED** (different from Botella 2022; included as T1 in vivo supplementary) |

---

## What is local now (updated 2026-06-25)

**Genomes:** `data/genomes/` — 7 genomes (FASTA+GFF): E. coli, B. subtilis, MTB, P. aeruginosa, S. aureus, C. crescentus, H. volcanii.

**Promoters:** `data/promoters/PPD/` — 129,149 promoters; windows built in `data/datasets/promoters/ALL.tsv` (28,098 rows); 80/20 split + census committed to `splits/promoters/`.

**RBS TSVs** (gitignored, in `data/rbs_database/raw/<SOURCE>/`):
- `ecoli_riboret_tis_sites.tsv` — 2,296 rows (Meydan 2019; BW25113 K-12 only, BL21 dropped)
- `ecoli_tetrp_tis_sites.tsv` — 249 rows (Nakahigashi 2016)
- `mtb_riboseq_tis_sites.tsv` — 3,569 rows (Sawyer 2021, 98% GFF-matched)
- `bsub_riboseq_tis_sites.tsv` — 4,176 rows (Lalanne 2018)
- `bsub_spore_tis_sites.tsv` — 4,332 rows (Iwańska 2024)
- `caulo_riboseq_tis_sites.tsv` — 3,885 rows (Schrader 2014)
- `saur_exsd_tis_sites.tsv` — 46 rows (Kohl 2026, novel sORFs only)
- `hvolc_riboseq_tis_sites.tsv` — 1,555 rows (Gelsinger 2020: 1,413 annotated + 142 novel/extension)

**Rho TSVs** (gitignored, in `data/rho_database/raw/<SOURCE>/`):
- `mtb_rhoduc_sites.tsv` — 303 RSRs (Botella 2017)
- `bsub_hselex_sites.tsv` — 4,789 rut peaks (B. subtilis H-SELEX, in vitro T2)
- `intrinsic_termite_sites.tsv` — 5,646 sites (TERMITe, E. coli + B. subtilis only)
- `rhotermpredict_sites.tsv` — 23,976 predictions (cross-check only, never headline)

**All experimental label sources now PULLED.** Ready to proceed to window-building.

## Next (Week 1, per research_plan §16)
- Manual downloads above (Peters 2012, Botella 2022)
- Build RBS 300-bp windows + same-context decoys → `data/datasets/rbs/`
- Run RBS census gate + 80/20 ≤60%-identity split
- Build Rho windows + census + split → `data/datasets/rho/`
