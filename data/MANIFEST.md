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
| *E. coli* Ribo-RET | Meydan 2019, GEO **GSE122129** | TIS table in journal SI → **MANUAL**; GEO signal **PIPELINE** |
| *E. coli* TetRP | Nakahigashi 2016, DDBJ **PRJDB2960** | **PIPELINE** (raw reads) |
| *E. coli* ΔaSD | Saito 2020, GEO **GSE135906** | wig signal **PULLED** (`GSE135906_62wigfiles.tar.gz`) → PIPELINE |
| *M. tuberculosis* Ribo-seq | Zhu 2021, ArrayExpress **E-MTAB-8835** | per-gene SD/UNSD/leaderless table in SI → **MANUAL** |
| *M. tuberculosis* dRNA-seq leaderless | Cortes 2013, SRA **SRP028740** | **PIPELINE** |
| *B. subtilis* Ribo-seq | Lalanne 2017, GEO **GSE95211** | **PIPELINE** |
| *B. subtilis* sporulation | Bhatt 2024, GEO **GSE249450** | **PIPELINE** |
| *S. aureus* extended-SD | Kohl 2026, *Nat Commun* 10.1038/s41467-026-69079-8 | accession in SI → **MANUAL** |
| *C. crescentus* | Schrader 2014, GEO **GSE54883** | **PIPELINE** |
| *H. volcanii* (archaea) | Gelsinger 2020, *NAR* 48:5201 | accession in SI → **MANUAL** |
| *P. aeruginosa* | PGAP-derived from GCF_000006765.1 GFF (T2) | derivable from pulled genome |

## 4. Rho-dependent terminators (§5d; full catalog `data/rho_database/sources.tsv`)

| Organism | Source / accession | Status |
|---|---|---|
| *E. coli* Term-seq 3′-ends | NAR 2018 46:6797, GEO **GSE109766** | per-position signal **PULLED** (`GSE109766_…counts_per_position.txt.gz`) → PIPELINE |
| *E. coli* BCM Rho | Peters 2012, *PNAS* (≈1,000 sites) | site table in SI → **MANUAL** |
| *M. tuberculosis* RhoDUC | Botella 2022, ArrayExpress **E-MTAB-11753** (1,385 sites) | raw BAMs (200–350 MB ea., skipped); site table in SI → **MANUAL** |
| *B. subtilis* H-SELEX (in vitro) | PMC12350095 (600 BsRho rut) | SI → **MANUAL** |
| Intrinsic decoy (TERMITe) | PMC12207403, Zenodo + GitHub (13 species) | **MANUAL** (Zenodo) |
| RhoTermPredict (cross-check only) | PMC6407284 | **MANUAL** (algorithm-derived; never headline, §12) |

---

## What is local now
`data/genomes/` (7 genomes, FASTA+GFF) · `data/promoters/PPD/` (29 species, 129,149 promoters) ·
`data/rbs_database/raw/` (E. coli ΔaSD wig signal) · `data/rho_database/raw/` (E. coli Rho
Term-seq per-position signal). All under `data/download/checksums.txt`.

## Next (Week 1, per research_plan §16)
Run the TIS/peak-calling pipeline on the PIPELINE rows to convert signal → element coordinates;
hand-fetch the MANUAL journal-SI tables (DOIs in the per-element `sources.tsv`); then build
300-bp windows + same-context decoys and the 80/20 ≤60%-identity split (§5b, §5c).
