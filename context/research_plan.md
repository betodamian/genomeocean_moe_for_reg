# Research Plan: Expert-Routed Annotation of Bacterial Regulatory Elements with GenomeOcean-MoE

**Author:** Beto Damian
**Status:** Draft v1 — pre-registration phase
**Builds directly on:** [`junhos_work.md`](junhos_work.md) (Hong, *GenomeOcean: Sparse Upcycling and Expert Specialization in Genomic MoE*, May 2026)

---

## 0. One-paragraph summary

Junho demonstrated that GenomeOcean-MoE develops an **expert-routing channel that is genome-invariant and function-aligned** (routing AMI = 0.15 vs hidden-state AMI = 0.80), that specific experts act as **causally verified** detectors of structural classes (e.g., the upcycled L7 E7 tRNA detector, DiD = +0.559; null in the scratch model), and that **routing fingerprints solve in-domain prokaryotic classification with 8× fewer features than embeddings**. This plan asks whether that same routing channel can *identify and annotate three bacterial regulatory elements* — **σ-dependent promoters, Shine-Dalgarno ribosome binding sites (RBS), and Rho-dependent terminators** — across phylogenetically diverse microbial genomes, beating named classical tools and named dense genomic language models, **and whether the advantage is causally attributable to expert specialization such that a model with no specialized experts cannot reproduce it.** The MoE-necessity claim is the load-bearing scientific contribution and is tested by construction (the primary detector consumes routing features that dense models do not possess) and by causal expert ablation (Junho's gold-standard method extended to regulatory classes).

---

## 1. Focused research question

> **Primary RQ.** Across phylogenetically diverse bacterial and archaeal genomes, do GenomeOcean-MoE's expert-routing fingerprints identify and annotate σ-dependent promoters, Shine-Dalgarno ribosome binding sites, and Rho-dependent terminators *more accurately and more transferably* than (a) classical element-specific tools and (b) dense genomic language models — **and is any advantage causally carried by specialized experts, such that an architecture with no experts cannot achieve the same result?**

This decomposes into three sub-questions that mirror and extend Junho's RQ1 (competitiveness) / RQ2 (specialization):

- **RQ-A — Detection & annotation accuracy.** Does a *routing-aware* detector beat the named baselines on (i) element presence/absence classification, (ii) sub-token-resolution localization, and (iii) sub-typing (sigma class; leadered vs leaderless; Rho-dependent vs intrinsic)?
- **RQ-B — MoE necessity (the core claim).** Is the gain (i) carried by the routing channel *beyond* anything an embedding can express, and (ii) destroyed by ablating the responsible experts? A dense model has neither a routing channel nor experts to ablate, so a positive result here is, by construction, unreachable without specialized experts.
- **RQ-C — Specialization & genome-invariant transfer.** Do dedicated experts emerge for each regulatory class, and do their detectors transfer across phyla (function-aligned, taxonomy-invariant) the way Junho's structural-class experts did — and does the *upcycled* model produce causal regulatory experts where the *scratch* model silently collapses?

### Falsifiable predictions (pre-registered before any test genome is touched)

| ID | Prediction | Falsified if |
|----|-----------|-------------|
| P1 | `routing_concat` (864-d) > `embedding_only` (768-d) on the same MoE backbone for all three elements, paired across leave-one-genome-out folds | No significant paired improvement (BH-FDR q < 0.05) |
| P2 | Best routing-aware detector > every named classical tool and every named dense LGM on held-out genomes (AUPRC, boundary-F1) | A dense LGM or classical tool matches/beats it within bootstrap CI |
| P3 | ≥1 expert per element class shows significant log₂ enrichment over the marginal-masked null; ablating it raises element-token cross-entropy with a Difference-in-Differences (DiD) ≫ 0 vs CDS control | No causal expert exists for any class (DiD CI overlaps 0) |
| P4 | Upcycled regulatory experts are causal (DiD ≫ 0); scratch counterparts are weaker or causally null | Scratch matches upcycled on causality |
| P5 | Routing detectors transfer across phyla better than dense embeddings (smaller held-out-phylum performance drop) | Routing transfer ≤ embedding transfer |
| P6 | Bacterial (in-domain) promoter routing recovers, unlike the **eukaryotic GUE promoter** task where Junho saw routing collapse (routing_only MCC 0.039 vs embedding 0.779) | Routing collapses on in-domain bacterial promoters too |

P6 is a deliberate honesty hook: Junho's *only* promoter result was **out-of-domain eukaryotic** and routing failed there. This plan tests whether **in-domain bacterial** promoters behave differently. We commit to reporting this regardless of outcome.

---

## 2. Why this is an MoE-necessary design (the central argument)

The requirement is that "a model with no specialized experts could not achieve the same result." We satisfy this on two independent grounds:

1. **Architectural exclusivity (by construction).** The *primary* detector's feature vector is the routing fingerprint: per-token router softmax over 12 layers × 8 experts (96-d, pooled) plus per-position expert assignments for localization, optionally concatenated with the residual stream. A dense model **emits no routing distribution** — these features do not exist for it. The dense models can therefore only be evaluated on the embedding-only sub-detector, which is exactly the control that P1 isolates. If `routing_concat` beats `embedding_only`, the surplus is, definitionally, information a dense model cannot supply.

2. **Causal exclusivity (by ablation).** We extend Junho's expert-masking protocol (mask expert logits to −10⁹, measure cross-entropy shift on held-out tokens, compute DiD vs a same-layer control expert) to the three regulatory classes. A significant DiD localizes the capability to a specific expert. **A dense network has no expert to mask**, so the mechanism cannot exist there; and per P4, an MoE *without* successful specialization (the scratch model) is predicted to fail the same test — proving it is *specialized experts*, not mere MoE capacity, that carry the result.

Together these make the MoE the necessary substrate: the detector cannot be instantiated without a routing channel, and the routing channel's regulatory signal is causally tied to identifiable experts that a non-expert model does not and cannot have.

---

## 3. Tools and models to beat (named)

### 3a. Classical, element-specific tools

| Element | Tools the MoE detector must beat | Role |
|---------|----------------------------------|------|
| **Promoters (σ-dependent)** | **BPROM** (Softberry σ70 predictor), **bTSSfinder** (multi-σ TSS/promoter), **G4PromFinder**, **PromoterHunter** (PHISITE), and an ML predictor **iPromoter-2L / MULTiPly** | Primary promoter baselines |
| **Ribosome binding sites (Shine-Dalgarno)** | **Prodigal** (built-in RBS/SD motif scoring), **Salis Lab RBS Calculator v2.1** (biophysical ΔG model), **RBSfinder**, **Free2Bind** (anti-SD base-pairing) | Primary RBS baselines |
| **Rho-dependent terminators** | **RhoTermPredict** (the dedicated Rho-dependent tool), with **TransTermHP**, **ARNold**, **RNIE**, **WebGeSTer** as *intrinsic*-terminator comparators (to test the harder Rho-vs-intrinsic discrimination) | Primary + discrimination baselines |

### 3b. Dense / non-MoE genomic language models (LGMs)

- **GenomeOcean dense family — GO-100M, GO-500M, GO-4B** (the direct dense ablation controls; same data lineage, no experts).
- **Nucleotide Transformer** — NT-2.5B-multispecies and NT-500M.
- **DNABERT-2**.
- **Evo 2** (Arc Institute; prokaryote-heavy training — the strongest in-domain LGM comparator), with **Evo 1** if compute allows.
- **ProkBERT** (prokaryote-specialized; ships promoter tasks — the toughest task-matched comparator).
- **HyenaDNA** and **Caduceus** as long-context / state-space comparators (secondary tier).

All LGMs are evaluated with the **identical frozen-backbone probing protocol** (Section 7) on the **identical held-out loci** to keep the comparison fair; classical tools are run at default *and* threshold-tuned operating points.

---

## 4. Datasets and ground truth (every label traceable and version-pinned)

"Completely verifiable" means every positive/negative label is traceable to a public, version-pinned source, and every number can be regenerated from released code + pinned accessions.

### 4a. Genome panel (diversity backbone — reuse Junho's panels for continuity)
- **High-confidence core (5):** *E. coli* K-12 MG1655, *B. subtilis* 168, *M. tuberculosis* H37Rv, *P. aeruginosa* PAO1, *S. aureus* — the exact genomes from Junho §4, GC 33–67%.
- **Diversity set (20):** Junho's 20 PGAP-annotated genomes (7 phyla + 2 archaea). Stratify all analyses by GC content and phylum.

### 4b. Element ground truth (source → confidence tier)

**Promoters / TSS**
- **RegulonDB** (E. coli K-12): experimentally supported promoters with TSS coordinate, σ factor, and evidence code. Keep only "Confirmed"/strong-evidence entries (primer extension, dRNA-seq, RNA-seq).
- **PRODORIC** and **DBTBS / SubtiWiki** (B. subtilis) for curated multi-σ promoters.
- **dRNA-seq TSS catalogs** for diversity: e.g., *H. pylori* (Sharma et al. 2010), *M. tuberculosis* (Cortes et al. 2013 / Shell et al.), *Synechocystis*, *Campylobacter jejuni*, *Salmonella*, condition-resolved *B. subtilis* (Nicolas et al. 2012).

**Ribosome binding sites (Shine-Dalgarno)**
- Truth anchored on **ribosome-profiling translation-initiation-site (TIS) maps** where available (E. coli Ribo-seq / antibiotic-arrest TIS profiling — Meydan/Weaver et al.); SD region defined as −20..−1 nt relative to confirmed start codons.
- SD strength scored by base-pairing to the organism-specific **anti-SD** (16S rRNA 3′ tail), giving a continuous biophysical target.
- **Leaderless transcripts** (no SD; abundant in Actinobacteria such as *M. tuberculosis*) included as a built-in contrastive negative — this is where genome diversity stresses the detector.

**Rho-dependent terminators**
- E. coli genome-wide Rho-dependent terminator maps from **NET-seq + bicyclomycin (BCM) readthrough** (Peters et al. 2012) and **Term-seq** (Dar et al. 2016 across multiple bacteria).
- **RhoTermPredict** curated positives.
- **Intrinsic terminators** (hairpin + poly-U; TransTermHP/RNIE-labeled) included as the *contrastive* class — Rho-vs-intrinsic discrimination is the discriminating test where motif tools are weakest.

### 4c. Negative sets and confound matching
- Negatives are GC- and length-matched intergenic/genic windows lacking the element, plus **hard decoys**: intrinsic terminators as decoys for Rho terminators; leaderless starts as decoys for RBS; non-promoter intergenic for promoters.
- Confounds explicitly controlled (regression covariates and matched sampling): **GC content, element/window length, genic vs intergenic context, strand, and BPE-boundary phase.** This mirrors Junho's discipline (his marginal-masked P(e) null, and his note that CDS JSD was near-zero purely because CDS dominates token mass).

---

## 5. Coordinate → token mapping (deterministic and released)

Regulatory elements are small relative to GenomeOcean's BPE resolution (~4.89 bp/token): a σ70 −10 box ≈ 6 bp ≈ 1–2 tokens, an SD motif ≈ 6–8 bp ≈ 1–2 tokens, a Rho *rut* site ≈ 60–100 nt ≈ 12–20 tokens.

Steps:
1. Tokenize each genome (both strands) with the **GO-4B BPE tokenizer**; persist a deterministic bp↔token index table per genome (released artifact).
2. Define each element's **core token span** (tokens overlapping the annotated motif) and **flank windows** (±N tokens).
3. Report the fraction of elements that are token-resolvable; for sub-token motifs, attribute the motif to the overlapping token(s) plus immediate neighbors and record the BPE phase as a covariate.
4. Localization metrics are reported in **bp** (de-tokenized) so resolution limits are transparent and comparable to classical tools that operate at bp resolution.

This table is the backbone of reproducibility: identical inputs → identical token spans → identical features.

---

## 6. Feature extraction

For **GenomeOcean-MoE** (upcycled pilot4 as primary; pilot3 and scratch as comparators), extract per locus:
- **Embedding** — 768-d residual stream (mean-pooled over the locus; also positional for localization).
- **Routing-only** — 12×8 router softmax, mean-pooled to **96-d** (Junho's fingerprint), plus **per-position** routing vectors for localization.
- **routing_concat** — 864-d (embedding ⊕ routing).
- **per-expert bins** — 768×k expert-conditioned features (Junho's fourth probe variant; ≥1000-token binning fallback).

For **dense LGMs** (GO-100M/500M/4B, NT, DNABERT-2, Evo 2, ProkBERT, HyenaDNA, Caduceus): **embedding-only** (each model's native hidden state), since no routing channel exists. This asymmetry *is* the experiment.

---

## 7. Detection / annotation models (frozen-backbone, held-out)

Following Junho's frozen-weights probing exactly so results are attributable to the representation, not to fine-tuning:
- Probes: **Logistic Regression** (L2, C=1.0) and **balanced XGBoost** (Junho's two-probe protocol).
- Cross-validation: **leave-one-genome-out (LOGO)**, plus a stricter **leave-one-phylum-out** for transfer (RQ-C). This prevents the taxonomic leakage Junho flagged (bulk CDS splits by lineage/GC).
- Tasks:
  - **(T1) Classification** — element present at locus? (per class, one-vs-negatives + one-vs-hard-decoy).
  - **(T2) Localization** — predict element position within window; report median bp error and boundary-F1.
  - **(T3) Sub-typing** — σ class (promoters), leadered vs leaderless (RBS), Rho-dependent vs intrinsic (terminators).
- Metrics: **MCC, macro-F1, AUPRC** (primary, given class imbalance), boundary-F1, bp localization error. Genome-level **bootstrap CIs** and paired tests across folds.

---

## 8. The MoE-necessity experiments (RQ-B core)

Three converging, independently sufficient lines of evidence:

### 8a. Incremental-value / channel ablation (tests P1, P2)
- On the same MoE backbone, compare `routing_concat` vs `embedding_only` vs `routing_only` vs `per_expert_bins`, paired across LOGO folds (Wilcoxon, BH-FDR). Junho's prior: `routing_concat` was *often better and never worse* in-domain; routing collapsed only on OOD eukaryotic tasks. P6 tests which regime bacterial promoters fall into.
- Then place the MoE `embedding_only` against every dense LGM's embedding to rule out "it's just a better backbone." The MoE claim survives only if **routing_concat > embedding_only AND the surplus features are dense-inaccessible** (true by construction).

### 8b. Expert-detector discovery + causal ablation (tests P3, the gold standard)
1. **Discovery:** for each regulatory class, compute per-(layer, expert) log₂ enrichment on element token-spans vs the **marginal-masked null** (Junho's P(e) ≥ 0.01 baseline), with **G-test + BH-FDR**. Rank detectors as Junho did (his Intergenic L6 E7 +1.75 detector is the most promising prior, since promoters/RBS/Rho-rut all live in or near intergenic/UTR space).
2. **Causality:** mask the candidate expert's logits to −10⁹ on held-out validation tokens; measure cross-entropy increase on **element tokens vs CDS control tokens** → **Difference-in-Differences**. Require DiD ≫ 0 on the target with a near-zero DiD on a same-layer control expert (Junho's L7 E7 tRNA test: DiD +0.559; control L7 E3 DiD −0.065).
3. **Functional consequence:** re-run the T1/T2 detectors with the expert ablated; require a measurable performance drop. **Dense models have no expert to mask — this entire sub-experiment is impossible for them.**

### 8c. Upcycled vs Scratch (tests P4 — specialization, not capacity)
- Repeat 8b for **MoE-Scratch**. Prediction (from Junho's silent-collapse finding: scratch had 14/96 experts at P(e) < 0.01 on OOD genomes, and a causally *null* tRNA detector, DiD +0.006): scratch regulatory experts are weaker or causally null. This isolates *successful specialization* (an upcycling outcome) as the active ingredient — an MoE alone is insufficient, a *specialized* MoE is necessary.

---

## 9. Genome-invariant transfer (RQ-C, tests P5)
- Train detectors on a subset of phyla, test on a held-out phylum. Compare the performance drop for routing features vs dense embeddings.
- Expectation from Junho's AMI contrast (routing 0.15 / function-aligned vs hidden-state 0.80 / taxonomy-aligned): routing detectors should transfer with a smaller drop because they are function-aligned and genome-invariant, whereas dense embeddings carry taxonomic confounds and should transfer worse. A smaller routing drop is a second, independent MoE-advantage axis.

---

## 10. Statistical rigor, controls, and honesty commitments
- **Nulls:** label permutation; marginal-masked routing null; GC/length-matched negatives.
- **Multiple testing:** BH-FDR across all (layer, expert, class, genome) cells, as Junho did.
- **Uncertainty:** genome-level bootstrap CIs; paired fold-level tests; pre-declared effect-size thresholds.
- **Confounds:** GC, length, genic/intergenic context, strand, BPE phase as covariates and in matched sampling.
- **Pre-registered failure criteria:** the falsification rows in §1.
- **Anti-overclaim commitments (in Junho's spirit — he retracted the "1,266 Pfam families" claim):** we will not report "an expert *for* promoters" without a causal DiD; we will not report cross-genome wins without leave-one-phylum-out; we will report P6 (in-domain promoter routing recovery or collapse) regardless of direction.

---

## 11. Deliverables
1. A **genome-annotation pipeline** emitting GFF3 tracks of predicted promoters / RBS / Rho-dependent terminators with calibrated scores, validated against held-out experimental tracks.
2. A **public benchmark** (genome panel + version-pinned labels + token-mapping tables + negative/decoy sets) so the comparison against BPROM/bTSSfinder/Prodigal/RBS Calculator/RhoTermPredict and NT/DNABERT-2/Evo 2/ProkBERT is reproducible by third parties.
3. A **figure/number regeneration suite**: one command rebuilds every reported value.

---

## 12. Reproducibility checklist ("completely verifiable")
- [ ] All datasets pinned by version/accession; download + checksum scripts committed.
- [ ] Deterministic bp↔token mapping tables released per genome.
- [ ] Frozen-backbone feature extraction with fixed seeds and committed configs.
- [ ] LOGO and leave-one-phylum-out splits committed as files (no on-the-fly randomness).
- [ ] Expert-ablation masks and DiD computation scripts committed; raw CE tables released.
- [ ] Baseline tool versions, command lines, and thresholds committed (BPROM, bTSSfinder, Prodigal, RBS Calculator v2.1, RhoTermPredict, etc.).
- [ ] LGM checkpoints and probing code committed (GO dense, NT, DNABERT-2, Evo 2, ProkBERT, …).
- [ ] Environment lockfile + hardware notes (Junho: 4× H200 train; single H200 / A40 eval).
- [ ] One-command regeneration of every figure and table.

---

## 13. Mapping back to Junho's findings (traceability)

| This plan uses | Junho's result it builds on |
|----------------|------------------------------|
| Routing fingerprint as primary feature (96-d) | Routing solved in-domain bio classification (drug-resistance 0.980) with 8× fewer features |
| `routing_concat` as primary detector | "Often better, never worse" than embedding in-domain |
| Expert-detector + causal DiD ablation | L7 E7 tRNA detector, DiD +0.559; control DiD −0.065 |
| Upcycled-vs-scratch specialization test | Scratch silent collapse (14/96 experts < 0.01); causally null tRNA detector (DiD +0.006) |
| Leave-one-phylum-out transfer | Routing AMI 0.15 (function) vs hidden-state AMI 0.80 (taxonomy) |
| In-domain promoter recovery hypothesis (P6) | Routing collapsed on **OOD eukaryotic** GUE promoter (routing_only 0.039 vs embedding 0.779) |
| Intergenic-anchored detector search | Intergenic detector L6 E7 (+1.75) — promoters/RBS/Rho-rut live in/near intergenic space |
| GC/phylum stratification + marginal-masked nulls | His 5-genome (GC 33–67%) and 20-genome PGAP panels; P(e) ≥ 0.01 null; G-test BH-FDR |
