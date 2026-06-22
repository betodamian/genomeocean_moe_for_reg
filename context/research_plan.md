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
| P3 | ≥1 expert per element class shows significant log₂ enrichment over the marginal-masked null; ablating it raises element-token cross-entropy with a Difference-in-Differences (DiD) ≫ 0 vs **matched intergenic non-element control** (not CDS) | No causal expert exists for any class (DiD CI overlaps 0) |
| P4 | Upcycled regulatory experts are causal (DiD ≫ 0); scratch counterparts are weaker or causally null | Scratch matches upcycled on causality |
| P5 | Routing detectors transfer across phyla better than dense embeddings (smaller held-out-phylum performance drop) | Routing transfer ≤ embedding transfer |
| P6 | Bacterial (in-domain) promoter routing recovers, unlike the **eukaryotic GUE promoter** task where Junho saw routing collapse (routing_only MCC 0.039 vs embedding 0.779) | Routing collapses on in-domain bacterial promoters too |
| P7 | **The win is regulatory, not structural.** Routing discriminates each element from its **same-context, same-position decoy** (e.g., promoter vs non-promoter intergenic; promoter vs RBS), and the advantage survives **partialling out** the known structural-class routing axis (Junho's intergenic/CDS detectors) | Discrimination collapses to chance once structural class/position is held fixed or partialled out — i.e., the model only re-detected "intergenic vs coding" |

P6 guards a known risk: Junho's *only* promoter result was **out-of-domain eukaryotic**, where routing collapsed, so the **in-domain bacterial** regime is tested separately. P7 guards the more dangerous risk: the routing channel already encodes intergenic-vs-CDS (the L6 E7 intergenic detector, +1.75), so a routing "win" could reflect that structural distinction rather than anything regulatory. Every primary result is therefore defined against **same-context, same-position decoys** (§4c–4d), and the causal test (§8b) controls against matched intergenic non-element tokens, not CDS. Both predictions are reported regardless of outcome.

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

### 4c. Negative sets — same-context, same-position decoys are *primary*

**The shortcut to avoid:** all three targets live in or beside intergenic/UTR space, and Junho showed the routing channel already encodes the intergenic-vs-CDS distinction (L6 E7, +1.75). So a classifier with *coding* negatives — or negatives that merely differ in length or position — can score high by re-detecting "this is intergenic," learning nothing element-specific. The negative design must remove every such shortcut so the only thing left to discriminate is the regulatory signal itself.

Two tiers of negatives, with the harder tier carrying the headline numbers:
- **Tier-2 (PRIMARY) — matched-context, matched-position hard decoys.** Each positive is contrasted with a negative drawn from the *same structural class* and the *same position relative to the gene*:
  - **Promoter** vs non-promoter upstream/intergenic windows at matched TSS-distance.
  - **RBS (Shine-Dalgarno)** vs **leaderless** gene starts (same position — immediately 5′ of a start codon — but no SD).
  - **Rho-dependent terminator** vs **intrinsic** terminators *and* non-terminating 3′/downstream windows.
  Because both classes share context, length, and position, neither the intergenic-ness nor the start-vs-end location is discriminative — only the element is.
- **Tier-1 (SANITY CHECK ONLY) — generic background.** GC- and length-matched intergenic/genic windows lacking the element. Reported transparently as the *contaminated, easy* setting and never as a primary result, precisely because it is gameable by the structural shortcut.

Confounds explicitly controlled (matched sampling **and** regression covariates): **GC content, element/window length, distance-to-gene-feature, strand, intergenic position, genic vs intergenic context, and BPE-boundary phase.** Null model: **permutation within matched strata** (shuffle labels inside GC/length/position bins) so a "significant" result cannot ride a confound. This extends Junho's discipline (marginal-masked P(e) null; his caution that CDS JSD was near-zero only because CDS dominates token mass).

### 4d. The multi-label "which element" head (structural confound cancels)

Beyond pairwise decoys, we train a single **multi-label / multi-class head** over `{promoter, RBS, Rho-terminator, intrinsic-terminator, plain-intergenic, CDS}`. Forcing the model to say *which* regulatory element (or none) makes the shared intergenic signal **common-mode** across the regulatory classes, so it cancels and cannot drive class separation.

- **Acid test — promoter vs RBS.** Both sit at gene 5′ ends, so they share context *and* position; nothing structural or positional distinguishes them. Routing that cleanly separates promoter from RBS is therefore evidence of genuinely regulatory, element-specific representation — the cleanest possible refutation of the "it only learned intergenic" objection.
- Multi-label (not just multi-class) because a single window can legitimately carry more than one label (e.g., overlapping promoter and 5′-UTR/RBS context); independent per-label heads avoid forcing false mutual exclusivity.

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
- Probes (two, with distinct jobs — Junho's two-probe protocol):
  - **Linear probe — primary.** Logistic Regression (L2, C=1.0). Carries all **MoE-necessity / interpretability** claims (RQ-B): a linear readout makes the clean, non-manufacturable statement "this element is *linearly* decodable from the router," so the model — not a powerful classifier — is doing the work.
  - **XGBoost (balanced) — secondary, two jobs.** (i) The **fair, expressive readout for the head-to-head vs the tools** (RQ-A), since the classical baselines (RhoTermPredict, the Salis RBS Calculator, the ML promoter tools) are themselves nonlinear — restricting GO-MoE to linear would self-handicap the comparison. (ii) A **nonlinear ceiling check**: if XGBoost succeeds where the linear probe fails, the signal is present but nonlinear (rescuing the result from a false negative); if the two are close, the signal is cleanly linear.
- Cross-validation: **leave-one-genome-out (LOGO)**, plus a stricter **leave-one-phylum-out** for transfer (RQ-C). This prevents the taxonomic leakage Junho flagged (bulk CDS splits by lineage/GC).
- Tasks:
  - **(T1) Classification** — element present at locus? Scored **primarily** against the Tier-2 matched-context, matched-position decoys (§4c); the one-vs-generic-background setting is reported only as a labeled sanity check.
  - **(T2) Localization** — predict element position within window; report median bp error and boundary-F1.
  - **(T3) Sub-typing & multi-label** — σ class (promoters), leadered vs leaderless (RBS), Rho-dependent vs intrinsic (terminators), plus the unified multi-label "which element (or none)" head from §4d, with **promoter-vs-RBS** as the structural-confound acid test.
- Metrics: **MCC, macro-F1, AUPRC** (primary, given class imbalance), boundary-F1, bp localization error. Genome-level **bootstrap CIs** and paired tests across folds. **Headline claims use Tier-2 decoys only;** Tier-1 numbers are shown beside them to expose the gap the shortcut would have inflated.

---

## 8. The MoE-necessity experiments (RQ-B core)

Three converging, independently sufficient lines of evidence:

### 8a. Incremental-value / channel ablation (tests P1, P2)
- On the same MoE backbone, compare `routing_concat` vs `embedding_only` vs `routing_only` vs `per_expert_bins`, paired across LOGO folds (Wilcoxon, BH-FDR). **This MoE-necessity comparison uses the linear probe** (§7) — the surplus must be *linearly* readable so the gain is attributable to the representation, not to a classifier finding a boundary. Junho's prior: `routing_concat` was *often better and never worse* in-domain; routing collapsed only on OOD eukaryotic tasks. P6 tests which regime bacterial promoters fall into.
- Then place the MoE `embedding_only` against every dense LGM's embedding to rule out "it's just a better backbone." The MoE claim survives only if **routing_concat > embedding_only AND the surplus features are dense-inaccessible** (true by construction).

### 8b. Expert-detector discovery + causal ablation (tests P3, the gold standard)
1. **Discovery:** for each regulatory class, compute per-(layer, expert) log₂ enrichment on element token-spans vs the **marginal-masked null** (Junho's P(e) ≥ 0.01 baseline), with **G-test + BH-FDR**. Junho's Intergenic L6 E7 (+1.75) is the relevant prior *and the relevant hazard* — promoters/RBS/Rho-rut all sit in intergenic/UTR space, so we must separate an element-specific expert from a generic intergenic expert (see step 2's control).
2. **Causality (control fixed):** mask the candidate expert's logits to −10⁹ on held-out validation tokens; measure cross-entropy increase as **element tokens vs *matched intergenic non-element* control tokens** (GC/length/position-matched) → **Difference-in-Differences**. The control is deliberately **not CDS**: a CDS control would let the intergenic-vs-CDS expert masquerade as an element detector. Require DiD ≫ 0 on the target with a near-zero DiD on a same-layer control expert (paralleling Junho's L7 E7 tRNA test: DiD +0.559; control L7 E3 DiD −0.065).
3. **Structural-axis partial-out:** project routing features onto, and remove, the directions that encode Junho's structural classes (intergenic/CDS/tRNA/rRNA). Require the element's enrichment and detector accuracy to **survive** this removal. If they vanish, the signal was the structural shortcut and we report it as such (P7).
4. **Functional consequence:** re-run the T1/T2/T3 detectors with the expert ablated; require a measurable drop on the **Tier-2 decoy** task specifically. **Dense models have no expert to mask — this entire sub-experiment is impossible for them.**

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
- **Confounds:** GC, length, distance-to-gene-feature, strand, intergenic position, genic/intergenic context, BPE phase as covariates and in matched sampling; permutation-within-strata null.
- **Structural-shortcut control:** the intergenic-vs-CDS distinction is a known routing axis (Junho's L6 E7), so it is treated as a confound to defeat, not a result. Headline numbers use Tier-2 same-context/same-position decoys (§4c–4d); the causal control is matched intergenic non-element tokens, not CDS (§8b); and every claimed advantage must survive partialling out the structural axis (P7).
- **Pre-registered failure criteria:** the falsification rows in §1.
- **Anti-overclaim commitments (in Junho's spirit — he retracted the "1,266 Pfam families" claim):** we will not report "an expert *for* promoters" without a causal DiD against an intergenic (not CDS) control; we will not report a detection win from the Tier-1 background setting; we will not report cross-genome wins without leave-one-phylum-out; we will report P6 (in-domain promoter routing recovery) and P7 (regulatory-not-structural) regardless of direction.

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
| Structural-shortcut control: intergenic (not CDS) ablation control + partial-out + same-context decoys (§4c–4d, §8b, P7) | The same L6 E7 intergenic detector is *also the confound* — routing already separates intergenic from CDS, so a naive win could be purely structural |
| GC/phylum stratification + marginal-masked nulls | His 5-genome (GC 33–67%) and 20-genome PGAP panels; P(e) ≥ 0.01 null; G-test BH-FDR |

---

## 14. Timeline (5-week execution plan)

Five weeks is aggressive for the full program, so the schedule front-loads the **5 high-confidence core genomes** and the two load-bearing questions (RQ-A detection, RQ-B MoE-necessity); the 20-genome diversity sweep (RQ-C) runs only as far as time allows and is the first thing to defer if a week slips.

| Week | Focus | Concrete deliverables | Plan refs |
|------|-------|------------------------|-----------|
| **1** | **Benchmark & setup** | Pre-register hypotheses (P1–P7). Pull the genome panel (5 core first), version-pin ground-truth labels (RegulonDB, Ribo-seq/anti-SD, Peters/Term-seq), build same-context decoys, and generate the deterministic bp↔token mapping. Lock LOGO / leave-one-phylum-out splits to files. | §0, §4, §5 |
| **2** | **Feature extraction** | Run GenomeOcean-MoE (upcycled + scratch) and dense baselines over all loci; extract embeddings, 96-d routing fingerprints, per-expert bins. De-risk end-to-end on E. coli with one pilot detector. | §6 |
| **3** | **Detection & baselines (RQ-A)** | Train per-element detectors (T1/T2/T3, LOGO CV); run every named classical tool (BPROM, bTSSfinder, Prodigal SD, RBS Calculator, RhoTermPredict) and dense LGM (NT, DNABERT-2, Evo 2, ProkBERT, GO-dense) on identical held-out loci; produce the head-to-head scoreboard. | §3, §7 |
| **4** | **MoE-necessity (RQ-B — core)** | Channel ablation (`routing_concat` vs `embedding_only`); expert-detector discovery (enrichment + G-test/BH-FDR); causal ablation (DiD vs intergenic control + structural partial-out); upcycled-vs-scratch; the promoter-vs-RBS acid test. | §4d, §8 |
| **5** | **Transfer, stats & write-up** | Leave-one-phylum-out transfer (as far as the 20 diverse genomes allow); bootstrap CIs + multiple-testing correction; assemble the annotation deliverable, benchmark, and one-command regeneration suite; write results, figures, slides. Buffer for slippage. | §9, §10, §11, §12 |

**Critical path & risk:** Weeks 1–2 (data + features) gate everything — any slip there cascades, so the pilot E. coli detector at the end of Week 2 is the go/no-go checkpoint. If behind by Week 4, protect the MoE-necessity experiments (the unique contribution) and trim the diversity sweep, not the other way around.
