# Research Plan: Expert-Routed Annotation of Bacterial Regulatory Elements with GenomeOcean-MoE

**Author:** Beto Damian
**Status:** Draft v2 — overhauled after the Jun 23 2026 project review
**Builds on:** [`junhos_work.md`](junhos_work.md) (Hong, *GenomeOcean: Sparse Upcycling and Expert Specialization in Genomic MoE*, May 2026)

---

## 0. Summary

Junho showed that GenomeOcean-MoE develops an **expert-routing channel that is genome-invariant and function-aligned** (routing AMI = 0.15 vs hidden-state AMI = 0.80), that specific experts act as **causally verified** detectors of structural classes (the upcycled L7 E7 tRNA detector, DiD = +0.559; null in scratch), and that **routing fingerprints solve in-domain prokaryotic classification with 8× fewer features than embeddings**. This plan asks whether that same routing channel can *identify and annotate* three bacterial regulatory elements — **σ-dependent promoters, Shine-Dalgarno ribosome binding sites (RBS), and Rho-dependent terminators** — across diverse microbial genomes, beating named classical tools and named dense genomic language models, **with the advantage causally attributable to expert specialization such that a no-expert model cannot reproduce it**, and — the central addition in v2 — **demonstrated on sequence-dissimilar data the model has not seen**, so the result is generalization, not memorization.

### What changed in v2 (Jun 23 2026 review)

Restructured around four team-aligned decisions plus Nic K's validation guidance:

- **Phase 0 "ceiling benchmark / smoke test" is a gate (§4).** Before building any annotator, verify the *frozen* embeddings + routing fingerprints discriminate the target elements at all. Go/no-go, with a diagnosis path if they don't. This reframes the current phase as proof-of-concept — "learn whether the signal exists and ask the right questions" — not a finished annotator (Zhong Wang).
- **Frozen-model annotation is the primary methodology (§7).** Features come from frozen embeddings + routing paths; no retraining on the critical path.
- **An independent, *sequence-dissimilar* validation set is now central (§5, Nic K).** Generalization is proven on elements held out by **sequence similarity (≤ ~60% identity)**, not only by taxonomy/GC. Two regimes: classify *seen* data (sanity) vs discover *unseen/novel* elements (the real test, the "Achilles heel").
- **An optional retraining branch (§10) is added for validation only** — to (a) correct the upcycling weight-rescale bug found in the code review and confirm conclusions are robust, and (b) retrain with the validation set *excluded from pretraining*, closing the pretraining-memorization question. Secondary to the frozen path; the code review concluded full retraining is likely unnecessary, and training from scratch is out of intern scope.

---

## 1. Focused research question

> **Primary RQ.** Across diverse bacterial/archaeal genomes — and on elements the model has **not seen** (low sequence similarity to training) — do GenomeOcean-MoE's expert-routing fingerprints identify and annotate σ-dependent promoters, Shine-Dalgarno RBS, and Rho-dependent terminators *more accurately and more generalizably* than (a) classical element-specific tools and (b) dense genomic language models — **and is any advantage causally carried by specialized experts, such that an architecture with no experts cannot achieve it?**

Sub-questions:

- **RQ-A — Detection & annotation accuracy.** Does a routing-aware detector beat the named baselines on presence/absence, localization, and sub-typing — scored against same-context decoys?
- **RQ-B — MoE necessity (core).** Is the gain carried by the routing channel beyond any embedding, and destroyed by ablating the responsible experts? A dense model has neither a routing channel nor experts to ablate.
- **RQ-C — Generalization, not memorization (Nic K's emphasis).** Does the detector hold up on **sequence-dissimilar, unseen** elements (≤60% identity), and on held-out phyla — rather than riding sequence similarity to training? This is the "Achilles heel" the team flagged: proving the model is not simply re-recognizing inputs.

### Falsifiable predictions (pre-registered before any validation element is touched)

| ID | Prediction | Falsified if |
|----|-----------|-------------|
| **P0 (gate)** | **Phase-0 smoke test:** on *seen* data, frozen `routing_concat` discriminates each element from its same-context decoy above a **GC-content-matched** chance baseline | No element clears the GC-matched baseline → Phase 0 fails → enter the diagnosis path (§4) instead of building the annotator |
| P1 | `routing_concat` (864-d) > `embedding_only` (768-d) on the same MoE backbone for all three elements, paired across folds | No significant paired improvement (BH-FDR q < 0.05) |
| P2 | Best routing-aware detector > every named classical tool and dense LGM on held-out genomes (AUPRC, boundary-F1) | A dense LGM or classical tool matches/beats it within bootstrap CI |
| P3 | ≥1 expert per element shows significant log₂ enrichment over the marginal-masked null; ablating it raises element-token cross-entropy with DiD ≫ 0 vs a **matched intergenic non-element** control (not CDS) | No causal expert exists for any class (DiD CI overlaps 0) |
| P4 | Upcycled regulatory experts are causal (DiD ≫ 0); scratch counterparts weaker or null | Scratch matches upcycled on causality |
| P5 | Routing detectors transfer across held-out phyla better than dense embeddings (smaller drop) | Routing transfer ≤ embedding transfer |
| P6 | In-domain bacterial promoter routing recovers, unlike the OOD eukaryotic GUE promoter task (routing_only MCC 0.039 vs embedding 0.779) | Routing collapses on in-domain bacterial promoters too |
| P7 | **Regulatory, not structural.** Routing discriminates each element from its same-context, same-position decoy, and the advantage survives partialling out the structural-class routing axis | Discrimination collapses once structural class/position is held fixed — the model only re-detected "intergenic vs coding" |
| **P8 (Nic K)** | **Generalization, not memorization.** On the **sequence-dissimilar (≤60% identity) unseen split**, the routing-aware detector stays above the GC-matched baseline | Performance collapses toward chance on dissimilar elements → the model was riding sequence similarity to training |
| **P9 (retraining branch, §10)** | The **bug-corrected** retrained model reproduces the expert-specialization and detection results, **and** a model pretrained with the validation set **excluded** preserves the unseen-split performance | Conclusions flip under the corrected model, or unseen-split performance depends on the validation data having been in pretraining (i.e., leakage) |

P0, P8, and P9 are the v2 additions: a go/no-go signal check, an explicit not-memorizing test, and a leakage/robustness control. We report all predictions regardless of direction.

---

## 2. Why this is an MoE-necessary design (the central argument)

"A model with no specialized experts could not achieve the same result," on two independent grounds:

1. **Architectural exclusivity (by construction).** The primary detector's feature vector is the routing fingerprint — per-token router softmax over 12 layers × 8 experts (96-d, pooled) plus per-position expert assignments — optionally concatenated with the residual stream. A dense model **emits no routing distribution**; these features do not exist for it. Dense models can only be evaluated on the embedding-only sub-detector (the P1 control). If `routing_concat` beats `embedding_only`, the surplus is information a dense model cannot supply.

2. **Causal exclusivity (by ablation).** We extend Junho's expert-masking protocol (mask logits to −10⁹, measure cross-entropy shift, compute DiD vs a control expert) to the regulatory classes. A dense network has no expert to mask, and per P4 a *non-specialized* MoE (scratch) is predicted to fail the same test — so it is *specialized* experts, not mere MoE capacity, that carry the result.

*Reviewer caveat (Zhong Wang, Rob Egan):* turning off **any** part of a network degrades performance, so the ablation claim is stated as a **differential** (DiD: element vs matched control), never as absolute degradation. A bare "switching the expert off breaks function" is not claimed; only "it breaks the *element's* function significantly more than a matched control's" (§8b).

---

## 3. Tools and models to beat (named)

### 3a. Classical, element-specific tools

| Element | Tools to beat | Role |
|---------|---------------|------|
| **Promoters (σ-dependent)** | **BPROM** (Softberry σ70), **bTSSfinder** (multi-σ), **G4PromFinder**, **PromoterHunter** (PHISITE), **iPromoter-2L / MULTiPly** (ML) | Primary promoter baselines |
| **RBS (Shine-Dalgarno)** | **Prodigal** RBS/SD scoring, **Salis RBS Calculator v2.1** (ΔG), **RBSfinder**, **Free2Bind** | Primary RBS baselines |
| **Rho-dependent terminators** | **RhoTermPredict**; with **TransTermHP**, **ARNold**, **RNIE**, **WebGeSTer** as intrinsic-terminator comparators | Primary + discrimination |

### 3b. Dense / non-MoE genomic language models

- **GenomeOcean dense — GO-100M, GO-500M, GO-4B** (direct dense controls; same data lineage, no experts).
- **Nucleotide Transformer** (NT-2.5B-multispecies, NT-500M), **DNABERT-2**, **Evo 2** (strongest in-domain), **ProkBERT** (task-matched), **HyenaDNA**/**Caduceus** (secondary).

All LGMs use the **identical frozen-backbone probing protocol** (§7) on **identical held-out loci**; classical tools run at default *and* threshold-tuned operating points; the chance baseline is **GC-content-matched** (§12), not naive shuffling.

---

## 4. Phase 0 — the ceiling benchmark / smoke test (the gate)

**Adopted by the team as a Phase-0 gate.** Before investing in a full annotator, prove the *frozen* model already contains discriminative signal for the targets. This de-risks everything and reframes the current phase as a proof-of-concept (Zhong Wang).

**The smoke test.** On **seen** data (elements from training-distribution genomes), extract `routing_concat` and ask a linear probe to discriminate each element from its **same-context decoy** (§5c). Success criterion = clearing a **GC-content-matched** chance baseline by a pre-set margin (P0). This is the minimum viable signal: if frozen features can't separate a promoter from a same-context non-promoter even on familiar data, nothing downstream will work.

**If Phase 0 passes** → proceed to the full detector (§7–§9) and the unseen-generalization tests (§5, §11).

**If Phase 0 fails** → switch from building to *diagnosing*, using a shared failure-mode tree (adapted from the team's ceiling-benchmark framework):

| Failure mode | Symptom | Response |
|---|---|---|
| 1 — confounded space | features separate by GC / taxonomy, not element | add confounder controls (§5c, §12); report the ceiling honestly |
| 2 — element not recovered | windows don't surface the element signal | change windowing / feature views (routing_only, per-expert bins, per-position) |
| 3 — short-window instability | sub-token / partial windows unstable | targeted representation repair only if needed (links to §10) |

**Deliverable of Phase 0:** a "tested path forward" — a prototype if the signal exists, a diagnosis if it does not. This is the Week-5 fallback floor even in the worst case.

---

## 5. Dataset design — training & validation (Nic K's focus)

The dataset is now the centerpiece, because the team's main concern is **proving generalization, not memorization**. Every label remains traceable to a public, version-pinned source.

### 5a. Formal train / validation split (80 / 20)
- A single, committed **80/20 train–validation split** built once and frozen to files.
- The split is made **by sequence-similarity cluster, not at random** (next point), so the 20% validation cannot be near-duplicates of the 80% training.

### 5b. Sequence-similarity holdout — the key rigor upgrade (Nic K)
- Cluster all element windows (and their flanking sequence) by **sequence identity** with **MMseqs2 / CD-HIT**; hold out **whole clusters** so the validation set contains elements with **≤ ~60% identity** to anything in training.
- Rationale (Nic K): taxonomic/GC diversity is *secondary*; what proves generalization is performance on **independent, low-similarity** sequence. A model that only succeeds on ≥90%-identical elements is memorizing.
- **Two evaluation regimes** map onto Nic K's seen-vs-unseen distinction:
  - **Regime A — seen / in-distribution (sanity + Phase 0).** Classify elements similar to training. Necessary, not sufficient.
  - **Regime B — unseen / novel (the real test, P8).** Detect elements at ≤60% identity to anything seen — i.e., *new sites not already known*. Headline generalization claims rest here. The "Achilles heel" is passing Regime B without leaking through Regime A.
- Report a **performance-vs-similarity curve** (accuracy as a function of max identity to training) so the drop-off is explicit.

### 5c. Window construction and negatives
- **300 bp windows** (≈ 60 tokens at ~4.89 bp/token) centered so the **target element and its functional context dominate the window** (Nic K: extract windows where the majority of the sequence is the target material), rather than a 6-bp motif lost in a long intergenic stretch.
- **Negative / decoy class = same-context lookalikes**, never empty or coding windows (so the model cannot win on "gene vs non-gene"):
  - **Promoter** vs non-promoter upstream/intergenic at matched TSS-distance.
  - **RBS** vs **leaderless** gene starts (same position, no SD).
  - **Rho terminator** vs **intrinsic** terminators and non-terminating 3′ windows.
- A two-tier scheme: **Tier-2 same-context, same-position decoys are PRIMARY** (carry headline numbers); a Tier-1 generic-background setting is a *labeled sanity check only*, since it is gameable by the intergenic shortcut.
- **Multi-label "which element (or none)" head** over `{promoter, RBS, Rho-term, intrinsic-term, intergenic, CDS}`: independent per-element heads (so co-occurring elements are handled, a window can be promoter *and* RBS), and the structural classes as their own labels so the shared intergenic signal cancels. The **promoter-vs-RBS** discrimination (same context *and* position) is the acid test that a win is regulatory, not structural.

### 5d. Genome panel and ground truth
- **High-confidence core (5):** *E. coli* K-12 MG1655, *B. subtilis* 168, *M. tuberculosis* H37Rv, *P. aeruginosa* PAO1, *S. aureus* (experimentally curated labels). **Diverse set (20):** PGAP-annotated genomes spanning 7 phyla + 2 archaea, from NCBI / top reference genomes. Diversity is now a *secondary* axis; the **sequence-similarity holdout is primary**.
- Ground-truth sources (experimental, never software-predicted):
  - **Promoters/TSS:** RegulonDB (strong-evidence only), PRODORIC, DBTBS/SubtiWiki, dRNA-seq TSS catalogs (H. pylori, M. tuberculosis, Synechocystis, Campylobacter, Salmonella, B. subtilis).
  - **RBS:** ribosome-profiling TIS maps; SD defined −20..−1 of confirmed starts; anti-SD pairing; leaderless transcripts as contrast.
  - **Rho terminators:** NET-seq + bicyclomycin readthrough (Peters et al. 2012), Term-seq (Dar et al. 2016), RhoTermPredict positives; intrinsic terminators as contrast. *Limited mostly to E. coli → run as a deep case study; carry cross-genome claims on promoters and RBS.*

### 5e. Pretraining-leakage caveat (motivates §10)
GO-MoE was pretrained on ~645 Gbp of metagenomic data, so a validation genome can overlap the **backbone's** pretraining even when it is sequence-dissimilar to the **probe's** training set. The sequence-similarity holdout controls leakage at the probe level; fully closing the question at the backbone level requires the optional retraining branch (§10).

---

## 6. Coordinate → token mapping (deterministic, released)

Elements are small vs ~4.89 bp/token: −10 box ≈ 6 bp ≈ 1–2 tokens, SD ≈ 6–8 bp ≈ 1–2 tokens, Rho *rut* ≈ 60–100 nt ≈ 12–20 tokens. Steps: tokenize both strands with the GO-4B BPE tokenizer and persist a bp↔token table per genome; define each element's core token span + flank to the 300 bp window; report the fraction token-resolvable and record **BPE phase** (clean vs split) as a covariate; report localization in **bp** (de-tokenized) for fair comparison with bp-resolution tools.

---

## 7. Feature extraction (frozen backbone — primary methodology)

Per the team decision, the annotation methodology uses **frozen** GO-MoE; no retraining on the critical path. For **GenomeOcean-MoE** (upcycled pilot4 primary; pilot3 and scratch as comparators), per window extract: **embedding** (768-d), **routing-only** (96-d pooled router softmax + per-position routing), **routing_concat** (864-d), **per-expert bins**. For **dense LGMs**: embedding-only (no routing channel exists). This asymmetry is the experiment.

---

## 8. Detection / annotation models (frozen-backbone, held-out)

- **Probes (two roles):**
  - **Linear probe (logistic regression) — primary.** Carries all MoE-necessity claims (RQ-B): a linear readout makes the non-manufacturable claim "the element is *linearly* decodable from the router," so the model — not a powerful classifier — does the work.
  - **XGBoost — secondary.** (i) The fair, expressive readout for the head-to-head vs the nonlinear tools (RQ-A); (ii) a nonlinear ceiling check (linear-fail + XGBoost-success = signal present but nonlinear).
- **Splits:** the §5 sequence-similarity train/val split (primary), plus **leave-one-genome-out** and **leave-one-phylum-out** for transfer (RQ-C).
- **Tasks:** (T1) classification vs Tier-2 decoys; (T2) localization (bp error, boundary-F1); (T3) sub-typing + the multi-label head (promoter-vs-RBS acid test).
- **Metrics:** MCC, macro-F1, **AUPRC** (primary, given class imbalance), boundary-F1, bp error, with **genome-level bootstrap CIs**. Headline numbers use Tier-2 decoys and the **unseen (Regime B)** split; seen/Tier-1 numbers are shown beside them to expose the gap.

---

## 9. The MoE-necessity experiments (RQ-B core)

**9a. Incremental value / channel ablation (P1, P2).** On the same backbone compare `routing_concat` vs `embedding_only` vs `routing_only` vs `per_expert_bins`, paired, using the **linear probe** so the surplus is linearly attributable to the representation. Then place MoE `embedding_only` against each dense LGM's embedding to rule out "just a better backbone." The MoE claim survives only if `routing_concat > embedding_only` and the surplus is dense-inaccessible (true by construction).

**9b. Expert-detector discovery + causal ablation (P3).** Per class, compute per-(layer, expert) log₂ enrichment vs the marginal-masked null (G-test + BH-FDR). Then mask the candidate expert and measure cross-entropy increase as **element vs matched intergenic non-element** control → **DiD** (control is deliberately *not* CDS, or the intergenic expert masquerades as an element detector). Require DiD ≫ 0 on target with ≈0 on a same-layer control expert. **Structural-axis partial-out:** remove the directions encoding Junho's structural classes and require the signal to survive (else report it as the structural shortcut, P7). Dense models have no expert to mask — this is impossible for them.

**9c. Upcycled vs scratch (P4).** Repeat 9b for MoE-Scratch; predict its regulatory experts are weaker/causally null (Junho: scratch silent collapse, null tRNA DiD +0.006) — isolating *successful specialization*, not MoE capacity, as the active ingredient.

---

## 10. Optional retraining branch (validation purposes only — secondary)

The frozen path (§7) is primary. This branch is a **stretch / collaborative** control (training is heavy — Junho used 4× H200, ~50k steps, ~52B tokens — and is beyond intern scope per Zhong Wang; run only with mentor/Junho support). The code review concluded **full retraining is likely unnecessary**, so this branch exists to *stress-test the conclusions*, not to produce the deliverable.

**10a. Bug-correction control (robustness).** Zhong Wang's code review found a **math error in the upcycling**: weights were not rescaled by 2× to preserve the original distribution after the 50% expert-dropout step. Validation loss still reached expected levels by ~50k steps, so the curriculum is effective — but to confirm conclusions are robust to the bug, re-run upcycling **with the 2× rescale corrected** and check that expert specialization (the enrichment + DiD experts of §9) and the detection results reproduce. *Pass:* conclusions hold under the corrected model (P9, part 1).

**10b. Leakage-free control (closing the memorization question).** To answer Nic K's "Achilles heel" at the backbone level: produce a model variant whose **pretraining excludes the §5 validation genomes / sequence clusters** (continue-pretrain or upcycle from a leakage-free checkpoint). Re-run Regime B (unseen, ≤60% identity). *Pass:* unseen-split performance is preserved even when the validation data was never in pretraining → the result is generalization, not pretraining memorization (P9, part 2). *Fail:* performance depended on having seen the data in pretraining → report the leakage honestly and scope the claim down.

**Scope/honesty:** this branch is explicitly optional and secondary. If compute or collaboration is unavailable, the sequence-similarity holdout (§5b) is the practical leakage control and the frozen path stands on its own; the retraining branch only *strengthens* the leakage argument, it is not required for a valid result.

---

## 11. Generalization & transfer (RQ-C)

- **Primary (Nic K):** the **sequence-dissimilar unseen split** (Regime B, §5b, P8) — performance vs max-identity-to-training curve, against the GC-matched baseline.
- **Secondary:** **leave-one-phylum-out** transfer (P5) — routing should drop less than dense embeddings on held-out phyla (Junho's AMI 0.15 vs 0.80), a second, independent MoE-advantage axis.

---

## 12. Statistical rigor, controls, honesty commitments

- **Chance baseline = GC-content-matched profiles** (Zhong Wang / Rob Egan), not naive shuffling — "better than random" must mean better than GC-matched random.
- **Nulls:** label permutation; marginal-masked routing null; permutation-within-matched-strata (GC/length/position/BPE-phase).
- **Confounds controlled:** GC, length, distance-to-feature, strand, intergenic position, genic/intergenic context, BPE phase — as covariates and in matched sampling.
- **Multiple testing:** BH-FDR across all (layer, expert, class, genome) cells.
- **Class imbalance:** AUPRC/MCC primary; balanced probes; per-class reporting (rare Rho-terminators and σ-subtypes never hidden in an aggregate).
- **Ablation stated as a differential** (DiD vs matched control), never absolute degradation (§2 caveat).
- **Structural shortcut treated as a confound to defeat**, not a result (Tier-2 decoys, intergenic-not-CDS control, partial-out).
- **Pre-registered failure criteria** = the P0–P9 falsification rows.
- **Anti-overclaim:** no "expert *for* promoters" without a causal DiD against an intergenic control; no detection win reported from the Tier-1 setting; no generalization claim without the unseen (Regime B) split; P6/P7/P8 reported regardless of direction.

---

## 13. Deliverables

1. **Phase-0 ceiling report** (go/no-go + diagnosis) — the guaranteed Week-5 floor.
2. If Phase 0 passes: a **retrieval-style annotation prototype** emitting GFF3 promoter/RBS/Rho tracks with calibrated scores **and an abstain option when confidence is low**, validated on the unseen split and against held-out experimental tracks.
3. A **public benchmark** (genome panel + version-pinned labels + 80/20 sequence-similarity split + token-mapping tables + decoys) reproducing the comparison vs BPROM/Prodigal/RhoTermPredict and NT/Evo 2/ProkBERT.
4. A **one-command figure/number regeneration suite**.

---

## 14. Reproducibility checklist

- [ ] Datasets pinned by version/accession; download + checksum scripts committed.
- [ ] **80/20 sequence-similarity split committed as files** (MMseqs2/CD-HIT params + cluster assignments); LOGO / leave-one-phylum-out splits committed.
- [ ] Deterministic bp↔token mapping tables released per genome.
- [ ] Frozen-backbone feature extraction with fixed seeds + committed configs.
- [ ] Phase-0 smoke-test script + GC-matched baseline generator committed.
- [ ] Expert-ablation masks + DiD scripts; raw CE tables released.
- [ ] Baseline tool versions/commands/thresholds committed.
- [ ] LGM checkpoints + probing code committed.
- [ ] (If run) retraining-branch configs: the 2× rescale fix and the leakage-free pretraining exclusion list committed.
- [ ] Environment lockfile + hardware notes.

---

## 15. Traceability

### To Junho's findings
| This plan uses | Junho's result |
|---|---|
| Routing fingerprint as primary feature | Routing solved in-domain classification with 8× fewer features |
| `routing_concat` primary detector | "Often better, never worse" than embedding in-domain |
| Expert-detector + causal DiD ablation | L7 E7 tRNA detector DiD +0.559; control −0.065 |
| Upcycled-vs-scratch test | Scratch silent collapse; null tRNA DiD +0.006 |
| Leave-one-phylum-out transfer | Routing AMI 0.15 (function) vs hidden 0.80 (taxonomy) |
| In-domain promoter recovery (P6) | Routing collapsed on OOD eukaryotic GUE promoter |
| Intergenic-anchored search + structural-shortcut control | Intergenic detector L6 E7 (+1.75) is both prior and confound |

### To the Jun 23 2026 review
| This plan reflects | Source |
|---|---|
| Phase-0 ceiling benchmark / smoke test gate (§4) | Team decision; ceiling-benchmark framework |
| Frozen-model annotation, no retraining on critical path (§7) | Team decision |
| Sequence-similarity holdout, ≤60% identity, seen-vs-unseen, 80/20 (§5) | Nic K |
| 300 bp windows, majority-target windows, NCBI sourcing (§5c–5d) | Nic K |
| GC-content-matched chance baseline; differential ablation claim (§12, §2) | Zhong Wang / Rob Egan |
| Retraining branch: 2× rescale bug fix + leakage-free pretraining (§10) | Code review (upcycling error) + Nic K (memorization) |
| Supervised verification as a valid first step | Rob Egan |
| Same-context (RBS-style) negatives, not intergenic-vs-coding (§5c) | Team / prior feedback |

---

## 16. Timeline (5-week execution plan)

Restructured to put the Phase-0 gate and the sequence-similarity split first. The retraining branch (§10) is **off the critical path** — pursued only if Phase 0 passes early and mentor/compute support exists.

| Week | Focus | Deliverables | Refs |
|------|-------|--------------|------|
| **1** | **Dataset + split** | Pull 5 core genomes + experimental labels; build 300 bp windows + same-context decoys; **cluster by sequence identity and commit the 80/20 ≤60% split**; lock all splits to files. | §5, §6 |
| **2** | **Phase 0 smoke test (gate)** | Frozen feature extraction; run the smoke test vs the GC-matched baseline on seen data (P0). Go/no-go; if no-go, run the §4 diagnosis tree. | §4, §7 |
| **3** | **Detection + baselines (RQ-A)** | Train linear-probe detectors; score head-to-head vs classical tools and dense LGMs on held-out genomes. | §3, §8, §9a |
| **4** | **MoE-necessity + generalization (RQ-B/C)** | Expert ablation (DiD, intergenic control, partial-out), upcycled-vs-scratch, promoter-vs-RBS acid test; **unseen (Regime B) generalization curve** (P8) + leave-one-phylum-out. | §9, §11 |
| **5** | **Stats, write-up, (optional) retraining** | GC-matched nulls + bootstrap CIs; annotation prototype or Phase-0 diagnosis; figures/slides. *If time + support:* launch §10 bug-fix / leakage-free controls. | §10, §12, §13 |

**Critical path & risk:** Weeks 1–2 gate everything; the Phase-0 smoke test is the explicit go/no-go. If behind by Week 4, protect the MoE-necessity experiments and the unseen-split generalization test (the two unique contributions), and defer the diversity sweep and the retraining branch. Even a Phase-0 *failure* yields a publishable diagnosis (§13.1), so there is a guaranteed deliverable.
