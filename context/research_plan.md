# Research Plan: Expert-Routed Annotation of Bacterial Regulatory Elements with GenomeOcean-MoE

**Author:** Beto Damian
**Status:** Draft v2.7 — overhauled after the Jun 23 2026 project review; Jun 24 2026 additions: data-sufficiency audit (§5f), verified database sourcing (§5d), curated RBS database (`data/rbs_database/`), RBS robustness audit, comparison-integrity / anti-circularity safeguards (§12), **element-strengthening pass** (RBS → 4 phyla + archaea; Rho → 2-phylum in vivo panel, `data/rho_database/`), and **trivial-specialization-trap control** (§9d / P10, Zhong Wang) — context-conditional MI given token-identity to separate function-aware routing (H₂) from motif-token hashing (H₁)
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
| P5 | Routing detectors transfer across held-out phyla better than dense embeddings (smaller drop). **Breadth differs by element (§5f):** promoters ≥7 phyla, RBS 4 phyla + archaea, **Rho a 2-phylum transfer (E. coli↔MTB)** | Routing transfer ≤ embedding transfer |
| P6 | In-domain bacterial promoter routing recovers, unlike the OOD eukaryotic GUE promoter task (routing_only MCC 0.039 vs embedding 0.779) | Routing collapses on in-domain bacterial promoters too |
| P7 | **Regulatory, not structural.** Routing discriminates each element from its same-context, same-position decoy, and the advantage survives partialling out the structural-class routing axis | Discrimination collapses once structural class/position is held fixed — the model only re-detected "intergenic vs coding" |
| **P8 (Nic K)** | **Generalization, not memorization.** On the **sequence-dissimilar (≤60% identity) unseen split** — *cross-genome* for all three elements (promoters ≥7 phyla, RBS 4 phyla + archaea, **Rho 2-phylum E. coli↔MTB**) plus intra-genome dissimilar splits (§5f) — the routing-aware detector stays above the GC-matched baseline | Performance collapses toward chance on dissimilar elements → the model was riding sequence similarity to training |
| **P9 (retraining branch, §10)** | The **bug-corrected** retrained model reproduces the expert-specialization and detection results, **and** a model pretrained with the validation set **excluded** preserves the unseen-split performance | Conclusions flip under the corrected model, or unseen-split performance depends on the validation data having been in pretraining (i.e., leakage) |
| **P10 (Zhong Wang — trivial-specialization trap)** | **Routing is context-conditional for regulatory classes, not token-ID hashing.** For motif/element tokens appearing in ≥2 functional contexts, routing carries class information *conditional on token identity*: **I(expert ; element_class \| token_id) > 0** above a **within-token label-shuffle null**, at the informative layers — i.e. the same BPE token (TATAAT, TTGACA, GGAGG, C-rich *rut*) routes differently when functional vs decoy | Conditional MI ≈ null (CI overlaps the within-token shuffle) → routing only re-encodes motif/token presence (H₁ "trivial specialization"); the "function-aware expert" claim collapses to motif-token detection and **RQ-B is not supported** (RQ-A detection still valid) |

P0, P8, and P9 are the v2 additions: a go/no-go signal check, an explicit not-memorizing test, and a leakage/robustness control. **P10 (v2.7, Zhong Wang)** imports the H₀/H₁/H₂ disambiguation that validated Junho's structural-expert specialization (junhos_work.md §4 follow-up) and applies it to the regulatory classes. We report all predictions regardless of direction.

---

## 2. Why this is an MoE-necessary design (the central argument)

"A model with no specialized experts could not achieve the same result," on two independent grounds:

1. **Architectural exclusivity (by construction).** The primary detector's feature vector is the routing fingerprint — per-token router softmax over 12 layers × 8 experts (96-d, pooled) plus per-position expert assignments — optionally concatenated with the residual stream. A dense model **emits no routing distribution**; these features do not exist for it. Dense models can only be evaluated on the embedding-only sub-detector (the P1 control). If `routing_concat` beats `embedding_only`, the surplus is information a dense model cannot supply.

2. **Causal exclusivity (by ablation).** We extend Junho's expert-masking protocol (mask logits to −10⁹, measure cross-entropy shift, compute DiD vs a control expert) to the regulatory classes. A dense network has no expert to mask, and per P4 a *non-specialized* MoE (scratch) is predicted to fail the same test — so it is *specialized* experts, not mere MoE capacity, that carry the result. **Necessary precondition (P10, §9d):** the specialization must be **context-conditional**, not token-ID hashing — for short-motif elements (TATAAT, GGAGG, *rut*) this is the load-bearing control, since otherwise "function-aware expert" reduces to "the motif-token expert." Junho's work passed exactly this test (I(expert ; class | token_id) = 0.243 bits at L7, 42× the within-token null; junhos_work.md §4 follow-up); the regulatory work must reproduce it.

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
  - **RBS** (SD-containing) vs **UNSD-leadered** starts **within the same organism** (same position, same 5′-UTR context, same GC; differ only in the SD motif) — this is the PRIMARY RBS decoy. *Leaderless* starts are a **secondary** decoy only, because experimentally-confirmed leaderless starts concentrate in GC-rich *M. tuberculosis* (65.5% GC) while SD positives concentrate in lower-GC genomes (E. coli 51%, B. subtilis 43.5%, S. aureus 33%) — so a cross-organism SD-vs-leaderless contrast is confounded with GC/taxonomy and is gutted by the GC-matched baseline (§5d, §12).
  - **Rho terminator** vs **intrinsic** terminators and non-terminating 3′ windows.
- A two-tier scheme: **Tier-2 same-context, same-position decoys are PRIMARY** (carry headline numbers); a Tier-1 generic-background setting is a *labeled sanity check only*, since it is gameable by the intergenic shortcut.
- **Multi-label "which element (or none)" head** over `{promoter, RBS, Rho-term, intrinsic-term, intergenic, CDS}`: independent per-element heads (so co-occurring elements are handled, a window can be promoter *and* RBS), and the structural classes as their own labels so the shared intergenic signal cancels. The **promoter-vs-RBS** discrimination (same context *and* position) is the acid test that a win is regulatory, not structural.

### 5d. Genome panel and ground truth
- **High-confidence core (5):** *E. coli* K-12 MG1655, *B. subtilis* 168, *M. tuberculosis* H37Rv, *P. aeruginosa* PAO1, *S. aureus* (experimentally curated labels). **Diverse set (20):** PGAP-annotated genomes spanning 7 phyla + 2 archaea, sequences from NCBI RefSeq. Diversity is now a *secondary* axis; the **sequence-similarity holdout is primary**.
- Ground-truth sources (experimental, never software-predicted):

  **σ-Promoters / TSS — primary source: PPD**
  - **PPD** (Prokaryotic Promoter Database, Lin group, http://lin-group.cn/database/ppd/): **129,148 manually curated, experimentally verified promoters across 63 prokaryotic species / 74 strains**. This is the primary download for promoter positives; it consolidates dRNA-seq catalogs, RegulonDB, and DBTBS into one flat-file download with per-entry evidence codes and organism IDs. Use PPD as the label source; filter to strong/experimental evidence only.
  - **RegulonDB v12** (http://regulondb.ccg.unam.mx/, downloadable flat files): 4,050 E. coli K-12 promoters with evidence-level flags — use as a high-fidelity cross-check and source for E. coli σ-factor sub-typing labels not in PPD.
  - **dRNA-seq TSS catalogs** (SRA/GEO, accessed per-paper): *M. tuberculosis* H37Rv (~4,000+ TSS, Cortes et al. 2013 PMC3898074), *Campylobacter jejuni* (Dugar et al. 2013 PMC3656092, 3 isolates), *Salmonella* Typhimurium (3,583 TSS, Kroger et al. 2013 PMC3963941), *Synechocystis* sp. (TSSNote 2022 PMC9745317), *Bacillus subtilis* BSGatlas, *H. pylori* (Sharma et al. 2010 — the method's origin). Use these to supplement PPD where PPD coverage is sparse and to verify cross-organism leaderless-TSS calls.
  - **PRODORIC** (https://www.prodoric.de/, 27 organisms): TFBS-centric; use for σ-factor binding site motifs and as a secondary cross-check for organisms in PPD.
  - *Excluded:* CDBProm (24M predicted promoters, software-generated — not experimental ground truth per §12).

  **RBS / Shine-Dalgarno — curated in-house database: `data/rbs_database/`**

  No public multi-organism experimental RBS database exists. A curated database was built from verified experimental sources (see [`data/rbs_database/README.md`](../data/rbs_database/README.md) for full schema, evidence hierarchy, curation rules, and per-organism fetch instructions). Summary:

  | Organism | Primary source | T1 TIS entries (est.) | Best leaderless-negative source |
  |---|---|---|---|
  | *E. coli* K-12 | Meydan et al. 2019 Ribo-RET (GSE122129); Nakahigashi 2016 TetRP (PRJDB2960); Saito 2020 ΔaSD Ribo-seq (GSE135906) | ~2,000–4,290 | UNSD-leadered from Meydan supp. table (~54% SD → ~46% non-SD) |
  | *M. tuberculosis* H37Rv | Zhu et al. 2021 Ribo-seq (E-MTAB-8835); Cortes 2013 dRNA-seq (SRP028740) | ~3,500 | **497 ribosome-profiling-confirmed leaderless + ~1,040 dRNA-seq-confirmed leaderless** — primary leaderless negative class |
  | *B. subtilis* 168 | Lalanne et al. 2017 Ribo-seq (GSE95211); Bhatt 2024 (GSE249450) | ~4,200 | rare leaderless (~6%); use UNSD-leadered as negative |
  | *S. aureus* NCTC 8325 | Kohl et al. 2026 Nat Commun (extended-SD Ribo-seq, PMID 41680142) | ~2,700 | very few leaderless; extended SD motifs — note species-specific motif in SD window |
  | *P. aeruginosa* PAO1 | PGAP-derived (Tier-2, no T1 Ribo-seq found as of Jun 2026) | ~5,570 (T2) | no leaderless call possible without TSS data |

  **Evidence hierarchy** (T1 = ribosome-profiling TIS arrest; T2 = dRNA-seq leaderless or PGAP-derived; predicted labels excluded per §12). Expected total: ~18,000–20,000 entries; ~11,000 SD positives; ~2,500–3,500 **UNSD-leadered** entries (the primary negative); ~1,800 leaderless (secondary negative). P. aeruginosa T2/PGAP-derived entries are excluded from headline Tier-2 results until upgraded.

  **Primary decoy = SD vs UNSD-leadered, within organism** (§5c): both classes are leadered, ribosome-bound, same position, same GC — they differ only in the SD motif, so a win is attributable to SD detection, not GC/taxonomy. UNSD exists in every Ribo-seq organism (E. coli ~46% of starts, MTB 678, B. subtilis ~16%, S. aureus low), so this decoy is **GC-clean and cross-genome-capable**. Leaderless is the secondary decoy only (GC-confounded, §5c).

  **For the head-to-head vs classical tools/gLMs**, the SD-vs-UNSD split (ΔG-defined) is *not* the headline task — it would be circular against the ΔG-based RBS baselines (Salis, Free2Bind). The headline RBS comparison is **experimentally-grounded start detection/localization** (Ribo-seq-confirmed TIS vs same-context non-start); SD-vs-UNSD is reported separately and drives the MoE-internal routing analysis (§9). See §12 anti-circularity.

  **Organism-specific SD labeling (required).** A single E. coli-tuned rule mislabels divergent SD systems — *S. aureus* uses **extended** start-codon-proximal SD motifs and its native starts are not decoded by *E. coli* ribosomes (Kohl et al. 2026). SD class is therefore assigned with **organism-specific anti-SD sequences** (the 16S rRNA 3′ tail of each genome) and **organism-specific ΔG thresholds**, not one global cutoff. Default ΔG ≤ −3.4 kcal/mol (E. coli/Salis) is used only where no organism-specific calibration exists, and the anti-SD sequence + threshold used per entry is recorded in the database `notes`.

  **RBS is the highest Phase-0 risk element (window/motif-size mismatch).** The SD core is 6–8 bp ≈ 1–2 tokens inside a 300 bp ≈ 60-token window (§6), so the **pooled** `routing_concat` can dilute the SD signal below detectability. RBS therefore pre-registers a **window-size + pooling sensitivity sweep** (e.g. 300/120/60 bp; mean-pool vs per-position vs per-expert-bin) and leads with the **per-position / per-expert-bin** routing views rather than the pooled fingerprint. If RBS fails Phase-0 P0, the §4 failure-mode tree (mode 2: element not recovered → change windowing/feature views) is entered for RBS specifically.

  **Phylum coverage for RBS = 4 bacterial phyla + an archaeal domain axis** (strengthened Jun 24 2026):
  - **Gammaproteobacteria** — *E. coli* (Ribo-RET/TetRP, §above).
  - **Actinobacteria** — *M. tuberculosis* (Ribo-seq; leaderless-rich).
  - **Firmicutes** — *B. subtilis* + *S. aureus* (extended-SD).
  - **Alphaproteobacteria** — ***Caulobacter crescentus*** (Schrader et al. 2014, GEO **GSE54883**: ribosome profiling + RNA-seq + global 5′-RACE + LC-MS proteomics; 3,235 CDS, start codons corrected for 12.8%, 94 new CDS). A 4th, phylogenetically distinct bacterial phylum with multi-method experimental TIS — adds a leave-one-phylum-out fold.
  - **Archaea (domain-transfer axis)** — ***Haloferax volcanii*** (NAR 2020, 48:5201: HHT ribosome profiling; 1,413 annotated TIS, >70% leaderless). Tests whether RBS routing generalizes across the **bacteria→archaea domain boundary** — a stronger generalization probe than any within-bacteria fold, and a second leaderless-rich organism complementing MTB. Reported as an optional domain-transfer result, consistent with Junho's diverse panel (7 phyla + 2 archaea).

  Leave-one-phylum-out (P5) for RBS now spans **4 bacterial phyla**; **S. aureus extended-SD remains the strongest within-bacteria divergent-motif test**, and *H. volcanii* adds an archaeal domain-transfer test. *P. aeruginosa* (Gammaproteobacteria) still contributes no experimental RBS truth (PGAP-derived T2 only) and is not counted in the experimental phylum tally.

  **Rho-dependent terminators — 2-phylum in-vivo cross-genome panel (strengthened Jun 24 2026; was E. coli-only)**

  Rho is no longer a single-genome case study: genome-wide *in vivo* experimental Rho-termination maps now exist in two phyla via **orthogonal methods**, plus an *in vitro* third organism. See [`data/rho_database/README.md`](../data/rho_database/README.md).

  | Organism / phylum | Source | Method | Sites | Evidence tier |
  |---|---|---|---|---|
  | *E. coli* / Gammaproteobacteria | **Peters et al. 2012** (PNAS 109:15584) | NET-seq + bicyclomycin (BCM) readthrough | ~1,000 | T1 in vivo |
  | *E. coli* / Gammaproteobacteria | **Term-seq 3′-end map 2018** (NAR 46:6797, GEO **GSE109766**) | Term-seq exact RNA 3′-ends, BCM-validated | 144 high-confidence, bp-resolution | T1 in vivo (concordance) |
  | *M. tuberculosis* / Actinobacteria | **Botella et al. 2022/2023** (PMID 37096044, ArrayExpress **E-MTAB-11753**) | Term-seq + **RhoDUC genetic Rho depletion** (Mtb Rho is BCM-resistant → orthogonal to E. coli) | **1,385** RD-TTS | T1 in vivo |
  | *B. subtilis* / Firmicutes | **H-SELEX rut map** (PMC12350095) | In vitro Rho-helicase SELEX (BsRho) | 600 BsRho rut sites | **T2 in vitro** (flagged) |

  - **High-confidence E. coli positives = concordance** of Peters 2012 (BCM readthrough) ∩ NAR 2018 (Term-seq exact 3′-ends) — two orthogonal in vivo signals, giving both presence/absence and **bp-resolution boundaries** for the localization task (T2 of §8).
  - **Cross-genome Rho is now real but modest:** E. coli ↔ *M. tuberculosis* = 2 phyla, 2 orthogonal in vivo methods (~1,000 and ~1,385 sites). This supports leave-one-genome-out and a cross-phylum unseen test for Rho (previously impossible). *B. subtilis* H-SELEX adds a 3rd phylum but is **in vitro** → reported with the evidence-type caveat, not co-mingled with in vivo positives.
  - **Intrinsic terminators as contrast class:** **TERMITe 2024** (Nucleic Acids Research, PMC12207403, GitHub + Zenodo): experimentally confirmed intrinsic terminator atlas across **13 bacterial species** (E. coli 691–957 sites, *B. subtilis* 635–1,214, *L. monocytogenes* 862, *E. faecalis* 796, 6 *Streptomyces* spp., *Z. mobilis*). Use the per-organism intrinsic set as the Tier-2 "same-context, same-position" decoy (Rho vs intrinsic at matched 3′ positions) — available in E. coli **and** B. subtilis, and Streptomyces gives an Actinobacterial intrinsic contrast near MTB.
  - **RhoTermPredict positives** (Gaviria-Cantin et al. 2019 PMC6407284, E. coli + *B. subtilis* + *Salmonella*) retained as supplementary confirmed sites — but note these are partly algorithm-derived, so used only as cross-checks, never as headline ground truth (§12 anti-circularity).
  - *Honest scope:* the in vivo Rho cross-genome panel is **2 phyla** (E. coli, M. tuberculosis) — genuine but modest versus the promoter panel's ≥7 phyla. The cross-phylum Rho generalization claim is therefore stated as a **2-organism transfer** (train E. coli → test MTB and vice versa), not a broad multi-phylum sweep; *B. subtilis* (in vitro) and the Streptomyces intrinsic contrast (TERMITe) are supporting, not headline.

### 5e. Pretraining-leakage caveat (motivates §10)
GO-MoE was pretrained on ~645 Gbp of metagenomic data, so a validation genome can overlap the **backbone's** pretraining even when it is sequence-dissimilar to the **probe's** training set. The sequence-similarity holdout controls leakage at the probe level; fully closing the question at the backbone level requires the optional retraining branch (§10).

### 5f. Data provenance and per-element sufficiency (verified Jun 24 2026 — the binding constraint)

**NCBI alone is not a label source, and was never one.** A direct check of the best-annotated bacterial genome on NCBI — *E. coli* K-12 RefSeq `NC_000913.3` — returns **zero `regulatory_class` features** (no `promoter`, no `ribosome_binding_site`, no `terminator`; the only "terminator" strings are protein product names). NCBI RefSeq/PGAP supplies only **(i) the genome sequence** and **(ii) structural gene coordinates** (CDS / tRNA / rRNA / ncRNA — the same scaffold Junho used; verified gene counts: E. coli 4,290 CDS, B. subtilis 4,237, M. tuberculosis 3,906, P. aeruginosa 5,572, S. aureus 2,767). In this project NCBI therefore plays exactly two roles: the **sequence substrate**, and the **structural scaffold** that defines CDS starts (for RBS windows), intergenic spans, and the same-context decoys (§5c). **Every regulatory label comes from the experimental databases below, not from NCBI.** A pipeline that expects to read promoter/RBS/terminator truth out of NCBI is reading empty features and must be corrected.

**Per-element data census (verified counts, Jun 24 2026):**

| Element | Experimental ground-truth source(s) | Volume | Organism / phylum breadth | Cross-genome unseen split (P5, P8) feasible? |
|---|---|---|---|---|
| **σ-promoters / TSS** | **PPD** (primary): **129,148** experimentally verified promoters across **63 species / 74 strains**; RegulonDB v12 (E. coli 4,050, cross-check + σ sub-typing); dRNA-seq catalogs (MTB ~4k+, Campylobacter, Salmonella 3.5k, Synechocystis, B. subtilis, H. pylori — supplement PPD) | **129k+ across ≥63 species / ≥7 phyla** | **Yes** — robust multi-phylum coverage; supports leave-one-phylum-out and ≤60% cross-genome unseen split |
| **Shine-Dalgarno RBS** | Curated DB `data/rbs_database/` (10 Ribo-seq/dRNA-seq sources); SD vs **UNSD-leadered** primary decoy (§5c) | ~11k+ SD positives; ~3k UNSD negatives; ~1.8k leaderless (secondary). **4 bacterial phyla + archaea** (Gamma/Actino/Firmicutes/Alpha + *H. volcanii*) | **Yes via UNSD** (GC-clean, every Ribo-seq organism); 4-phylum leave-one-out + archaeal domain-transfer. S. aureus extended-SD = strongest within-bacteria test |
| **Rho-dependent terminators** | Curated DB `data/rho_database/`: Peters 2012 BCM + NAR 2018 Term-seq (E. coli in vivo), Botella 2022 RhoDUC depletion (MTB in vivo, E-MTAB-11753), B. subtilis H-SELEX (in vitro); TERMITe intrinsic decoy | E. coli ~1,000 + MTB ~1,385 in vivo; B. subtilis 600 in vitro | **Yes — modest:** 2-phylum in vivo (E. coli ↔ MTB, orthogonal methods) cross-genome + intra-genome ≤60%-identity splits |

**Consequences for the methodology (binding):**

1. **All three elements now carry cross-genome generalization claims, at differing breadth** (P5, P8): promoters ≥7 phyla (broad); RBS 4 bacterial phyla + archaeal domain-transfer; **Rho a modest but genuine 2-phylum in vivo panel** (E. coli ↔ M. tuberculosis, orthogonal methods — BCM chemical inhibition vs RhoDUC genetic depletion), plus intra-genome ≤60%-identity splits in each. Rho's cross-phylum claim is stated as a **2-organism transfer**, not a broad sweep; *B. subtilis* (in vitro H-SELEX) and Streptomyces intrinsic contrasts are supporting. This **upgrades** the earlier "Rho = E. coli case study" scoping (Jun 24 strengthening).
2. The RBS Tier-2 negative is **SD vs UNSD-leadered within organism** (GC-clean, §5c), not SD-vs-leaderless (GC-confounded — leaderless concentrates in 65.5%-GC MTB while SD positives sit at 33–51% GC). SD labels use **organism-specific anti-SD + ΔG thresholds** (S. aureus extended-SD would be mislabeled by an E. coli rule). RBS is the **highest Phase-0 risk** element (6–8 bp SD diluted in a 60-token window) → pre-registered window-size/pooling sweep. Full audit in [`data/rbs_database/README.md`](../data/rbs_database/README.md).
3. **CDBProm** (~24M predicted promoters across ~6,000 organisms) is **software-predicted** and therefore **excluded from ground truth** per §12; it may serve only as a weak-label or baseline comparator, never as a positive label.

**Week-1 data-census gate (new — precedes the Phase-0 gate).** Before clustering, count — per element and per regime — the number of **independent ≤60%-identity clusters** that survive MMseqs2/CD-HIT deduplication. Pre-register a minimum (target: ≥30 held-out clusters spanning ≥2 phyla for a *cross-genome* unseen claim; ≥30 clusters for an *intra-genome* claim). Any element/regime below threshold is formally demoted to "case study, no generalization claim" **before** Phase 0, so a headline generalization number is never reported on a class that lacks the data to support it. This makes the §5e leakage argument and P8 contingent on a counted, committed census rather than an assumption.

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

**9d. Token-identity disambiguation — the trivial-specialization trap (P10, Zhong Wang).** The single most dangerous confound for the *biological-semantics* claim, and one the §12 GC-baseline and §9b structural partial-out do **not** catch. Because promoters/RBS/Rho are defined by **short conserved motifs that map to specific BPE tokens** (−10 TATAAT, −35 TTGACA, SD GGAGG, C-rich *rut*), "routing detects the element" can be 100% "routing detects the motif token" — the same trap Junho's work had to rule out for its structural experts (junhos_work.md §4 follow-up). We replicate that exact disambiguation against three hypotheses: **H₀** routing is load-balanced noise (≈ uniform), **H₁** routing is token-ID / surface-composition hashing (context adds nothing), **H₂** routing is context-conditional (function-aware). Procedure:

- Identify **element-defining tokens that also occur in non-element contexts** (e.g. GGAGG inside a CDS, TATAAT in a non-promoter window) — the analog of Junho's 1,248 multi-class token IDs.
- Estimate the **context-conditional mutual information I(expert ; element_class | token_id)** per layer, with a **within-token label-shuffle null** (shuffle class labels among instances of the *same* token ID). H₂ requires conditional MI significantly above the within-token null at the informative layers; H₁ predicts conditional MI ≈ null; H₀ predicts near-zero unconditional MI too (check P(expert | class) is non-uniform beyond the marginal-masked null of §9b).
- **Token-controlled DiD:** in the §9b ablation, the matched control must **hold token-identity fixed** (compare element-context vs non-element-context instances of the *same* token), so a positive DiD reflects context-conditional routing, not removal of a token's only router path.

**Scope (important, not a blanket disqualifier).** Token-ID/motif detection is **legitimate for RQ-A** (the detection/annotation head-to-head — BPROM also keys on the motif; a routing detector that finds TATAAT is a fair, useful detector). The trap bites **only RQ-B**: the claims in §2 that experts encode regulatory *function*. So 9d gates the biological-semantics / MoE-necessity narrative, **not** the detection comparison. If 9d fails (H₁), we report the regulatory experts honestly as motif-token detectors and retain the RQ-A detection results.

---

## 10. Optional retraining branch (validation purposes only — secondary)

The frozen path (§7) is primary. This branch is a **stretch / collaborative** control (training is heavy — Junho used 4× H200, ~50k steps, ~52B tokens — and is beyond intern scope per Zhong Wang; run only with mentor/Junho support). The code review concluded **full retraining is likely unnecessary**, so this branch exists to *stress-test the conclusions*, not to produce the deliverable.

**10a. Bug-correction control (robustness).** Zhong Wang's code review found a **math error in the upcycling**: weights were not rescaled by 2× to preserve the original distribution after the 50% expert-dropout step. Validation loss still reached expected levels by ~50k steps, so the curriculum is effective — but to confirm conclusions are robust to the bug, re-run upcycling **with the 2× rescale corrected** and check that expert specialization (the enrichment + DiD experts of §9) and the detection results reproduce. *Pass:* conclusions hold under the corrected model (P9, part 1).

**10b. Leakage-free control (closing the memorization question).** To answer Nic K's "Achilles heel" at the backbone level: produce a model variant whose **pretraining excludes the §5 validation genomes / sequence clusters** (continue-pretrain or upcycle from a leakage-free checkpoint). Re-run Regime B (unseen, ≤60% identity). *Pass:* unseen-split performance is preserved even when the validation data was never in pretraining → the result is generalization, not pretraining memorization (P9, part 2). *Fail:* performance depended on having seen the data in pretraining → report the leakage honestly and scope the claim down.

**Scope/honesty:** this branch is explicitly optional and secondary. If compute or collaboration is unavailable, the sequence-similarity holdout (§5b) is the practical leakage control and the frozen path stands on its own; the retraining branch only *strengthens* the leakage argument, it is not required for a valid result.

---

## 11. Generalization & transfer (RQ-C)

- **Primary (Nic K):** the **sequence-dissimilar unseen split** (Regime B, §5b, P8) — performance vs max-identity-to-training curve, against the GC-matched baseline. **Cross-genome for all three elements**, at differing breadth (§5f): promoters ≥7 phyla; RBS 4 bacterial phyla **+ a bacteria→archaea domain-transfer test** (*H. volcanii*); Rho a **2-phylum transfer (E. coli ↔ M. tuberculosis)** via orthogonal in vivo methods. Every element also gets an intra-genome ≤60%-identity split.
- **Secondary:** **leave-one-phylum-out** transfer (P5) — routing should drop less than dense embeddings on held-out phyla (Junho's AMI 0.15 vs 0.80), a second, independent MoE-advantage axis. Now available for **all three elements** (promoters broad, RBS 4-phylum, Rho 2-phylum E. coli↔MTB) — Rho's is stated as a 2-organism transfer, not a broad sweep.

---

## 12. Statistical rigor, controls, honesty commitments

- **Chance baseline = GC-content-matched profiles** (Zhong Wang / Rob Egan), not naive shuffling — "better than random" must mean better than GC-matched random.
- **Nulls:** label permutation; marginal-masked routing null; permutation-within-matched-strata (GC/length/position/BPE-phase).
- **Confounds controlled:** GC, length, distance-to-feature, strand, intergenic position, genic/intergenic context, BPE phase, **and token-identity / surface k-mer composition** (the H₁ "trivial-specialization" confound, §9d) — as covariates and in matched sampling; token-identity is additionally controlled by the *conditional*-MI test and the token-fixed DiD (§9d), since matched sampling alone cannot remove a confound that *is* the element's defining motif.
- **Multiple testing:** BH-FDR across all (layer, expert, class, genome) cells.
- **Class imbalance:** AUPRC/MCC primary; balanced probes; per-class reporting (rare Rho-terminators and σ-subtypes never hidden in an aggregate).
- **Ablation stated as a differential** (DiD vs matched control), never absolute degradation (§2 caveat).
- **Structural shortcut treated as a confound to defeat**, not a result (Tier-2 decoys, intergenic-not-CDS control, partial-out).
- **Pre-registered failure criteria** = the P0–P9 falsification rows.
- **Anti-overclaim:** no "expert *for* promoters" without a causal DiD against an intergenic control; no detection win reported from the Tier-1 setting; no generalization claim without the unseen (Regime B) split; P6/P7/P8 reported regardless of direction.
- **Comparison integrity / anti-circularity (added Jun 24 2026 — protects the head-to-head that is the project's main goal).** A baseline may never be benchmarked against labels it helped define:
  - **(a) RBS vs energy-based tools.** The SD-vs-UNSD split is defined by a ΔG rule (§5d), and two named RBS baselines — **Salis RBS Calculator v2.1** and **Free2Bind** — *are* ΔG / anti-SD-pairing models, so scoring them against ΔG-defined labels is circular. The **headline RBS head-to-head is therefore experimentally-grounded translation-start detection / localization** (predict the ribosome-profiling-confirmed TIS vs same-context non-translated start windows), which is ΔG-free and on which all tools (Prodigal start-calling, Salis, GO-MoE, dense gLMs) compete equally. The ΔG-defined **SD-vs-UNSD subtyping is reported separately**, with energy-based baselines marked as partially circular (agreement sanity-check, not a benchmark win). The ΔG split stays valid for the MoE-internal routing analysis (§9), where no external ΔG baseline is involved.
  - **(b) Promoters vs promoter-prediction tools.** PPD aggregates from many upstream sources; any entry traceable to a tool under benchmark (BPROM, bTSSfinder, iPromoter-2L/MULTiPly, G4PromFinder) is excluded. Only experimental-evidence entries are kept (dRNA-seq / TSS-mapping, RegulonDB strong-evidence) — so promoter tools are scored against experimental truth, never against predictions they (or their training lineage) produced.
  - **(c) Identical loci / identical protocol.** All dense gLMs and classical tools are evaluated on the **same held-out windows** under the §7 frozen-probing protocol (gLMs) or at default + tuned thresholds (classical), so any GO-MoE margin is attributable to the model, not to an easier evaluation set.

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
| **1** | **Dataset + split + data census** | Pull 5 core genome **sequences from NCBI** and **regulatory labels from the experimental DBs** (RegulonDB/PRODORIC/dRNA-seq/NET-seq — *not* NCBI, §5f); build 300 bp windows + same-context decoys; **run the per-element data-census gate (§5f) and demote any class below the cluster threshold to case-study status before Phase 0**; **cluster by sequence identity and commit the 80/20 ≤60% split**; lock all splits to files. | §5, §6 |
| **2** | **Phase 0 smoke test (gate)** | Frozen feature extraction; run the smoke test vs the GC-matched baseline on seen data (P0). Go/no-go; if no-go, run the §4 diagnosis tree. | §4, §7 |
| **3** | **Detection + baselines (RQ-A)** | Train linear-probe detectors; score head-to-head vs classical tools and dense LGMs on held-out genomes. | §3, §8, §9a |
| **4** | **MoE-necessity + generalization (RQ-B/C)** | Expert ablation (DiD, intergenic control, partial-out), upcycled-vs-scratch, promoter-vs-RBS acid test; **unseen (Regime B) generalization curve** (P8) + leave-one-phylum-out. | §9, §11 |
| **5** | **Stats, write-up, (optional) retraining** | GC-matched nulls + bootstrap CIs; annotation prototype or Phase-0 diagnosis; figures/slides. *If time + support:* launch §10 bug-fix / leakage-free controls. | §10, §12, §13 |

**Critical path & risk:** Weeks 1–2 gate everything; the Phase-0 smoke test is the explicit go/no-go. If behind by Week 4, protect the MoE-necessity experiments and the unseen-split generalization test (the two unique contributions), and defer the diversity sweep and the retraining branch. Even a Phase-0 *failure* yields a publishable diagnosis (§13.1), so there is a guaranteed deliverable.
