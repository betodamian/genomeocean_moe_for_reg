# Research Plan — Plain Language Summary

A simplified, step-by-step explanation of the full research plan in [`research_plan.md`](research_plan.md), written for a general audience.

---

## What This Research Is Trying to Do (The Goal)

Bacteria have DNA that includes not just genes, but also **control switches** — short stretches of DNA that tell the cell *when* and *how much* to read each gene. These are called **regulatory elements**. Finding them automatically is important for understanding how bacteria work, how they resist antibiotics, etc.

This project asks: **Can a new AI model (GenomeOcean-MoE) find these control switches better than existing tools?**

---

## What's Special About This AI Model

Most AI models for DNA are "dense" — every part of the model processes every input. GenomeOcean-MoE uses a **Mixture-of-Experts (MoE)** design: it has 8 specialized sub-networks ("experts") and for each chunk of DNA, only 2 of the 8 experts activate. Think of it like a hospital with specialists — a tRNA chunk goes to the "structural RNA specialist," an intergenic region goes to the "non-coding specialist," etc.

The model was trained on **645 billion base pairs** of bacterial DNA from metagenomics (environmental samples), without being told anything about regulatory elements.

---

## Previous Work That This Builds On (Junho's Findings)

Junho Hong (May 2026) showed that:

1. The MoE model **naturally developed expert specialists** without supervision — one expert became a tRNA detector (~13× more likely to fire on tRNA than random), proven by showing that silencing it hurt tRNA prediction far more than anything else.
2. The **routing channel** (which expert fires) encodes **biological function**, while the **hidden states** (the main representations) encode **taxonomy** (what species it is). These are two separate information channels.
3. The routing fingerprint (a 96-number summary of which experts fired) could classify bacterial gene types **with 8× fewer features** than the full model output.
4. A follow-up confirmed this is **genuinely context-aware** — the same DNA word routes differently depending on whether it's in a functional vs. non-functional context. It's not just memorizing which words look like regulatory elements.

---

## The Three Control Switches Being Studied

| Element | What it does |
|---|---|
| **σ-dependent Promoters** | Where bacteria start reading a gene — the "on switch" |
| **Shine-Dalgarno RBS (Ribosome Binding Sites)** | A short sequence that tells ribosomes where to start making protein |
| **Rho-dependent Terminators** | Where bacteria stop reading — the "off switch" using a special protein called Rho |

---

## The Core Questions

1. **Can the model detect these elements at all?** (Does its routing channel carry relevant signal?)
2. **Does it beat existing specialist tools and other DNA AI models?**
3. **Is the advantage actually due to the expert specialists** — or would any AI do as well?
4. **Does it work on DNA sequences the model has never seen before** — or is it just memorizing?

---

## The Data Going In

**Where the labels (ground truth) come from** — this is all experimental data from published research, never predictions:

- **Promoters:** 129,148 experimentally verified promoters from the Prokaryotic Promoter Database (PPD), covering 63+ bacterial species across 7+ phyla.
- **RBS:** A hand-curated database built from ribosome-profiling experiments in 5 organisms (*E. coli*, *M. tuberculosis*, *B. subtilis*, *S. aureus*, *Caulobacter*) + archaea (*Haloferax volcanii*). ~11,000 positive examples.
- **Rho terminators:** Experimental maps from *E. coli* (~1,000 sites, two independent methods) and *M. tuberculosis* (~1,385 sites, a different method). Importantly, *M. tuberculosis* Rho is resistant to the drug that blocks *E. coli* Rho, so these are truly independent confirmations.

**Critically:** NCBI (the public genome database) has **zero useful regulatory labels** — the plan explicitly calls this out. All labels come from the specialized experimental databases above.

---

## How the Data Is Split (Train vs. Test)

This is the key rigor upgrade. Instead of randomly splitting data, the team:

1. Groups all sequences by **similarity to each other** (using a clustering algorithm).
2. Puts 80% in training, 20% in testing — but the 20% **must be <60% identical** to anything in the 80%. So the test set is genuinely novel sequence.

This creates two test regimes:
- **Regime A (sanity check):** Test on sequences similar to training. "Can it do things it's seen before?"
- **Regime B (the real test):** Test on truly dissimilar sequences. "Does it generalize, or is it memorizing?"

---

## What the Model Outputs

From each 300-base-pair window of DNA, the frozen model extracts:
- **Embeddings** (768 numbers — the standard AI representation)
- **Routing fingerprint** (96 numbers — a summary of which experts fired at each layer)
- **Combined** (864 numbers)

Then a simple logistic regression (a very basic classifier) uses these numbers to predict: *is this a promoter? An RBS? A Rho terminator? Or none of the above?*

---

## What Gets Measured (Metrics)

| Metric | What it measures |
|---|---|
| **AUPRC** (primary) | How well it detects elements when they're rare — a more honest measure than accuracy |
| **MCC** | Overall classification quality accounting for imbalanced classes |
| **Boundary-F1 / bp error** | How precisely it locates the element (not just "yes there's a promoter," but *where* exactly) |
| **DiD (Difference-in-Differences)** | Whether silencing a specific expert hurts detection of an element *more* than it hurts a control — proving that expert is causally responsible |
| **Conditional Mutual Information** | Whether routing is context-aware (function-sensitive) vs. just word-sensitive |

Everything is computed with confidence intervals and multiple-testing correction.

---

## The Key Experiments (In Order)

**Week 1 — Build and count the data.** Assemble all experimental labels, build the 300bp windows, cluster by similarity, commit the 80/20 split. Count how many independent clusters exist per element — if there aren't enough, demote that element to "case study" before running anything.

**Week 2 — Smoke test (the gate).** On *seen* data only: can the frozen model's routing fingerprint distinguish a real promoter from a fake one, better than a GC-content-matched random baseline? If yes → proceed. If no → switch to diagnosing *why* (bad window size? wrong feature view?) instead of building the annotator.

**Week 3 — Detection vs. baselines.** Train the linear probe and compare head-to-head against: classical bioinformatics tools (BPROM for promoters, Prodigal for RBS, RhoTermPredict for terminators) and other DNA AI models (Nucleotide Transformer, DNABERT-2, Evo 2, ProkBERT). All evaluated on the same held-out windows.

**Week 4 — MoE-necessity + generalization.** Silence specific experts and measure the damage (DiD). Compare the specialized upcycled model vs. a model trained from scratch (which Junho showed develops weaker specialists). Test on the truly novel sequences (Regime B). Plot performance vs. sequence similarity to training.

**Week 5 — Stats and write-up.** Add bootstrap confidence intervals, GC-matched baselines. Build the annotation prototype if Phase 0 passed. Optionally: fix a known math bug in the upcycling code and rerun to make sure results hold.

---

## The Pre-Registered Predictions (What Would Prove It Wrong)

The team wrote down 11 specific, falsifiable predictions before touching the validation data. Key ones:

- The routing fingerprint beats GC-matched chance on *seen* data (gate — if this fails, stop and diagnose)
- Routing + embedding beats embedding alone, consistently
- The MoE beats all classical tools and other AI models on novel sequences
- At least one expert per element is causally responsible (silencing it specifically hurts that element)
- The same expert routing patterns appear on sequences the model never saw (generalization, not memorization)
- Routing is context-aware for these elements — the same DNA word routes differently when functional vs. not (ruling out "it's just detecting the motif word")

All results are reported regardless of direction — if a prediction is wrong, that's reported honestly.

---

## What Gets Delivered

1. A **go/no-go ceiling report** from the smoke test (guaranteed deliverable even if everything else fails)
2. If successful: a **working annotation tool** that outputs predicted locations of promoters, RBS, and terminators — with confidence scores and an "I'm not sure" option for low-confidence calls
3. A **public benchmark** with all data, splits, and labels so anyone can reproduce the comparison
4. Scripts to regenerate every number in one command

---

## Bottom Line

The project is testing whether a specialized AI — one that has developed internal "expert" networks tuned to different parts of bacterial genomes — can find gene control switches more accurately and more generally than any existing tool, and whether the advantage is genuinely because of those specialized experts rather than just model size or memorization.
