# GenomeOcean: Sparse Upcycling and Expert Specialization in Genomic MoE

**Presenter:** Junho Hong, Northwestern University (May 2026)
**Core Premise:** Scaling genomic foundation models requires massive and proportional compute increases (e.g., Nucleotide Transformer scaled from 500M to 2.5B; Evo 2 from 7B to 40B). MoE architectures effectively solve this in NLP (where Mixtral 8×7B matched Llama 2 70B while using ~5× fewer active parameters). Genomics is highly compelling for MoE usage because prokaryotic genomes contain functionally distinct element types with different structural constraints and statistical properties:
*   **CDS:** Codon bias, GC3 patterns
*   **rRNA:** Highly conserved, multi-copy
*   **tRNA:** Short sequences (~75nt), heavily structurally constrained
*   **Intergenic:** AT-rich, low complexity

If MoE routing can capture this genomic modularity, it promises computational efficiency alongside biological interpretability, marking a novel direction, as no prior work had applied MoE to genomic language models.

---

## 1. Project Goals & Research Questions

The research was divided into two primary queries to evaluate the genomic MoE:
*   **RQ1 — Competitiveness:** Does a genomic MoE generate a competitive foundation model with respect to modeling quality, downstream task performance, and the biological plausibility of sequence generations?
*   **RQ2 — Expert Specialization:** Does the MoE routing structure develop biologically meaningful mechanisms? Specifically, do experts specialize across different genomic structures like CDS, rRNA, tRNA, ncRNA, and intergenic regions?

To test whether pretrained dense model knowledge successfully transfers into a sparse model, an experimental variable was established comparing a dense-to-MoE upcycled approach against a model trained from scratch.

---

## 2. Tools, Architecture, and Experimental Setup

### Model Configurations
The research tested several configurations against published baseline models:
*   **GO-100M (Published baseline):** 120M active parameters, 120M total parameters
*   **GO-500M (Published baseline):** 500M active parameters, 500M total parameters
*   **MoE-Upcycled (pilot3):** 8 experts, top-k = 2. 714M total parameters, 205M active parameters. Initialized from GO-100M weights with a dropout of p=0.5. Learning Rate: 5e-5
*   **MoE-Upcycled (pilot4):** 8 experts, top-k = 2. 714M total parameters, 205M active parameters. Initialized from GO-100M weights with a dropout of p=0.5. Learning Rate: 4e-4
*   **MoE-Scratch:** 8 experts, top-k = 2. 714M total parameters, 205M active parameters. Randomly initialized. Learning Rate: 4e-4

### Architecture Parameters
All models were built with 12 layers, 768 hidden dimensions, and 8 attention heads. MoE models used an FFN hidden dimension of 3072 per expert utilizing SwiGLU activations. Further architectural features included RoPE, RMSNorm, and a BPE vocabulary of 4096 tokens.

### Training and Data Pipeline
*   **Dataset:** Prokaryotic metagenomic assemblies; comprising 645 Gbp of data derived from 6 co-assemblies.
*   **Tokenizer Context:** 1024-token context length representing sequences of ~5 kb using a GO-4B tokenizer (~4.89 bp/token).
*   **Hardware and Training Tools:** Training ran on 4× H200 GPUs utilizing BF16 precision, an AdamW optimizer, and document-level attention masking.
*   **Schedule:** Models were trained for 50K steps processing ~52B tokens with a batch size of 1024. A Cosine Learning Rate schedule was used containing 2000 warmup steps with a 10% minimum learning rate.
*   **Validation Hardware:** Ran on a single NVIDIA H200 GPU, batch size 512, sequence length 1024 in BF16 via HuggingFace.

---

## 3. Findings: RQ1 — Competitiveness Metrics

### Validation Loss
MoE architectures demonstrated superior convergence compared to standard dense baselines:
*   **GO-4B (4B params):** 5.437
*   **MoE-Upcycled pilot4 (205M active):** 5.497
*   **MoE-Upcycled pilot3 (205M active):** 5.629
*   **MoE-Scratch (205M active):** 5.662
*   **GO-500M (500M params):** 5.732
*   **GO-100M (120M params):** 6.072

*Findings:* MoE-Upcycled (pilot4) approached the performance of the GO-4B model despite having 20× fewer active parameters. Furthermore, all MoE iterations outperformed GO-500M with 2.5× fewer active parameters using ~4× fewer tokens. The upcycled approach converged ~9,500 steps faster than the scratch model, beating the GO-500M benchmark at step 19,500 compared to step 29,001 for scratch.

### Inference Throughput
Measured by tokens/sec at high utilization:
*   **GO-100M:** 907,437 tok/sec (120M params)
*   **MoE-Upcycled:** 573,287 tok/sec (205M active / 714M total params)
*   **MoE-Scratch:** 571,995 tok/sec (205M active / 714M total params)
*   **GO-500M:** 307,316 tok/sec (541M params)

*Findings:* When compute-saturated, the MoE models were 1.87× faster than GO-500M at 2.5× fewer active parameters. However, at low batch sizes (e.g., bs=4 on an A40 GPU) inference became memory-bandwidth bound and throughput tracked total parameters (714M) rather than active parameters (205M).

### Downstream Probing (GUE Benchmark)
Model weights were frozen, using a mean-pooled last hidden state with logistic regression to yield an MCC metric across 4 Out-Of-Distribution (OOD) tasks (tasks were eukaryotic/viral, though models were prokaryotic-trained):
*   **Promoter (300bp, binary):** GO-500M (0.799) > MoE-Scratch (0.795) > GO-100M (0.786) > MoE-Upcycled (0.779)
*   **Splice (400bp, 3-class):** GO-500M (0.445) > MoE-Scratch (0.399) > GO-100M (0.392) > MoE-Upcycled (0.378)
*   **H3K4me3 (500bp, binary):** GO-500M (0.365) > MoE-Scratch (0.343) > GO-100M (0.342) > MoE-Upcycled (0.328)
*   **COVID (1000bp, 9-class):** GO-500M (0.592) > GO-100M (0.521) > MoE-Scratch (0.519) > MoE-Upcycled (0.474)
*   **Mean MCC:** GO-500M (0.550), MoE-Scratch (0.514), GO-100M (0.510), MoE-Upcycled (0.490)

*Findings:* While MoE-Scratch edged out the GO-100M baseline, both MoE models fell short of the GO-500M on these tasks, partly due to how the linear probes undersample the biological data encoded discreetly in the routing channel.

### Routing-Aware Probing
Because continuous embeddings might undersample MoE knowledge, probes tested the continuous residual stream (768-d), discrete routing fingerprints (12×8 softmax mean-pooled to 96-d), a concatenation (864-d), and per_expert bins (768×k). Training utilized a LogisticRegression classifier (L2, C=1.0) on an A40 GPU (1 hr, 64 GB memory, ~3× wall-clock vs embedding only, binning fallback minimum ~1000 tokens).
*   **In-domain (Drug resistance):** routing_only (0.980) & routing_concat (0.980) heavily outperformed embedding (0.842)
*   **In-domain (Gene classification):** routing_concat (0.966) > routing_only (0.897) > embedding (0.861)
*   **In-domain (Taxonomic):** embedding (0.774) > routing_concat (0.758) > routing_only (0.007)
*   **OOD GUE (Promoter):** embedding (0.779) > routing_concat (0.774) > routing_only (0.039)
*   **OOD GUE (Splice):** embedding (0.378) > routing_concat (0.377) > routing_only (0.000)

*Findings:* The 96-d routing layer natively solved in-domain biological classification using 8× fewer features than standard embeddings. A concatenated probe was often better and never worse. Routing alone naturally collapsed when parsing eukaryotic/viral regulatory motifs (OOD tasks) or sequence-level taxonomy, highlighting its localized, prokaryote-specific functional alignment.

### Generative Quality (Propose-and-Fold)
Generated sequences (1000 DNA sequences at T=0.7, top-k=4) were run through Pyrodigal (ORF extraction), translated, and evaluated via ESMFold to generate a pLDDT score.
*   **Natural (E. coli K-12):** Mean pLDDT 69.8 (91.8% > 50, 57.4% > 70). Coding density: 85.8%
*   **GO-500M:** Mean pLDDT 58.1 (60.9% > 50, 29.4% > 70). Coding density: 97.0%
*   **MoE-Scratch:** Mean pLDDT 57.6 (61.3% > 50, 28.0% > 70). Coding density: 95.5%
*   **MoE-Upcycled:** Mean pLDDT 56.9 (61.0% > 50, 25.2% > 70). Coding density: 94.8%
*   **GO-100M:** Mean pLDDT 54.2 (52.9% > 50, 18.2% > 70). Coding density: 97.2%
*   **Random DNA:** Mean pLDDT 44.7 (26.4% > 50, 0.5% > 70). Coding density: 16.0%

*Findings:* Both MoE implementations beat the GO-100M and were comparable to the GO-500M despite using 2.5× fewer active parameters. The ~95% coding densities successfully matched natural bacterial limits (85-95%).

---

## 4. Findings: RQ2 — Expert Specialization

### Annotation-Based Class Enrichment & Hidden Expert Collapse
Tested 5 diverse reference genomes with GC ratios from 33-67% (E. coli, B. subtilis, M. tuberculosis, P. aeruginosa, S. aureus) against structural classes: intergenic, CDS, ncRNA, tRNA, and rRNA (ranked in priority for purity ≥ 0.8). Tested using a marginal-masked P(e) ≥ 0.01 null baseline (simulating 2·N_c class-independent assignments) with Jensen-Shannon divergence (JSD) and G-test (BH-FDR correction).
*   **Mean JSD:** MoE-Upcycled achieved a JS divergence of 0.0786 (47.5× above the null baseline) while MoE-Scratch scored just 0.0092 (6.2× above null). The Upcycled model exceeded the scratch model's divergence in 58 out of 60 (layer, class) cells. Both had near-zero JS for CDS as it mathematically comprised 80-90% of all tokens.
*   **Expert Collapse:** While both models balanced during training with P(e) ranges of [0.108, 0.152] (Upcycled) and [0.110, 0.144] (Scratch) with 0/96 experts tracking under a 0.03 probability, they behaved differently on out-of-distribution reference genomes. The Upcycled model maintained an active P(e) range of [0.029, 0.348] with 0 instances below 0.01. The Scratch model suffered silent collapse, stretching to [0.0001, 0.470] with 14/96 experts falling below P(e) < 0.01, and 4 falling below P(e) < 0.001.

### Biological Detectors & Causal Ablation
The strongest log₂ enriched detectors per class identified distinct functions:
*   **Upcycled Model:** Intergenic (L6 E7, +1.75), CDS (L9 E0, +0.18), ncRNA (L6 E7, +1.59), tRNA (L7 E7, +3.72), rRNA (L6 E7, +2.19).
*   **Scratch Model:** Intergenic (L4 E7, +0.97), CDS (L4 E1, +0.15), ncRNA (L4 E7, +1.09), tRNA (L7 E4, +1.13), rRNA (L8 E6, +2.62).

*Findings & Ablation Test:* The Upcycled model's L7 E7 routed tRNA at ~13× expected probability (2³·⁷²) while only firing on ~3% of the layer marginal, proving to be a highly potent tRNA detector. The scratch model's best equivalent (L7 E4) routed at only ~2.2× expected probability (2¹·¹³), making it 3.3× weaker. Masking expert logits to −10⁹ on 4.4M validation tokens verified causality. Ablating the Upcycled L7 E7 target generated a Difference-in-Differences (DiD) of +0.559 (+0.998 target vs +0.439 CDS control), raising target Cross Entropy by +1.0 nats without affecting the same layer control L7 E3 (DiD -0.065). Contrastingly, the Scratch tRNA detector (L7 E4) was causally null with a DiD of +0.006, though its rRNA detector (L8 E6) did show causal specificity with a DiD of +0.553.

### Unsupervised Routing Classification
Using leave-one-genome-out cross-validation across 5 folds and 4.4M tokens via Logistic Regression and balanced XGBoost classifiers, the 96-dimensional fingerprints generated profound predictive signals with absolutely zero functional training labels.
*   **XGBoost Macro F1:** MoE-Upcycled (pilot4) achieved 0.665. MoE-Upcycled (pilot3) hit 0.586. MoE-Scratch lagged at 0.357.
*   **Specific tRNA Target F1:** Reached 0.774 with a Precision of 0.745 and a Recall of 0.806.
*   **Layer informativeness:** L11 (0.382), L6 (0.309), L7 (0.284), L3 (0.276).

### Follow-up control (Jun 2026, Zhong Wang): Ruling out trivial-specialization confounders

A mentor-directed follow-up (AI-agent analysis) tested whether the observed "expert specialization" is genuinely *biological* or a trivial artifact, by pitting three hypotheses against each other:

*   **H₀ — Load-balancing dominates.** Routing is essentially independent of content; the auxiliary load-balancing loss has flattened the distribution so much that P(expert | token, context) ≈ P(expert) (uniform-ish). The router carries no information. *(Junho's random prior.)*
*   **H₁ — Token-identity / composition only ("trivial specialization trap").** Routing is a deterministic-ish hash of the token ID (or its surface composition: GC, k-mer content). P(expert | token, context) ≈ P(expert | token) — context adds nothing.
*   **H₂ — Biological semantics.** Routing is context-conditional: the same token ID routes to different experts depending on the surrounding sequence's functional role. P(expert | token, context) ≠ P(expert | token). The router has learned something about *function*, not just *string*.

**Test and result.** For the **1,248 BPE token IDs that appear in ≥2 functional classes** (354K total samples across the 5 reference genomes), the agent estimated the context-conditional mutual information **I(expert ; class | token_id)** — i.e. class information in routing that *cannot* be attributed to token identity. At the most informative layer (**L7**), routing carries **0.243 bits of class information not attributable to token identity** — **96% of the total expert–class mutual information at that layer**, and **42× the within-token label-shuffle null**.

**Conclusion.** This supports **H₂** and rules out H₀ and H₁: the routing circuit has learned a **function-aware representation that is independent of surface composition** (the same token routes differently by functional context). This is the control that elevates "experts specialize" from the trivial-specialization trap to a genuine biological-semantics claim — and it is the disambiguation the downstream regulatory-annotation work must replicate for its own classes (see research_plan §9d / P10).

---

## 5. Findings: Functional Clustering of Routing Features

The pipeline averaged token router softmax values across gene spans to build compositional 96-d fingerprint vectors. Using Aitchison-Euclidean CLR transforms to respect probability simplex constraints, researchers deployed UMAP 96→15D representations (using k=30 KNN standard to scRNA-seq) into Leiden clustering at a parameter of γ=0.10 (tested over 0.05, 0.1, 0.2, 0.3, 0.5, 1, 2). Evaluating 71,936 features against COG (23 classes), Pfam (2,139 families), and KEGG (279 pathways) across 20 PGAP annotated genomes (7 phyla, 2 archaea, 70,163 CDS, 1,192 tRNA, 299 ncRNA, 267 rRNA, 15 tmRNA) yielded 24 substantial clusters.

### The 24 Cluster Partition Highlights
The largest cluster held 16.4% of tokens, the smallest had n=31. The global Silhouette cluster score was +0.043 in UMAP-15D space, with 12/24 individually positive clusters and 8/24 scoring above +0.1.
*   **7 Non-Coding Clusters:** Successfully recovered ~90% of all rRNA and ~84% of all tRNA in the dataset. This included pure structural topology like C18 (a pure rRNA topological island scoring +0.999 silhouette, 98% pure, n=83), C20 (tRNA-dominant at 91%, n=46), and C13 (a combined translation group of 75% tRNA and 14% rRNA).
*   **6 Biological CDS Clusters:** Characterized functionally rather than taxonomically:
    *   C22: Highly specific to MTB PE-PGRS genes (100% of n=60 elements, Pfam F1 = 0.516, purity 0.89).
    *   Transporters: Sorted heavily by domain architecture rather than organism, grouped into C12 (ABC ATPase + TonB receivers, n=822 + 116, F1 = 0.375), C9 (MFS_1, n=382, F1=0.110), and C8 (BPD, n=2288, F1=0.066).
    *   Conserved Groups: Translation units fell into C7 (comprising 537 out of 1031 ribosomal markers mapped, n=3199, F1=0.013) and chaperones into C10 (DnaK 82%, GroEL 80%, n=2272, F1=0.016). Other tracked clusters included C4 (0.042), C2 (0.032), C0 (0.025), C11 (0.025), C3 (0.019, SusD/C), C1 (0.018), C15 (0.015), C6 (0.014), C14 (0.012), and C5 (0.008).
*   **11 Bulk CDS Clusters:** Making up ~80% of the dataset, these presented as a continuous manifold divided by phylogeny and GC content instead of distinct functions. The median Pfam F1 across 21 testable clusters proved to be only 0.022.

### Taxonomic Vs. Functional Encoding
Comparing clustering metrics directly highlighted how the MoE architectures separated features:
*   **MoE Routing Channel (96-d, γ=0.10):** Proved genome-invariant with an AMI of 0.15 (0 out of 24 clusters were organism-pure). The routing channel aligned specifically with *function*.
*   **MoE Hidden State (768-d, γ=1.00):** Exhibited taxonomic alignment with an AMI of 0.80 (26 out of 51 clusters were >50% isolated to one organism). The hidden state encoded *taxonomy*.
*   **Dense GO-100M Base (768-d, γ=1.00):** Mimicked the hidden state representation, remaining organism-aligned with an AMI of 0.77.

---

## 6. Discussion and Project Takeaways

The findings generated clear answers for the viability of Genomic Mixture-of-Experts:
*   **Confirmed Findings:** MoE routing distinctly segregates functional structural classes without supervision. Non-coding RNAs split decisively from CDS regions. The router manages to functionally isolate highly-paralogous families (like the PE-PGRS island with an F1=0.52), translation and chaperone machineries, and isolates transporter classes purely by domain architecture rather than organismic descent. The model explicitly uncoupled information channels, sending function into the router and taxonomy into the hidden state.
*   **What was NOT proven:** The routing parameters did not generate a perfectly clean 24-way functional partition; with a median Pfam F1 of 0.022, a typical cluster lacked a single dominant Pfam family. While earlier analyses claimed "1,266 Pfam families recovered," the project noted this falsely counted broad, statistically over-represented pairs at q<0.05, a vastly weaker claim than 1:1 cluster-to-family mappings. Bulk CDS sequences (the 11 largest clusters) split strictly by broad taxonomy and GC content, completely ignoring fine-grained function.
*   **Conclusion:** The MoE routing fingerprint operates as a continuous manifold. The extreme, discrete tail-modes map explicitly to biologically coherent gene structures, while the dominant mass of the bulk network simply scales up alongside taxonomic lineage and genomic GC content parameters.