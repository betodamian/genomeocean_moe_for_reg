# Curated Rho-dependent Terminator Database

**Created:** 2026-06-24 (Rho strengthening pass)
**Status:** Sources verified; cross-genome panel established. Raw data to be fetched Week 1.

This directory lifts Rho-dependent terminators from a **single-genome E. coli case study**
to a **2-phylum in vivo cross-genome panel** using genuine genome-wide experimental data
mapped by **orthogonal methods**. It does not relax the experimental-evidence standard
(plan §12), does not change the project goal (GO-MoE vs classical tools / dense gLMs), and
states the cross-phylum claim at its true, modest breadth.

---

## Why this strengthening is rigorous (not an overclaim)

The earlier scoping ("Rho = E. coli only") rested on Peters 2012 being the sole genome-wide
Rho map. A verified literature check (Jun 24 2026) found two more in vivo-relevant sources:

1. **A second, orthogonal E. coli map** (Term-seq exact 3′-ends, GSE109766) that is
   *concordant* with Peters 2012 (BCM) — enabling **concordance-filtered high-confidence
   positives** and **bp-resolution boundaries** (the localization task, plan §8 T2).
2. **A second phylum** — *M. tuberculosis* (Botella 2022) — mapped by **genetic Rho
   depletion (RhoDUC)** rather than bicyclomycin. This matters: Mtb Rho is
   **BCM-resistant**, so the E. coli chemical method physically cannot be used; the genetic
   method is a genuinely independent in vivo readout, not a re-run of the same assay.

So the cross-genome claim rests on **two phyla mapped by two different in vivo methods**,
with ~1,000 (E. coli) and ~1,385 (MTB) sites — comparable depth. That is a real, if modest,
generalization panel. *B. subtilis* H-SELEX adds a third phylum but is **in vitro**, so it is
carried with an explicit evidence-type caveat and never merged into the in vivo positive set.

---

## Sources (see `sources.tsv` for the machine-readable table)

| Organism / phylum | Source | Method | Sites | Tier |
|---|---|---|---|---|
| *E. coli* / Gammaproteobacteria | Peters et al. 2012 | NET-seq + bicyclomycin readthrough | ~1,000 (region) | T1 in vivo |
| *E. coli* / Gammaproteobacteria | Term-seq 3′-end map 2018 (GSE109766) | Term-seq exact 3′-ends, BCM-validated | 144 (bp-resolution) | T1 in vivo |
| *M. tuberculosis* / Actinobacteria | Botella et al. 2022 (E-MTAB-11753, PMID 37096044) | Term-seq + RhoDUC genetic depletion | 1,385 | T1 in vivo |
| *B. subtilis* / Firmicutes | H-SELEX rut map (PMC12350095) | In vitro Rho-helicase SELEX (BsRho) | 600 | **T2 in vitro** |
| 13 species (E. coli, B. subtilis, Streptomyces…) | TERMITe 2024 | Term-seq intrinsic atlas | 165–1,214/sp | **decoy** (not Rho positives) |
| E. coli / B. subtilis / Salmonella | RhoTermPredict (PMC6407284) | algorithm-derived | small | **cross-check only** |

---

## Evidence hierarchy (consistent with plan §12)

| Tier | Evidence | Use |
|---|---|---|
| **T1 in vivo** | BCM readthrough; Term-seq 3′-ends; genetic Rho depletion | Headline positives |
| **T2 in vitro** | H-SELEX rut activation | Supporting only, evidence-type flagged |
| **decoy** | TERMITe intrinsic terminators | Same-context, same-position negative (Rho vs intrinsic) |
| **cross-check only** | RhoTermPredict (algorithm-derived) | Candidate cross-reference; never headline truth |

---

## How the panel is used in the methodology

- **High-confidence E. coli positives:** Peters 2012 (BCM) ∩ Term-seq 2018 (exact 3′-ends).
  Concordance raises label reliability; the Term-seq map supplies bp-resolution boundaries
  for localization scoring.
- **Cross-genome / leave-one-genome-out (P5, P8):** train on E. coli → test on MTB and the
  reverse. Two orthogonal in vivo methods make a "win" attributable to the Rho signal, not
  to a shared assay artifact. Stated as a **2-organism transfer**, not a broad sweep.
- **Intra-genome ≤60% identity split (P8):** within each of E. coli and MTB separately —
  *rut* sites are C-rich, unstructured, weakly conserved, so a sequence-dissimilar holdout
  is a genuine not-memorizing test even within one genome.
- **Tier-2 decoy (the regulatory-not-structural acid test, P7):** Rho vs **intrinsic**
  terminator (TERMITe) at matched 3′ positions — same context, same position, differ only in
  termination mechanism. Available in E. coli and B. subtilis; Streptomyces intrinsic gives an
  Actinobacterial contrast near MTB.
- **Anti-circularity (plan §12):** RhoTermPredict and any algorithm-derived calls are never
  headline ground truth; the head-to-head vs the classical Rho tool (RhoTermPredict) is scored
  against the *experimental* in vivo maps, not against RhoTermPredict's own predictions.

---

## Honest scope statement

The in vivo Rho cross-genome panel is **2 phyla** (E. coli, M. tuberculosis) — genuine and
method-orthogonal, but far thinner than the promoter panel (≥7 phyla) or the RBS panel
(4 phyla + archaea). Rho's cross-phylum generalization is therefore reported as a **2-organism
transfer with intra-genome dissimilar splits**, with per-element results never hidden inside
an aggregate (plan §12). This is a strict upgrade over the prior "E. coli case study only"
scoping, achieved without weakening rigor.

---

## Fetch instructions (Week 1)

```bash
# 1. E. coli Term-seq 3'-ends (bp-resolution, primary boundary source)
#    GEO GSE109766  (NAR 2018 46:6797). Interactive browser at Weizmann Institute.

# 2. M. tuberculosis Rho (2nd phylum, genetic depletion)
#    ArrayExpress E-MTAB-11753  (Botella 2022, PMID 37096044). 1,385 RD-TTS table in supp.

# 3. E. coli BCM Rho regions (Peters 2012)
#    Retrieve supplementary site table from the paper; verify GEO accession at fetch time.

# 4. B. subtilis H-SELEX rut (in vitro, supporting)
#    From PMC12350095 supplementary; flag every entry as T2_invitro in `notes`.

# 5. Intrinsic-terminator decoy (TERMITe)
#    GitHub + Zenodo (PMC12207403). Per-organism intrinsic terminator coordinates.
```
