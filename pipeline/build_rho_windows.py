#!/usr/bin/env python3
"""
Build Rho-terminator 300-bp windows + same-context decoys (research_plan §5c) for
all Rho label sources in data/rho_database/raw/.

Window design
-------------
300 bp, termination site centred at position 150 (0-indexed).
  - Upstream  (window[0:150])  = transcription direction 5' context:
      Rho positives: C-rich rut site where Rho loads
      Intrinsic decoys: hairpin + T-stretch
  - Centre   (window[150])    = termination 3'-end (or best approximation)
  - Downstream (window[150:300]) = sequence after termination

Centering strategy per source (different region sizes)
------------------------------------------------------
  "start"    : site_start == site_end (bp-resolution TTS) — centre directly
  "end"      : large RSR regions — use site_end (3' edge, closest to termination)
  "midpoint" : rut peaks or intrinsic terminators — use (start+end)//2

Labels
------
  positive_rho_t1      : in-vivo Rho termination sites (headline; T1)
                          E. coli  : Term-seq ∩ BCM concordance (151 bp-resolution sites)
                          M. tuberculosis : RhoDUC2 exact TTS (439 sites)
  positive_rho_t1_rsr  : MTB Rho-Sensitive Regions (Botella 2017, 303 sites; T1 in vivo
                          but lower resolution — site_end centering; supplementary)
  positive_rho_t2      : in-vitro rut sites (T2; not co-mingled with T1 in census)
                          E. coli  : EcRho H-SELEX rut peaks (4,189 sites)
                          B. subtilis : BsRho H-SELEX rut peaks (600 sites)
  decoy_intrinsic      : TERMITe intrinsic terminators (Tier-2 primary decoy)
                          E. coli  : 1,648 sites   B. subtilis : 3,998 sites
  decoy_intergenic     : random intergenic windows (Tier-1 sanity check, all organisms)

Output : data/datasets/rho/ALL.tsv
Run    : .venv/bin/python pipeline/build_rho_windows.py
"""
import csv, os, random
from collections import defaultdict

random.seed(20260625)
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RHO  = os.path.join(ROOT, "data/rho_database/raw")
GEN  = os.path.join(ROOT, "data/genomes")
OUT  = os.path.join(ROOT, "data/datasets/rho")
WIN  = 300
HALF = WIN // 2   # centre at index 150 (0-indexed)
MASK = 250        # Tier-1 decoys must be >= MASK bp from any termination site
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")

def rc(s): return s.translate(COMP)[::-1]

# ── source configuration ──────────────────────────────────────────────────────
# Each entry: (source_id, tsv_file, label, centering, filter_organism_or_None)
# centering: "start" | "end" | "midpoint"
# filter_organism: if not None, only load rows where organism == this value
SOURCES = {
    "ecoli_K12_MG1655": {
        "genome_dir": "ecoli_K12_MG1655",
        "phylum":     "Gammaproteobacteria",
        "positives": [
            ("ECOLI_BCM_RHO", "ecoli_rho_termseq_sites.tsv",
             "positive_rho_t1", "start", None),
            ("BSUB_HSELEX", "bsub_hselex_sites.tsv",
             "positive_rho_t2", "midpoint", "ecoli_K12_MG1655"),
        ],
        "decoys_intrinsic": [
            ("INTRINSIC_TERMITE", "intrinsic_termite_sites.tsv",
             "midpoint", "ecoli_K12_MG1655"),
        ],
    },
    "mtuberculosis_H37Rv": {
        "genome_dir": "mtuberculosis_H37Rv",
        "phylum":     "Actinobacteria",
        "positives": [
            ("MTB_RHODUC2", "mtb_rhoduc2_sites.tsv",
             "positive_rho_t1", "start", None),
            ("MTB_RHODUC",  "mtb_rhoduc_sites.tsv",
             "positive_rho_t1_rsr", "end", None),
        ],
        "decoys_intrinsic": [],   # no MTB intrinsic set → intergenic only
    },
    "bsubtilis_168": {
        "genome_dir": "bsubtilis_168",
        "phylum":     "Firmicutes",
        "positives": [
            ("BSUB_HSELEX", "bsub_hselex_sites.tsv",
             "positive_rho_t2", "midpoint", "bsubtilis_168"),
        ],
        "decoys_intrinsic": [
            ("INTRINSIC_TERMITE", "intrinsic_termite_sites.tsv",
             "midpoint", "bsubtilis_168"),
        ],
    },
}

# ── genome helpers (identical to build_rbs_windows.py) ───────────────────────

def load_fasta(path):
    seqs, cid, buf = {}, None, []
    with open(path) as fh:
        for ln in fh:
            if ln.startswith(">"):
                if cid: seqs[cid] = "".join(buf).upper()
                cid = ln[1:].split()[0]; buf = []
            else:
                buf.append(ln.strip())
    if cid: seqs[cid] = "".join(buf).upper()
    return seqs

def genome_paths(gen_dir):
    d = os.path.join(GEN, gen_dir)
    fna = next(os.path.join(d, f) for f in os.listdir(d) if f.endswith(".fna"))
    gff = next(os.path.join(d, f) for f in os.listdir(d) if f.endswith(".gff"))
    return fna, gff

def gene_spans(gff_path):
    spans = defaultdict(list)
    with open(gff_path) as fh:
        for ln in fh:
            if ln.startswith("#"): continue
            f = ln.rstrip("\n").split("\t")
            if len(f) < 8 or f[2] != "gene": continue
            spans[f[0]].append((int(f[3]) - 1, int(f[4])))
    for c in spans: spans[c].sort()
    return spans


# ── window extraction ─────────────────────────────────────────────────────────

def site_center(row, centering):
    """Compute 1-based genomic centre for the window."""
    s, e = int(row["site_start"]), int(row["site_end"])
    if centering == "start":   return s
    if centering == "end":     return e
    return (s + e) // 2          # midpoint

def extract_window(seq, center, strand):
    """300-bp window; centre at index 150 (0-indexed). Returns string or None."""
    t = center - 1   # 0-indexed
    if strand == "+":
        s, e = t - HALF, t + HALF
        if s < 0 or e > len(seq): return None
        return seq[s:e]
    else:
        s, e = t - HALF + 1, t + HALF + 1
        if s < 0 or e > len(seq): return None
        return rc(seq[s:e])

def find_contig(center, contigs):
    """Return first contig where center falls within bounds."""
    for cid, seq in contigs.items():
        if 1 <= center <= len(seq):
            return cid
    return None


# ── load TSV ─────────────────────────────────────────────────────────────────

def load_tsv(source_id, tsv_file, filter_org=None):
    path = os.path.join(RHO, source_id, tsv_file)
    if not os.path.exists(path):
        return []
    with open(path) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    if filter_org:
        rows = [r for r in rows if r.get("organism") == filter_org]
    return rows


# ── per-organism builder ──────────────────────────────────────────────────────

def build_org(org, cfg, writer):
    fna, gff = genome_paths(cfg["genome_dir"])
    contigs  = load_fasta(fna)
    spans    = gene_spans(gff)
    phylum   = cfg["phylum"]

    # ── positives ────────────────────────────────────────────────────────────
    rows_pos = []
    seen     = {}          # (center, strand) -> source_id
    ctr_0idx_by_contig = defaultdict(set)   # for decoy masking (0-indexed)
    n_miss   = 0

    for source_id, tsv_file, label, centering, filter_org in cfg["positives"]:
        for r in load_tsv(source_id, tsv_file, filter_org):
            try:
                ctr = site_center(r, centering)
                st  = r["strand"]
            except (ValueError, KeyError):
                continue
            if st not in ("+", "-"): continue
            key = (ctr, st)
            if key in seen: continue

            cid = find_contig(ctr, contigs)
            if cid is None: n_miss += 1; continue
            win = extract_window(contigs[cid], ctr, st)
            if win is None or len(win) != WIN or "N" in win:
                n_miss += 1; continue

            seen[key] = source_id
            wid = f"{org}_rho_{label}_{len(rows_pos)}"
            rows_pos.append((wid, org, phylum, label, source_id,
                             cid, ctr, st, win))
            ctr_0idx_by_contig[cid].add(ctr - 1)

    # ── Tier-2 intrinsic decoys ───────────────────────────────────────────────
    rows_t2  = []
    seen_dec = {}

    for source_id, tsv_file, centering, filter_org in cfg["decoys_intrinsic"]:
        for r in load_tsv(source_id, tsv_file, filter_org):
            try:
                ctr = site_center(r, centering)
                st  = r["strand"]
            except (ValueError, KeyError):
                continue
            if st not in ("+", "-"): continue
            key = (ctr, st)
            if key in seen_dec: continue

            cid = find_contig(ctr, contigs)
            if cid is None: continue
            win = extract_window(contigs[cid], ctr, st)
            if win is None or len(win) != WIN or "N" in win: continue

            seen_dec[key] = source_id
            wid = f"{org}_rho_decoy_intrinsic_{len(rows_t2)}"
            rows_t2.append((wid, org, phylum, "decoy_intrinsic", source_id,
                            cid, ctr, st, win))
            ctr_0idx_by_contig[cid].add(ctr - 1)   # mask intrinsic sites too

    # ── Tier-1 intergenic decoys ──────────────────────────────────────────────
    rows_t1  = []
    contigs_with_spans = [c for c in contigs if c in spans] or list(contigs)
    # Generate one intergenic decoy per positive (T1 + T1_rsr + T2 combined)
    target   = len(rows_pos)
    attempts = 0
    while len(rows_t1) < target and attempts < target * 100:
        attempts += 1
        cid = random.choice(contigs_with_spans)
        seq = contigs[cid]; L = len(seq)
        if L < WIN: continue
        pos = random.randint(HALF, L - HALF - 1)   # 0-indexed
        if any(s <= pos < e for s, e in spans.get(cid, [])): continue
        if any(abs(pos - t) < MASK for t in ctr_0idx_by_contig[cid]): continue
        strand = random.choice("+-")
        win = extract_window(seq, pos + 1, strand)
        if win is None or len(win) != WIN or "N" in win: continue
        wid = f"{org}_rho_decoy_intergenic_{len(rows_t1)}"
        rows_t1.append((wid, org, phylum, "decoy_intergenic", "DECOY",
                        cid, pos + 1, strand, win))

    for row in rows_pos + rows_t2 + rows_t1:
        writer.writerow(row)

    n_t1pos  = sum(1 for r in rows_pos if "t1" in r[3])
    n_t2pos  = sum(1 for r in rows_pos if "t2" in r[3])
    print(f"  {org:22s}  positives={len(rows_pos):5d} "
          f"(T1={n_t1pos} T2={n_t2pos})  "
          f"miss={n_miss}  "
          f"decoy_intrinsic={len(rows_t2)}  "
          f"decoy_intergenic={len(rows_t1)}")
    return len(rows_pos), len(rows_t2), len(rows_t1)


def main():
    os.makedirs(OUT, exist_ok=True)
    out_path = os.path.join(OUT, "ALL.tsv")
    fields   = ["id", "organism", "phylum", "label", "source_id",
                "contig", "center_pos", "strand", "window_seq"]
    tp = td_t2 = td_t1 = 0
    print("Building Rho windows ...")
    with open(out_path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(fields)
        for org, cfg in SOURCES.items():
            np_, nt2, nt1 = build_org(org, cfg, w)
            tp += np_; td_t2 += nt2; td_t1 += nt1
    print(f"\nTOTAL  positives={tp:,}  "
          f"decoy_intrinsic(Tier-2)={td_t2:,}  "
          f"decoy_intergenic(Tier-1)={td_t1:,}  -> {out_path}")


if __name__ == "__main__":
    main()
