#!/usr/bin/env python3
"""
Build RBS 300-bp windows + same-context decoys (research_plan §5c) for all
experimental TIS sources in data/rbs_database/raw/.

Window design
-------------
300 bp, TIS centred at position 150 (0-indexed), so window[150:153] = start codon.
- Positive_SD      : Ribo-seq-confirmed TIS with upstream SD motif (primary positive).
- Positive_UNSD    : Ribo-seq-confirmed TIS without SD motif — Tier-2 decoy for the
                     SD-detection task; secondary positive for TIS-detection (§5c).
- Decoy_intergenic : Intergenic 300-bp window >= MASK bp from any TIS — Tier-1 sanity
                     check only (gameable by the intergenic shortcut; §5c).

SD labelling (v1)
-----------------
Search the 22 nt upstream of each start codon (window[128:150]) for SD-core motifs.
Uses organism-specific patterns and position windows (S. aureus extended SD at −12..−20).
Any hit ≥ 4 nt → SD; no hit → UNSD (includes leaderless, which has no upstream SD).

Liftover (ECOLI_RIBORET / ECOLI_TETRP)
---------------------------------------
These TSVs carry BW25113 coordinates. We first try the position directly in MG1655
(most genes are syntenic); on failure we look up gene_name in the GFF. Positions that
fail both checks are dropped and counted as liftover_miss.

Dedup
-----
Multiple sources per organism (BSUB_RIBOSEQ+BSUB_SPORE, SAUR_RIBORET+SAUR_EXSD) are
merged and deduplicated by (tis_pos, strand); first source encountered wins.

Output: data/datasets/rbs/ALL.tsv
Run:    .venv/bin/python pipeline/build_rbs_windows.py
"""
import csv, os, random, re
from collections import defaultdict

random.seed(20260625)
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RBS  = os.path.join(ROOT, "data/rbs_database/raw")
GEN  = os.path.join(ROOT, "data/genomes")
OUT  = os.path.join(ROOT, "data/datasets/rbs")
WIN  = 300
HALF = WIN // 2      # TIS at index 150 (0-indexed)
MASK = 250           # Tier-1 decoys must be >= MASK bp from any TIS
STARTS = {"ATG", "GTG", "TTG"}
COMP   = str.maketrans("ACGTNacgtn", "TGCANtgcan")

def rc(s): return s.translate(COMP)[::-1]

# Organism → (genome_dir, phylum, [(source_id, tsv_file, coord_space)])
# coord_space: "native" | "bw25113"
ORGS = {
    "ecoli_K12_MG1655": (
        "ecoli_K12_MG1655", "Gammaproteobacteria",
        [("ECOLI_RIBORET", "ecoli_riboret_tis_sites.tsv", "bw25113"),
         ("ECOLI_TETRP",   "ecoli_tetrp_tis_sites.tsv",   "bw25113")],
    ),
    "ecoli_BL21_DE3": (
        "ecoli_BL21_DE3", "Gammaproteobacteria",
        [("ECOLI_BL21_RIBORET", "ecoli_bl21_riboret_tis_sites.tsv", "native")],
    ),
    "bsubtilis_168": (
        "bsubtilis_168", "Firmicutes",
        [("BSUB_RIBOSEQ", "bsub_riboseq_tis_sites.tsv", "native"),
         ("BSUB_SPORE",   "bsub_spore_tis_sites.tsv",   "native")],
    ),
    "mtuberculosis_H37Rv": (
        "mtuberculosis_H37Rv", "Actinobacteria",
        [("MTB_RIBOSEQ", "mtb_riboseq_tis_sites.tsv", "native")],
    ),
    "ccrescentus_NA1000": (
        "ccrescentus_NA1000", "Alphaproteobacteria",
        [("CAULO_RIBOSEQ", "caulo_riboseq_tis_sites.tsv", "native")],
    ),
    "saureus_HG001": (
        "saureus_HG001", "Firmicutes",
        [("SAUR_RIBORET", "saur_riboret_tis_sites.tsv", "native"),
         ("SAUR_EXSD",    "saur_exsd_tis_sites.tsv",    "native")],
    ),
    "hvolcanii_DS2": (
        "hvolcanii_DS2", "Archaea",
        [("HVOLC_RIBOSEQ", "hvolc_riboseq_tis_sites.tsv", "native")],
    ),
}

# SD search per organism:
#   (upstream_start_in_win, upstream_end_in_win, set_of_patterns)
# window[upstream_start:upstream_end] is searched for any SD pattern.
# Default: -3 to -22 (window[128:150]); S. aureus extended SD reaches to -20.
# Pattern set: SD core motifs ≥4 nt matching the 16S 3' anti-SD complement.
_SD_DEFAULT = {"GGAGG", "AGGAG", "GAGGA", "GAGG", "AGGA", "GGAG", "AGGG"}
_SD_SAUR    = _SD_DEFAULT | {"AAGGAGG", "AAGGAG", "AAGGA", "AAGG"}
SD_CFG = {
    "default":       (128, 150, _SD_DEFAULT),
    "saureus_HG001": (128, 150, _SD_SAUR),
}

def sd_label(win, org):
    cfg = SD_CFG.get(org, SD_CFG["default"])
    upstream = win[cfg[0]:cfg[1]]
    patterns = cfg[2]
    if any(p in upstream for p in patterns):
        return "SD"
    return "UNSD"


# ── Genome helpers ────────────────────────────────────────────────────────────

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
    """contig -> sorted list of (start0, end0_exclusive) for gene features."""
    spans = defaultdict(list)
    with open(gff_path) as fh:
        for ln in fh:
            if ln.startswith("#"): continue
            f = ln.rstrip("\n").split("\t")
            if len(f) < 8 or f[2] != "gene": continue
            spans[f[0]].append((int(f[3]) - 1, int(f[4])))
    for c in spans: spans[c].sort()
    return spans

def gff_gene_tis(gff_path):
    """name/locus_tag -> (tis_pos_1based, strand) from CDS features (for BW25113 liftover)."""
    result = {}
    with open(gff_path) as fh:
        for ln in fh:
            if ln.startswith("#") or "\t" not in ln: continue
            cols = ln.rstrip("\n").split("\t")
            if len(cols) < 9 or cols[2] != "CDS": continue
            s, e, st = int(cols[3]), int(cols[4]), cols[6]
            tis = s if st == "+" else e
            for key in ("Name", "gene", "locus_tag"):
                m = re.search(rf'{key}=([^;]+)', cols[8])
                if m:
                    result.setdefault(m.group(1), (tis, st))
    return result


# ── Window extraction ─────────────────────────────────────────────────────────

def codon_at(seq, tis_pos, strand):
    """3-nt start codon at tis_pos (1-based) in the given strand direction."""
    t = tis_pos - 1
    if strand == "+":
        return seq[t:t + 3] if t + 3 <= len(seq) else ""
    return rc(seq[t - 2:t + 1]) if t >= 2 else ""

def extract_window(seq, tis_pos, strand):
    """300-bp window, TIS at index 150 (0-indexed). Returns window string or None."""
    t = tis_pos - 1
    if strand == "+":
        s, e = t - HALF, t + HALF
        if s < 0 or e > len(seq): return None
        return seq[s:e]
    else:
        s, e = t - HALF + 1, t + HALF + 1
        if s < 0 or e > len(seq): return None
        return rc(seq[s:e])

def find_contig(tis_pos, strand, contigs):
    """Return first contig name that has a valid start codon at tis_pos, or None."""
    for cid, seq in contigs.items():
        if 1 <= tis_pos <= len(seq) and codon_at(seq, tis_pos, strand) in STARTS:
            return cid
    return None


# ── Per-organism builder ──────────────────────────────────────────────────────

def load_tsv(source_id, tsv_file):
    path = os.path.join(RBS, source_id, tsv_file)
    if not os.path.exists(path):
        return []
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))

def build_org(org, gen_dir, phylum, sources, writer):
    fna, gff = genome_paths(gen_dir)
    contigs = load_fasta(fna)
    spans   = gene_spans(gff)

    needs_liftover = any(cs == "bw25113" for _, _, cs in sources)
    gff_tis = gff_gene_tis(gff) if needs_liftover else {}

    # ── collect and deduplicate TIS ──────────────────────────────────────────
    seen  = {}   # (tis_pos, strand) -> source_id
    tis_list = []
    n_loaded = liftover_miss = 0

    for source_id, tsv_file, coord_space in sources:
        rows = load_tsv(source_id, tsv_file)
        n_loaded += len(rows)
        for r in rows:
            try:
                p = int(r["tis_pos"])
            except (ValueError, TypeError, KeyError):
                continue
            st = r.get("strand", "")
            if st not in ("+", "-"):
                continue

            if coord_space == "bw25113":
                # Try direct position in MG1655
                cid_direct = find_contig(p, st, contigs)
                if cid_direct is None:
                    # Fallback: GFF gene-name lookup
                    gname = r.get("gene_name", "")
                    if gname and gname in gff_tis:
                        gff_p, gff_st = gff_tis[gname]
                        if gff_st == st and find_contig(gff_p, st, contigs):
                            p = gff_p
                        else:
                            liftover_miss += 1; continue
                    else:
                        liftover_miss += 1; continue

            key = (p, st)
            if key not in seen:
                seen[key] = source_id
                tis_list.append((p, st, source_id))

    # ── build positive windows ────────────────────────────────────────────────
    rows_pos = []
    tis_0idx_by_contig = defaultdict(set)   # contig -> set of 0-indexed TIS positions
    n_miss = 0

    for tis_pos, strand, source_id in tis_list:
        cid = find_contig(tis_pos, strand, contigs)
        if cid is None:
            n_miss += 1; continue
        win = extract_window(contigs[cid], tis_pos, strand)
        if win is None or len(win) != WIN:
            n_miss += 1; continue
        if win[HALF:HALF + 3] not in STARTS or "N" in win:
            n_miss += 1; continue

        sc  = win[HALF:HALF + 3]
        slab = sd_label(win, org)
        label = f"positive_{slab}"
        wid = f"{org}_rbs_{label}_{len(rows_pos)}"
        rows_pos.append((wid, org, phylum, label, source_id,
                         cid, tis_pos, strand, sc, win))
        tis_0idx_by_contig[cid].add(tis_pos - 1)

    # ── Tier-1 intergenic decoys (sanity check) ───────────────────────────────
    rows_dec = []
    contigs_with_spans = [c for c in contigs if c in spans] or list(contigs)
    target   = len(rows_pos)
    attempts = 0
    while len(rows_dec) < target and attempts < target * 100:
        attempts += 1
        cid = random.choice(contigs_with_spans)
        seq = contigs[cid]; L = len(seq)
        if L < WIN: continue
        pos = random.randint(HALF, L - HALF - 1)   # 0-indexed window centre
        if any(s <= pos < e for s, e in spans.get(cid, [])): continue
        if any(abs(pos - t) < MASK for t in tis_0idx_by_contig[cid]): continue
        strand = random.choice("+-")
        win = extract_window(seq, pos + 1, strand)
        if win is None or len(win) != WIN or "N" in win: continue
        wid = f"{org}_rbs_decoy_intergenic_{len(rows_dec)}"
        rows_dec.append((wid, org, phylum, "decoy_intergenic", "DECOY",
                         cid, pos + 1, strand, win[HALF:HALF + 3], win))

    for row in rows_pos + rows_dec:
        writer.writerow(row)

    n_sd   = sum(1 for r in rows_pos if r[3] == "positive_SD")
    n_unsd = len(rows_pos) - n_sd
    print(f"  {org:22s}  loaded={n_loaded:5d}  unique={len(tis_list):5d}  "
          f"windows={len(rows_pos):5d} (SD={n_sd} UNSD={n_unsd})  "
          f"miss={n_miss}  lift_miss={liftover_miss}  "
          f"Tier1_decoys={len(rows_dec)}")
    return len(rows_pos), len(rows_dec), n_sd, n_unsd


def main():
    os.makedirs(OUT, exist_ok=True)
    out_path = os.path.join(OUT, "ALL.tsv")
    fields   = ["id", "organism", "phylum", "label", "source_id",
                "contig", "tis_pos", "strand", "start_codon", "window_seq"]

    total_pos = total_dec = total_sd = total_unsd = 0
    print("Building RBS windows ...")
    with open(out_path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(fields)
        for org, (gen_dir, phylum, sources) in ORGS.items():
            np, nd, nsd, nunsd = build_org(org, gen_dir, phylum, sources, w)
            total_pos += np; total_dec += nd
            total_sd  += nsd; total_unsd += nunsd

    print(f"\nTOTAL  positives={total_pos:,} (SD={total_sd:,}  UNSD={total_unsd:,})"
          f"  Tier1_decoys={total_dec:,}  -> {out_path}")


if __name__ == "__main__":
    main()
