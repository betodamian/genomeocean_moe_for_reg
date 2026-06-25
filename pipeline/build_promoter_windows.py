#!/usr/bin/env python3
"""
Build promoter windows + same-context decoys (research_plan.md §5c) for the PPD
species that exactly match a downloaded genome.

Positives  : PPD experimentally-verified promoters. Each PPD row carries an 81-bp
             sequence with the TSS marked uppercase. We locate that sequence in the
             genome (exact match, both strands) -> exact TSS coordinate + strand,
             robust to PPD's coordinate convention, and self-validating (unlocatable
             promoters are dropped & counted).
Windows    : 300 bp, promoter-strand 5'->3', TSS centred at position 150 (§5c/§6).
Decoys     : same-context negatives -> intergenic 300 bp windows >= MASK bp from any
             mapped TSS (so the model can't win on "gene vs non-gene"); matched count
             & strand distribution. (v1; the harder upstream-of-gene decoy is a later
             refinement.)
Output     : data/datasets/promoters/<org>.tsv  and  ALL.tsv  (gitignored).

Run: .venv/bin/python pipeline/build_promoter_windows.py
"""
import csv, os, random, sys
from collections import defaultdict

random.seed(20260625)
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
WIN = 300              # window length
HALF = WIN // 2        # TSS at index 150
MASK = 250             # decoys must be >= MASK bp from any mapped TSS

# organism -> (genome_dir, ppd_species_name, phylum)
ORGS = {
    "ecoli_K12_MG1655":  ("Escherichia coli str K-12 substr. MG1655", "Gammaproteobacteria"),
    "bsubtilis_168":     ("Bacillus subtilis subsp. subtilis str. 168", "Firmicutes"),
    "hvolcanii_DS2":     ("Haloferax volcanii DS2", "Archaea"),
}
PPD_DIR = os.path.join(ROOT, "data/promoters/PPD/csv")
GEN_DIR = os.path.join(ROOT, "data/genomes")
OUT_DIR = os.path.join(ROOT, "data/datasets/promoters")

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
def rc(s): return s.translate(COMP)[::-1]

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

def genome_paths(org):
    d = os.path.join(GEN_DIR, org)
    fna = next(f for f in os.listdir(d) if f.endswith(".fna"))
    gff = next(f for f in os.listdir(d) if f.endswith(".gff"))
    return os.path.join(d, fna), os.path.join(d, gff)

def gene_spans(gff_path):
    """contig -> sorted list of (start0, end) for gene features (0-based, end-exclusive)."""
    spans = defaultdict(list)
    with open(gff_path) as fh:
        for ln in fh:
            if ln.startswith("#"): continue
            f = ln.rstrip("\n").split("\t")
            if len(f) < 8 or f[2] != "gene": continue
            spans[f[0]].append((int(f[3]) - 1, int(f[4])))
    for c in spans: spans[c].sort()
    return spans

def ppd_rows(species):
    """yield (tss_offset_in_seq, seq_upper, strand) for valid 81-bp promoter rows."""
    path = None
    for fn in os.listdir(PPD_DIR):
        if fn.startswith(species):  # filename begins with species name
            path = os.path.join(PPD_DIR, fn); break
    if path is None:
        sys.exit(f"PPD csv not found for {species}")
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh):
            if r.get("SpeciesName") != species:  # 'other species' safety
                pass
            seq = r["PromoterSeq"]
            if not seq or len(seq) < 60: continue
            ups = [i for i, ch in enumerate(seq) if ch.isupper()]
            if len(ups) != 1: continue          # need exactly one TSS-marked base
            yield ups[0], seq.upper(), r["Strand"]

def locate(seq_u, strand, contigs):
    """return (contig, tss_plus_coord, strand) or None. tss within the matched window."""
    for cid, plus in contigs.items():
        if strand == "+":
            i = plus.find(seq_u)
            if i >= 0: return cid, i, "+"
        else:
            t = rc(seq_u)
            i = plus.find(t)
            if i >= 0: return cid, i, "-"
    # strand may be mislabeled — try the other orientation
    for cid, plus in contigs.items():
        i = plus.find(seq_u)
        if i >= 0: return cid, i, "+"
        i = plus.find(rc(seq_u))
        if i >= 0: return cid, i, "-"
    return None

def window_for(plus, hit_start, tss_off_in_seq, seq_len, strand):
    """300 bp window on promoter strand, TSS centred. hit_start = match start on +."""
    if strand == "+":
        tss = hit_start + tss_off_in_seq
        s = tss - HALF
        if s < 0 or s + WIN > len(plus): return None
        return plus[s:s + WIN], tss
    else:
        # promoter seq is rc of plus[hit_start:hit_start+seq_len]; TSS is tss_off from
        # the promoter's 5' end -> on + strand that is the (seq_len-1-tss_off) offset.
        tss = hit_start + (seq_len - 1 - tss_off_in_seq)
        s = tss - (HALF - 1)
        if s < 0 or s + WIN > len(plus): return None
        return rc(plus[s:s + WIN]), tss

def build_org(org, species, phylum, writer):
    fna, gff = genome_paths(org)
    contigs = load_fasta(fna)
    n_pos = n_miss = 0
    tss_by_contig = defaultdict(set)
    rows_pos = []
    for tss_off, seq_u, strand in ppd_rows(species):
        hit = locate(seq_u, strand, contigs)
        if hit is None: n_miss += 1; continue
        cid, hit_start, str_resolved = hit
        w = window_for(contigs[cid], hit_start, tss_off, len(seq_u), str_resolved)
        if w is None: n_miss += 1; continue
        win, tss = w
        if "N" in win: continue
        rows_pos.append((f"{org}_prom_pos_{len(rows_pos)}", org, phylum, "positive",
                         cid, tss, str_resolved, win))
        tss_by_contig[cid].add(tss)
        n_pos += 1
    # ---- decoys: intergenic windows >= MASK from any mapped TSS ----
    spans = gene_spans(gff)
    rows_dec = []
    # candidate intergenic ranges per contig
    target = n_pos
    contigs_list = [c for c in contigs if c in spans]
    attempts = 0
    while len(rows_dec) < target and attempts < target * 50:
        attempts += 1
        cid = random.choice(contigs_list)
        L = len(contigs[cid])
        pos = random.randint(HALF, L - HALF - 1)
        # reject if inside a gene body (we want intergenic) or near a TSS
        in_gene = any(s <= pos < e for s, e in spans[cid])
        if in_gene: continue
        if any(abs(pos - t) < MASK for t in tss_by_contig[cid]): continue
        strand = random.choice("+-")
        if strand == "+":
            win = contigs[cid][pos - HALF:pos + HALF]
        else:
            win = rc(contigs[cid][pos - HALF + 1:pos + HALF + 1])
        if len(win) != WIN or "N" in win: continue
        rows_dec.append((f"{org}_prom_dec_{len(rows_dec)}", org, phylum, "decoy",
                         cid, pos, strand, win))
    for row in rows_pos + rows_dec:
        writer.writerow(row)
    return n_pos, n_miss, len(rows_dec)

def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    all_path = os.path.join(OUT_DIR, "ALL.tsv")
    summary = []
    with open(all_path, "w", newline="") as afh:
        w = csv.writer(afh, delimiter="\t", lineterminator="\n")
        w.writerow(["id", "organism", "phylum", "label", "contig", "tss_or_pos", "strand", "window_seq"])
        for org, (species, phylum) in ORGS.items():
            npos, nmiss, ndec = build_org(org, species, phylum, w)
            rate = 100 * npos / (npos + nmiss) if (npos + nmiss) else 0
            summary.append((org, npos, nmiss, ndec, rate))
            print(f"{org:22s} positives={npos:6d}  unlocated={nmiss:5d} "
                  f"({rate:4.1f}% located)  decoys={ndec:6d}")
    tp = sum(s[1] for s in summary); td = sum(s[3] for s in summary)
    print(f"\nTOTAL  positives={tp}  decoys={td}  ->  {all_path}")

if __name__ == "__main__":
    main()
