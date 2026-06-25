#!/usr/bin/env python3
"""
Auto-fetch the PROCESSED RBS/Rho label tables (journal supplementary files) per
data/{rbs,rho}_database/sources.tsv.

Route (https, no anti-bot wall): PMID/DOI --(PMC idconv)--> PMCID
  --> EuropePMC supplementaryFiles endpoint, which hosts the OA supplementary files
      itself (the PMC /bin/ links redirect to publisher CDNs that serve a CAPTCHA;
      EuropePMC does not). Returns a zip; nested zips are unpacked one level.

Open-access papers => real tables pulled. Papers not in the EuropePMC OA subset
(Elsevier author-manuscripts, not-in-PMC, RegulonDB SPA) => NEED-MANUAL + a link.

Downloads -> data/{rbs,rho}_database/raw/<source_id>/   (gitignored)
Report    -> printed + data/label_fetch_report.md
Run: .venv/bin/python pipeline/fetch_labels.py
"""
import io, json, os, re, subprocess, time, urllib.parse, zipfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RBS = os.path.join(ROOT, "data/rbs_database/raw")
RHO = os.path.join(ROOT, "data/rho_database/raw")
KEEP = (".xlsx", ".xls", ".csv", ".tsv", ".txt", ".gff", ".gff3", ".bed", ".gz")
UA = "Mozilla/5.0 (academic data fetch; research@example.org)"
EPMC = "https://www.ebi.ac.uk/europepmc/webservices/rest/{}/supplementaryFiles"

# (source_id, element, identifier, manual_link)
SOURCES = [
    # ECOLI_RIBORET: Meydan 2019 Mol Cell (PMID 30904393, PMC7115971) — NOT OA in EuropePMC;
    #   download Table S1–S6 from https://doi.org/10.1016/j.molcel.2019.02.017 (institutional access)
    #   or PMC author-manuscript https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7115971/
    ("ECOLI_RIBORET", "rbs", "pmid:30904393",
     "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7115971/ "
     "(Meydan 2019 Mol Cell; Table S1-S6 at PMC or via institutional access DOI 10.1016/j.molcel.2019.02.017)"),
    # ECOLI_NOMRNA: Saito 2020 eLife 55002 (PMCID PMC7043885) — OA but EuropePMC only has figures.
    #   Raw data at GEO GSE135906 (WIG / FASTQ). For processed TIS table: download from
    #   https://elifesciences.org/articles/55002 (Source data links within figures).
    ("ECOLI_NOMRNA",  "rbs", "pmcid:PMC7043885",
     "https://elifesciences.org/articles/55002 "
     "(Saito 2020 eLife; data in GEO GSE135906 — WIG files; or source data attached to figures on eLife)"),
    ("ECOLI_TETRP",   "rbs", "pmid:27013550",
     "https://academic.oup.com/dnaresearch/article/23/3/193/2469929 (Nakahigashi 2016; Supp Table 2-5)"),
    # MTB_RIBOSEQ: Sawyer 2021 Cell Rep (PMID 33535039, PMCID PMC7856553) — OA; fetched automatically.
    ("MTB_RIBOSEQ",   "rbs", "pmcid:PMC7856553",
     "https://doi.org/10.1016/j.celrep.2021.108695 "
     "(Sawyer et al. 2021 Cell Reports; Table S = mmc2/mmc3/mmc4.xlsx)"),
    # BSUB_RIBOSEQ: Lalanne 2018 Cell (PMID 29606352, PMCID PMC5978003) — NOT OA in EuropePMC;
    #   download supplementary tables from https://doi.org/10.1016/j.cell.2018.03.007 (institutional access).
    ("BSUB_RIBOSEQ",  "rbs", "pmid:29606352",
     "https://doi.org/10.1016/j.cell.2018.03.007 "
     "(Lalanne et al. 2018 Cell; B. subtilis TIS data in supplementary tables; institutional access required)"),
    # BSUB_SPORE: Iwanska 2024 Nat Commun (PMID 39169056, PMCID PMC11339384) — OA; fetched automatically.
    ("BSUB_SPORE",    "rbs", "pmcid:PMC11339384",
     "https://www.nature.com/articles/s41467-024-51654-6 "
     "(Iwanska et al. 2024 Nat Commun; B. subtilis sporulation translation; MOESM xlsx files)"),
    ("CAULO_RIBOSEQ", "rbs", "doi:10.1371/journal.pgen.1004463",
     "https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1004463 (Schrader 2014; S18/S23 datasets)"),
    ("SAUR_EXSD",     "rbs", "pmid:41680142",
     "https://www.nature.com/articles/s41467-026-69079-8 (Kohl 2026; Supplementary Data)"),
    # HVOLC_RIBOSEQ: Gelsinger 2020 NAR (PMID 32382758, PMCID PMC7261190) — OA; fetched automatically.
    #   H. volcanii DS2 ribosome profiling; 1,413 annotated TIS + novel/extension TIS.
    ("HVOLC_RIBOSEQ", "rbs", "pmcid:PMC7261190",
     "https://academic.oup.com/nar/article/48/10/5201/5831329 "
     "(Gelsinger 2020 NAR; Ribo_MS_TableS1_final.xlsx in nested supp zip)"),
    # MTB_RHODUC: Botella 2017 Nat Commun (PMID 28348398, PMCID PMC5379054) — OA; fetched automatically.
    #   Genome-wide Rho-dependent termination sites in M. tuberculosis (xlsx S2–S9 + s10.txt).
    ("MTB_RHODUC",    "rho", "pmcid:PMC5379054",
     "https://www.nature.com/articles/ncomms14731 "
     "(Botella 2017 Nat Commun; MTB Rho termination sites in supplementary xlsx S2-S9)"),
    ("BSUB_HSELEX",   "rho", "pmcid:PMC12350095",
     "https://pmc.ncbi.nlm.nih.gov/articles/PMC12350095/ (B. subtilis H-SELEX rut sites)"),
    ("INTRINSIC_TERMITE", "rho", "pmcid:PMC12207403",
     "https://academic.oup.com/nar/article/53/12/gkaf553 (TERMITe intrinsic atlas)"),
    ("RHOTERMPREDICT","rho", "pmcid:PMC6407284",
     "https://pmc.ncbi.nlm.nih.gov/articles/PMC6407284/ (RhoTermPredict; cross-check only)"),
]

def curl(url, out=None, timeout=120):
    cmd = ["curl", "-s", "-L", "-m", str(timeout), "-A", UA]
    cmd += (["-o", out] if out else [])
    cmd.append(url)
    r = subprocess.run(cmd, capture_output=(out is None))
    return r.stdout.decode("utf-8", "replace") if out is None else (r.returncode == 0)

def idconv(ident):
    kind, val = ident.split(":", 1)
    if kind == "pmcid":
        return val
    q = urllib.parse.quote(val)
    url = (f"https://pmc.ncbi.nlm.nih.gov/tools/idconv/api/v1/articles/"
           f"?ids={q}&format=json&tool=goreg&email=research@example.org")
    try:
        d = json.loads(curl(url, timeout=40) or "{}")
        r = (d.get("records") or [{}])[0]
        return r.get("pmcid") if r.get("status") != "error" else None
    except Exception:
        return None

def extract_zip(zbytes, dest, got, depth=0):
    try:
        zf = zipfile.ZipFile(io.BytesIO(zbytes))
    except zipfile.BadZipFile:
        return
    for name in zf.namelist():
        base = os.path.basename(name)
        if not base:
            continue
        data = zf.read(name)
        if base.lower().endswith(".zip") and depth < 3:
            extract_zip(data, dest, got, depth + 1)
        elif base.lower().endswith(KEEP):
            with open(os.path.join(dest, base), "wb") as fh:
                fh.write(data)
            got.append((base, len(data)))

def fetch(sid, element, ident, manual):
    dest = os.path.join(RBS if element == "rbs" else RHO, sid)
    pmcid = idconv(ident)
    if not pmcid:
        return ("NOT-IN-PMC", None, [], manual)
    raw = subprocess.run(["curl", "-s", "-L", "-m", "300", "-A", UA, EPMC.format(pmcid)],
                         capture_output=True).stdout
    if len(raw) < 1000 or raw[:2] != b"PK":       # empty XML or not a zip => no OA supp
        return ("NO-OA-SUPP", pmcid, [], manual)
    os.makedirs(dest, exist_ok=True)
    got = []
    extract_zip(raw, dest, got)
    return ("FETCHED" if got else "ZIP-NO-TABLES", pmcid, got, manual)

def main():
    for d in (RBS, RHO):  # clear prior bogus attempts (keep root-level GEO signal files)
        for sub in os.listdir(d) if os.path.isdir(d) else []:
            p = os.path.join(d, sub)
            if os.path.isdir(p):
                for f in os.listdir(p): os.remove(os.path.join(p, f))
    rep = ["# RBS/Rho label auto-fetch report\n\n",
           "Route: EuropePMC OA supplementary-files endpoint. "
           "FETCHED = real supplementary table(s) pulled.\n\n"]
    fetched, manual = [], []
    for sid, element, ident, link in SOURCES:
        status, pmcid, got, mlink = fetch(sid, element, ident, link)
        time.sleep(1.0)
        print(f"[{element}] {sid:16s} {status:13s} pmc={pmcid}")
        for fn, sz in got: print(f"        -> {sz:>11,} {fn}")
        if status == "FETCHED":
            fetched.append(sid)
            rep.append(f"- **{sid}** ({element}) — FETCHED from {pmcid}: "
                       + ", ".join(f"`{fn}` ({sz:,}B)" for fn, sz in got) + "\n")
        else:
            manual.append((sid, element, status, mlink))
    rep.append("\n## Need manual fetch (open the link, save the supplementary table into "
               "`data/<el>_database/raw/<source_id>/`)\n\n")
    for sid, element, status, link in manual:
        rep.append(f"- **{sid}** ({element}) — _{status}_ — {link}\n")
    with open(os.path.join(ROOT, "data/label_fetch_report.md"), "w") as fh:
        fh.write("".join(rep))
    print(f"\nFETCHED {len(fetched)} | NEED-MANUAL {len(manual)}  ->  data/label_fetch_report.md")

if __name__ == "__main__":
    main()
