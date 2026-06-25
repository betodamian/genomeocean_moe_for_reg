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
    ("ECOLI_RIBORET", "rbs", "pmid:30904393", "https://www.cell.com/molecular-cell/fulltext/S1097-2765(19)30119-3 (Meydan 2019; Mol Cell SI tables S1-S6)"),
    ("ECOLI_NOMRNA",  "rbs", "pmid:32424076", "https://elifesciences.org/articles/55498/figures (Saito 2020; eLife supp files)"),
    ("ECOLI_TETRP",   "rbs", "pmid:27013550", "https://academic.oup.com/dnaresearch/article/23/3/193/2469929 (Nakahigashi 2016; Supp Table 2-5)"),
    ("MTB_RIBOSEQ",   "rbs", "pmid:33513356", "https://www.cell.com/cell-reports/fulltext/S2211-1247(20)31641-9 (Zhu 2021; Cell Rep Table S1-S2)"),
    ("BSUB_RIBOSEQ",  "rbs", "pmid:29144454", "https://www.cell.com/cell-reports/fulltext/S2211-1247(17)31523-3 (Lalanne 2017; Cell Rep supp)"),
    ("BSUB_SPORE",    "rbs", "pmid:39179838", "https://www.nature.com/articles/s41467-024-51654-6#Sec (Bhatt 2024; Supplementary Data)"),
    ("CAULO_RIBOSEQ", "rbs", "doi:10.1371/journal.pgen.1004463", "https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1004463 (Schrader 2014; S18/S23 datasets)"),
    ("SAUR_EXSD",     "rbs", "pmid:41680142", "https://www.nature.com/articles/s41467-026-69079-8#Sec (Kohl 2026; Supplementary Data)"),
    ("MTB_RHODUC",    "rho", "pmid:37096044", "https://www.cell.com/molecular-cell/fulltext/S1097-2765(23)00224-6 (Botella 2022; Table S = 1,385 RD-TTS)"),
    ("BSUB_HSELEX",   "rho", "pmcid:PMC12350095", "https://pmc.ncbi.nlm.nih.gov/articles/PMC12350095/ (B. subtilis H-SELEX rut sites)"),
    ("INTRINSIC_TERMITE", "rho", "pmcid:PMC12207403", "https://academic.oup.com/nar/article/53/12/gkaf553 (TERMITe intrinsic atlas)"),
    ("RHOTERMPREDICT","rho", "pmcid:PMC6407284", "https://pmc.ncbi.nlm.nih.gov/articles/PMC6407284/ (RhoTermPredict; cross-check only)"),
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
