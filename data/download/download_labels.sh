#!/usr/bin/env bash
# Download regulatory-element LABEL data (processed supplementary tables) per research_plan.md.
# Scope: directly-accessible processed tables from GEO/ArrayExpress supplementary files.
#   - Skips raw reads (.sra/.bam/.bigwig/.bw) and files > MAXSIZE (raw data needs the
#     Week-1 TIS/peak-calling pipeline, not a blind bulk pull).
#   - Sources, accessions, and citations: see data/MANIFEST.md.
# Idempotent; logs to stdout. Items requiring manual/journal fetch are listed at the end.
set -uo pipefail
cd "$(dirname "$0")/../.."   # repo root
MAXSIZE=$((60*1024*1024))    # 60 MB per-file guard

geo_suppl_dir () {            # GSE123456 -> https://ftp.ncbi.nlm.nih.gov/geo/series/GSE123nnn/GSE123456/suppl/
  local id="$1" num prefix
  num="${id#GSE}"; prefix="${num%???}"
  echo "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE${prefix}nnn/${id}/suppl/"
}

pull_geo () {                # $1=GSE id  $2=dest dir  $3=description
  local id="$1" dest="$2" desc="$3" dir files f url
  echo ">>> $id  ($desc)"
  mkdir -p "$dest"
  dir="$(geo_suppl_dir "$id")"
  files="$(curl -s -m 40 "$dir" | grep -oiE 'href="[^"]+"' | sed -E 's/href="//;s/"//' \
           | grep -viE '^/|^https?:' | grep -iE '\.(txt|tsv|csv|gff|gff3|bed|xlsx|gz)$' \
           | grep -viE '\.(bam|sra|bw|bigwig|wig)\.gz$' || true)"
  if [ -z "$files" ]; then echo "    (no small supplementary tables listed; see MANIFEST)"; return; fi
  while read -r f; do
    [ -z "$f" ] && continue
    url="${dir}${f}"
    sz="$(curl -s -m 30 -I -L "$url" | tr -d '\r' | awk -F': ' 'tolower($1)=="content-length"{print $2}' | tail -1)"
    if [ -n "${sz:-}" ] && [ "$sz" -gt "$MAXSIZE" ]; then
      echo "    SKIP (>60MB, raw — fetch in Week-1 pipeline): $f"; continue
    fi
    echo "    get: $f"
    curl -s -m 300 -L -o "$dest/$f" "$url"
  done <<< "$files"
}

# ---------- RBS / TIS sources (research_plan.md §5d, data/rbs_database/sources.tsv) ----------
pull_geo GSE122129 data/rbs_database/raw "E. coli Ribo-RET TIS (Meydan 2019)"
pull_geo GSE135906 data/rbs_database/raw "E. coli ribosome profiling, daSD (Saito 2020)"
pull_geo GSE95211  data/rbs_database/raw "B. subtilis Ribo-seq (Lalanne 2017)"
pull_geo GSE249450 data/rbs_database/raw "B. subtilis sporulation Ribo-seq (Bhatt 2024)"
pull_geo GSE54883  data/rbs_database/raw "Caulobacter crescentus Ribo-seq+5'-RACE (Schrader 2014)"

# ---------- Rho sources (research_plan.md §5d, data/rho_database/sources.tsv) ----------
pull_geo GSE109766 data/rho_database/raw "E. coli Rho Term-seq 3'-ends (NAR 2018)"

# ---------- ArrayExpress (MTB) — listed via BioStudies API ----------
for acc in E-MTAB-8835 E-MTAB-11753; do
  echo ">>> $acc  (ArrayExpress — MTB; $( [ "$acc" = E-MTAB-8835 ] && echo 'Ribo-seq Zhu 2021' || echo 'Rho RhoDUC Botella 2022' ))"
  curl -s -m 40 "https://www.ebi.ac.uk/biostudies/api/v1/studies/$acc" \
    | python -c "import sys,json
try: d=json.load(sys.stdin)
except Exception: print('    (API listing unavailable — manual fetch, see MANIFEST)'); sys.exit()
def walk(o):
    if isinstance(o,dict):
        if 'path' in o: print('    file:',o.get('path'),'size',o.get('size'))
        for v in o.values(): walk(v)
    elif isinstance(o,list):
        for v in o: walk(v)
walk(d)" 2>&1 | head -20
done

echo ""
echo "================ MANUAL / WEEK-1-PIPELINE FETCH (not auto-pullable) ================"
echo " RegulonDB v12 (E. coli promoters/TUs/terminators, sigma-subtyping): SPA site;"
echo "   fetch release files from https://regulondb.ccg.unam.mx (Downloads) — see MANIFEST."
echo " E. coli TetRP TIS (Nakahigashi 2016): DDBJ BioProject PRJDB2960 (raw reads)."
echo " S. aureus extended-SD TIS (Kohl 2026): GEO/SRA in paper data-availability (verify)."
echo " M. tuberculosis dRNA-seq leaderless (Cortes 2013): SRA SRP028740."
echo " TERMITe intrinsic terminator atlas: Zenodo + GitHub (per PMC12207403)."
echo " Peters 2012 E. coli Rho (BCM): PNAS supplementary tables."
echo " B. subtilis H-SELEX rut sites: PMC12350095 supplementary."
echo " Journal supplementary tables (per-gene TIS class) for Meydan/Zhu/Botella:"
echo "   download from each paper's SI (DOIs in data/*_database/sources.tsv)."
echo "===================================================================================="
