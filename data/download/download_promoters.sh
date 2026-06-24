#!/usr/bin/env bash
# Download promoter LABEL data per research_plan.md §5d (anti-circularity §12).
# Primary source: PPD — Prokaryotic Promoter Database (experimentally verified).
#   Liu et al. 2021, J Mol Biol 433:166860; http://lin-group.cn/database/ppd/
#   129,148 experimental promoters across 63 species / 74 strains.
# RegulonDB v12 is a SPA (no stable bulk URL) — fetched manually; see MANIFEST.md.
set -euo pipefail
cd "$(dirname "$0")/../.."   # repo root
OUT=data/promoters
mkdir -p "$OUT"

echo ">>> PPD (Prokaryotic Promoter Database) — experimental promoters"
curl -s -m 180 -L -o "$OUT/PPD_csv.zip" "http://lin-group.cn/database/ppd/todownload/csv.zip"
python -c "import zipfile; zipfile.ZipFile('$OUT/PPD_csv.zip').extractall('$OUT/PPD')"
n=$(for f in "$OUT"/PPD/csv/*.csv; do echo $(( $(wc -l < "$f") - 1 )); done | paste -sd+ - | bc 2>/dev/null \
    || awk 'END{print "(count: run wc)"}')
echo "    extracted $(ls "$OUT"/PPD/csv/*.csv | wc -l) species CSVs; ~$n promoters total"
echo ""
echo "NOTE (per §12 anti-circularity): before use, filter PPD entries to experimental"
echo "evidence and exclude any traceable to a benchmarked tool (BPROM/bTSSfinder/iPromoter)."
echo "RegulonDB v12 (E. coli sigma-subtyping cross-check): fetch release from"
echo "https://regulondb.ccg.unam.mx (Downloads tab; SPA — not scriptable). See MANIFEST.md."
