#!/usr/bin/env bash
# Download all genome sequences + annotations for the GenomeOcean-MoE regulatory project.
# Source: NCBI Datasets API v2 (https://www.ncbi.nlm.nih.gov/datasets/).
# Genomes pinned by RefSeq assembly accession (research_plan.md §5d, §5f).
# Each genome: GENOME_FASTA (.fna) + GENOME_GFF (.gff) extracted to data/genomes/<label>/.
# Idempotent: re-running re-downloads and overwrites. Records SHA256 in checksums.txt.
set -euo pipefail

cd "$(dirname "$0")/../.."   # repo root
OUT=data/genomes
TMP=data/download/_tmp
mkdir -p "$OUT" "$TMP"
CHECK=data/download/checksums.txt
: > "$CHECK"

# accession <tab> short_label  (role per research_plan.md)
read -r -d '' GENOMES <<'EOF' || true
GCF_000005845.2	ecoli_K12_MG1655
GCF_000009045.1	bsubtilis_168
GCF_000195955.2	mtuberculosis_H37Rv
GCF_000006765.1	paeruginosa_PAO1
GCF_000013425.1	saureus_NCTC8325
GCF_000022005.1	ccrescentus_NA1000
GCF_000025685.1	hvolcanii_DS2
EOF

API="https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession"

while IFS=$'\t' read -r ACC LABEL; do
  [ -z "${ACC:-}" ] && continue
  echo ">>> $ACC ($LABEL)"
  ZIP="$TMP/$ACC.zip"
  curl -s -L -o "$ZIP" \
    "$API/$ACC/download?include_annotation_type=GENOME_FASTA&include_annotation_type=GENOME_GFF"
  DEST="$OUT/$LABEL"
  mkdir -p "$DEST"
  python - "$ZIP" "$DEST" "$ACC" <<'PY'
import sys, zipfile, os, shutil
zippath, dest, acc = sys.argv[1], sys.argv[2], sys.argv[3]
z = zipfile.ZipFile(zippath)
for name in z.namelist():
    if name.endswith(".fna") or name.endswith(".gff"):
        data = z.read(name)
        base = os.path.basename(name)
        if base == "genomic.gff":
            base = f"{acc}_genomic.gff"
        with open(os.path.join(dest, base), "wb") as fh:
            fh.write(data)
        print("   extracted", base, f"({len(data):,} bytes)")
PY
  rm -f "$ZIP"
done <<< "$GENOMES"

# checksums
echo ">>> writing SHA256 checksums"
( cd "$OUT" && find . -type f \( -name '*.fna' -o -name '*.gff' \) -print0 \
  | sort -z | xargs -0 sha256sum ) > "$CHECK"
rm -rf "$TMP"
echo "DONE. Files in $OUT ; checksums in $CHECK"
