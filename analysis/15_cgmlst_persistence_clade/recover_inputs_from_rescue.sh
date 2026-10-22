#!/usr/bin/env bash
set -euo pipefail

# Recover only the files needed by the cgMLST persistence/clade analysis
# from the verified manuscript-rescue GitHub Release.
#
# This avoids extracting the entire rescued project, but the 18 split
# archive parts (~33 GB total) still have to be downloaded once.

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
INPUT_DIR="$HERE/inputs"
WORK_DIR="${RESCUE_WORK_DIR:-$HOME/Ecoli_epi_rescue_minimal}"
ASSET_DIR="$WORK_DIR/release_assets"
EXTRACT_DIR="$WORK_DIR/extracted"

REPO="ajulojays/Ecoli_epi_manuscript_backup"
TAG="manuscript-rescue-2026-10-09"
PREFIX="Ecoli_epi_MANUSCRIPT_RESCUE_20261009.tar.zst.part-"

mkdir -p "$INPUT_DIR" "$ASSET_DIR" "$EXTRACT_DIR"

for cmd in gh zstd tar sha256sum; do
    command -v "$cmd" >/dev/null 2>&1 || {
        echo "[ERROR] Missing required command: $cmd" >&2
        echo "Ubuntu/Debian: sudo apt update && sudo apt install -y zstd" >&2
        echo "Install/authenticate GitHub CLI separately if gh is missing." >&2
        exit 2
    }
done

echo "============================================================"
echo " Minimal recovery from manuscript-rescue release"
echo " Repo: $REPO"
echo " Tag:  $TAG"
echo " Work: $WORK_DIR"
echo "============================================================"

echo "[1/5] Confirming GitHub authentication/private-repo access..."
gh repo view "$REPO" >/dev/null

echo "[2/5] Downloading rescue checksum manifests and split archive parts..."
gh release download "$TAG" \
  --repo "$REPO" \
  --dir "$ASSET_DIR" \
  --pattern 'SHA256SUMS_parts.txt' \
  --pattern 'Ecoli_epi_MANUSCRIPT_RESCUE_20261009.tar.zst.part-*' \
  --clobber

cd "$ASSET_DIR"

n_parts=$(find . -maxdepth 1 -type f -name "${PREFIX}*" | wc -l)
echo "[INFO] Archive parts found: $n_parts"
if [[ "$n_parts" -ne 18 ]]; then
    echo "[ERROR] Expected 18 archive parts, found $n_parts." >&2
    exit 3
fi

echo "[3/5] Verifying all split-part SHA256 checksums..."
sha256sum -c SHA256SUMS_parts.txt

echo "[4/5] Streaming the verified archive and extracting only needed files..."
rm -rf "$EXTRACT_DIR"
mkdir -p "$EXTRACT_DIR"

# --no-anchored lets basename-like patterns match regardless of the rescue root.
# All four files were part of the original project; the first two are required
# for computation, while clade/tree files improve reproducibility/visualization.
cat ${PREFIX}* | zstd -dc | tar -x -C "$EXTRACT_DIR" \
  --wildcards --no-anchored \
  '*target_metadata_704_preharvest_non_dairy.csv' \
  '*cgMLST95.tsv' \
  '*all_preharvest_704_clades.csv' \
  '*all_preharvest_704_core_phylo.treefile'

echo "[5/5] Copying recovered files into module inputs/..."

find_one () {
    local name="$1"
    find "$EXTRACT_DIR" -type f -name "$name" -print -quit
}

META="$(find_one target_metadata_704_preharvest_non_dairy.csv || true)"
CG="$(find_one cgMLST95.tsv || true)"
CLADES="$(find_one all_preharvest_704_clades.csv || true)"
TREE="$(find_one all_preharvest_704_core_phylo.treefile || true)"

[[ -n "$META" ]] || { echo "[ERROR] target metadata not found after extraction." >&2; exit 4; }
[[ -n "$CG" ]] || { echo "[ERROR] cgMLST95.tsv not found after extraction." >&2; exit 4; }

cp -f "$META" "$INPUT_DIR/metadata_704.csv"
cp -f "$CG" "$INPUT_DIR/cgMLST95.tsv"

if [[ -n "$CLADES" ]]; then
    cp -f "$CLADES" "$INPUT_DIR/all_preharvest_704_clades.csv"
fi
if [[ -n "$TREE" ]]; then
    cp -f "$TREE" "$INPUT_DIR/all_preharvest_704_core_phylo.treefile"
fi

echo
echo "[SUCCESS] Recovered analysis inputs:"
ls -lh "$INPUT_DIR"/metadata_704.csv "$INPUT_DIR"/cgMLST95.tsv
[[ -f "$INPUT_DIR/all_preharvest_704_clades.csv" ]] && ls -lh "$INPUT_DIR/all_preharvest_704_clades.csv"
[[ -f "$INPUT_DIR/all_preharvest_704_core_phylo.treefile" ]] && ls -lh "$INPUT_DIR/all_preharvest_704_core_phylo.treefile"

echo
echo "You may now run:"
echo "  bash run.sh"
echo
echo "To reclaim disk space AFTER successful recovery, remove:"
echo "  $WORK_DIR"
