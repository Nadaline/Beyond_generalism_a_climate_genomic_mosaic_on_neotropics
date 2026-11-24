#!/usr/bin/env bash
# Copy earlier plot exports from the private research repository into the
# public archive without replacing the six figures selected for the R7 thesis.
set -euo pipefail

SOURCE="https://github.com/Nadaline/Beyond-Generalism-staging.git"
DEST="https://github.com/Nadaline/Beyond_generalism_a_climate_genomic_mosaic_on_neotropics.git"
TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT

git clone --depth 1 --branch release/ameiva-official-2026 "$SOURCE" "$TMP/source"
git clone --depth 1 "$DEST" "$TMP/public"

mkdir -p "$TMP/public/figures"
rsync -a --ignore-existing "$TMP/source/figures/" "$TMP/public/figures/"

cd "$TMP/public"
git add figures/
if git diff --cached --quiet; then
  echo "The earlier plot exports are already present."
  exit 0
fi

git commit -m "Archive earlier climate-association figure exports"
git push origin main
