#!/usr/bin/env bash
# Copy plot exports from the separate research repository while retaining existing figures.
set -euo pipefail

# 1. Define the source, destination and temporary workspace ----
SOURCE="https://github.com/Nadaline/Beyond-Generalism-staging.git"
DEST="https://github.com/Nadaline/Beyond_generalism_a_climate_genomic_mosaic_on_neotropics.git"
TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT

# 2. Clone both repositories ----
git clone --depth 1 --branch release/ameiva-official-2026 "$SOURCE" "$TMP/source"
git clone --depth 1 "$DEST" "$TMP/public"

# 3. Copy figures absent from the public repository ----
mkdir -p "$TMP/public/figures"
rsync -a --ignore-existing "$TMP/source/figures/" "$TMP/public/figures/"

# 4. Commit and publish newly copied figures ----
cd "$TMP/public"
git add figures/
if git diff --cached --quiet; then
  echo "The earlier plot exports are already present."
  exit 0
fi

git commit -m "Archive earlier climate-association figure exports"
git push origin main
