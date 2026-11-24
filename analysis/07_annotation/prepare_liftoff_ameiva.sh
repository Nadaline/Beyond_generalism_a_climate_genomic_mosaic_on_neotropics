#!/usr/bin/env bash
set -euo pipefail

# Prepara a transferencia da anotacao entre as versoes do assembly ao qual o
# VCF foi realmente alinhado. Os arquivos grandes permanecem no HD externo.

BASE="/media/nadaline/NADALINE II/The-Edge-of-Adaptation/genome_annotation"
NCBI="$BASE/ncbi_packages"
WORK="$BASE/work"

mkdir -p "$NCBI" "$WORK" "$BASE/results"

AMEIVA_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/046/270/845/GCA_046270845.1_ASM4627084v1/GCA_046270845.1_ASM4627084v1_genomic.fna.gz"
NEW_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/023/333/525/GCA_023333525.2_rAspTig1.1.pri/GCA_023333525.2_rAspTig1.1.pri_genomic.fna.gz"
GFF_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/023/333/525/GCA_023333525.2_rAspTig1.1.pri/GCA_023333525.2_rAspTig1.1.pri_genomic.gff.gz"

curl -L --fail --continue-at - "$AMEIVA_URL" -o "$NCBI/GCA_046270845.1_genomic.fna.gz"
curl -L --fail --continue-at - "$NEW_URL" -o "$NCBI/GCA_023333525.2_genomic.fna.gz"
curl -L --fail --continue-at - "$GFF_URL" -o "$NCBI/GCA_023333525.2_genomic.gff.gz"

gzip -cd "$NCBI/GCA_046270845.1_genomic.fna.gz" > "$WORK/GCA_046270845.1_genomic.fna"
gzip -cd "$NCBI/GCA_023333525.2_genomic.fna.gz" > "$WORK/GCA_023333525.2_genomic.fna"
gzip -cd "$NCBI/GCA_023333525.2_genomic.gff.gz" > "$WORK/GCA_023333525.2_genomic.gff3"

if ! command -v liftoff >/dev/null 2>&1; then
  echo "Liftoff nao encontrado. Instale-o antes de continuar." >&2
  echo "Sugestao: conda create -p /tmp/liftoff_env -c conda-forge -c bioconda liftoff" >&2
  exit 1
fi

liftoff \
  -g "$WORK/GCA_023333525.2_genomic.gff3" \
  -o "$WORK/Aspidoscelis_on_Ameiva_GCA_046270845.1.gff3" \
  -u "$WORK/unmapped_features.txt" \
  -dir "$WORK/liftoff_intermediate" \
  -p 4 \
  "$WORK/GCA_046270845.1_genomic.fna" \
  "$WORK/GCA_023333525.2_genomic.fna"

echo "Liftoff concluido: $WORK/Aspidoscelis_on_Ameiva_GCA_046270845.1.gff3"
