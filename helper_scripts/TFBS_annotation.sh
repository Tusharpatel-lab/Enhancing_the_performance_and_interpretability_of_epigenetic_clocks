#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# TFBS annotation pipeline
# Description: Download JASPAR motifs for hg38 and intersect with combined human blood TF ChIP-seq profiles. 
#
# Requirements:
#   - genomes_GRCh38_v106_chr.info: Chromosome lengths for making bins. (not provided)
#   - combined_TF.bed: Merged TF ChIP-seq profiles from cistrome (provided as ../Data/combined_TF.bed). 
#
# Others:
#   - helper::lapply (from bashbone, More: https://github.com/Hoffmann-Lab/bashbone)
###############################################################################

# -----------------------
# Configurable paths
# -----------------------

source <path/of/installation>/latest/bashbone/activate.sh -c true

GENOME_INFO="../Data/genomes_GRCh38_v106_chr.info"
JASPAR_URL="http://hgdownload.soe.ucsc.edu/gbdb/hg38/jaspar/JASPAR2024.bb"

OUT_JASPAR="../Data/JASPAR_whole_genome_TFs.bed"
COMBINED_TF="../Data/combined_TF.bed"

THREADS=128
WINDOW_SIZE=10000

# -----------------------
# Download JASPAR motifs
# -----------------------

dl() {
  read -r chr sta sto < /dev/stdin
  bigBedToBed "$JASPAR_URL" \
    -chrom="$chr" \
    -start="$sta" \
    -end="$sto" \
    stdout
}
export -f dl

bedtools makewindows \
  -w "$WINDOW_SIZE" \
  -g "$GENOME_INFO" \
| helper::lapply -l 1 -d 1 -t "$THREADS" -k -c dl \
> "$OUT_JASPAR"

# -----------------------
# Intersect with TF ChIP-seq profiles
# -----------------------

awk 'FNR==NR {keys[$4]; next} $7 in keys' \
  "$COMBINED_TF" \
  "$OUT_JASPAR" \
> ../Data/filtered_file1.bed

awk 'FNR==NR {keys[$7]; next} $4 in keys' \
  "$OUT_JASPAR" \
  "$COMBINED_TF" \
> ../Data/filtered_file2.bed

bedtools intersect \
  -a ../Data/filtered_file1.bed \
  -b ../Data/filtered_file2.bed \
> ../Data/intersection.bed

sort -k1,1 -k2,2n ../Data/intersection.bed \
> ../Data/sorted_intersection.bed

bedtools merge \
  -i ../Data/sorted_intersection.bed \
  -c 7 \
  -o distinct \
> ../Data/merged_intersection.bed