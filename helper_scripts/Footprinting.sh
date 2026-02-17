#!/usr/bin/env bash
set -euo pipefail

# -----------------------------------------------------------------------------
# ATAC-seq Postprocessing and Footprinting pipeline
# Description: Processes ATAC-seq BAM files, performs size selection, blacklist filtering, peak calling, Tn5 bias correction, footprint scoring, and differential analysis for biological conditions.
# 
# Requirements:
#   - hg38-blacklist.bed: The ChIP-seq blacklisted regions from ENCODE_v3. (not provided)
#   - GRCh38.fa: The hg38 fasta genome sequence (not provided)
#   - gaining_met.bed/losing_met.bed: Genomic regions bed files as outputs from Analysis.Rmd (Provided in ../Data/)
# -----------------------------------------------------------------------------

# Activate bashbone environment (More: https://github.com/Hoffmann-Lab/bashbone)

source <path/of/installation>/latest/bashbone/activate.sh -c true

# -----------------------------------------------------------------------------
# Step 1: Size selection and blacklist filtering
# -----------------------------------------------------------------------------
for bam in ../Data/ATACseq/*.bam; do
    echo "[INFO] Processing BAM: $bam"
    # Size selection
    alignment::postprocess -t 32 -j sizeselect -x 0:120 "$bam"
    sizeselect_bam="${bam%.bam}.sizeselected.bam"
    # Remove blacklisted regions
    alignment::postprocess -t 32 -j blacklist -x ../Data/hg38-blacklist.bed "$sizeselect_bam"
done

# -----------------------------------------------------------------------------
# Step 2: Peak calling with MACS2
# -----------------------------------------------------------------------------
mkdir -p ../Data/ATACseq/macs2_output
macs2 callpeak \
    -t ../Data/ATACseq/*.sizeselected.blacklisted.bam \
    -f BAMPE \
    -g hs \
    -n Combined \
    -q 0.01 \
    --nomodel \
    --outdir ../Data/ATACseq/macs2_output

# -----------------------------------------------------------------------------
# Step 3: Correcting for Tn5 transposase bias (TOBIAS ATACorrect)
# -----------------------------------------------------------------------------
mkdir -p ../Data/ATACseq/ATACorrect_test2
for bam in ../Data/ATACseq/*.sizeselected.blacklisted.bam; do
    echo "[INFO] Correcting Tn5 bias for: $bam"
    TOBIAS ATACorrect \
        --bam "$bam" \
        --genome ../Data/GRCh38.fa \
        --peaks ../Data/ATACseq/macs2_output/Combined_peaks.narrowPeak \
        --blacklist ../Data/hg38-blacklist.bed \
        --outdir ../Data/ATACseq/ATACorrect_test2 \
        --cores 32
done

# -----------------------------------------------------------------------------
# Step 4: Compute footprint scores
# -----------------------------------------------------------------------------
mkdir -p ../Data/ATACseq/TOBIAS_score2
for corrected_bw in ../Data/ATACseq/ATACorrect_test2/*_corrected.bw; do
    sample=$(basename "$corrected_bw" "_corrected.bw")
    echo "[INFO] Calculating footprint scores for: $sample"
    TOBIAS FootprintScores \
        --signal "$corrected_bw" \
        --regions ../Data/ATACseq/macs2_output/Combined_peaks.narrowPeak \
        --output ../Data/ATACseq/TOBIAS_score2/"$sample".bw \
        --cores 32
done

# -----------------------------------------------------------------------------
# Step 5: Aggregate replicates by condition
# -----------------------------------------------------------------------------
bigwigAverage \
    -b ../Data/ATACseq/TOBIAS_score2/SRR17466788.bw ../Data/ATACseq/TOBIAS_score2/SRR17466798.bw ../Data/ATACseq/TOBIAS_score2/SRR17466717.bw ../Data/ATACseq/TOBIAS_score2/SRR17466739.bw \
    -o ../Data/ATACseq/TOBIAS_score2/Y.bw -v -p max

bigwigAverage \
    -b ../Data/ATACseq/TOBIAS_score2/SRR17466791.bw ../Data/ATACseq/TOBIAS_score2/SRR17466779.bw ../Data/ATACseq/TOBIAS_score2/SRR17466729.bw ../Data/ATACseq/TOBIAS_score2/SRR17466768.bw \
    -o ../Data/ATACseq/TOBIAS_score2/O.bw -v -p max

# -----------------------------------------------------------------------------
# Step 6: Differential signal (log2 fold change)
# -----------------------------------------------------------------------------
bigwigCompare \
    -b1 ../Data/ATACseq/TOBIAS_score2/O.bw \
    -b2 ../Data/ATACseq/TOBIAS_score2/Y.bw \
    --operation log2 \
    --binSize 1 \
    --skipZeroOverZero \
    --numberOfProcessors max/2 \
    -o ../Data/ATACseq/TOBIAS_score2/diff.bw

# -----------------------------------------------------------------------------
# Step 7: Signal in TFBS regions gaining/losing methylation
# -----------------------------------------------------------------------------
computeMatrix scale-regions \
    -S ../Data/ATACseq/TOBIAS_score2/diff.bw \
    -R ../Data/gaining_met.bed ../Data/losing_met.bed \
    -a 5000 -b 5000 \
    --missingDataAsZero \
    --binSize 1 \
    --numberOfProcessors max/2 \
    -o ../Data/ATACseq/TOBIAS_score2/Footprinting_diff.gz

# -----------------------------------------------------------------------------
# Step 8: Plotting footprinting signal
# -----------------------------------------------------------------------------
plotProfile \
    -m ../Data/ATACseq/TOBIAS_score2/Footprinting_diff.gz \
    -o ../Data/ATACseq/TOBIAS_score2/Footprinting_diff.pdf \
    --regionsLabel MetGain MetLoose \
    --averageType mean \
    --dpi 100 \
    --plotType lines \
    --startLabel TF \
    --endLabel BS

echo "[INFO] Pipeline completed successfully."
