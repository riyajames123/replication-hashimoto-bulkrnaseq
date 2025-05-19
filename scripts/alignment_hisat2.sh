#!/bin/bash
#SBATCH --job-name=hisat2_cleaned
#SBATCH --output=/scratch/james.ri/trans_proj/logs/hisat2_cleaned_%j.out
#SBATCH --error=/scratch/james.ri/trans_proj/logs/hisat2_cleaned_%j.err
#SBATCH --partition=courses
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G

# Activate conda environment
source ~/.bashrc
conda activate hisat2_env

# Define paths
WORKDIR="/scratch/james.ri/trans_proj"
REF_DIR="$WORKDIR/ref"
INPUT_DIR="$WORKDIR/ribodetector_cleaned"
OUTPUT_DIR="$WORKDIR/hisat2_cleaned_output"
INDEX_PREFIX="$REF_DIR/hisat2_index/GRCh38_genome"
LOGDIR="$WORKDIR/logs"
SUMMARY_REPORT="$LOGDIR/hisat2_alignment_summary_cleaned.tsv"

# Create output directories if they don't exist
mkdir -p "$OUTPUT_DIR"
mkdir -p "$LOGDIR"

# Write header for summary report
echo -e "Sample\tTotal_Reads\tAligned_Reads\tAlignment_Rate" > "$SUMMARY_REPORT"

# Loop through all forward reads (assuming matching reverse reads exist)
for fwd_read in "$INPUT_DIR"/*_cleaned_1.fq; do
    sample_name=$(basename "$fwd_read" _cleaned_1.fq)
    rev_read="$INPUT_DIR/${sample_name}_cleaned_2.fq"

    echo "Processing $sample_name..."
    hisat2 \
        -p 16 \
        --dta \
        --rna-strandness RF \
        -x "$INDEX_PREFIX" \
        -1 "$fwd_read" \
        -2 "$rev_read" \
        -S "$OUTPUT_DIR/${sample_name}.sam" \
        2> "$LOGDIR/${sample_name}_hisat2_cleaned.log"

    # Convert SAM to BAM
    samtools view -@ 16 -bS "$OUTPUT_DIR/${sample_name}.sam" | samtools sort -@ 16 -o "$OUTPUT_DIR/${sample_name}.sorted.bam"
    samtools index "$OUTPUT_DIR/${sample_name}.sorted.bam"
    rm "$OUTPUT_DIR/${sample_name}.sam"

    # Extract alignment summary from log
    total_reads=$(grep "pairs of reads in total" "$LOGDIR/${sample_name}_hisat2_cleaned.log" | awk '{print $1}')
    aligned_reads=$(grep "aligned concordantly 0 times" "$LOGDIR/${sample_name}_hisat2_cleaned.log" | awk '{print $1}')
    aligned_percent=$(grep "overall alignment rate" "$LOGDIR/${sample_name}_hisat2_cleaned.log" | awk '{print $1}')

    # Append to summary report
    echo -e "${sample_name}\t${total_reads}\t${aligned_reads}\t${aligned_percent}" >> "$SUMMARY_REPORT"
done

echo "All alignments complete."
