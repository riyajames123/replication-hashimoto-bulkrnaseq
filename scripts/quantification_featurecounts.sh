#!/bin/bash
#SBATCH --job-name=featurecounts
#SBATCH --output=/scratch/james.ri/trans_proj/logs/featurecounts_%j.out
#SBATCH --error=/scratch/james.ri/trans_proj/logs/featurecounts_%j.err
#SBATCH --partition=courses
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G

# Activate conda environment
source ~/.bashrc
conda activate subread_env  

# Paths
WORKDIR="/scratch/james.ri/trans_proj"
GTF="$WORKDIR/ref/gencode.v46.annotation.gtf"
BAM_DIR="$WORKDIR/hisat2_cleaned_output"
OUTFILE="$WORKDIR/featurecounts_new_output/gene_counts.tsv"
LOGDIR="$WORKDIR/logs"

# Create output directory if it doesn't exist
mkdir -p "$(dirname "$OUTFILE")"
mkdir -p "$LOGDIR"

# Run featureCounts
featureCounts \
  -T 4 \
  -a "$GTF" \
  -o "$OUTFILE" \
  -g gene_id \
  -t exon \
  -p \
  -s 2 \
  "$BAM_DIR"/*.sorted.bam \
  2> "$LOGDIR/featurecounts_summary.log"

echo "featureCounts quantification complete. Output: $OUTFILE"
