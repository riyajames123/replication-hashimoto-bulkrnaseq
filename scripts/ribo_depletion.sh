#!/bin/bash
#SBATCH --job-name=ribodetector_array
#SBATCH --partition=courses
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=16:00:00
#SBATCH --array=0-9  # Adjust this to match number of samples - 1
#SBATCH --output=/scratch/james.ri/trans_proj/logs/ribodetector_%A_%a.out
#SBATCH --error=/scratch/james.ri/trans_proj/logs/ribodetector_%A_%a.err

# Load environment
source ~/.bashrc
conda activate ribodetector_env

# Set paths
INPUT_DIR="/scratch/forbes.ai/trans_final_data/final_project/fastp_no_overrep/trimmed_fastp"
OUTPUT_DIR="/scratch/james.ri/trans_proj/ribodetector_cleaned"
LOG_DIR="/scratch/james.ri/trans_proj/logs"
SUMMARY="$LOG_DIR/ribodetector_summary.tsv"
THREADS=16
READ_LENGTH=141

# Create output dirs if needed
mkdir -p "$OUTPUT_DIR"
mkdir -p "$LOG_DIR"

# Get list of samples
SAMPLES=($(ls $INPUT_DIR/*_f1.fastq.gz | xargs -n 1 basename | sed 's/_f1.fastq.gz//'))
SAMPLE_NAME=${SAMPLES[$SLURM_ARRAY_TASK_ID]}
FWD="$INPUT_DIR/${SAMPLE_NAME}_f1.fastq.gz"
REV="$INPUT_DIR/${SAMPLE_NAME}_r2.fastq.gz"
OUT1="$OUTPUT_DIR/${SAMPLE_NAME}_cleaned_1.fq"
OUT2="$OUTPUT_DIR/${SAMPLE_NAME}_cleaned_2.fq"

# Skip if already done
if [[ -f "$OUT1" && -f "$OUT2" ]]; then
    echo "[$(date)] Skipping $SAMPLE_NAME — already processed" >> "$LOG_DIR/ribodetector_job.log"
    exit 0
fi

# Run RiboDetector
ribodetector_cpu \
    -t $THREADS \
    -l $READ_LENGTH \
    -i "$FWD" "$REV" \
    -e rrna \
    --chunk_size 256 \
    -o "$OUT1" "$OUT2" \
    >> "$LOG_DIR/${SAMPLE_NAME}_ribodetector.out" \
    2>> "$LOG_DIR/${SAMPLE_NAME}_ribodetector.err"

# Count and log if successful
if [[ -f "$OUT1" ]]; then
    orig_reads=$(zcat "$FWD" | wc -l)
    orig_reads=$((orig_reads / 4))
    cleaned_reads=$(grep -c "^+$" "$OUT1")
    removed=$((orig_reads - cleaned_reads))
    percent_removed=$(awk -v r="$removed" -v t="$orig_reads" 'BEGIN { printf "%.2f", (r / t) * 100 }')

    # Ensure summary file has a header (if task 0 and file missing)
    if [[ "$SLURM_ARRAY_TASK_ID" -eq 0 && ! -f "$SUMMARY" ]]; then
        echo -e "Sample\tOriginal_Reads\tCleaned_Reads\tRemoved\tPercent_Removed" > "$SUMMARY"
    fi

    echo -e "${SAMPLE_NAME}\t${orig_reads}\t${cleaned_reads}\t${removed}\t${percent_removed}%" >> "$SUMMARY"
    echo "[$(date)] Finished $SAMPLE_NAME — Removed ${removed} reads (${percent_removed}%)" >> "$LOG_DIR/ribodetector_job.log"
else
    echo "[$(date)] WARNING: $SAMPLE_NAME failed or output missing." >> "$LOG_DIR/ribodetector_job.log"
    exit 1
fi
