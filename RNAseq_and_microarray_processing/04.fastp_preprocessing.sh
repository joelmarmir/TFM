#!/bin/bash

############################################################################
##### Quality control and pre-processing of FASTQs (in parallell jobs) #####
############################################################################

# -----------------------------------------
# -------- Setup the job in SLURM ---------
# -----------------------------------------

#SBATCH --job-name=fastp_preprocess
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=80
#SBATCH --mem=500G
#SBATCH --partition=bigmem
#SBATCH --output=/users/genomics/jmartinez/fastp_preprocess_2_log.out
#SBATCH --nodelist=node17

# --------------------------------------------
# ------ Use FASTP to preprocess FASTQs -------
# ---------------------------------------------

# Load fastp module
module load fastp

# Directories
BASE_FASTQ_DIR="/users/genomics/jmartinez/data/01_fastqs/00_temp/"
BASE_OUTPUT_DIR="/users/genomics/jmartinez/data/02_preprocessed_fastqs/00_temp/"

# Define which studies to process
STUDY_FOLDERS+=(
    "Conerly_2016"
)

start_time=$(date +%s)

# Function to process FASTQ files
process_fastq() {
    local FILE="$1"
    local OUTPUT_STUDY_DIR="$2"

    # Case 1: Skip _2.fastq.gz files as they will be processed with their pairs
    if [[ "$FILE" =~ _2.fastq.gz$ ]]; then
        return

    # Case 2: Process _1.fastq.gz files with their pairs   
    elif [[ "$FILE" =~ _1.fastq.gz$ ]]; then
        PAIR_1="$FILE"
        PAIR_2="${FILE/_1.fastq.gz/_2.fastq.gz}"

        if [[ -f "$PAIR_2" ]]; then
            fastp -i "$PAIR_1" -I "$PAIR_2" \
                  -o "${OUTPUT_STUDY_DIR}/$(basename "$PAIR_1" | sed 's/_1.fastq.gz/_1_trimmed.fastq.gz/')" \
                  -O "${OUTPUT_STUDY_DIR}/$(basename "$PAIR_2" | sed 's/_2.fastq.gz/_2_trimmed.fastq.gz/')" \
                  -h "${OUTPUT_STUDY_DIR}/$(basename "$PAIR_1" | sed 's/_1.fastq.gz/_fastp_report_paired.html/')" \
                  -j "${OUTPUT_STUDY_DIR}/$(basename "$PAIR_1" | sed 's/_1.fastq.gz/_fastp_report_paired.json/')" \
                  --thread 16 --detect_adapter_for_pe
        else
            echo "Warning: Pair file for ${PAIR_1} not found. Skipping..."
        fi
    # Case 3: Process single-end FASTQ files
    else
        fastp -i "$FILE" -o "${OUTPUT_STUDY_DIR}/$(basename "$FILE" | sed 's/.fastq.gz/_trimmed.fastq.gz/')" \
              -h "${OUTPUT_STUDY_DIR}/$(basename "$FILE" | sed 's/.fastq.gz/_fastp_report_single.html/')" \
              -j "${OUTPUT_STUDY_DIR}/$(basename "$FILE" | sed 's/.fastq.gz/_fastp_report_single.json/')" \
              --thread 16
    fi
}

export -f process_fastq

# Set maximum number of parallel jobs
MAX_JOBS=5

# Process all studies and files in parallel
for STUDY_NAME in "${STUDY_FOLDERS[@]}"; do
    STUDY_SUBFOLDER="${BASE_FASTQ_DIR}/${STUDY_NAME}/"
    OUTPUT_STUDY_DIR="${BASE_OUTPUT_DIR}/${STUDY_NAME}"
    mkdir -p "$OUTPUT_STUDY_DIR"
    echo "Processing study: $STUDY_NAME"

    # Use process substitution to avoid subshell issues
    while IFS= read -r FILE; do
        process_fastq "$FILE" "$OUTPUT_STUDY_DIR" &
        # Enforce job limit dynamically
        while [[ $(jobs -r | wc -l) -ge $MAX_JOBS ]]; do
            wait -n  # Wait for any job to finish
        done
    done < <(find "$STUDY_SUBFOLDER" -maxdepth 1 -name "*.fastq.gz")  # Process substitution
done

# Wait for ALL jobs to finish before calculating time
wait

end_time=$(date +%s)
elapsed_time=$((end_time - start_time))

echo "FASTQ processing completed!"
echo "Total time taken: $elapsed_time seconds"