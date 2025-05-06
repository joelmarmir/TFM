#!/bin/bash

#####################################################################
##### STAR: align reads to reference genome (in parallell jobs) #####
#####################################################################

# -----------------------------------------
# -------- Setup the job in SLURM ---------
# -----------------------------------------

#SBATCH --job-name=STAR_alingment
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=90
#SBATCH --mem=500G
#SBATCH --partition=bigmem
#SBATCH --output=/users/genomics/jmartinez/STAR_alingment_6_log.out
#SBATCH --nodelist=node17

# ------------------------
# ------ Run STAR  -------
# ------------------------

# Load fastp module
module load STAR/2.7.8a-GCC-10.2.0

# Directories
input_dir="/users/genomics/jmartinez/data/02_preprocessed_fastqs/00_temp"
output_dir="/users/genomics/jmartinez/data/03_STAR_files/00_temp"
human_index_dir="/users/genomics/jmartinez/data/00_reference_genomes/human/indices"
mouse_index_dir="/users/genomics/jmartinez/data/00_reference_genomes/mouse/indices"

STUDY_FOLDERS=(
    "Aydin_2019"
    "Casey_2018"
    # Complete with the desired studies
)

start_time=$(date +%s)

# Function to process FASTQ files
STAR_align() {
    local FILE="$1"
    local OUTPUT_STUDY_DIR="$2"

    # Determine genome indices based on the file name
    if [[ "$FILE" == *"hg38"* ]]; then
        genome_indices="$human_index_dir"
    elif [[ "$FILE" == *"mm10"* ]]; then
        genome_indices="$mouse_index_dir"
    else
        echo "Warning: Genome reference not found in file name. Skipping $FILE"
        return
    fi

    # Case 1: Skip _2.fastq.gz files as they will be processed with their pairs
    if [[ "$FILE" =~ _2_trimmed.fastq.gz$ ]]; then
        return

    # Case 2: Process _1.fastq.gz files with their pairs   
    elif [[ "$FILE" =~ _1_trimmed.fastq.gz$ ]]; then
        PAIR_1="$FILE"
        PAIR_2="${FILE/_1_trimmed.fastq.gz/_2_trimmed.fastq.gz}"

        if [[ -f "$PAIR_2" ]]; then
            STAR \
              --runThreadN 18 \
              --genomeDir "${genome_indices}" \
              --readFilesIn "$PAIR_1" "$PAIR_2" \
              --readFilesCommand zcat \
              --outFileNamePrefix "${OUTPUT_STUDY_DIR}/$(basename "$PAIR_1" | sed 's/_1_trimmed.fastq.gz/_paired_mapped_/')" \
              --outSAMtype BAM SortedByCoordinate

        else
            echo "Warning: Pair file for ${PAIR_1} not found. Skipping..."
        fi
    # Case 3: Process single-end FASTQ files
    else
        STAR \
          --runThreadN 18 \
          --genomeDir "${genome_indices}" \
          --readFilesIn "$FILE" \
          --readFilesCommand zcat \
          --outFileNamePrefix "${OUTPUT_STUDY_DIR}/$(basename "$FILE" | sed 's/_trimmed.fastq.gz/_single_mapped_/')" \
          --outSAMtype BAM SortedByCoordinate
    fi
}

export -f STAR_align

# Set maximum number of parallel jobs
MAX_JOBS=5

# Process all studies and files in parallel
for STUDY_NAME in "${STUDY_FOLDERS[@]}"; do
    STUDY_SUBFOLDER="${input_dir}/${STUDY_NAME}/"
    OUTPUT_STUDY_DIR="${output_dir}/${STUDY_NAME}"
    mkdir -p "$OUTPUT_STUDY_DIR"
    echo "Processing study: $STUDY_NAME"

    # Use process substitution to avoid subshell issues
    while IFS= read -r FILE; do
        STAR_align "$FILE" "$OUTPUT_STUDY_DIR" &
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