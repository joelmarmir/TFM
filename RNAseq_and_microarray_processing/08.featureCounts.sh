#!/bin/bash

#######################################################################
##### Extract the counts of the mapped reads (in parrallell jobs) #####
#######################################################################

# -----------------------------------------
# -------- Setup the job in SLURM ---------
# -----------------------------------------

#SBATCH --job-name=featureCounts
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=90
#SBATCH --mem=500G
#SBATCH --partition=bigmem
#SBATCH --output=/users/genomics/jmartinez/featureCounts_3_log.out
#SBATCH --nodelist=node17

# ---------------------------------
# ------ Run FEATURECOUNTS  -------
# ---------------------------------

# Load fastp module
module load Subread/2.0.3

# Directories
input_dir="/users/genomics/jmartinez/data/03_STAR_files/"
output_dir="/users/genomics/jmartinez/data/04_counts"
human_gtf="/users/genomics/jmartinez/data/00_reference_genomes/human/raw/gencode.v47.primary_assembly.annotation.gtf"
mouse_gtf="/users/genomics/jmartinez/data/00_reference_genomes/mouse/raw/gencode.vM36.primary_assembly.annotation.gtf"

STUDY_FOLDERS=(
    "Sabo_2014"
    "See_2022"
    "Vainorius_2023"
    "Walz_2014"
    "Wang_2023"
    "Wapinski_2013"
    "Weber_2016"
    "Woods_2023"
    "Yu_2023"
)

start_time=$(date +%s)

# Function to process FASTQ files
read_count() {
    local FILE="$1"
    local OUTPUT_STUDY_DIR="$2"

    # Determine genome indices based on the file name
    if [[ "$FILE" == *"hg38"* ]]; then
        genome_annot="$human_gtf"
    elif [[ "$FILE" == *"mm10"* ]]; then
        genome_annot="$mouse_gtf"
    else
        echo "Warning: Genome reference not found in file name. Skipping $FILE"
        return
    fi

    # Case 1: count single-end BAM files
    if [[ "$FILE" == *"single"* ]]; then
        featureCounts \
            -T 15 \
            -a "${genome_annot}" \
            -o "${OUTPUT_STUDY_DIR}/$(basename "$FILE" | sed 's/_mapped_/_counts_/').txt" \
            "$FILE"
    # Case 2: count paired-end BAM files    
    elif [[ "$FILE" == *"paired"* ]]; then
        featureCounts \
            -T 15 \
            -p \
            --countReadPairs \
            -B \
            -a "${genome_annot}" \
            -o "${OUTPUT_STUDY_DIR}/$(basename "$FILE" | sed 's/_mapped_/_counts_/').txt" \
            "$FILE"
    # Error handling
    else
        echo "Warning: File not found. Skipping $FILE"
    fi
}
export -f read_count

# Set maximum number of parallel jobs
MAX_JOBS=6

# Process all studies and files in parallel
for STUDY_NAME in "${STUDY_FOLDERS[@]}"; do
    STUDY_SUBFOLDER="${input_dir}/${STUDY_NAME}/"
    OUTPUT_STUDY_DIR="${output_dir}/${STUDY_NAME}"
    mkdir -p "$OUTPUT_STUDY_DIR"
    echo "Processing study: $STUDY_NAME"

    # Use process substitution to avoid subshell issues
    while IFS= read -r FILE; do
        read_count "$FILE" "$OUTPUT_STUDY_DIR" &
        # Enforce job limit dynamically
        while [[ $(jobs -r | wc -l) -ge $MAX_JOBS ]]; do
            wait -n  # Wait for any job to finish
        done
    done < <(find "$STUDY_SUBFOLDER" -maxdepth 1 -name "*.out.bam")  # Process substitution
done

# Wait for ALL jobs to finish before calculating time
wait

end_time=$(date +%s)
elapsed_time=$((end_time - start_time))

echo "Counts completed!"
echo "Total time taken: $elapsed_time seconds"