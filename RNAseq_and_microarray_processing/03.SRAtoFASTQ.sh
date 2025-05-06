#!/bin/bash

#################################################################
##### Extract FASTQ files from SRAs IDs (in parallell jobs) #####
#################################################################

# -----------------------------------------
# -------- Setup the job in SLURM ---------
# -----------------------------------------

#SBATCH --nodelist=node06
#SBATCH --partition=normal
#SBATCH --cpus-per-task=10
#SBATCH --mem=100G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=SRAtoFASTQ
#SBATCH --output=/users/genomics/jmartinez/SRRtoFASTQ_test_log.out

# ---------------------------------
# ------ Get FASTQ from SRA -------
# ---------------------------------

# Load SRA toolkit module
module load SRA-Toolkit/3.0

# Define input file and output directory
srr_list="/users/genomics/jmartinez/temp_data/study_srr_list_test.txt"
output_dir="/users/genomics/jmartinez/data/01_fastqs/00_temp"
mkdir -p "$output_dir"

# Function to process each SRR
process_srr() {
    local SRR="$1"
    local study_folder="$2"
    echo "Processing: $SRR (Study: $study_folder)"
    
    # Extract fastqs from SRA
    fasterq-dump --split-3 --outdir "$output_dir/$study_folder" "$SRR"
    gzip "$output_dir/$study_folder"/*.fastq
}

# Read SRR IDs and their corresponding study folder names from srr_list
while IFS= read -r line; do
    # Check if the line is a study folder name
    if [[ "$line" =~ ^\> ]]; then
        study_folder="${line#>}"
        mkdir -p "$output_dir/$study_folder"
    elif [[ "$line" =~ ^SRR ]]; then
        SRR="$line"
        process_srr "$SRR" "$study_folder" &
        
        # Limit the number of background processes to the number of available CPU cores
        while [[ $(jobs -r -p | wc -l) -ge $SLURM_CPUS_PER_TASK ]]; do
            wait -n
        done
    fi
done < "$srr_list"

# Wait for all background processes to finish
wait

echo "FASTQ extraction completed!"
