#!/bin/bash

########################################################################################
##### Differential expression analysis (in parallell jobs) - part 1: sending jobs  #####
########################################################################################

# -----------------------------------------
# -------- Setup the job in SLURM ---------
# -----------------------------------------

#SBATCH --job-name=Deseq2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --mem=240G
#SBATCH --partition=long
#SBATCH --output=/users/genomics/jmartinez/Deseq2_6_log.out
#SBATCH --nodelist=node14

# -------------------------------
# ------ Run DE analysis  -------
# -------------------------------

# Load the R module if needed
module load R

# Base directories
phenodata_base_dir="/users/genomics/jmartinez/data/05_tables"
counts_base_dir="/users/genomics/jmartinez/data/04_counts"
results_dir="/users/genomics/jmartinez/data/06_log2fc"

# List of study names to process
STUDY_FOLDERS=(
    "Casey_2018"
    "Costa_2022"
    "Croci_2017"
    "Jaenicke_2016"
    "Joung_2023"
    "Jung_2017"
    "Kress_2016"
    "Letourneau_2014"
    "Lorenzin_2016"
    "Manandhar_2016"
    "Matsuda_2016"
    "Park_2017"
    "Pataskar_2016"
    "Pereira_2024"
    "Sabo_2014"
    "See_2022"
    "Smith_2016"
    "Vainorius_2023"
    "Wang_2023"
    "Wapinski_2013"
    "Woods_2023"
)


# Function to process a single phenodata file
process_phenodata_file() {
    local phenodata_file="$1"
    local study_counts_dir="$2"
    local study_results_dir="$3"
    
    # Get the base name (without the .csv extension)
    base_name=$(basename "$phenodata_file" .csv)
    
    # Replace "phenotab" with "DESeq2_results" in the file name
    result_base_name="${base_name/phenotab/DESeq2_results}"
    
    # Define the output file path using the study-specific results directory
    output_file="${study_results_dir}/${result_base_name}.csv"
    
    echo "Using phenodata file: $phenodata_file"
    echo "Counts directory: $study_counts_dir"
    echo "Output file: $output_file"
    
    # Execute the R script with arguments:
    # 1. Phenodata CSV path, 2. Counts directory, 3. Output file path
    Rscript /users/genomics/jmartinez/scripts/differential_expression_analysis_v3.r "$phenodata_file" "$study_counts_dir" "$output_file"
}

# Iterate over each study in the list
for study in "${STUDY_FOLDERS[@]}"; do
    echo "Processing study: $study"
    
    study_phenodata_dir="${phenodata_base_dir}/${study}"
    study_counts_dir="${counts_base_dir}/${study}"
    study_results_dir="${results_dir}/${study}"
    
    # Create study-specific results directory if it doesn't exist
    mkdir -p "$study_results_dir"
    
    # Check if the counts directory exists
    if [ ! -d "$study_counts_dir" ]; then
        echo "Counts directory not found: $study_counts_dir. Skipping study: $study"
        continue
    fi

    # Process each CSV file in the study's phenodata directory in parallel
    for phenodata_file in "$study_phenodata_dir"/*.csv; do
        if [ ! -f "$phenodata_file" ]; then
            echo "No CSV files found in $study_phenodata_dir. Skipping..."
            continue
        fi
        
        # Process the phenodata file in the background
        process_phenodata_file "$phenodata_file" "$study_counts_dir" "$study_results_dir" &
        
        # Limit the number of background jobs to the number of CPUs allocated
        while [ "$(jobs -r | wc -l)" -ge 20 ]; do
            wait -n
        done
    done
done

# Wait for all background jobs to finish
wait

echo "All studies processed."