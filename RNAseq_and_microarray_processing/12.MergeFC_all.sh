#!/bin/bash

####################################
##### Merge FoldChange Tables  #####
####################################

# -----------------------------------------
# -------- Setup the job in SLURM ---------
# -----------------------------------------

#SBATCH --job-name=MergedFC
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=20G
#SBATCH --partition=normal
#SBATCH --output=/users/genomics/jmartinez/MergedFC_with_micro_1_log.out
#SBATCH --nodelist=node04

# ------------------------------------------------
# -------- Merge FoldChanges per species ---------
# ------------------------------------------------

MAIN_DIR="/users/genomics/jmartinez/data/06_fc"

STUDY_FOLDERS_MOUSE=(
    "Vainorius_2023" "Casey_2018" "Aydin_2019" "Pataskar_2016"
    "Fong_2015" "Fong_2012" "Sabo_2014" "dePretis_2017" "Li_2019"
    "Croci_2017" "Kress_2016" "Matsuda_2016" "Pereira_2024" "Wapinski_2013"
    "Weber_2016" "Conerly_2016" "Lee_2020" "Chalamasetty_2014"
    "Liu_2014"
)

STUDY_FOLDERS_HUMAN=(
    "Jung_2017" "Lorenzin_2016" "See_2022" "Walz_2014" "Jaenicke_2016"
    "Manandhar_2016" "Smith_2016" "Li_2023" "Woods_2023" "Park_2017"
    "Wang_2023" "Yu_2023" "Joung_2023" "Barfeld_2017"
)

MOUSE_OUT="/users/genomics/jmartinez/data/07_for_results/fc_all_mouse_with_micro.csv"
HUMAN_OUT="/users/genomics/jmartinez/data/07_for_results/fc_all_human_with_micro.csv"

# Create headers
echo "File,ENSEMBL,GeneSymbol,log2FoldChange,pvalue,padj" > "$MOUSE_OUT"
echo "File,ENSEMBL,GeneSymbol,log2FoldChange,pvalue,padj" > "$HUMAN_OUT"

# Function to process files
process_file() {
    local file="$1"
    local output="$2"
    local study_type="$3"

    echo "  → Processing file: $(basename "$file")"

    # Use FPAT to parse quoted CSV fields and remove quotes
    tail -n +2 "$file" | awk -v FPAT='([^,]*)|("([^"]|"")*")' \
        -v output="$output" \
        -v study_type="$study_type" \
        -v filename="$(basename "$file")" '
        BEGIN { OFS="," }
        {
            # Remove quotes from every field
            for (i=1; i<=NF; i++) {
                gsub(/^"|"$/, "", $i);
            }

            # Check if the file is RNAseq or microarray
            if ($1 ~ /^ENS/) {
                # RNAseq processed file
                ensembl = $1;
                gene_symbol = $2;
                log2fc = $3;
                pvalue = $6;
                padj = $7;
            } else {
                # Microarray processed file
                ensembl = "NA";  # Add NA for ENSEMBL
                gene_symbol = $1;
                log2fc = $2;
                pvalue = $5;
                padj = $6;
            }

            # Output the merged data with the file name as the first column
            print filename, ensembl, gene_symbol, log2fc, pvalue, padj >> output;
        }'
}

# Process mouse studies
for study in "${STUDY_FOLDERS_MOUSE[@]}"; do
    study_path="$MAIN_DIR/$study"
    if [[ -d "$study_path" ]]; then
        echo "Processing mouse study: $study"
        for file in "$study_path"/*.csv; do
            if [[ -f "$file" ]]; then
                process_file "$file" "$MOUSE_OUT" "mouse"
            fi
        done
    else
        echo "Skipping mouse study: $study (not found)" >&2
    fi
done

# Process human studies
for study in "${STUDY_FOLDERS_HUMAN[@]}"; do
    study_path="$MAIN_DIR/$study"
    if [[ -d "$study_path" ]]; then
        echo "Processing human study: $study"
        for file in "$study_path"/*.csv; do
            if [[ -f "$file" ]]; then
                process_file "$file" "$HUMAN_OUT" "human"
            fi
        done
    else
        echo "Skipping human study: $study (not found)" >&2
    fi
done

echo "Processing complete! Outputs: $MOUSE_OUT and $HUMAN_OUT"