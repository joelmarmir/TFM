#!/bin/bash

######################################################
#### Create a list with the SRAs grouped by study ####
######################################################

# Set input and output file
input="/users/genomics/jmartinez/temp_data/expression_datasets_csv_with_sra_v01.csv"
output="/users/genomics/jmartinez/temp_data/study_srr_list.txt"

# Process the input file, skipping the header (tail -n +2), and use awk to extract and format the data
tail -n +2 "$input" | awk -F',' '  # Use a comma as the field separator
{
    # Extract the study name from the second column
    study = $2
    gsub(/^[ \t"]+|[ \t"]+$/, "", study)  # Remove leading and trailing whitespace or quotes from the study name

    # Loop through columns 8 and 9 (assumed to contain SRAs IDs)
    for (col = 8; col <= 9; col++) {
        srr = $col  # Extract the value from the current column
        gsub(/^[ \t"]+|[ \t"]+$/, "", srr)  # Remove leading and trailing whitespace or quotes from the SRR ID

        # Check if the value matches the SRR format (e.g., SRR123456)
        if (srr ~ /^SRR[0-9]+$/) {
            # Append the SRR ID to the corresponding study in the `studies` array
            studies[study] = studies[study] srr "\n"

            # If the study has not been seen before, add it to the `order` array to preserve order
            if (!(study in seen)) {
                order[++n] = study  # Add the study to the order array
                seen[study] = 1    # Mark the study as seen
            }
        }
    }
}
END {
    # After processing all lines, output the studies and their SRR IDs in the desired format
    for (i = 1; i <= n; i++) {  # Iterate through the studies in the order they were encountered
        study = order[i]  # Get the study name
        printf(">%s\n%s", study, studies[study])  # Print the study name and its associated SRR IDs
    }
}' > "$output"  # Redirect the output to the specified file