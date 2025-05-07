############################################
##### Differential Expression Analysis #####
############################################

# Load libraries
library(DESeq2)
library(dplyr)
library(readr)
library(pheatmap)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(clusterProfiler)
library(ggplot2)


#--------------------------------------
#------------ Functions ---------------
#--------------------------------------

############## Data preparation ##############

# Define file paths
# Read command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: differential_expression_analysis_v3.r <phenodata_path> <counts_dir> <output_file>")
}

# Extract arguments
phenodata_path <- args[1]
counts_dir <- args[2]
output_file <- args[3]


# Read the phenotype table
phenodata <- read.csv(phenodata_path, header = TRUE, stringsAsFactors = TRUE)

# Ensure Replicate is a factor
phenodata$Replicate <- factor(phenodata$Replicate)

# Extract SRR codes (assuming they are in a specific column, adjust as needed)
# Let's assume the column is named "SRR_code"
srr_codes <- phenodata$SRA_id

# Initialize an empty list to store data
count_data_list <- list()

# Iterate over SRR codes and find corresponding count files
for (srr in srr_codes) {
    print(paste("Processing SRR code:", srr))  # Troubleshooting print
    
    count_file <- list.files(counts_dir, pattern = paste0(srr, ".*txt$"), full.names = TRUE)
    
    if (length(count_file) == 1) {  # Ensure a single match
        print(paste("Found count file:", count_file))  # Troubleshooting print
        
        count_table <- read.delim(count_file, header = TRUE, stringsAsFactors = FALSE, skip = 1)
        
        # Extract relevant columns: Geneid and counts (assumed to be column 7)
        count_data <- count_table[, c(1, 7)]
        colnames(count_data) <- c("Geneid", srr)
        
        # Store data
        count_data_list[[srr]] <- count_data
    } else {
        print(paste("Count file not found or multiple files found for SRR code:", srr))  # Troubleshooting print
    }
}

# Merge all count tables by Geneid
counts <- Reduce(function(x, y) full_join(x, y, by = "Geneid"), count_data_list)

# Set Geneid as row names
rownames(counts) <- counts$Geneid
counts <- counts[, -1]  # Remove Geneid column

# Generate deseq object with counts and condition information
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = phenodata,
                              design = ~Condition)


############## Deseq Analysis ##############

# Filter genes with less than 10 reads total
keep <- rowSums(counts(dds)) >= 10
dds_modified <- dds[keep,]

# Set uninduced as reference level
dds_modified$Condition <- relevel(dds_modified$Condition, ref = "uninduced")

# Manual dispersion setting (CRITICAL STEP FOR 1 REPLICATE)
dds_modified <- estimateSizeFactors(dds_modified)
# dds_modified <- estimateDispersionsGeneEst(dds_modified)  # Will warn but proceed

# Set dispersion to an arbitrary value (0.1 is example - adjust based on prior knowledge)
dispersions(dds_modified) <- 0.1

# Continue with Wald test
dds_modified <- nbinomWaldTest(dds_modified)
resultsNames(dds_modified)
res <- results(dds_modified, name = "Condition_induced_vs_uninduced")

# Extract log2FoldChange and ENSEMBL identifier
res_df <- as.data.frame(res)
res_df$ENSEMBL <- rownames(res_df)
res_df$ENSEMBL_short <- gsub("\\..*", "",row.names(res_df))

# Determine the appropriate annotation database
if ("hg38" %in% phenodata$Genome) {
  annotation_db <- org.Hs.eg.db
} else if ("mm10" %in% phenodata$Genome) {
  annotation_db <- org.Mm.eg.db
} else {
  stop("Unsupported genome version. Please use 'hg38' or 'mm10'.")
}

# Map ENSEMBL IDs to gene symbols
gene_symbols <- mapIds(annotation_db, keys = res_df$ENSEMBL_short, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")

# Add gene symbols to the results data frames
res_df$GeneSymbol <- gene_symbols

# Reorder columns to place GeneSymbol first
res_df <- res_df[, c("ENSEMBL", "GeneSymbol", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj", "baseMean")]

# Remove rownames
rownames(res_df) <- NULL
  
# Save results to file
write.csv(res_df, file = output_file, row.names = FALSE)