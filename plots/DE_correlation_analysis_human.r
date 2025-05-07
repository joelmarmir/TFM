######################################################
#### Correlation analysis of gene expression data ####
######################################################

# Input: csv
# Output: pdf

# Load libraries
library(tidyverse)
library(pheatmap)
library(ggplot2)
library(dplyr)
library(dendextend)

# Load data
df <- read.csv("/users/genomics/jmartinez/data/07_for_results/fc_all_human_with_micro3_github.csv", header = TRUE)

# Select relevant columns
df_selected <- df %>% dplyr::select(File, GeneSymbol, log2FoldChange)

# Clean study names by keeping only the part before the 4th underscore (if present)
df_selected <- df_selected %>%
  mutate(File = sub("^(([^_]+_){3}[^_]+).*", "\\1", File))


# Pivot the data: Rows = Genes, Columns = Studies (Files)
# I added "values_fn = mean", because some genes have multiple entries in the same study
# This is likely the case for RNAseq data, where multiple transcripts are reported for the same gene
# For microarray data I already calculated the mean for each gene
df_wide <- df_selected %>%
  pivot_wider(names_from = File, values_from = log2FoldChange, values_fn = mean)

head(df_wide)
nrow(df_wide)

# Remove rows with "NA" in the column GeneSymbol
df_wide <- df_wide %>% filter(!is.na(GeneSymbol))


# Check for genes without any NA values
nrow(df_wide)

# Convert to matrix (removing ENSEMBL column)
expr_matrix <- df_wide %>%
  column_to_rownames(var = "GeneSymbol") %>%
  as.matrix()

# Keep only the last timepoint from each study
studies_not_in_list <- c("Matsuda_2016_Neurod1_microglia_Day2", "Aydin_2019_Neurog2_eb_12h", "Aydin_2019_Ascl1_eb_12h",
                         "Joung_2023_Ascl1_hesc_day7", "Barfeld_2017_Myc_lncap_5h", "Chalamasetty_2014_Msgn1_eb_12h",
                         "Chalamasetty_2014_Msgn1_eb_24h")

# Remove early timepoints samples from the expression matrices
expr_matrix <- expr_matrix[, !sapply(colnames(expr_matrix), function(col) any(grepl(paste0("^", studies_not_in_list, collapse = "|"), col)))]
colnames(expr_matrix)

# NEW: Z-score normalization within each study (column)
#expr_matrix_scaled <- scale(expr_matrix)  # Normalize raw matrix
#expr_matrix_clean_scaled <- scale(expr_matrix_clean)  # Normalize cleaned matrix

# Compute correlation matrix across studies
cor_matrix <- cor(expr_matrix, use = "pairwise.complete.obs", method = "pearson")

# Highlight studies without replicates
studies_no_replicates <- c("Conerly_2016", "Fong_2015", "Li_2019", "Liu_2023", "Walz_2014", "Weber_2016", "Yu_2023")

# Add asterisks to studies without replicates that match the beginning of colnames
colnames(cor_matrix) <- ifelse(sapply(colnames(cor_matrix), function(col) any(grepl(paste0("^", studies_no_replicates, collapse = "|"), col))),
                paste0(colnames(cor_matrix), "*"),
                colnames(cor_matrix))

rownames(cor_matrix) <- ifelse(sapply(rownames(cor_matrix), function(row) any(grepl(paste0("^", studies_no_replicates, collapse = "|"), row))),
                paste0(rownames(cor_matrix), "*"),
                rownames(cor_matrix))

#########################################
########## Extract cell types ###########
#########################################

# Extract cell types from study names (between the 3rd _ and 4th _), ignoring trailing '*'
extract_cell_type <- function(name) {
  parts <- unlist(strsplit(name, "_"))
  if (length(parts) >= 4) {
    cell_type <- parts[4]
    return(gsub("\\*", "", cell_type))  # Remove any trailing "*"
  } else {
    return(NA)
  }
}


# Create a data frame for column annotations
cell_types <- sapply(colnames(cor_matrix), extract_cell_type)
annotation_col <- data.frame(CellType = cell_types)
rownames(annotation_col) <- colnames(cor_matrix)

# Define colors for cell types
unique_cell_types <- unique(cell_types)
cell_type_colors <- setNames(colorRampPalette(c("blue", "green", "red", "purple", "orange"))(length(unique_cell_types)), unique_cell_types)
annotation_colors <- list(CellType = cell_type_colors)


################################################
############# Extract bhlh clades ##############
################################################

# Load bHLH clades
bhlh_clades <- readLines("/users/genomics/jmartinez/temp_data/bhlh_class.txt")

# Parse clades into a named list
clades <- list()
current_clade <- NULL
for (line in bhlh_clades) {
  if (startsWith(line, ">")) {
    current_clade <- sub(">", "", line)  # Extract clade name and convert to lowercase
    clades[[current_clade]] <- c()
  } else if (nzchar(line)) {
    clades[[current_clade]] <- c(clades[[current_clade]], tolower(line))  # Add factors to the clade
  }
}

# Extract transcription factors from column names (between 3rd _ and 4th _)
extract_tf <- function(name) {
  parts <- unlist(strsplit(name, "_"))
  if (length(parts) >= 4) {
    return(tolower(parts[3]))  # Extract the 3rd part and convert to lowercase
  } else {
    return(NA)  # Return NA if the format is unexpected
  }
}

# Classify transcription factors into clades
tfs <- sapply(colnames(cor_matrix), extract_tf)
tf_clades <- sapply(tfs, function(tf) {
  clade <- names(clades)[sapply(clades, function(factors) tf %in% factors)]
  if (length(clade) > 0) {
    return(clade[1])  # Return the first matching clade
  } else {
    return("unknown")  # Default to "unknown" if no match is found
  }
})

# Create a data frame for column annotations
annotation_col <- data.frame(
  CellType = cell_types,  # From previous code
  Clade = tf_clades
)
rownames(annotation_col) <- colnames(cor_matrix)

# Define colors for clades
unique_clades <- unique(tf_clades)
clade_colors <- setNames(colorRampPalette(c("yellow", "cyan", "magenta", "brown", "orange", "gray"))(length(unique_clades)), unique_clades)
annotation_colors$Clade <- clade_colors

##################################################
############ Add bhlh specific factor ############
##################################################

# Convert bhlh type list to lowercase
specific_bhlh_type <- tolower(c("NEUROD1", "NEUROD2", "NEUROG2", "ATOH1", "TWIST2", "TWIST1", "OLIG2", "MSGN1", "MESP1", 
                                 "MYOD1", "MYF5", "ASCL1", "ASCL2", "TCF4", "SIM2", "MYC", "TFEB", "MYCN", "MITF", "HEY1", "HEY2", "TFE3"))
specific_bhlh_colors <- c("dodgerblue1", "dodgerblue3", "dodgerblue4", "steelblue4", "olivedrab1", "olivedrab3", "turquoise", 
                          "coral", "indianred", "indianred1", "tomato3", "purple1", "purple3", "darkorange1", 
                          "orange2", "gold1", "gold3", "khaki1", "lightgoldenrod2", "plum1", "plum3", "springgreen2")

# Map the specific bhlh types to the column names (case-insensitive match)
specific_bhlh_type_mapping <- sapply(tolower(colnames(cor_matrix)), function(name) {
  match <- specific_bhlh_type[sapply(specific_bhlh_type, function(cell) grepl(cell, name))]
  if (length(match) > 0) {
    return(match[1])  # Return the first match
  } else {
    return("unknown")  # Default to "unknown" if no match is found
  }
})

# Add the specific bhlh type to the annotation_col data frame
annotation_col$SpecificBHLH <- specific_bhlh_type_mapping

# Define colors for the specific cell types
specific_bhlh_type_colors <- setNames(specific_bhlh_colors, specific_bhlh_type)
specific_bhlh_type_colors["unknown"] <- "gray"  # Add a default color for "unknown"

# Filter specific bHLH types to include only those present in the heatmap
present_bhlh_types <- unique(annotation_col$SpecificBHLH)
filtered_bhlh_type_colors <- specific_bhlh_type_colors[names(specific_bhlh_type_colors) %in% present_bhlh_types]

# Add the filtered bHLH type colors to the annotation_colors list
annotation_colors$SpecificBHLH <- filtered_bhlh_type_colors

# Reorder the columns in annotation_col
annotation_col <- annotation_col[, c("SpecificBHLH", "CellType", "Clade")]

############### Create a dendogram  ######################

# Convert correlation matrix to a distance matrix: I use 1 - cor_matrix, so that the values
# with a higher correlation, are closer to 0, and cluster together
dist_mat <- as.dist(1 - cor_matrix)

# Hierarchical clustering
hc <- hclust(dist_mat, method = "ward.D2")

# Convert to dendrogram
dend <- as.dendrogram(hc)

# Check dend
pdf("/users/genomics/jmartinez/dend_plot2.pdf")
plot(dend, main = "Column Dendrogram", horiz = FALSE)
dev.off()

# Use the same dendrogram for rows and columns (since cor_matrix is symmetric)
pdf("/users/genomics/jmartinez/heatmap71_v5_all_pearson_human.pdf", width = 10, height = 10)
pheatmap(
  cor_matrix,
  cluster_cols = hc,
  cluster_rows = hc,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  breaks = seq(-1, 1, length.out = 101),
  main = "Study Similarity Heatmap (Rotated)",
  annotation_col = annotation_col,
  annotation_colors = annotation_colors,
  fontsize_row = 10,
  fontsize_col = 10,
  border_color = NA
)
dev.off()
