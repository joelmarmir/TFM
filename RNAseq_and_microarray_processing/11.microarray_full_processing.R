#####################################################
#### Microarray Differential Expression Pipeline ####
#####################################################

#----------------------------------------------------------
#---- Data pre-processing and ExpressionSet generation ----
#----------------------------------------------------------

# Load required libraries
library(GEOquery)
library(limma)
library(Biobase)
library(affy)
library(mouse4302.db)
library(dplyr)
library(ggplot)
library(ggrepel)

#----------------------------------------------------------
#------------------ Read and Normalize CEL Files ----------
#----------------------------------------------------------

# SET WD
setwd("/users/genomics/jmartinez/data/08_cells/Fong_2012/GSE34907_RAW/")

# Read CEL files and exclude the .CEL.gz part of the names
rawData <- ReadAffy(filenames = list.files(pattern = "*.CEL.gz"))

# Normalize data
normData <- affy::rma(rawData)

# Extract exprs
exprs_matrix <- exprs(normData)

# Extract the GSMs needed
selected_GSMs <- grep("^GSM857409|^GSM857410|^GSM857411|^GSM857412|^GSM857413|^GSM857414", colnames(exprs_matrix), value = TRUE)
exprs_matrix <- exprs_matrix[, selected_GSMs]

head(exprs_matrix)
tail(exprs_matrix)

#----------------------------------------------------------
#-------------- Assign Experimental Conditions ------------
#----------------------------------------------------------

# Add condition column to as pData
condition <- factor(c(
  "uninduced",
  "uninduced",
  "uninduced",
  "induced",
  "induced",
  "induced"),
  levels = c("uninduced", "induced"))
 
# Extract components
phenoData <- data.frame(row.names = colnames(exprs_matrix), condition = condition)

#----------------------------------------------------------
#------------------ ExpressionSet generation --------------
#----------------------------------------------------------

# Annotate genes with PROBE IDs
probe_ids <- rownames(exprs_matrix)

# Annotation with mouse4302.db
mapped_genes <- AnnotationDbi::select(
  mouse4302.db,
  keys = probe_ids,
  keytype = "PROBEID",
  columns = c("PROBEID", "SYMBOL")
)

# Remove duplicates
mapped_genes <- mapped_genes %>%
  group_by(PROBEID) %>%
  slice(1) %>%  # Keeps the first occurrence
  ungroup()

# Set PROBE ID as row names
feature_data <- as.data.frame(mapped_genes)
rownames(feature_data) <- feature_data$PROBEID
feature_data$PROBEID <- NULL 

# Rename columns
colnames(feature_data) <- c("GeneSymbol")
head(feature_data)

# Creation of eset
eset <- ExpressionSet(assayData = exprs_matrix,
                      phenoData = AnnotatedDataFrame(phenoData),
                      featureData = AnnotatedDataFrame(feature_data))



#----------------------------------------------------------
#--------- Collapse probes to gene symbols (avereps) -------
#----------------------------------------------------------

# Extract gene symbols from featureData
gene_symbols <- fData(eset)$GeneSymbol

# Filter probes without gene symbols
keep <- !is.na(gene_symbols) & gene_symbols != ""
exprs_filtered <- eset[keep, ]
gene_symbols <- gene_symbols[keep]

# Collapse probes using avereps
exprs_collapsed <- avereps(exprs(exprs_filtered), ID = gene_symbols)

# Create new featureData
new_fdata <- data.frame(
  GeneSymbol = rownames(exprs_collapsed),
  row.names = rownames(exprs_collapsed)
)

# Create averaged ExpressionSet
eset_avg <- ExpressionSet(
  assayData = exprs_collapsed,
  phenoData = AnnotatedDataFrame(phenoData),
  featureData = AnnotatedDataFrame(new_fdata)
)


#----------------------------------------------------------
#------------------ Calculate contrasts -------------------
#----------------------------------------------------------

# Design matrix
design <- model.matrix(~0 + pData(eset_avg)$condition)
colnames(design) <- levels(pData(eset_avg)$condition)

# Contrast matrix
contrast_matrix <- makeContrasts(
  ind_min_unind = induced - uninduced,
  levels = design
)

# Fit linear model
fit <- lmFit(eset_avg, design)
fit_main <- contrasts.fit(fit, contrast_matrix)
fit_main <- eBayes(fit_main)

# Get results
toptable <- topTable(fit_main, 
                     coef = 1, 
                     number = Inf,
                     sort.by = "logFC")


#----------------------------------------------------------
#--------------------- Data cleaning ----------------------
#----------------------------------------------------------

# Remove rownames
rownames(toptable) <- NULL

# Order results by logFC
toptable <- toptable[order(toptable$logFC, decreasing = TRUE), ]

# View number of occurrences
nrow(toptable)
head(toptable)

#----------------------------------------------------------
#--------------------- Save results -----------------------
#----------------------------------------------------------

write.csv(toptable,
          file = "/users/genomics/jmartinez/data/06_fc/Fong_2012/Fong_2012_Neurod2_mef_microarray_results1.csv",
          row.names = FALSE)

#----------------------------------------------------------
#---------------- Volcano plot (optional) -----------------
#----------------------------------------------------------

# Save the plot as a PDF
pdf("/users/genomics/jmartinez/volcano_plot7.pdf", width = 8, height = 6)

# Create the volcano plot
ggplot(toptable, aes(x = logFC, y = -log10(P.Value))) +
  geom_point(aes(color = adj.P.Val < 0.05)) +
  scale_color_manual(values = c("black", "red")) +
  geom_vline(xintercept = c(-2, 2), linetype = "dashed", color = "blue") +
  labs(title = "Volcano Plot", x = "Log Fold Change", y = "-Log10 P-Value") +
  theme_minimal() +
  geom_text_repel(data = subset(toptable, abs(logFC) > 2), 
                  aes(label = GeneSymbol), 
                  size = 3)

# Close the PDF device
dev.off()


