###########################################
##### Extract SRA codes from GSM codes ####
###########################################

# Load libraries
library(rentrez)

# Set input file path
input_path <- "/users/genomics/jmartinez/temp_data/expression_datasets.csv"
# Set output file path
output_path <- "/users/genomics/jmartinez/temp_data/expression_datasets_csv_with_sra.csv"

#------------------------
#---- Load the data -----
#------------------------

# Load the data
data <- read.csv(input_path,
                  header = TRUE,
                  sep = ";",
                  na.strings = "")
as.data.frame(data)
head(data)

# Rename column one as due to formatting issues it contains " X..."
colnames(data)[1] <- "Induced_TF"
head(data)

# Check for missing values (it is normal to have some NA values, due to how data is structured)
sum(is.na(data))

#------------------------------
#---- GSM to SRA function -----
#------------------------------

# Function to retrieve SRA information from a GSM code
get_sra_from_gsm <- function(gsm_code) {
    # Return NA if the GSM code is missing
    if (is.na(gsm_code)) {
        return(NA)
    }
    # Search for the GSM code in the SRA database
    search_result <- entrez_search(db = "sra", term = gsm_code)
  
    # If there are results, fetch the runinfo
    if (search_result$count > 0) {
        sra_id <- search_result$ids[1]
        runinfo <- entrez_fetch(db = "sra", id = sra_id, rettype = "runinfo", retmode = "text")
    
        # Extract the SRR code from the runinfo
        runinfo_lines <- strsplit(runinfo, "\n")[[1]]
        srr_code <- substr(runinfo_lines[2:length(runinfo_lines)], 1, 11)
        srr_code <- sub(",.*$", "", srr_code)
        # If more than one SRR code is detected, concatenate them with underscores
        if (length(srr_code) > 1) {
            srr_code <- paste(srr_code, collapse = "_")}
        return(srr_code)
  } else {
        return(NA)
  }
}


# Example usage
gsm_code <- "GSM2232829"
srr_code <- get_sra_from_gsm(gsm_code)
print(srr_code)

#---------------------------------
#--- Apply the function to csv ---
#---------------------------------

# Apply the function to the data (take in mind that it will take a some time)
data$Uninduced_SRA <- sapply(data$Uninduced, get_sra_from_gsm)
data$Induced_SRA <- sapply(data$Induced, get_sra_from_gsm)
head(data)

# Save the data
write.csv(data, output_path, row.names = FALSE)
