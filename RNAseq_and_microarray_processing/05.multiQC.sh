#!/bin/bash

#############################################################
###### Summarize quality checks after fastp processing ######
#############################################################

# -----------------
# -- Use MultiQC --
# -----------------

# Load Modules
module load Miniconda3/20240927
multiqc --version

# Run multiqc to gather qc metrics
multiqc /users/genomics/jmartinez/data/02_preprocessed_fastqs/00_temp/ -o /users/genomics/jmartinez/data/02_preprocessed_fastqs/00_temp/z_multiqc_report