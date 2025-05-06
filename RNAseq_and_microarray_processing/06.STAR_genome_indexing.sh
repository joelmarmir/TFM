#!/bin/bash

#################################################
##### STAR: genome indexing (mouse + human) #####
#################################################

# -----------------------------------------
# -------- Setup the job in SLURM ---------
# -----------------------------------------

#SBATCH --job-name=genome_indexing
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=250G
#SBATCH --partition=long
#SBATCH --nodelist=node07
#SBATCH --output=/users/genomics/jmartinez/genome_indexing_log.out

# ------------------------------------------------------------
# -- Run STAR on the latest human and mouse genome versions --
# ------------------------------------------------------------

# Module load
module load STAR/2.7.8a-GCC-10.2.0

# Set timer
start_time=$(date +%s)

# STAR genome indexing (HUMAN)
STAR \
  --runThreadN 16 \
  --runMode genomeGenerate \
  --genomeDir /users/genomics/jmartinez/data/00_reference_genomes/human/indices \
  --genomeFastaFiles /users/genomics/jmartinez/data/00_reference_genomes/human/raw/GRCh38.primary_assembly.genome.fa \
  --sjdbGTFfile /users/genomics/jmartinez/data/00_reference_genomes/human/raw/gencode.v47.primary_assembly.annotation.gtf \
  --sjdbOverhang 100

STAR genome indexing (MOUSE)
STAR \
  --runThreadN 16 \
  --runMode genomeGenerate \
  --genomeDir /users/genomics/jmartinez/data/00_reference_genomes/mouse/indices \
  --genomeFastaFiles /users/genomics/jmartinez/data/00_reference_genomes/mouse/raw/GRCm39.primary_assembly.genome.fa \
  --sjdbGTFfile /users/genomics/jmartinez/data/00_reference_genomes/mouse/raw/gencode.vM36.primary_assembly.annotation.gtf \
  --sjdbOverhang 100

# Stop timer
echo "Finished at: $(date)"
end_time=$(date +%s)
elapsed_time=$((end_time - start_time))
echo "Total time taken: $elapsed_time seconds"