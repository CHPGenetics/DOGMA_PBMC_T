#!/bin/bash
#SBATCH --job-name=ArchR
#SBATCH --output=/ix1/wchen/Shiyue/Projects/2023_07_DOGMA_Revision/SCARlink/Code/Slurm_Out/%j.out
#SBATCH --error=/ix1/wchen/Shiyue/Projects/2023_07_DOGMA_Revision/SCARlink/Code/Slurm_Out/%j.err
#SBATCH --cluster=smp
#SBATCH --partition=smp
#SBATCH --cpus-per-task=64
#SBATCH --mem=768000
#SBATCH --time=20:00:00

module load gcc/8.2.0 r/4.0.0
# Rscript /ix1/wchen/Shiyue/Projects/2023_07_DOGMA_Revision/SCARlink/Code/01_00_ArchR_Proj.R
Rscript /ix1/wchen/Shiyue/Projects/2023_07_DOGMA_Revision/SCARlink/Code/01_PreProcessing.R
# Rscript /ix1/wchen/Shiyue/Projects/2023_07_DOGMA_Revision/SCARlink/Code/X_test.R

