#!/bin/bash
#SBATCH --job-name=SCARlink
#SBATCH --output=/ix1/rduerr/shared/rduerr_wchen/Shiyue/Chen_Files/2023_07_DOGMA_Revision/SCARlink/Code/Slurm_Out/%j.out
#SBATCH --error=/ix1/rduerr/shared/rduerr_wchen/Shiyue/Chen_Files/2023_07_DOGMA_Revision/SCARlink/Code/Slurm_Out/%j.err
#SBATCH --cpus-per-task=16
#SBATCH --cluster=smp
#SBATCH --partition=smp
#SBATCH --mem=512000
#SBATCH --time=2:00:00

SCARlink=/ix1/wchen/Shiyue/Biosoft/SCARlink/
WORK_PATH=/ix1/rduerr/shared/rduerr_wchen/Shiyue/Chen_Files/2023_07_DOGMA_Revision/SCARlink/Result/

source activate /ihome/wchen/sht175/.conda/envs/scarlink-env

python ${SCARlink}scarlink/preprocessing/create_h5_files.py \
--scrna ${WORK_PATH}New_50k/dogma_scrna_subset.rds \
--scatac ${WORK_PATH}New_50k/DOGMA_ATAC_Subset -o ${WORK_PATH}New_50k -nc 16
  