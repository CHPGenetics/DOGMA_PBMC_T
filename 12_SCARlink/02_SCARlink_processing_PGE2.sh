#!/bin/bash
#SBATCH --job-name=SCARlink
#SBATCH --output=/ix1/wchen/Shiyue/Projects/2023_07_DOGMA_Revision/SCARlink/Code/Slurm_Out/%j.out
#SBATCH --error=/ix1/wchen/Shiyue/Projects/2023_07_DOGMA_Revision/SCARlink/Code/Slurm_Out/%j.err
#SBATCH --cpus-per-task=16
#SBATCH --cluster=smp
#SBATCH --partition=smp
#SBATCH --mem=512000
#SBATCH --time=2:00:00

SCARlink=/ix1/wchen/Shiyue/Biosoft/SCARlink/
WORK_PATH=/ix1/wchen/Shiyue/Projects/2023_07_DOGMA_Revision/SCARlink/Result/

source activate /ihome/wchen/sht175/.conda/envs/scarlink-env

python ${SCARlink}scarlink/preprocessing/create_h5_files.py \
--scrna ${WORK_PATH}Act_IL1B_IL23_PGE2/dogma_scrna_subset.rds \
--scatac ${WORK_PATH}Act_IL1B_IL23_PGE2/DOGMA_ATAC_Subset -o ${WORK_PATH}Act_IL1B_IL23_PGE2 -nc 16
