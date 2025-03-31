#!/bin/bash
#SBATCH --job-name=SCARlink
#SBATCH --output=/ix1/wchen/Shiyue/Projects/2023_07_DOGMA_Revision/SCARlink/Code/Slurm_Out/%j.out
#SBATCH --error=/ix1/wchen/Shiyue/Projects/2023_07_DOGMA_Revision/SCARlink/Code/Slurm_Out/%j.err
#SBATCH --cpus-per-task=16
#SBATCH --mem=512000
#SBATCH --time=2:00:00

module load gcc/8.2.0
module load python/ondemand-jupyter-python3.8

SCARlink=/ix1/wchen/Shiyue/Biosoft/SCARlink/
WORK_PATH=/ix1/wchen/Shiyue/Projects/2023_07_DOGMA_Revision/SCARlink/Result/

conda activate scarlink-env

module load gcc/12.2.0 r/4.3.0

python ${SCARlink}scarlink/preprocessing/create_h5_files.py --scrna ${WORK_PATH}dogma_scrna_subset.rds \
--scatac ${WORK_PATH}DOGMA_ATAC_Subset -o ${WORK_PATH}processing -nc 16


