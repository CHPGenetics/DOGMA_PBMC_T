#!/bin/bash
#SBATCH --job-name=SCARlink
#SBATCH --output=/ix1/rduerr/shared/rduerr_wchen/Shiyue/Chen_Files/2023_07_DOGMA_Revision/SCARlink/Code/Slurm_Out/%j.out
#SBATCH --error=/ix1/rduerr/shared/rduerr_wchen/Shiyue/Chen_Files/2023_07_DOGMA_Revision/SCARlink/Code/Slurm_Out/%j.err
#SBATCH --cpus-per-task=64
#SBATCH --mem=768000
#SBATCH --time=20:00:00

source activate /ihome/wchen/sht175/.conda/envs/scarlink-env

SCARlink=/ix1/wchen/Shiyue/Biosoft/SCARlink/
WORK_PATH=/ix1/rduerr/shared/rduerr_wchen/Shiyue/Chen_Files/2023_07_DOGMA_Revision/SCARlink/Result/

export XLA_FLAGS=--xla_gpu_cuda_data_dir=/ihome/wchen/sht175/.conda/envs/scarlink-env/lib/
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/ihome/wchen/sht175/.conda/envs/scarlink-env/lib/
export TF_CPP_MIN_LOG_LEVEL=2

python ${SCARlink}scarlink/run_scarlink_tiles.py \
-o ${WORK_PATH}All -c celltype_updated

