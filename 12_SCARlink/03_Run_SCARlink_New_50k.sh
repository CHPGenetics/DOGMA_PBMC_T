#!/bin/bash
#SBATCH --job-name=SCARlink
#SBATCH --output=/ix1/rduerr/shared/rduerr_wchen/Shiyue/Chen_Files/2023_07_DOGMA_Revision/SCARlink/Code/Slurm_Out/%j.out
#SBATCH --error=/ix1/rduerr/shared/rduerr_wchen/Shiyue/Chen_Files/2023_07_DOGMA_Revision/SCARlink/Code/Slurm_Out/%j.err
#SBATCH --array=1-100
#SBATCH --cpus-per-task=2
#SBATCH --cluster=smp
#SBATCH --partition=smp
#SBATCH --mem=64000
#SBATCH --time=200:00:00

SCARlink=/ix1/wchen/Shiyue/Biosoft/SCARlink/
WORK_PATH=/ix1/rduerr/shared/rduerr_wchen/Shiyue/Chen_Files/2023_07_DOGMA_Revision/SCARlink//Result/

source activate /ihome/wchen/sht175/.conda/envs/scarlink-env

export XLA_FLAGS=--xla_gpu_cuda_data_dir=/ihome/wchen/sht175/.conda/envs/scarlink-env/lib/
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/ihome/wchen/sht175/.conda/envs/scarlink-env/lib/
export TF_CPP_MIN_LOG_LEVEL=2

python ${SCARlink}scarlink/run_scarlink.py -o ${WORK_PATH}New_50k \
-g hg38 -c celltype_updated \
-p $SLURM_ARRAY_TASK_ID -n 100


