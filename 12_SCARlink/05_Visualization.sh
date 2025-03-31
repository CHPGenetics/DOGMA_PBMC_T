#!/bin/bash
#SBATCH --job-name=SCARlink
#SBATCH --output=/ix1/wchen/Shiyue/Projects/2023_07_DOGMA_Revision/SCARlink/Code/Slurm_Out/%j.out
#SBATCH --error=/ix1/wchen/Shiyue/Projects/2023_07_DOGMA_Revision/SCARlink/Code/Slurm_Out/%j.err
#SBATCH --cluster=smp
#SBATCH --partition=smp
#SBATCH --cpus-per-task=10
#SBATCH --mem=68000
#SBATCH --time=1:00:00

conda activate scarlink-env

module load gcc/8.2.0
module load python/ondemand-jupyter-python3.8
module load r/4.0.0

SCARlink=/ix1/wchen/Shiyue/Biosoft/SCARlink/
WORK_PATH=/ix1/wchen/Shiyue/Projects/2023_07_DOGMA_Revision/SCARlink/Result/

python ${SCARlink}scarlink/run_scarlink_visualization.py \
-o ${WORK_PATH}processing -c celltype_updated --genes CCL20





