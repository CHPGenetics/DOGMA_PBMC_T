#!/bin/bash
#SBATCH --job-name=SCARlink
#SBATCH --output=/ix1/wchen/Shiyue/Projects/2023_07_DOGMA_Revision/SCARlink/Code/Slurm_Out/%j.out
#SBATCH --error=/ix1/wchen/Shiyue/Projects/2023_07_DOGMA_Revision/SCARlink/Code/Slurm_Out/%j.err
#SBATCH --ntasks-per-node=1
#SBATCH --gres=gpu:1
#SBATCH --cluster=gpu
#SBATCH --partition=a100
#SBATCH --mem=112000
#SBATCH --time=2:00:00

# SCARlink=/ix1/wchen/Shiyue/Biosoft/SCARlink/
WORK_PATH=/ix1/wchen/Shiyue/Projects/2023_07_DOGMA_Revision/SCARlink/Result/


source activate /ix1/wchen/xiangyu/conda_env/scarlink/

module load gcc/12.2.0 r/4.4.0
# module load gcc/8.2.0 python/ondemand-jupyter-python3.8

export XLA_FLAGS=--xla_gpu_cuda_data_dir=/ihome/wchen/xiy231/.conda/envs/scarlink/lib/
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/ihome/wchen/xiy231/.conda/envs/scarlink/lib/
# export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/ihome/wchen/xiy231/.local/lib/python3.8/site-packages/tensorrt_libs/

# python ${SCARlink}scarlink/run_scarlink.py -o ${WORK_PATH}processing -g hg38
scarlink -o ${WORK_PATH}pbmc_all_out_10k -g hg38 -c celltype_updated \
--gene_list ${WORK_PATH}pbmc_all_out_10k/gene_subset.txt

