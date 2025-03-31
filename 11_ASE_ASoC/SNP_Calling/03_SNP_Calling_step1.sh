#!/bin/bash
#SBATCH --job-name=SNPcall
#SBATCH --output=/ix1/wchen/Shiyue/Projects/2024_05_Allelic_Imbalance_Method/Code/Slurm_Out/%j.out
#SBATCH --error=/ix1/wchen/Shiyue/Projects/2024_05_Allelic_Imbalance_Method/Code/Slurm_Out/%j.err
#SBATCH --cpus-per-task=3
#SBATCH --mem=92480
#SBATCH --time=200:00:00
#SBATCH --array=1-25%25

Monopogen=/ix1/wchen/Shiyue/Biosoft/Monopogen/
REF_PATH=/ix1/wchen/Shiyue/References/Fasta/hg38/hg38.fa
KGP_PATH=/ix1/rduerr/shared/rduerr_wchen/Shiyue/Reference/KGP_Phased/
WORK_PATH=/ix1/wchen/Shiyue/Projects/2024_05_Allelic_Imbalance_Method/Result/03_SNP/

SAMPLE_NAME=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${WORK_PATH}Sample.txt)

module load python/ondemand-jupyter-python3.10
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${Monopogen}/apps

# Germline SNV calling
## Prepare
python  ${Monopogen}src/Monopogen.py preProcess \
-a ${Monopogen}apps \
-b ${WORK_PATH}${SAMPLE_NAME}/bam.lst \
-o ${WORK_PATH}${SAMPLE_NAME}

