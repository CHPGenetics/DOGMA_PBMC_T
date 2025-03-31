#!/bin/bash
#SBATCH --job-name=SNPcall
#SBATCH --output=/ix1/wchen/Shiyue/Projects/2024_05_Allelic_Imbalance_Method/Code/Slurm_Out/%j.out
#SBATCH --error=/ix1/wchen/Shiyue/Projects/2024_05_Allelic_Imbalance_Method/Code/Slurm_Out/%j.err
#SBATCH --cpus-per-task=3
#SBATCH --mem=10480
#SBATCH --time=300:00:00
#SBATCH --array=1-264%264

Monopogen=/ix1/wchen/Shiyue/Biosoft/Monopogen/
REF_PATH=/ix1/wchen/Shiyue/References/Fasta/hg38/hg38.fa
KGP_PATH=/ix1/rduerr/shared/rduerr_wchen/Shiyue/Reference/KGP_Phased/
WORK_PATH=/ix1/wchen/Shiyue/Projects/2024_05_Allelic_Imbalance_Method/Result/03_SNP/

SAMPLE_INDEX=$(( (${SLURM_ARRAY_TASK_ID} - 1) / 22 + 1 ))
CHR_ID=$(( (${SLURM_ARRAY_TASK_ID} - 1) % 22 + 1 ))
SAMPLE_NAME=$(sed -n "${SAMPLE_INDEX}p" ${WORK_PATH}Sample_1.txt)

module load python/ondemand-jupyter-python3.10
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${Monopogen}/apps

# Germline SNV calling
cd ${WORK_PATH}${SAMPLE_NAME}

python ${Monopogen}src/Monopogen.py germline \
-a ${Monopogen}apps -r ${Monopogen}Region/region_chr${CHR_ID}.lst \
-s varScan -p ${KGP_PATH} \
-g ${REF_PATH} -o ${WORK_PATH}${SAMPLE_NAME}
