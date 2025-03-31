#!/bin/bash
#SBATCH --job-name=VCF_GT
#SBATCH --output=/ix1/wchen/Shiyue/Projects/2024_05_Allelic_Imbalance_Method/Code/Slurm_Out/%j.out
#SBATCH --error=/ix1/wchen/Shiyue/Projects/2024_05_Allelic_Imbalance_Method/Code/Slurm_Out/%j.err
#SBATCH --cluster=smp
#SBATCH --partition=smp
#SBATCH --cpus-per-task=2
#SBATCH --mem=70480
#SBATCH --time=24:00:00
#SBATCH --array=1-25%25

WORK_PATH=/ix1/wchen/Shiyue/Projects/2024_05_Allelic_Imbalance_Method/Result/03_SNP/
SAMPLE_NAME=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${WORK_PATH}Sample_test.txt)
bcftools=/ix1/wchen/Shiyue/Biosoft/bcftools/bcftools/bcftools
export BCFTOOLS_PLUGINS=/ix1/wchen/Shiyue/Biosoft/bcftools/bcftools/plugins

${bcftools} +tag2tag ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_filtered.vcf.gz -- -r --PL-to-GT | \
${bcftools} view -Oz -o ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_genotyped.vcf.gz


