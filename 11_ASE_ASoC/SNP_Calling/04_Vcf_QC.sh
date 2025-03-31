#!/bin/bash
#SBATCH --job-name=VCF
#SBATCH --output=/ix1/wchen/Shiyue/Projects/2024_05_Allelic_Imbalance_Method/Code/Slurm_Out/%j.out
#SBATCH --error=/ix1/wchen/Shiyue/Projects/2024_05_Allelic_Imbalance_Method/Code/Slurm_Out/%j.err
#SBATCH --cluster=smp
#SBATCH --partition=smp
#SBATCH --cpus-per-task=2
#SBATCH --mem=70480
#SBATCH --time=24:00:00
#SBATCH --array=1-25%25

module unload gcc/8.2.0
module unload htslib/1.9
module unload bcftools/1.15.1

module load gcc/8.2.0
module load htslib/1.9
module load bcftools/1.15.1

WORK_PATH=/ix1/wchen/Shiyue/Projects/2024_05_Allelic_Imbalance_Method/Result/03_SNP/
SAMPLE_NAME=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${WORK_PATH}Sample_test.txt)

## Concat
for chr in `seq 1 22`;do
file_post_list+=" ${WORK_PATH}${SAMPLE_NAME}/germline/chr${chr}.gl.vcf.gz"
done

for chr in `seq 1 22`; do
    bcftools index -f ${WORK_PATH}${SAMPLE_NAME}/germline/chr${chr}.gl.vcf.gz
done

bcftools concat -Oz -o ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_raw.vcf.gz ${file_post_list}

gunzip ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_raw.vcf.gz
awk 'BEGIN {FS="\t"; OFS="\t"} {if ($1 ~ /^#/ || NF == 11) print $0}' ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_raw.vcf > ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_fixed.vcf
bgzip ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_fixed.vcf
## Basic Statistics
bcftools stats ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_fixed.vcf.gz > ${WORK_PATH}${SAMPLE_NAME}/germline/vcf_stats_raw.txt

# QC
## Filtering Low Quality Variants
bcftools view -e 'INFO/DP<10' ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_fixed.vcf.gz \
-Oz -o ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_filtered.vcf.gz

rm ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_raw.vcf
rm ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_fixed.vcf.gz

# 
# bcftools view -v snps ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_RNA_10.vcf.gz | grep -vc "^#"
# bcftools view ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_RNA.vcf.gz | tail -n 20
# bcftools view ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_RNA_30.vcf.gz | tail -n 5
# bcftools view ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_ATAC.vcf.gz | tail -n 20
# bcftools view ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_ATAC_30.vcf.gz | tail -n 5
# bcftools view ${WORK_PATH}${SAMPLE_NAME}/germline/chr20.gl.vcf.gz | tail -n 5


