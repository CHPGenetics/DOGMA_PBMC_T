#!/bin/bash
#SBATCH --job-name=Genotype
#SBATCH --output=/ix1/wchen/Shiyue/Projects/2024_05_Allelic_Imbalance_Method/Code_New/Slurm_Out/%j.out
#SBATCH --error=/ix1/wchen/Shiyue/Projects/2024_05_Allelic_Imbalance_Method/Code_New/Slurm_Out/%j.err
#SBATCH --cluster=smp
#SBATCH --partition=smp
#SBATCH --cpus-per-task=2
#SBATCH --mem=70480
#SBATCH --time=20:00:00
#SBATCH --array=1-25%25

module load gcc/8.2.0
module load htslib/1.9
module load bcftools/1.15.1
module load plink/1.90b6.7

WORK_PATH=/ix1/wchen/Shiyue/Projects/2024_05_Allelic_Imbalance_Method/Result/03_SNP/
SAMPLE_NAME=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${WORK_PATH}Sample_test.txt)

# Rename vcf file
# Combine RNA and ATAC (Merged)
for depth in 10 20 30 40; do

## Step 1: Identify conflicting variants
bcftools index ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_ATAC_${depth}.vcf.gz
bcftools index ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_RNA_${depth}.vcf.gz

bcftools isec -p ${WORK_PATH}${SAMPLE_NAME}/germline/isec_output \
-Oz ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_ATAC_${depth}.vcf.gz \
${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_RNA_${depth}.vcf.gz

bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t[%GT]\n' \
${WORK_PATH}${SAMPLE_NAME}/germline/isec_output/0000.vcf.gz\
> ${WORK_PATH}${SAMPLE_NAME}/germline/isec_output/0000_GT.txt

bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t[%GT]\n' \
${WORK_PATH}${SAMPLE_NAME}/germline/isec_output/0001.vcf.gz\
> ${WORK_PATH}${SAMPLE_NAME}/germline/isec_output/0001_GT.txt

awk 'NR==FNR {a[$1"\t"$2]=$0; next} {if ($1"\t"$2 in a && $6 != a[$1"\t"$2]) print $1"\t"$2}' \
${WORK_PATH}${SAMPLE_NAME}/germline/isec_output/0000_GT.txt \
${WORK_PATH}${SAMPLE_NAME}/germline/isec_output/0001_GT.txt\
> ${WORK_PATH}${SAMPLE_NAME}/germline/conflicting_variants.txt

if [ -s "${WORK_PATH}${SAMPLE_NAME}/germline/conflicting_variants.txt" ]; then
   
## Step 2: Remove conflicting variants
bcftools view -T ^${WORK_PATH}${SAMPLE_NAME}/germline/conflicting_variants.txt \
-o ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_ATAC_filtered_${depth}.vcf.gz \
-Oz ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_ATAC_${depth}.vcf.gz

bcftools view -T ^${WORK_PATH}${SAMPLE_NAME}/germline/conflicting_variants.txt \
-o ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_RNA_filtered_${depth}.vcf.gz \
-Oz ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_RNA_${depth}.vcf.gz

else

cp ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_RNA_${depth}.vcf.gz \
${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_RNA_filtered_${depth}.vcf.gz

cp ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_ATAC_${depth}.vcf.gz \
${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_ATAC_filtered_${depth}.vcf.gz

fi

## Step 3: Merge filtered files
bcftools index ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_ATAC_filtered_${depth}.vcf.gz
bcftools index ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_RNA_filtered_${depth}.vcf.gz

bcftools concat -a -o ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_merged_${depth}.vcf.gz \
-Oz ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_RNA_filtered_${depth}.vcf.gz \
${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_ATAC_filtered_${depth}.vcf.gz

bcftools index ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_merged_${depth}.vcf.gz

# rm ${WORK_PATH}${SAMPLE_NAME}/germline/conflicting_variants.txt
# rm -r ${WORK_PATH}${SAMPLE_NAME}/germline/isec_output
rm ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_ATAC_filtered_${depth}.vcf.gz
rm ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_RNA_filtered_${depth}.vcf.gz
rm ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_ATAC_filtered_${depth}.vcf.gz.csi
rm ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_RNA_filtered_${depth}.vcf.gz.csi

done

