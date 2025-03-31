#!/bin/bash
#SBATCH --job-name=Genotype
#SBATCH --output=/ix1/wchen/Shiyue/Projects/2024_05_Allelic_Imbalance_Method/Code_New/Slurm_Out/%j.out
#SBATCH --error=/ix1/wchen/Shiyue/Projects/2024_05_Allelic_Imbalance_Method/Code_New/Slurm_Out/%j.err
#SBATCH --cluster=smp
#SBATCH --partition=smp
#SBATCH --cpus-per-task=10
#SBATCH --mem=102480
#SBATCH --time=20:00:00

module unload gcc/8.2.0
module unload htslib/1.9
module unload bcftools/1.15.1
module unload plink/1.90b6.7

module load gcc/8.2.0
module load htslib/1.9
module load bcftools/1.15.1
module load plink/1.90b6.7

WORK_PATH=/ix1/wchen/Shiyue/Projects/2024_05_Allelic_Imbalance_Method/Result/03_SNP/
SAMPLES_FILE=${WORK_PATH}Sample_test.txt

for depth in 10 20 30 40; do

for omics in RNA ATAC merged; do

# Combine vcf of each sample (RNA/ATAC/Merge)
## Merge VCF files
VCF_LIST=${WORK_PATH}Combine/${omics}_${depth}_vcf.txt
rm ${VCF_LIST}
while read SAMPLE; do
    echo "${WORK_PATH}${SAMPLE}/germline/${SAMPLE}_${omics}_${depth}.vcf.gz" >> ${VCF_LIST}
done < ${SAMPLES_FILE}

bcftools merge -l ${VCF_LIST} \
-Oz -o ${WORK_PATH}Combine/Merged_${omics}_${depth}.vcf.gz

# Convert to Plink format (filter the duplicated records)
gunzip -c ${WORK_PATH}Combine/Merged_${omics}_${depth}.vcf.gz \
> ${WORK_PATH}Combine/Merged_${omics}_${depth}.vcf

plink --vcf ${WORK_PATH}Combine/Merged_${omics}_${depth}.vcf --maf 0.001 \
--make-bed --out ${WORK_PATH}Combine/Merged_${omics}_${depth}

rm ${WORK_PATH}Combine/Merged_${omics}_${depth}.vcf

done

done



