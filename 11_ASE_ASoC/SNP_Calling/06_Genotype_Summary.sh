#!/bin/bash
#SBATCH --job-name=Genotype
#SBATCH --output=/ix1/wchen/Shiyue/Projects/2024_05_Allelic_Imbalance_Method/Code/Slurm_Out/%j.out
#SBATCH --error=/ix1/wchen/Shiyue/Projects/2024_05_Allelic_Imbalance_Method/Code/Slurm_Out/%j.err
#SBATCH --cluster=smp
#SBATCH --partition=smp
#SBATCH --cpus-per-task=2
#SBATCH --mem=70480
#SBATCH --time=20:00:00
#SBATCH --array=1-25%25

module unload gcc/8.2.0
module unload htslib/1.9
module unload bcftools/1.15.1
module unload plink/1.90b6.7

module load gcc/8.2.0
module load htslib/1.9
module load bcftools/1.15.1
module load plink/1.90b6.7

WORK_PATH=/ix1/wchen/Shiyue/Projects/2024_05_Allelic_Imbalance_Method/Result/03_SNP/
SAMPLE_NAME=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${WORK_PATH}Sample_test.txt)

# Seperate Genotypes from ATAC and RNA
bcftools view -s ${SAMPLE_NAME}_ATAC -Oz -o ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_ATAC.vcf.gz \
${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_genotyped.vcf.gz
bcftools reheader -s <(echo "${SAMPLE_NAME}_ATAC ${SAMPLE_NAME}" | tr " " "\t") \
-o ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_ATAC_rename.vcf.gz \
${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_ATAC.vcf.gz

## Filter ./.
bcftools view -i 'GT!="./."' ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_ATAC_rename.vcf.gz \
-o ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_ATAC_filtered.vcf.gz
rm ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_ATAC.vcf.gz
rm ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_ATAC_rename.vcf.gz

bcftools view -s ${SAMPLE_NAME}_GEX -Oz -o ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_RNA.vcf.gz \
${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_genotyped.vcf.gz
bcftools reheader -s <(echo "${SAMPLE_NAME}_GEX ${SAMPLE_NAME}" | tr " " "\t") \
-o ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_RNA_rename.vcf.gz \
${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_RNA.vcf.gz

## Filter ./.
bcftools view -i 'GT!="./."' ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_RNA_rename.vcf.gz \
-o ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_RNA_filtered.vcf.gz
rm ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_RNA.vcf.gz
rm ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_RNA_rename.vcf.gz

# Filtering Low Quality Variants
## ATAC
bcftools view -i 'FORMAT/DP>=10' ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_ATAC_filtered.vcf.gz \
-Oz -o ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_ATAC_10.vcf.gz

bcftools view -i 'FORMAT/DP>=20' ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_ATAC_filtered.vcf.gz \
-Oz -o ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_ATAC_20.vcf.gz

bcftools view -i 'FORMAT/DP>=30' ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_ATAC_filtered.vcf.gz \
-Oz -o ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_ATAC_30.vcf.gz

bcftools view -i 'FORMAT/DP>=40' ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_ATAC_filtered.vcf.gz \
-Oz -o ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_ATAC_40.vcf.gz

## RNA
bcftools view -i 'FORMAT/DP>=10' ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_RNA_filtered.vcf.gz \
-Oz -o ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_RNA_10.vcf.gz

bcftools view -i 'FORMAT/DP>=20' ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_RNA_filtered.vcf.gz \
-Oz -o ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_RNA_20.vcf.gz

bcftools view -i 'FORMAT/DP>=30' ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_RNA_filtered.vcf.gz \
-Oz -o ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_RNA_30.vcf.gz

bcftools view -i 'FORMAT/DP>=40' ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_RNA_filtered.vcf.gz \
-Oz -o ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_RNA_40.vcf.gz

################# STOP HERE
# ## Plink
# ## Unzip vcf.gz
# gunzip -c ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_RNA_10.vcf.gz > ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_RNA_10.vcf
# gunzip -c ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_RNA_20.vcf.gz > ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_RNA_20.vcf
# gunzip -c ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_RNA_30.vcf.gz > ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_RNA_30.vcf
# gunzip -c ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_RNA_40.vcf.gz > ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_RNA_40.vcf
# 
# gunzip -c ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_ATAC_10.vcf.gz > ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_ATAC_10.vcf
# gunzip -c ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_ATAC_20.vcf.gz > ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_ATAC_20.vcf
# gunzip -c ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_ATAC_30.vcf.gz > ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_ATAC_30.vcf
# gunzip -c ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_ATAC_40.vcf.gz > ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_ATAC_40.vcf
# 
# ## Filter MAF
# ### 0
# plink --vcf ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_RNA_10.vcf \
# --make-bed --out ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_RNA_10
# plink --vcf ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_RNA_20.vcf \
# --make-bed --out ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_RNA_20
# plink --vcf ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_RNA_30.vcf \
# --make-bed --out ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_RNA_30
# plink --vcf ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_RNA_40.vcf \
# --make-bed --out ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_RNA_40
# 
# plink --vcf ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_ATAC_10.vcf \
# --make-bed --out ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_ATAC_10
# plink --vcf ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_ATAC_20.vcf \
# --make-bed --out ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_ATAC_20
# plink --vcf ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_ATAC_30.vcf \
# --make-bed --out ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_ATAC_30
# plink --vcf ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_ATAC_40.vcf \
# --make-bed --out ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_ATAC_40

# ### 0.01
# Because only have 1 sample, we cannot filter MAF like this
# plink --vcf ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_RNA_10.vcf --maf 0.01 \
# --make-bed --out ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_RNA_10_maf1
# plink --vcf ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_RNA_20.vcf --maf 0.01 \
# --make-bed --out ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_RNA_20_maf1
# plink --vcf ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_RNA_30.vcf --maf 0.01 \
# --make-bed --out ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_RNA_30_maf1
# plink --vcf ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_RNA_40.vcf --maf 0.01 \
# --make-bed --out ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_RNA_40_maf1
# 
# plink --vcf ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_ATAC_10.vcf --maf 0.01 \
# --make-bed --out ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_ATAC_10_maf1
# plink --vcf ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_ATAC_20.vcf --maf 0.01 \
# --make-bed --out ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_ATAC_20_maf1
# plink --vcf ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_ATAC_30.vcf --maf 0.01 \
# --make-bed --out ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_ATAC_30_maf1
# plink --vcf ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_ATAC_40.vcf --maf 0.01 \
# --make-bed --out ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_ATAC_40_maf1
# 
# ### 0.05
# plink --vcf ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_RNA_10.vcf --maf 0.05 \
# --make-bed --out ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_RNA_10_maf5
# plink --vcf ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_RNA_20.vcf --maf 0.05 \
# --make-bed --out ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_RNA_20_maf5
# plink --vcf ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_RNA_30.vcf --maf 0.05 \
# --make-bed --out ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_RNA_30_maf5
# plink --vcf ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_RNA_40.vcf --maf 0.05 \
# --make-bed --out ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_RNA_40_maf5
# 
# plink --vcf ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_ATAC_10.vcf --maf 0.05 \
# --make-bed --out ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_ATAC_10_maf5
# plink --vcf ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_ATAC_20.vcf --maf 0.05 \
# --make-bed --out ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_ATAC_20_maf5
# plink --vcf ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_ATAC_30.vcf --maf 0.05 \
# --make-bed --out ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_ATAC_30_maf5
# plink --vcf ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_ATAC_40.vcf --maf 0.05 \
# --make-bed --out ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_ATAC_40_maf5

# # Remove files
# rm ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_RNA_10.vcf
# rm ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_RNA_20.vcf
# rm ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_RNA_30.vcf
# rm ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_RNA_40.vcf
# 
# rm ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_ATAC_10.vcf
# rm ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_ATAC_20.vcf
# rm ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_ATAC_30.vcf
# rm ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_ATAC_40.vcf
# 
# rm ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_RNA_10.vcf.gz
# rm ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_RNA_20.vcf.gz
# rm ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_RNA_30.vcf.gz
# rm ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_RNA_40.vcf.gz
# 
# rm ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_ATAC_10.vcf.gz
# rm ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_ATAC_20.vcf.gz
# rm ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_ATAC_30.vcf.gz
# rm ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_ATAC_40.vcf.gz

# ## Basic Statistics
# bcftools stats ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_filtered.vcf.gz > ${WORK_PATH}${SAMPLE_NAME}/germline/vcf_stats.txt
# bcftools stats ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_RNA_10.vcf.gz > ${WORK_PATH}${SAMPLE_NAME}/germline/vcf_stats1.txt
# bcftools stats ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_ATAC_10.vcf.gz > ${WORK_PATH}${SAMPLE_NAME}/germline/vcf_stats2.txt
# 
# ## View
# bcftools view ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_filtered.vcf.gz | head -n 5
# bcftools view ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}.vcf.gz | tail -n 5
# 
# ## Rename
# bcftools reheader -s ${WORK_PATH}${SAMPLE_NAME}/Sample_rename.txt \
# -o ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}.vcf.gz ${WORK_PATH}${SAMPLE_NAME}/germline/${SAMPLE_NAME}_genotyped.vcf.gz



