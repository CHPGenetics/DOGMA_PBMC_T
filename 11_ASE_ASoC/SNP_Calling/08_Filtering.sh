#!/bin/bash
#SBATCH --job-name=Filtering
#SBATCH --output=/ix1/wchen/Shiyue/Projects/2024_05_Allelic_Imbalance_Method/Code_New/Slurm_Out/%j.out
#SBATCH --error=/ix1/wchen/Shiyue/Projects/2024_05_Allelic_Imbalance_Method/Code_New/Slurm_Out/%j.err
#SBATCH --cluster=smp
#SBATCH --partition=smp
#SBATCH --cpus-per-task=2
#SBATCH --mem=70480
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
KGP_SNP=/ix1/wchen/Shiyue/References/1KGP/hg38/

for depth in 10 20 30 40; do

for omics in RNA ATAC merged; do

# 1KGP SNPs
plink --bfile ${WORK_PATH}Combine/Final/${omics}_${depth} \
--extract ${KGP_SNP}snp_position.txt --make-bed \
--out ${WORK_PATH}Combine/Final/${omics}_${depth}_KGP

# 1KGP MAF 1% SNPs
plink --bfile ${WORK_PATH}Combine/Final/${omics}_${depth} \
--extract ${KGP_SNP}snp_position_0.01.txt --make-bed \
--out ${WORK_PATH}Combine/Final/${omics}_${depth}_KGP_0.01

# 1KGP MAF 5% SNPs
plink --bfile ${WORK_PATH}Combine/Final/${omics}_${depth} \
--extract ${KGP_SNP}snp_position_0.05.txt --make-bed \
--out ${WORK_PATH}Combine/Final/${omics}_${depth}_KGP_0.05

# Filter Geno in EVA-PR Genotype
plink --bfile ${WORK_PATH}Combine/Final/${omics}_${depth} \
--geno 0.5 --make-bed \
--out ${WORK_PATH}Combine/Final/${omics}_${depth}_Geno

# 1KGP SNPs
plink --bfile ${WORK_PATH}Combine/Final/${omics}_${depth}_Geno \
--extract ${KGP_SNP}snp_position.txt --make-bed \
--out ${WORK_PATH}Combine/Final/${omics}_${depth}_Geno_KGP

done

done


