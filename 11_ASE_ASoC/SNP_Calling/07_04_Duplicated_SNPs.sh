#!/bin/bash
#SBATCH --job-name=Duplicated
#SBATCH --output=/ix1/wchen/Shiyue/Projects/2024_05_Allelic_Imbalance_Method/Code_New/Slurm_Out/%j.out
#SBATCH --error=/ix1/wchen/Shiyue/Projects/2024_05_Allelic_Imbalance_Method/Code_New/Slurm_Out/%j.err
#SBATCH --cluster=smp
#SBATCH --partition=smp
#SBATCH --cpus-per-task=12
#SBATCH --mem=170480
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

for depth in 10 20 30 40; do

for omics in RNA ATAC merged; do

plink --file ${WORK_PATH}Combine/Final/Merged_${omics}_${depth} \
--make-bed --out ${WORK_PATH}Combine/Final/${omics}_${depth}

done

done

