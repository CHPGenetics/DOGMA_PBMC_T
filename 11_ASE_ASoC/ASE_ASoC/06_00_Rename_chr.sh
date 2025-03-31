#!/bin/bash
#SBATCH --job-name=WASP
#SBATCH --output=/ix1/wchen/Shiyue/Projects/2024_05_Allelic_Imbalance_Method/ASE/Code/Slurm_Out/%j.out
#SBATCH --error=/ix1/wchen/Shiyue/Projects/2024_05_Allelic_Imbalance_Method/ASE/Code/Slurm_Out/%j.err
#SBATCH --cluster=smp
#SBATCH --partition=smp
#SBATCH --cpus-per-task=16
#SBATCH --mem=102480
#SBATCH --time=2:00:00
#SBATCH --array=1-25%25

# Parameters
vcf_path=/ix1/wchen/Shiyue/Projects/2024_05_Allelic_Imbalance_Method/Result/06_Separate_Samples/

# Module load
module unload gatk/4.5.0.0
module unload bcftools/1.15.1
module unload samtools/1.9
module unload htslib/1.9
# module spider gatk
module load gatk/4.5.0.0
module load bcftools/1.15.1
module load gcc/8.2.0
module load samtools/1.9
module load htslib/1.9

WORK_PATH=/ix1/wchen/Shiyue/Projects/2024_05_Allelic_Imbalance_Method/Result/03_SNP/
SAMPLE_NAME=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${WORK_PATH}Sample_test.txt)

bcftools annotate ${vcf_path}${SAMPLE_NAME}.vcf.gz --rename-chrs ${vcf_path}chr_rename.txt \
-Oz -o ${vcf_path}${SAMPLE_NAME}_updated.vcf.gz

bcftools index -t ${vcf_path}${SAMPLE_NAME}_updated.vcf.gz
