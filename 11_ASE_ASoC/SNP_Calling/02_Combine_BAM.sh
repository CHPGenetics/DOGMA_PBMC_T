#!/bin/bash
#SBATCH --job-name=BAM
#SBATCH --output=/ix1/wchen/Shiyue/Projects/2024_05_Allelic_Imbalance_Method/Code/Slurm_Out/%j.out
#SBATCH --error=/ix1/wchen/Shiyue/Projects/2024_05_Allelic_Imbalance_Method/Code/Slurm_Out/%j.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=102480
#SBATCH --time=24:00:00
#SBATCH --array=1-15%8

# Parameters
data_path=/ix1/wchen/Shiyue/Projects/2024_05_Allelic_Imbalance_Method/Result/01_SNPs_Calling/
work_path=/ix1/wchen/Shiyue/Projects/2024_05_Allelic_Imbalance_Method/Result/02_Combine_BAM/

# Module load
module load gcc/8.2.0
module load samtools/1.9

# Duerr_20210419
DOGMA1=Duerr_20210419_DOGMAseq_DIG
DOGMA2=Duerr_20210419_DOGMAseq_PFA_LLL
## SB770848
sample=SB770848
### ATAC
samtools merge -f ${work_path}${sample}_ATAC_raw.bam ${data_path}${DOGMA1}/${sample}/ATAC.bam \
${data_path}${DOGMA2}/${sample}/ATAC.bam
samtools sort -o ${work_path}${sample}_ATAC.bam ${work_path}${sample}_ATAC_raw.bam
samtools index ${work_path}${sample}_ATAC.bam
rm ${work_path}${sample}_ATAC_raw.bam
### GEX
samtools merge -f ${work_path}${sample}_ATAC_raw.bam ${data_path}${DOGMA1}/${sample}/ATAC.bam \
${data_path}${DOGMA2}/${sample}/ATAC.bam
samtools sort -o ${work_path}${sample}_ATAC.bam ${work_path}${sample}_ATAC_raw.bam
samtools index ${work_path}${sample}_ATAC.bam
rm ${work_path}${sample}_ATAC_raw.bam




