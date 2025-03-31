#!/bin/bash
#SBATCH --job-name=BAM
#SBATCH --output=/ix1/wchen/Shiyue/Projects/2024_05_Allelic_Imbalance_Method/ASE/Code/Slurm_Out/%j.out
#SBATCH --error=/ix1/wchen/Shiyue/Projects/2024_05_Allelic_Imbalance_Method/ASE/Code/Slurm_Out/%j.err
#SBATCH --cluster=smp
#SBATCH --partition=smp
#SBATCH --cpus-per-task=16
#SBATCH --mem=102480
#SBATCH --time=24:00:00
#SBATCH --array=1-25%13

# Parameters
work_path=/ix1/wchen/Shiyue/Projects/2024_05_Allelic_Imbalance_Method/ASE/Result/01_Barcodes/
bam_path=/ix1/wchen/Shiyue/Projects/2024_05_Allelic_Imbalance_Method/Result/02_Combine_BAM/
SAMPLE_NAME=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${work_path}Sample.txt)

# Module load
module unload gcc/8.2.0
module unload samtools/1.9
module unload subset-bam/1.1.0
module unload bedtools/2.31.0

module load gcc/8.2.0
module load samtools/1.9
module load subset-bam/1.1.0
module load bedtools/2.31.0

for condition in Act_IL1B_IL23 Act_IL1B_IL23_PGE2 Act_IL1B_IL23_TGFB Act_IL1B_IL23_PGE2_TGFB; do

for celltype in ${work_path}${SAMPLE_NAME}/${condition}/*; do

celltype=$(basename "$celltype")

# Separate cell-type specific BAM files
## GEX
export TMPDIR=${work_path}${SAMPLE_NAME}/${condition}/${celltype}/
subset-bam --cores 10 --bam ${bam_path}${SAMPLE_NAME}_GEX.bam \
--cell-barcodes ${work_path}${SAMPLE_NAME}/${condition}/${celltype}/RNA/Barcodes.txt \
--out-bam ${work_path}${SAMPLE_NAME}/${condition}/${celltype}/GEX.bam
samtools index ${work_path}${SAMPLE_NAME}/${condition}/${celltype}/GEX.bam

## ATAC
subset-bam --cores 10 --bam ${bam_path}${SAMPLE_NAME}_ATAC.bam \
--cell-barcodes ${work_path}${SAMPLE_NAME}/${condition}/${celltype}/RNA/Barcodes.txt \
--out-bam ${work_path}${SAMPLE_NAME}/${condition}/${celltype}/ATAC.bam
samtools index ${work_path}${SAMPLE_NAME}/${condition}/${celltype}/ATAC.bam


done

done


