#!/bin/bash
#SBATCH --job-name=BAM
#SBATCH --output=/ix1/wchen/Shiyue/Projects/2024_05_Allelic_Imbalance_Method/Code/Slurm_Out/%j.out
#SBATCH --error=/ix1/wchen/Shiyue/Projects/2024_05_Allelic_Imbalance_Method/Code/Slurm_Out/%j.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=102480
#SBATCH --time=24:00:00
#SBATCH --array=1-3%3

# Parameters
data_path=/ix1/rduerr/shared/rduerr_wchen/Shiyue/02_Allelic_Imbalance_Method_2024/Data/
work_path=/ix1/wchen/Shiyue/Projects/2024_05_Allelic_Imbalance_Method/Result/01_SNPs_Calling/
# bam_path=/ix1/rduerr/shared/rduerr_wchen/raw/Cell_Ranger_output/cellranger-arc_count_2.0.1/
bam_path=/ix1/rduerr/shared/rduerr_wchen/raw/Cell_Ranger_output/cellranger-arc_count_2.0.2/

# Module load
module load gcc/8.2.0
module load samtools/1.9
module load subset-bam/1.1.0

DOGMA_NAME=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${data_path}samples_add.txt)
DOGMA_GEX=${bam_path}${DOGMA_NAME}_arc/gex_possorted_bam.bam
DOGMA_ATAC=${bam_path}${DOGMA_NAME}_arc/atac_possorted_bam.bam

cd ${work_path}${DOGMA_NAME}
samples=$(ls -d */ | sed 's#/##')
for BAM_SAMPLE in ${samples}; do

BAM_FILE=${work_path}${DOGMA_NAME}/${BAM_SAMPLE}/
BARCODES_GEX=${BAM_FILE}Barcodes_GEX.txt
BARCODES_ATAC=${BAM_FILE}Barcodes_ATAC.txt

# GEX
export TMPDIR=${BAM_FILE}
subset-bam --cores 10 --bam ${DOGMA_GEX} --cell-barcodes ${BARCODES_GEX} --out-bam ${BAM_FILE}GEX.bam
samtools index ${BAM_FILE}GEX.bam

# ATAC
subset-bam --cores 10 --bam ${DOGMA_ATAC} --cell-barcodes ${BARCODES_GEX} --out-bam ${BAM_FILE}ATAC.bam
samtools index ${BAM_FILE}ATAC.bam

echo "Filtered BAM file saved as ${BAM_FILE}"
done
