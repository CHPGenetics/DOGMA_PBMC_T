#!/bin/bash
#SBATCH --job-name=GATK
#SBATCH --output=/ix1/wchen/Shiyue/Projects/2024_05_Allelic_Imbalance_Method/ASE/Code/Slurm_Out/%j.out
#SBATCH --error=/ix1/wchen/Shiyue/Projects/2024_05_Allelic_Imbalance_Method/ASE/Code/Slurm_Out/%j.err
#SBATCH --cluster=smp
#SBATCH --partition=smp
#SBATCH --cpus-per-task=16
#SBATCH --mem=102480
#SBATCH --time=48:00:00
#SBATCH --array=1-25%25

# Parameters
vcf_path=/ix1/wchen/Shiyue/Projects/2024_05_Allelic_Imbalance_Method/Result/06_Separate_Samples/
ref_vcf_path=/ix1/wchen/Shiyue/References/Fasta/Homo_sapiens_assembly38/Homo_sapiens_assembly38.fasta
ref_vcf_path_updated=/ix1/wchen/Shiyue/References/Fasta/Homo_sapiens_assembly38/Homo_sapiens_assembly38_22.fasta
ref_path=/ix1/wchen/Shiyue/References/Fasta/hg38/hg38.fa

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

WORK_PATH=/ix1/wchen/Shiyue/Projects/2024_05_Allelic_Imbalance_Method/Result/
SAMPLE_NAME=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${WORK_PATH}03_SNP/Sample_test.txt)

condition="Act_IL1B_IL23_PGE2_TGFB"

work_path=/ix1/wchen/Shiyue/Projects/2024_05_Allelic_Imbalance_Method/ASE/Result/03_GATK_Allelic_Reads_ATAC/${condition}/
bam_path=/ix1/wchen/Shiyue/Projects/2024_05_Allelic_Imbalance_Method/ASE/Result/01_Barcodes/${SAMPLE_NAME}/${condition}/

for dir in ${bam_path}*; do
if [ -d "$dir" ]; then

celltype="$(basename "$dir")"

if [ -e ${bam_path}${celltype}/ATAC.bam ]; then

samtools sort -o ${bam_path}${celltype}/ATAC.sort.bam ${bam_path}${celltype}/ATAC.bam
samtools index ${bam_path}${celltype}/ATAC.sort.bam

rm ${bam_path}${celltype}/ATAC.bam
rm ${bam_path}${celltype}/ATAC.bam.bai

fi

samtools addreplacerg \
-r 'ID:1' \
-r 'PL:Illumina' \
-r 'SM:'${SAMPLE_NAME} \
-r 'LB:Library1' \
-o ${work_path}${SAMPLE_NAME}_${condition}_${celltype}_ATAC.bam \
${bam_path}${celltype}/ATAC.sort.bam

samtools index ${work_path}${SAMPLE_NAME}_${condition}_${celltype}_ATAC.bam

# Reads Count
gatk ASEReadCounter \
--disable-read-filter NotDuplicateReadFilter \
--variant ${vcf_path}${SAMPLE_NAME}_updated2.vcf.gz \
--output-format CSV \
--reference ${ref_path} \
--input ${work_path}${SAMPLE_NAME}_${condition}_${celltype}_ATAC.bam \
--output ${work_path}${SAMPLE_NAME}_${condition}_${celltype}.gatkAlleleReadCounts.csv

rm ${work_path}${SAMPLE_NAME}_${condition}_${celltype}_ATAC.bam
rm ${work_path}${SAMPLE_NAME}_${condition}_${celltype}_ATAC.bam.bai

fi
done
