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

condition="Act_IL1B_IL23_PGE2"

work_path=/ix1/wchen/Shiyue/Projects/2024_05_Allelic_Imbalance_Method/ASE/Result/03_GATK_Allelic_Reads_RNA/${condition}/
bam_path=/ix1/wchen/Shiyue/Projects/2024_05_Allelic_Imbalance_Method/ASE/Result/01_Barcodes/${SAMPLE_NAME}/${condition}/


gunzip ${vcf_path}${SAMPLE_NAME}_updated.vcf.gz
sed '/_Y_/d; /_K_/d; /_W_/d; /_R_/d' ${vcf_path}${SAMPLE_NAME}_updated.vcf > ${vcf_path}${SAMPLE_NAME}_selected.vcf
bgzip -c ${vcf_path}${SAMPLE_NAME}_selected.vcf > ${vcf_path}${SAMPLE_NAME}_selected.vcf.gz
bcftools index -t ${vcf_path}${SAMPLE_NAME}_selected.vcf.gz

gatk SelectVariants \
--reference ${ref_vcf_path_updated} \
--variant ${vcf_path}${SAMPLE_NAME}_selected.vcf.gz \
--restrict-alleles-to BIALLELIC \
--select 'vc.getHetCount()==1' \
--select-type-to-include SNP \
-O ${work_path}${SAMPLE_NAME}_dna_variants.selected.vcf.bgz

bcftools norm --rm-dup all ${work_path}${SAMPLE_NAME}_dna_variants.selected.vcf.bgz | bgzip > ${vcf_path}${SAMPLE_NAME}_updated2.vcf.gz
bcftools index -t ${vcf_path}${SAMPLE_NAME}_updated2.vcf.gz

rm ${work_path}${SAMPLE_NAME}_dna_variants.selected.vcf.bgz
rm ${work_path}${SAMPLE_NAME}_dna_variants.selected.vcf.bgz.tbi

for dir in ${bam_path}*; do
if [ -d "$dir" ]; then

celltype="$(basename "$dir")"

if [ -e ${bam_path}${celltype}/GEX.bam ]; then

samtools sort -o ${bam_path}${celltype}/GEX.sort.bam ${bam_path}${celltype}/GEX.bam
samtools index ${bam_path}${celltype}/GEX.sort.bam

rm ${bam_path}${celltype}/GEX.bam
rm ${bam_path}${celltype}/GEX.bam.bai

fi

samtools addreplacerg \
-r 'ID:1' \
-r 'PL:Illumina' \
-r 'SM:'${SAMPLE_NAME} \
-r 'LB:Library1' \
-o ${work_path}${SAMPLE_NAME}_${condition}_${celltype}.bam \
${bam_path}${celltype}/GEX.sort.bam

samtools index ${work_path}${SAMPLE_NAME}_${condition}_${celltype}.bam

# Reads Count
gatk ASEReadCounter \
--disable-read-filter NotDuplicateReadFilter \
--variant ${vcf_path}${SAMPLE_NAME}_updated2.vcf.gz \
--output-format CSV \
--reference ${ref_path} \
--input ${work_path}${SAMPLE_NAME}_${condition}_${celltype}.bam \
--output ${work_path}${SAMPLE_NAME}_${condition}_${celltype}.gatkAlleleReadCounts.csv

rm ${work_path}${SAMPLE_NAME}_${condition}_${celltype}.bam
rm ${work_path}${SAMPLE_NAME}_${condition}_${celltype}.bam.bai

fi
done


