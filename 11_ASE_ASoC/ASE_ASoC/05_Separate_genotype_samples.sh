#!/bin/bash
#SBATCH --job-name=Genotype
#SBATCH --output=/ix1/wchen/Shiyue/Projects/2024_05_Allelic_Imbalance_Method/ASE/Code/Slurm_Out/%j.out
#SBATCH --error=/ix1/wchen/Shiyue/Projects/2024_05_Allelic_Imbalance_Method/ASE/Code/Slurm_Out/%j.err
#SBATCH --cluster=smp
#SBATCH --partition=smp
#SBATCH --cpus-per-task=16
#SBATCH --mem=102480
#SBATCH --time=24:00:00

# Parameters
data_path=/ix1/wchen/Shiyue/Projects/2024_05_Allelic_Imbalance_Method/Result/05_Phasing_Heterozygous/
work_path=/ix1/wchen/Shiyue/Projects/2024_05_Allelic_Imbalance_Method/Result/06_Separate_Samples/

module unload gcc/8.2.0
module unload htslib/1.9
module unload bcftools/1.15.1

module load gcc/8.2.0
module load htslib/1.9
module load bcftools/1.15.1

sample_ids=$(bcftools query -l "${data_path}Phased_hg38.vcf.gz")

for sample_id in $sample_ids; do
    output_file="${work_path}${sample_id}.vcf.gz"

    bcftools view "${data_path}Phased_hg38.vcf.gz" -s "$sample_id" -Oz -o "$output_file"
    bcftools index "$output_file"
done
