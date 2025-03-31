#!/bin/bash
#SBATCH --job-name=Phasing
#SBATCH --output=/ix1/wchen/Shiyue/Projects/2024_05_Allelic_Imbalance_Method/ASE/Code/Slurm_Out/%j.out
#SBATCH --error=/ix1/wchen/Shiyue/Projects/2024_05_Allelic_Imbalance_Method/ASE/Code/Slurm_Out/%j.err
#SBATCH --cluster=smp
#SBATCH --partition=smp
#SBATCH --cpus-per-task=16
#SBATCH --mem=102480
#SBATCH --time=24:00:00

# Parameters
ref_path=/ix1/wchen/Shiyue/References/Phasing/
data_path=/ix1/wchen/Shiyue/Projects/2024_05_Allelic_Imbalance_Method/Result/03_SNP/Combine/Final/
work_path=/ix1/wchen/Shiyue/Projects/2024_05_Allelic_Imbalance_Method/Result/05_Phasing_Heterozygous/

module unload gcc/8.2.0
module unload htslib/1.9
module unload shapeit4/4.1.3
module unload bcftools/1.15.1
module unload plink/1.90b6.7

module load gcc/8.2.0
module load htslib/1.9
module load shapeit4/4.1.3
module load bcftools/1.15.1
module load plink/1.90b6.7

# Separate 22 chrom
for chr in `seq 1 22`;do

plink --bfile ${data_path}merged_10_Geno --chr ${chr} --recode vcf-iid --out ${work_path}chr${chr}_hg38
bgzip -c ${work_path}chr${chr}_hg38.vcf > ${work_path}chr${chr}_hg38.vcf.gz
tabix -p vcf ${work_path}chr${chr}_hg38.vcf.gz

echo chr${chr}

done

# Phasing
for chr in `seq 1 22`
do

shapeit4 --thread 20 --input ${work_path}chr${chr}_hg38.vcf.gz \
--map ${ref_path}chr${chr}.b38.gmap --region ${chr} \
--output ${work_path}Phased_chr${chr}_hg38.vcf.gz

done

# Concat
for i in `seq 1 22`
do
   file_list+=" ${work_path}Phased_chr${i}_hg38.vcf.gz"
done

bcftools concat -Oz -o ${work_path}Phased_hg38.vcf.gz ${file_list}
