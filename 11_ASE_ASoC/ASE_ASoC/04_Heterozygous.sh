#!/bin/bash
#SBATCH --job-name=Heterozygous
#SBATCH --output=/ix1/wchen/Shiyue/Projects/2024_05_Allelic_Imbalance_Method/ASE/Code/Slurm_Out/%j.out
#SBATCH --error=/ix1/wchen/Shiyue/Projects/2024_05_Allelic_Imbalance_Method/ASE/Code/Slurm_Out/%j.err
#SBATCH --cluster=smp
#SBATCH --partition=smp
#SBATCH --cpus-per-task=16
#SBATCH --mem=102480
#SBATCH --time=24:00:00

# Parameters
work_path=/ix1/wchen/Shiyue/Projects/2024_05_Allelic_Imbalance_Method/Result/05_Phasing_Heterozygous/

module unload gcc/8.2.0
module unload htslib/1.9
module unload bcftools/1.15.1

module load gcc/8.2.0
module load htslib/1.9
module load bcftools/1.15.1

bcftools query -Hf '%ID\t%CHROM\t%POS\t%REF\t%ALT\n' \
"${work_path}Phased_hg38.vcf.gz" \
| awk '/^#/ {print "ID CHROM POS REF ALT"; next} {print}' \
| gzip -c > "${work_path}variants_info.txt.gz"

samples=$(bcftools query -l ${work_path}Phased_hg38.vcf.gz | awk '{printf("%s\t", $0)}')
echo -e "ID\t$samples" > ${work_path}header.txt

{
    cat ${work_path}header.txt;
    bcftools query -f '%ID[\t%GT]\n' ${work_path}Phased_hg38.vcf.gz
} | gzip -f > ${work_path}variants_genotype.txt.gz
