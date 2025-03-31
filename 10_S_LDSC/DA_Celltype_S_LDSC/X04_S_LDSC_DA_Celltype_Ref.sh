#!/bin/bash
#SBATCH --job-name=ldsc_analysis
#SBATCH --output=/ix1/wchen/Shiyue/Projects/2023_07_DOGMA_Revision/Results/02_S_LDSC/Output/Pseudo_Celltype_DA/sbatch_out/%j.out
#SBATCH --error=/ix1/wchen/Shiyue/Projects/2023_07_DOGMA_Revision/Results/02_S_LDSC/Output/Pseudo_Celltype_DA/sbatch_out/%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=3
#SBATCH --array=1-440%20

cell_types=("CD4+_Regulatory_(Resting)" "CD4+_Memory_(Activated)_-_Other" "CD8+_Memory_(Resting)" "CD4+_Naive_(Activated)" "CD8+_Naive_(Resting)" "CD4+_Naive_(Resting)" "CD8+_Memory_(Activated)" "MAITs_(Resting)" "CD4+_Regulatory_(Activated)" "CD4+_Memory_(Resting)_-_Other" "CD8+_Naive_(Activated)" "CD4+_Memory_(Resting)_-_Th1" "MAITs_(Activated)" "Gamma_Delta" "CD4+_Memory_(Activated)_-_Th17"  "CD4+_Memory_(Activated)_-_Tfh"  "CD4+_Memory_(Activated)_-_Th1" "CD4+_Memory_(Resting)_-_Th17" "CD8+_Regulatory" "CD4+_Memory_(Resting)_-_Tfh")

cell_type_index=$(( ($SLURM_ARRAY_TASK_ID - 1) / 22 ))
chr_index=$(( ($SLURM_ARRAY_TASK_ID - 1) % 22 + 1 ))

celltype=${cell_types[$cell_type_index]}
chr=$chr_index

module load python/anaconda2.7-5.2.0
source activate /ix1/wchen/Shiyue/Biosoft/Python_Env/LDSC

PROJ_PATH=/ix1/wchen/Shiyue/Projects/2023_07_DOGMA_Revision/
REF=/ix1/wchen/Shiyue/References/LDSC/

#software
ldsc=/ix1/wchen/Shiyue/Biosoft/LDSC/ldsc.py
mk_annot=/ix1/wchen/Shiyue/Biosoft/LDSC/make_annot.py
munge_sumstats=/ix1/wchen/Shiyue/Biosoft/LDSC/munge_sumstats.py

#data path
input_path=${PROJ_PATH}Results/02_S_LDSC/Input/
hmp_path=${REF}hapmap3_snps/
ldsc_outpath=${PROJ_PATH}Results/02_S_LDSC/Output/Annotation_Peak3/
pf_path=${REF}1000G_Phrase3/1000G_EUR_Phase3_plink/
gs_file=${input_path}Peak3_list.txt
gc_file=${input_path}Peak3_coord_all.txt
b_ref=${pf_path}1000G.EUR.QC.${chr}
hmp_file=${hmp_path}hm.${chr}.snp

# make annot
python2 ${mk_annot} \
    --gene-set-file ${gs_file} \
    --gene-coord-file ${gc_file} \
    --windowsize 100000 \
    --bimfile ${b_ref}.bim \
    --annot-file ${ldsc_outpath}${celltype}_${chr}.gz

# ldsc
python2 ${ldsc} \
    --l2 --bfile ${b_ref} \
    --ld-wind-cm 1 \
    --annot ${ldsc_outpath}${celltype}_${chr}.gz \
    --thin-annot \
    --out ${ldsc_outpath}${celltype}_${chr} \
    --print-snps ${hmp_file}
    