#!/bin/bash
module load python/anaconda2.7-5.2.0
source activate /ix1/wchen/Shiyue/Biosoft/Python_Env/LDSC

PROJ_PATH=/ix1/wchen/Shiyue/Projects/2023_07_DOGMA_Revision/
REF=/ix1/wchen/Shiyue/References/LDSC/

#data path
input_path=${PROJ_PATH}Results/02_S_LDSC/Input/
hmp_path=${REF}hapmap3_snps/
ldsc_outpath=${PROJ_PATH}Results/02_S_LDSC/Output/Annotation_Peak3/
pf_path=${REF}1000G_Phrase3/1000G_EUR_Phase3_plink/

#software
ldsc=/ix1/wchen/Shiyue/Biosoft/LDSC/ldsc.py
mk_annot=/ix1/wchen/Shiyue/Biosoft/LDSC/make_annot.py
munge_sumstats=/ix1/wchen/Shiyue/Biosoft/LDSC/munge_sumstats.py

###ldsc reference
for chr in `seq 1 22`
do

gs_file=${input_path}Peak3_list.txt
gc_file=${input_path}Peak3_coord_all.txt
b_ref=${pf_path}1000G.EUR.QC.${chr}
hmp_file=${hmp_path}hm.${chr}.snp

#make annot
python2 ${mk_annot} \
--gene-set-file ${gs_file} \
--gene-coord-file ${gc_file} \
--windowsize 100000 \
--bimfile ${b_ref}.bim \
--annot-file ${ldsc_outpath}ref_${chr}.gz

#ldsc
python2 ${ldsc} \
--l2 --bfile ${b_ref} \
--ld-wind-cm 1 \
--annot ${ldsc_outpath}ref_${chr}.gz \
--thin-annot \
--out ${ldsc_outpath}ref_${chr} \
--print-snps ${hmp_file}
done

###ldsc
ldsc_outpath=${PROJ_PATH}Results/02_S_LDSC/Output/Pseudo_Celltype_DA/

for celltype in "CD4+_Regulatory_(Resting)" "CD4+_Memory_(Activated)_-_Other" "CD8+_Memory_(Resting)" "CD4+_Naive_(Activated)" "CD8+_Naive_(Resting)" "CD4+_Naive_(Resting)" "CD8+_Memory_(Activated)" "MAITs_(Resting)" "CD4+_Regulatory_(Activated)" "CD4+_Memory_(Resting)_-_Other" "CD8+_Naive_(Activated)" "CD4+_Memory_(Resting)_-_Th1" "MAITs_(Activated)" "Gamma_Delta" "CD4+_Memory_(Activated)_-_Th17"  "CD4+_Memory_(Activated)_-_Tfh"  "CD4+_Memory_(Activated)_-_Th1" "CD4+_Memory_(Resting)_-_Th17" "CD8+_Regulatory" "CD4+_Memory_(Resting)_-_Tfh" 
do
for chr in `seq 1 22`
do

gs_file=${input_path}Pseudo_Celltype_DA/Peak_list_${celltype}.txt
gc_file=${input_path}Peak3_coord_all.txt
b_ref=${pf_path}1000G.EUR.QC.${chr}
hmp_file=${hmp_path}hm.${chr}.snp

#make annot
python2 ${mk_annot} \
--gene-set-file ${gs_file} \
--gene-coord-file ${gc_file} \
--windowsize 100000 \
--bimfile ${b_ref}.bim \
--annot-file ${ldsc_outpath}${celltype}_${chr}.gz

#ldsc
python2 ${ldsc} \
--l2 --bfile ${b_ref} \
--ld-wind-cm 1 \
--annot ${ldsc_outpath}${celltype}_${chr}.gz \
--thin-annot \
--out ${ldsc_outpath}${celltype}_${chr} \
--print-snps ${hmp_file}

done
done

#GWAS 
###LDSC-cts
baseline_path=${REF}1000G_Phrase3/1000G_EUR_Phase3_baseline/
weight_path=${REF}weights_hm3_no_hla/
ldsc_outpath=${PROJ_PATH}Results/02_S_LDSC/Output/Pseudo_Celltype_DA/Result/

for pheno in IBD CD UC
do

python2 ${ldsc} \
--h2-cts ${ldsc_outpath}${pheno}.sumstats.gz \
--ref-ld-chr ${baseline_path}baseline. \
--out ${ldsc_outpath}${pheno}_S_LDSC_Result \
--ref-ld-chr-cts ${ldsc_outpath}ldcts \
--w-ld-chr ${weight_path}weights.

done

conda deactivate