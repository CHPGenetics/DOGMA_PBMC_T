#!/bin/bash
module load python/anaconda2.7-5.2.0
source activate /ix1/wchen/Shiyue/Biosoft/Python_Env/LDSC

PROJ_PATH=/ix1/wchen/Shiyue/Projects/2023_07_DOGMA_Revision/
REF=/ix1/wchen/Shiyue/References/LDSC/
pheno=IBD

#data path
cutoff="1E-5"
input_path=${PROJ_PATH}Results/02_S_LDSC/Input/
hmp_path=${REF}hapmap3_snps/
ldsc_outpath=${PROJ_PATH}Results/03_New_Enrichment/01_Celltype_DA/
pf_path=${REF}1000G_Phrase3/1000G_EUR_Phase3_plink/
baseline_path=${REF}1000G_Phrase3/1000G_EUR_Phase3_baseline/
weight_path=${REF}weights_hm3_no_hla/
frq_path=${REF}1000G_Phrase3/1000G_Phase3_frq/1000G.EUR.QC.

#software
ldsc=/ix1/wchen/Shiyue/Biosoft/LDSC/ldsc.py
munge_sumstats=/ix1/wchen/Shiyue/Biosoft/LDSC/munge_sumstats.py
mk_annot=/ix1/wchen/Shiyue/Biosoft/LDSC/make_annot.py

# Reference annotation
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
--annot-file ${ldsc_outpath}Ref/ref_${chr}.annot.gz

#ldsc
python2 ${ldsc} \
--l2 --bfile ${b_ref} \
--ld-wind-cm 1 \
--thin-annot \
--yes-really \
--print-snps ${hmp_file} \
--annot ${ldsc_outpath}Ref/ref_${chr}.annot.gz \
--out ${ldsc_outpath}Ref/ref_${chr}

done

#### Start from here
for pheno in ${pheno}
do
for chr in `seq 1 22`
do

b_ref=${pf_path}1000G.EUR.QC.${chr}
hmp_file=${hmp_path}hm.${chr}.snp

#LD score
python2 ${ldsc} \
--l2 --bfile ${b_ref} \
--ld-wind-cm 1 \
--yes-really \
--print-snps ${hmp_file} \
--annot ${ldsc_outpath}${pheno}/${cutoff}/S_LDSC/Anno_Merged_chr${chr}.annot.gz \
--out ${ldsc_outpath}${pheno}/${cutoff}/S_LDSC/Anno_Merged_chr${chr}

done
done

#Partitioned heritability
python2 ${ldsc} \
--h2 ${ldsc_outpath}${pheno}.sumstats.gz \
--ref-ld-chr ${ldsc_outpath}${pheno}/${cutoff}/S_LDSC/Anno_Merged_chr,${ldsc_outpath}Ref/ref_ \
--frqfile-chr ${frq_path} \
--w-ld-chr ${weight_path}weights. \
--overlap-annot \
--print-coefficients \
--print-delete-vals \
--out ${ldsc_outpath}Results/${pheno}_S_LDSC_Result

conda deactivate


