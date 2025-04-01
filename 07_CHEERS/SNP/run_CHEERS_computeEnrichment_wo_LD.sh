for i in $(ls ~/RWorkSpace/CITE-seq/Duerr/HTC/DOGMA_analysis/output/SNP/CHEERS/Disease | grep -v sh);
do python ~/RWorkSpace/CITE-seq/Duerr/HTC/DOGMA_analysis/data/SNP/CHEERS/CHEERS_computeEnrichment.py --input Merged_counts_normToMax_quantileNorm_euclideanNorm.txt --snp_list ~/RWorkSpace/CITE-seq/Duerr/HTC/DOGMA_analysis/output/SNP/CHEERS/Disease_SNP/${i}/SNP_list.txt --trait $i --outdir Merged_analysis_wo_LD/;
done

for i in $(ls ~/RWorkSpace/CITE-seq/Duerr/HTC/DOGMA_analysis/output/SNP/CHEERS/IBD | grep -v sh);
do python ~/RWorkSpace/CITE-seq/Duerr/HTC/DOGMA_analysis/data/SNP/CHEERS/CHEERS_computeEnrichment.py --input Merged_counts_normToMax_quantileNorm_euclideanNorm.txt --snp_list ~/RWorkSpace/CITE-seq/Duerr/HTC/DOGMA_analysis/output/SNP/CHEERS/IBD_SNP/${i}/SNP_list.txt --trait $i --outdir Merged_analysis_wo_LD/;
done

for i in $(ls ~/RWorkSpace/CITE-seq/Duerr/HTC/DOGMA_analysis/output/SNP/CHEERS/Disease | grep -v sh);
do python ~/RWorkSpace/CITE-seq/Duerr/HTC/DOGMA_analysis/data/SNP/CHEERS/CHEERS_computeEnrichment.py --input Individual_counts_normToMax_quantileNorm_euclideanNorm.txt --snp_list ~/RWorkSpace/CITE-seq/Duerr/HTC/DOGMA_analysis/output/SNP/CHEERS/Disease_SNP/${i}/SNP_list.txt --trait $i --outdir Individual_analysis_wo_LD/;
done

for i in $(ls ~/RWorkSpace/CITE-seq/Duerr/HTC/DOGMA_analysis/output/SNP/CHEERS/IBD | grep -v sh);
do python ~/RWorkSpace/CITE-seq/Duerr/HTC/DOGMA_analysis/data/SNP/CHEERS/CHEERS_computeEnrichment.py --input Individual_counts_normToMax_quantileNorm_euclideanNorm.txt --snp_list ~/RWorkSpace/CITE-seq/Duerr/HTC/DOGMA_analysis/output/SNP/CHEERS/IBD_SNP/${i}/SNP_list.txt --trait $i --outdir Individual_analysis_wo_LD/;
done