for i in $(ls ~/RWorkSpace/CITE-seq/Duerr/HTC/DOGMA_analysis/output/SNP/CHEERS/Disease | grep -v sh);
do python ~/RWorkSpace/CITE-seq/Duerr/HTC/DOGMA_analysis/data/SNP/CHEERS/CHEERS_computeEnrichment.py --input Merged_counts_normToMax_quantileNorm_euclideanNorm.txt --ld ~/RWorkSpace/CITE-seq/Duerr/HTC/DOGMA_analysis/output/SNP/CHEERS/Disease/${i} --trait $i --outdir Merged_analysis/;
done

for i in $(ls ~/RWorkSpace/CITE-seq/Duerr/HTC/DOGMA_analysis/output/SNP/CHEERS/IBD | grep -v sh);
do python ~/RWorkSpace/CITE-seq/Duerr/HTC/DOGMA_analysis/data/SNP/CHEERS/CHEERS_computeEnrichment.py --input Merged_counts_normToMax_quantileNorm_euclideanNorm.txt --ld ~/RWorkSpace/CITE-seq/Duerr/HTC/DOGMA_analysis/output/SNP/CHEERS/IBD/${i} --trait $i --outdir Merged_analysis/;
done

for i in $(ls ~/RWorkSpace/CITE-seq/Duerr/HTC/DOGMA_analysis/output/SNP/CHEERS/Disease | grep -v sh);
do python ~/RWorkSpace/CITE-seq/Duerr/HTC/DOGMA_analysis/data/SNP/CHEERS/CHEERS_computeEnrichment.py --input Individual_counts_normToMax_quantileNorm_euclideanNorm.txt --ld ~/RWorkSpace/CITE-seq/Duerr/HTC/DOGMA_analysis/output/SNP/CHEERS/Disease/${i} --trait $i --outdir Individual_analysis/;
done

for i in $(ls ~/RWorkSpace/CITE-seq/Duerr/HTC/DOGMA_analysis/output/SNP/CHEERS/IBD | grep -v sh);
do python ~/RWorkSpace/CITE-seq/Duerr/HTC/DOGMA_analysis/data/SNP/CHEERS/CHEERS_computeEnrichment.py --input Individual_counts_normToMax_quantileNorm_euclideanNorm.txt --ld ~/RWorkSpace/CITE-seq/Duerr/HTC/DOGMA_analysis/output/SNP/CHEERS/IBD/${i} --trait $i --outdir Individual_analysis/;
done