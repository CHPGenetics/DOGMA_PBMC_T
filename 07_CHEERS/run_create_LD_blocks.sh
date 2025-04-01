for i in $(ls -I *.sh);
do cd $i;
mkdir ../IBD/${i}
python ~/RWorkSpace/DOGMA-seq/PBMC/data/SNP/CHEERS/create_LD_blocks.py SNP.txt ../IBD/${i} ~/RWorkSpace/DOGMA-seq/PBMC/data/SNP/LD_GRCh38;
cd ../;
done