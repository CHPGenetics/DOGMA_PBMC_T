#!/bin/bash
#
#SBATCH -N 1
#SBATCH -t 6-00:00 # Runtime in D-HH:MM
#SBATCH -J R
#SBATCH --cpus-per-task=8 # Request that ncpus be allocated per process.
#SBATCH --mem-per-cpu=60G
#SBATCH --cluster=htc

#module load gcc/8.2.0 r/4.1.0
source ~/.bashrc
conda activate python2.7

ls ~/RWorkSpace/DOGMA-seq/PBMC/output/SNP/CHEERS/IBD | grep -v All | xargs -I var -n 1 -P 2 python ~/RWorkSpace/DOGMA-seq/PBMC/data/SNP/CHEERS/CHEERS_computeEnrichment.py --input Individual_condition_counts_normToMax_quantileNorm_euclideanNorm.txt --ld ~/RWorkSpace/DOGMA-seq/PBMC/output/SNP/CHEERS/IBD/var --trait var --outdir Individual_condition_analysis/
