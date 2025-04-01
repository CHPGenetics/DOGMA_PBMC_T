#!/bin/bash
#
#SBATCH -N 1
#SBATCH -t 3-00:00 # Runtime in D-HH:MM
#SBATCH -J R
#SBATCH --cpus-per-task=16 # Request that ncpus be allocated per process.
#SBATCH --mem-per-cpu=16G
#SBATCH --cluster=htc

#module load gcc/8.2.0 r/4.1.0
source ~/.bashrc
conda activate python2.7

python ~/RWorkSpace/DOGMA-seq/PBMC/data/SNP/CHEERS/CHEERS_normalize.py --input Merged/ --prefix Merged --outdir ./
python ~/RWorkSpace/DOGMA-seq/PBMC/data/SNP/CHEERS/CHEERS_normalize.py --input Individual/ --prefix Individual --outdir ./
python ~/RWorkSpace/DOGMA-seq/PBMC/data/SNP/CHEERS/CHEERS_normalize.py --input Merged_condition/ --prefix Merged_condition --outdir ./
