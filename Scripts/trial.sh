#!/bin/bash
#SBATCH -D /gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/PoolSNP
#SBATCH --mail-type=END
#SBATCH --mail-user=jnascimento@geomar.de
#SBATCH --partition=cluster
#SBATCH --nodes=1
#SBATCH --tasks-per-node=14
#SBATCH --cpus-per-task=1
#SBATCH --mem=100G
#SBATCH --qos=long
#SBATCH --time=150:00:00
#SBATCH --job-name=PoolSNP-trial-excludedsamples
#SBATCH --output=poolsnptrial-nofilters.out
#SBATCH --error=poolsnptrial-nofilters.err

bash /gxfs_home/geomar/smomw573/software/PoolSNP/PoolSNP.sh \
mpileup=/gxfs_work1/geomar/smomw573/seasonal_adaptation/analysis/variants/trial.mpileup  \
output=/gxfs_work1/geomar/smomw573/seasonal_adaptation/analysis/PoolSNP/blabla \
reference=/gxfs_work1/geomar/smomw573/seasonal_adaptation/genome/Eaffinis.Baltic.PseudoRef.Mar22.fasta \
names=EA_2009_T1,EA_2009_T2,EA_2009_T3,EA_2009_T4,EA_2011_T1,EA_2011_T2,EA_2015_T1,EA_2015_T2,EA_2015_T3,EA_2015_T4,EA_2022_T1,EA_2022_T2,EA_2022_T3,2022_T1,EA_2022_T2,EA_2022_ ,EA_2022_T4 \
min-cov=50 \
max-cov=0.99 \
miss-frac=0 \
base-quality 15 \
jobs=14 \
badsites=1 \
allsites=0

