#!/bin/bash
#SBATCH -D /gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/baypass/trial3
#SBATCH --mail-type=END
#SBATCH --mail-user=jnascimento@geomar.de
#SBATCH --partition=cluster
#SBATCH --nodes=1
#SBATCH --tasks-per-node=18
#SBATCH --cpus-per-task=1
#SBATCH --mem=50G
#SBATCH --time=48:00:00
#SBATCH --job-name=baypasshaplo100
#SBATCH --output=baypass3.out
#SBATCH --error=baypass3.err

inputdir=/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/variants
outputdir=/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/baypass/trial3
baypassdir=/gxfs_home/geomar/smomw573/software/baypass_public-master/sources

$baypassdir/g_baypass \
-gfile $inputdir/genobaypass \
-efile $inputdir/cov-baypass.txt \
-poolsizefile $inputdir/haploid-size-baypass-txt \
-outprefix $outputdir/final_correct_haplo -nthreads 18
