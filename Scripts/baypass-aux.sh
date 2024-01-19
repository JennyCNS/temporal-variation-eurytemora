#!/bin/bash
#SBATCH -D /gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/
#SBATCH --mail-type=END
#SBATCH --mail-user=jnascimento@geomar.de
#SBATCH --partition=base
#SBATCH --nodes=5
#SBATCH --tasks-per-node=2
#SBATCH --cpus-per-task=4
#SBATCH --mem=50G
#SBATCH --time=22:00:00
#SBATCH --job-name=bayenv-aux
#SBATCH --output=baypass-aux.out
#SBATCH --error=baypass-aux.err


inputdir=/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/variants
outputdir=/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/baypass
baypassdir=/gxfs_home/geomar/smomw573/software/baypass_public-master/sources

$baypassdir/g_baypass \
-gfile $inputdir/genobaypass \
-efile $inputdir/cov-baypass.txt.constrast \
-omegafile $outputdir/core-model_mat_omega.out \
-poolsizefile $inputdir/haploid-size-baypass-txt \
-d0yij 20 \
-outprefix $outputdir/aux-model -nthreads 10 -auxmodel 
