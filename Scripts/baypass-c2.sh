#!/bin/bash
#SBATCH -D /gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/
#SBATCH --mail-type=END
#SBATCH --mail-user=jnascimento@geomar.de
#SBATCH --partition=base
#SBATCH --nodes=5
#SBATCH --tasks-per-node=2
#SBATCH --cpus-per-task=4
#SBATCH --mem=50G
#SBATCH --time=24:00:00
#SBATCH --job-name=bp-c2
#SBATCH --output=baypass-c2.out
#SBATCH --error=baypass-c2.err


inputdir=/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/variants
outputdir=/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/baypass
baypassdir=/gxfs_home/geomar/smomw573/software/baypass_public-master/sources

$baypassdir/g_baypass \
-gfile $inputdir/genobaypass \
-efile $inputdir/cov-baypass.txt \
-omegafile $outputdir/core-model_mat_omega.out \
-poolsizefile $inputdir/haploid-size-baypass-txt \
-d0yij 20 \
-outprefix $outputdir/c2-model -nthreads 10 -contrastfile $inputdir/cov-baypass.txt
