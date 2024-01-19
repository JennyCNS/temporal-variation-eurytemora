#!/bin/bash
#SBATCH -D /gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/baypass/
#SBATCH --mail-type=END
#SBATCH --mail-user=jnascimento@geomar.de
#SBATCH --partition=cluster
#SBATCH --nodes=5
#SBATCH --tasks-per-node=2
#SBATCH --cpus-per-task=4
#SBATCH --mem=50G
#SBATCH --time=20:00:00
#SBATCH --job-name=bayenv-conv
#SBATCH --output=baypass.out
#SBATCH --error=baypass.err

inputdir=/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/variants
outputdir=/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/baypass
baypassdir=/gxfs_home/geomar/smomw573/software/baypass_public-master/sources

$baypassdir/g_baypass \
-gfile $inputdir/genobaypass \
-efile $inputdir/cov-baypass.txt \
-poolsizefile $inputdir/haploid-size-baypass-txt \
-outprefix $outputdir/correct-haplo -nthreads 10
