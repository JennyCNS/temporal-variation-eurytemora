#!/bin/bash
#SBATCH -D /gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/
#SBATCH --mail-type=END
#SBATCH --mail-user=jnascimento@geomar.de
#SBATCH --partition=base
#SBATCH --nodes=5
#SBATCH --tasks-per-node=2
#SBATCH --cpus-per-task=4
#SBATCH --mem=50G
#SBATCH --time=80:00:00
#SBATCH --job-name=bayenv-conv
#SBATCH --output=baypass-core-2.out
#SBATCH --error=baypass-core-2.err

inputdir=/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/variants
outputdir=/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/baypass
baypassdir=/gxfs_home/geomar/smomw573/software/baypass_public-master/sources

$baypassdir/g_baypass \
-gfile $inputdir/genobaypass \
-poolsizefile $inputdir/haploid-size-baypass-txt \
-d0yij 10 \
-npilot \
-outprefix $outputdir/core-model-2 -nthreads 10
