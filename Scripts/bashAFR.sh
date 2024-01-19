#!/bin/bash
#SBATCH -D /gxfs_work1/geomar/smomw573/seasonal_adaptation/analysis/variants
#SBATCH --mail-type=END
#SBATCH --mail-user=jnascimento@geomar.de
#SBATCH --partition=cluster
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=4:00:00
#SBATCH --job-name=AFmatrix
#SBATCH --output=afmatrix.out
#SBATCH --error=afmatrix.err


module load R/4.1.1
module load gcc/10.2.0

Rscript allele-freq-table.R

