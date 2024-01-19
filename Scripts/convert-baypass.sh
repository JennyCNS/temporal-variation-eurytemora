#!/bin/bash
#SBATCH -D /gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/
#SBATCH --mail-type=END
#SBATCH --mail-user=jnascimento@geomar.de
#SBATCH --partition=cluster
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=47:00:00
#SBATCH --job-name=trim
#SBATCH --output=baypass_convert.out
#SBATCH --error=baypass_convert.err


python Toolbox/reshaper_baypass.py /gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/variants/final_20230810_subset.vcf /gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/variants/popmap2.txt /gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/variants/final_20230810_subset.baypass

