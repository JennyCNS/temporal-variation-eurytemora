#!/bin/bash
#SBATCH -D /gxfs_home/geomar/smomw573/software/PGDSpider_2.1.1.5
#SBATCH --mail-type=END
#SBATCH --mail-user=jnascimento@geomar.de
#SBATCH --partition=cluster
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=42G
#SBATCH --time=40:00:00
#SBATCH --job-name=bayenv-conv
#SBATCH --output=bayenv_convert.out
#SBATCH --error=bayenv_convert.err

source activate java

java -Xmx33000m -Xms1024m -jar PGDSpider2-cli.jar -inputfile /gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/variants/final_20230810.vcf.gz -inputformat VCF -outputfile /gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/variants/final_20230810_trimmed.bayenv -outputformat BAYENV -spid /gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/variants/template_VCF_BAYENV.spid
