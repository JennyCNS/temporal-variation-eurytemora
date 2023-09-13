#creating final vcf

cd ~/work/seasonal_adaptation/analysis/PoolSNP

#filter bialelic positions




#!/bin/bash
#SBATCH -D /gxfs_home/geomar/smomw573/software/PoolSNP
#SBATCH --mail-type=END
#SBATCH --mail-user=jnascimento@geomar.de
#SBATCH --partition=cluster
#SBATCH --nodes=1
#SBATCH --tasks-per-node=5
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=10:00:00
#SBATCH --job-name=biallelic
#SBATCH --output=biallelic.out
#SBATCH --error=biallelic.err
module load vcftools/0.1.14
VCF=/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/PoolSNP/final-maf0.01-mincov50-mincount10-allsites0-excluding20072011.vcf.gz 
OUTPUT=/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/PoolSNP


vcftools --gzvcf $VCF --min-alleles 2 --max-alleles 2 --recode --out $OUT/finalfile_output_biallelic




#####


module load python/3.8.4

#run python script on vcf with filtered biallelic positions

python3 fix.py


zgrep "^#" finalfile_output_biallelic.recode.vcf.gz > header.txt

cat header.txt noheader2.txt > finalfile.vcf

#tabix the vcf
bgzip finalfile.vcf.gz
tabix -p vcf finalfile.vcf.gz
#add in header ##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">



#not needed but for future
#make bed from final vcf
module load bcftools/1.10.2
bcftools view -v snps input.vcf.gz -o output_snps.vcf.gz
bcftools query -f '%CHROM\t%POS0\t%POS1\n' finalfile.vcf > finalfile.bed

#add in header ##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">