#!/bin/bash
#SBATCH -D /gxfs_home/geomar/smomw573/software/PoolSNP
#SBATCH --mail-type=END
#SBATCH --mail-user=jnascimento@geomar.de
#SBATCH --partition=cluster
#SBATCH --nodes=1
#SBATCH --tasks-per-node=17
#SBATCH --cpus-per-task=1
#SBATCH --mem=100G
#SBATCH --qos=long
#SBATCH --time=150:00:00
#SBATCH --job-name=PoolSNP-trial5
#SBATCH --output=poolsnptrial5.out
#SBATCH --error=poolsnptrial5.err



bash /gxfs_home/geomar/smomw573/software/PoolSNP/PoolSNP.sh \
mpileup=/gxfs_work1/geomar/smomw573/seasonal_adaptation/analysis/PoolSNP/all.mpileup \
output=/gxfs_work1/geomar/smomw573/seasonal_adaptation/analysis/PoolSNP/output_maf0.01mincov50trial5 \
reference=/gxfs_work1/geomar/smomw573/seasonal_adaptation/genome/Eaffinis.Baltic.PseudoRef.Mar22.fasta \
names=EA_2007_T1,EA_2007_T2,EA_2009_T1,EA_2009_T2,EA_2009_T3,EA_2009_T4,EA_2011_T1,EA_2011_T2,EA_2011_T3,EA_2015_T1,EA_2015_T2,EA_2015_T3,EA_2015_T4,EA_2022_T1,EA_2022_T2,EA_2022_T3,EA_2022_T4 \
min-cov=50 \
max-cov=0.95 \
min-count=20 \
min-freq=0.01 \
miss-frac=0.1 \
base-quality 15 \
jobs=17 \
badsites=1 \
allsites=0

#trial 1 wrong sample name input
#trial 2 did not have scripts directory in the Poolsnp directory
#trial 3 ran together with nohup so missed header
#trial 4 fixed all of the above but ran with all sites = 1
#trial 5 = trial 4 but no sites = 0
#trial 6 = trial 5 but different working direcotry (softare) and all sites = 1
#trial 7 = trial 6 but min count 10


#so we will use trial 7 as a guide of which samples to keep or not because min count choice of 10 is more appropriate to our dataset and more precise.
# first remove SNPs from 7



module load bcftools/1.10.2
bcftools view -v snps output_maf0.01mincov50trial7.vcf.gz -o allsamples_snps.vcf.gz

###output
#output vcf trial5 -only snps min count 20
#21377
#output vcf trial6
#150450144 all variants



#output vcf trial7 -only snps min count 10
#13 882 122

#output vcf removing bad samples, min cov 50 maf 0.01 min count 5
#2 424 281 snps

#out-maf0.01-mincov50-mincount10-allsites-excluding20072011.vcf.gz
#150422214

#Min count 5 2424281
#Min count 10 2424259



###20.07.2023
#after running the inital SNP call with these filters we noticed that both 2007 samples and 2011_t3 have a lot of missing data.
#we decided to remove this samples (every sample missing more than 95% of all variants)

#made new mpileup file
#reran poolsnp without these samples
bamfiles.txt 
/gxfs_work1/geomar/smomw504/seasonal_adaptation/analysis/aligned/mkdups/EA_2009_T1.bam
/gxfs_work1/geomar/smomw504/seasonal_adaptation/analysis/aligned/mkdups/EA_2009_T2.bam
/gxfs_work1/geomar/smomw504/seasonal_adaptation/analysis/aligned/mkdups/EA_2009_T3.bam
/gxfs_work1/geomar/smomw504/seasonal_adaptation/analysis/aligned/mkdups/EA_2009_T4.bam
/gxfs_work1/geomar/smomw504/seasonal_adaptation/analysis/aligned/mkdups/EA_2011_T1.bam
/gxfs_work1/geomar/smomw504/seasonal_adaptation/analysis/aligned/mkdups/EA_2011_T2.bam
/gxfs_work1/geomar/smomw504/seasonal_adaptation/analysis/aligned/mkdups/EA_2015_T1.bam
/gxfs_work1/geomar/smomw504/seasonal_adaptation/analysis/aligned/mkdups/EA_2015_T2.bam
/gxfs_work1/geomar/smomw504/seasonal_adaptation/analysis/aligned/mkdups/EA_2015_T3.bam
/gxfs_work1/geomar/smomw504/seasonal_adaptation/analysis/aligned/mkdups/EA_2015_T4.bam
/gxfs_work1/geomar/smomw504/seasonal_adaptation/analysis/aligned/mkdups/EA_2022_T1.bam
/gxfs_work1/geomar/smomw504/seasonal_adaptation/analysis/aligned/mkdups/EA_2022_T2.bam
/gxfs_work1/geomar/smomw504/seasonal_adaptation/analysis/aligned/mkdups/EA_2022_T3.bam
/gxfs_work1/geomar/smomw504/seasonal_adaptation/analysis/aligned/mkdups/EA_2022_T4.bam


#21.07.2023
#check snps in the output files

module load bcftools/1.10.2
bcftools view -v snps input.vcf.gz -o output_snps.vcf.gz