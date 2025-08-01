#new analysis COPES - Eurytemora affinis new genome


# bash scripts
1.align_new.sh # aligning samples to new genome, RG tags, sam to bam
2.markdup # mark duplicates in bam files
4.stats-mkdups.sh #alignment stats with samtools
3.stats-beforedups.sh
5. to-mpileup.sh #creating mpileup file
6.poolsnp-all-samples.sh #all samples, min cov 50, min count 10 (93 snps)
6.2.poolsnp-all-samples.sh #min cov 10 instead of 50 (205 snps)
6.3.poolsnp-all-samples.sh # min count 5 instead of 10 (93 snps again)
6.4.poolsnp-no2007.sh #min count 5, min cov 30 (maf 0.01) excluding 2007 samples

# R Scripts


#date 19.11.24

I will rerun the analyses with the new genome given by Carol Lee's postdoc

```bash
#index the genome
#!/bin/bash
#SBATCH -D /gxfs_home/geomar/smomw573/work/seasonal_adaptation/log_files/array_jobs
#SBATCH --mail-type=END
#SBATCH --mail-user=jnascimento@geomar.de
#SBATCH --partition=base
#SBATCH --cpus-per-task=4
#SBATCH --mem=60G
#SBATCH --time=05:00:00
#SBATCH --job-name=bwa
#SBATCH --output=%x.%A.%a.out
#SBATCH --error=%x.%A.%a.err

source ~/miniconda3/bin/activate bwa


indir=/gxfs_work/geomar/smomw573/smomw573/seasonal_adaptation/raw_data/trimmed
my_bwa=~/miniconda3/bin/bwa-mem2
$my_bwa index /gxfs_work/geomar/smomw573/smomw573/seasonal_adaptation/new_genome/Purple_Eaffinis_genome_assembly.fa
```
script
1.align_new.sh

2.markdup 
mark duplicates in the bam file

script
3.stats-beforedups.sh

4.stats-mkdups.sh


#ok now that I have all the reads mapped and duplicates removed, I need to create the mpileup file with all samples

```bash 
find /gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/aligned_new/mkdups/ -name "*.bam" > bam-files-new.txt 
```

5. to-mpileup.sh
```bash
#!/bin/bash
#SBATCH -D /gxfs_home/geomar/smomw573/work/seasonal_adaptation/log_files/
#SBATCH --mail-type=END
#SBATCH --mail-user=jnascimento@geomar.de
#SBATCH --partition=base
#SBATCH --cpus-per-task=1
#SBATCH --mem=60G
#SBATCH --time=10:00:00
#SBATCH --job-name=to_mpileup
#SBATCH --output=%x.%A.out
#SBATCH --error=%x.%A.err

source ~/miniconda3/bin/activate bwa

cd $WORK/smomw573/seasonal_adaptation/analysis/variants

samtools mpileup -q 15 -Q 0 -d 8000 -R -A -B 
-f $WORK/smomw573/seasonal_adaptation/new_genome/Purple_Eaffinis_genome_assembly.fa 
-b $WORK/smomw573/seasonal_adaptation/analysis/variants/bam-files-new.txt  
-o all-new.mpileup

````

ok, now I need to call the SNPs.. 
I will first do this to all my samples. SNPs will be called with grenedalf
first I need to install GNU parallel
628687434 variants in mpileup
```bash
conda create --name parallel python=3.9 conda-forge::parallel
```
6.poolsnp-all-samples.sh
```bash
#!/bin/bash
#SBATCH -D /gxfs_home/geomar/smomw573/work/seasonal_adaptation/log_files/
#SBATCH --mail-type=END
#SBATCH --mail-user=jnascimento@geomar.de
#SBATCH --partition=base
#SBATCH --nodes=1
#SBATCH --tasks-per-node=17
#SBATCH --cpus-per-task=1
#SBATCH --mem=100G
#SBATCH --qos=long
#SBATCH --time=100:00:00
#SBATCH --job-name=PoolSNP-allsamples
#SBATCH --output=%x.%A.out
#SBATCH --error=%x.%A.err


source ~/miniconda3/bin/activate parallel

bash /gxfs_home/geomar/smomw573/software/PoolSNP/PoolSNP.sh \
mpileup=$WORK/smomw573/seasonal_adaptation/analysis/variants/all-new.mpileup \
output=$WORK/smomw573/seasonal_adaptation/analysis/PoolSNP/snps-all-samples-new.vcf  \
reference=$WORK/smomw573/seasonal_adaptation/new_genome/Purple_Eaffinis_genome_assembly.fa \
names=EA_2007_T1,EA_2007_T2,EA_2009_T1,EA_2009_T2,EA_2009_T3,EA_2009_T4,EA_2011_T1,EA_2011_T2,EA_2011_T3,EA_2015_T1,EA_2015_T2,EA_2015_T3,EA_2015_T4,EA_2022_T1,EA_2022_T2,EA_2022_T3,EA_2022_T4 \
min-cov=50 \
max-cov=0.95 \
min-count=10 \
min-freq=0.01 \
miss-frac=0.1 \
base-quality 15 \
jobs=17 \
badsites=1 \
allsites=0

```

*job was ready in ~ 14 hours. So I will use 25 next time
I only got 93 SNPs out, so will rerun everything.
I will try min-cov of 10 instead of 50
6.2.poolsnp-all-samples.sh
#
```bash
#!/bin/bash
#SBATCH -D /gxfs_home/geomar/smomw573/work/seasonal_adaptation/log_files/
#SBATCH --mail-type=END
#SBATCH --mail-user=jnascimento@geomar.de
#SBATCH --partition=base
#SBATCH --nodes=1
#SBATCH --tasks-per-node=17
#SBATCH --cpus-per-task=1
#SBATCH --mem=100G
#SBATCH --qos=long
#SBATCH --time=100:00:00
#SBATCH --job-name=PoolSNP-allsamples-2
#SBATCH --output=%x.%A.out
#SBATCH --error=%x.%A.err


source ~/miniconda3/bin/activate parallel

bash /gxfs_home/geomar/smomw573/software/PoolSNP/PoolSNP.sh \
mpileup=$WORK/smomw573/seasonal_adaptation/analysis/variants/all-new.mpileup \
output=$WORK/smomw573/seasonal_adaptation/analysis/PoolSNP/snps-all-samples-new-2.vcf  \
reference=$WORK/smomw573/seasonal_adaptation/new_genome/Purple_Eaffinis_genome_assembly.fa \
names=EA_2007_T1,EA_2007_T2,EA_2009_T1,EA_2009_T2,EA_2009_T3,EA_2009_T4,EA_2011_T1,EA_2011_T2,EA_2011_T3,EA_2015_T1,EA_2015_T2,EA_2015_T3,EA_2015_T4,EA_2022_T1,EA_2022_T2,EA_2022_T3,EA_2022_T4 \
min-cov=10 \
max-cov=0.95 \
min-count=10 \
min-freq=0.01 \
miss-frac=0.1 \
base-quality 15 \
jobs=17 \
badsites=1 \
allsites=0

```
and also will min-count=5 and min cov=50
6.3.poolsnp-all-samples.sh
```bash
#!/bin/bash
#SBATCH -D /gxfs_home/geomar/smomw573/work/seasonal_adaptation/log_files/
#SBATCH --mail-type=END
#SBATCH --mail-user=jnascimento@geomar.de
#SBATCH --partition=base
#SBATCH --nodes=1
#SBATCH --tasks-per-node=17
#SBATCH --cpus-per-task=1
#SBATCH --mem=100G
#SBATCH --qos=long
#SBATCH --time=100:00:00
#SBATCH --job-name=PoolSNP-allsamples-3
#SBATCH --output=%x.%A.out
#SBATCH --error=%x.%A.err


source ~/miniconda3/bin/activate parallel

bash /gxfs_home/geomar/smomw573/software/PoolSNP/PoolSNP.sh \
mpileup=$WORK/smomw573/seasonal_adaptation/analysis/variants/all-new.mpileup \
output=$WORK/smomw573/seasonal_adaptation/analysis/PoolSNP/snps-all-samples-new-3.vcf  \
reference=$WORK/smomw573/seasonal_adaptation/new_genome/Purple_Eaffinis_genome_assembly.fa \
names=EA_2007_T1,EA_2007_T2,EA_2009_T1,EA_2009_T2,EA_2009_T3,EA_2009_T4,EA_2011_T1,EA_2011_T2,EA_2011_T3,EA_2015_T1,EA_2015_T2,EA_2015_T3,EA_2015_T4,EA_2022_T1,EA_2022_T2,EA_2022_T3,EA_2022_T4 \
min-cov=50 \
max-cov=0.95 \
min-count=5 \
min-freq=0.01 \
miss-frac=0.1 \
base-quality 15 \
jobs=17 \
badsites=1 \
allsites=0

```
I will also try to produce the coverage file to add as an input on max-cov
```bash
conda activate bwa
samtools depth -a -f $WORK/smomw573/seasonal_adaptation/analysis/variants/bam-files-new.txt > $WORK/smomw573/seasonal_adaptation/analysis/variants/coverage-per-sample.txt
```
#ok I will do a run without the 2007 samples as they are too bad (very low alignment rates)
for this I need to create a new mpileup file...

5.2 to-mpileup-no2007.sh
```bash
#!/bin/bash
#SBATCH -D /gxfs_home/geomar/smomw573/work/seasonal_adaptation/log_files/
#SBATCH --mail-type=END
#SBATCH --mail-user=jnascimento@geomar.de
#SBATCH --partition=base
#SBATCH --cpus-per-task=1
#SBATCH --mem=60G
#SBATCH --time=24:00:00
#SBATCH --job-name=to_mpileup2
#SBATCH --output=%x.%A.out
#SBATCH --error=%x.%A.err

source ~/miniconda3/bin/activate bwa

cd $WORK/smomw573/seasonal_adaptation/analysis/variants

samtools mpileup -q 15 -Q 0 -d 8000 -R -A -B \
        -f $WORK/smomw573/seasonal_adaptation/new_genome/Purple_Eaffinis_genome_assembly.fa \
        -b $WORK/smomw573/seasonal_adaptation/analysis/variants/bam-files-new-no2007.txt \
        -o all-new-no2007.mpileup

````


6.4.poolsnp-no2007.sh
```bash
#!/bin/bash
#SBATCH -D /gxfs_home/geomar/smomw573/work/seasonal_adaptation/log_files/
#SBATCH --mail-type=END
#SBATCH --mail-user=jnascimento@geomar.de
#SBATCH --partition=base
#SBATCH --nodes=1
#SBATCH --tasks-per-node=17
#SBATCH --cpus-per-task=1
#SBATCH --mem=100G
#SBATCH --qos=long
#SBATCH --time=100:00:00
#SBATCH --job-name=PoolSNP-excluding2007
#SBATCH --output=%x.%A.out
#SBATCH --error=%x.%A.err


source ~/miniconda3/bin/activate parallel

bash /gxfs_home/geomar/smomw573/software/PoolSNP/PoolSNP.sh \
mpileup=$WORK/smomw573/seasonal_adaptation/analysis/variants/all-new.mpileup \
output=$WORK/smomw573/seasonal_adaptation/analysis/PoolSNP/snps-excluding2007-new.vcf  \
reference=$WORK/smomw573/seasonal_adaptation/new_genome/Purple_Eaffinis_genome_assembly.fa \
names=EA_2009_T1,EA_2009_T2,EA_2009_T3,EA_2009_T4,EA_2011_T1,EA_2011_T2,EA_2011_T3,EA_2015_T1,EA_2015_T2,EA_2015_T3,EA_2015_T4,EA_2022_T1,EA_2022_T2,EA_2022_T3,EA_2022_T4 \
min-cov=50 \
max-cov=0.95 \
min-count=5 \
min-freq=0.01 \
miss-frac=0.1 \
base-quality 15 \
jobs=17 \
badsites=1 \
allsites=0

```
#149696 snps

# check missing data per individual
```bash
# Get all sample names as a comma-separated list
samples=$(bcftools query -l combined.vcf.gz | tr '\n' ',' | sed 's/,$//')

# Query genotypes for all samples, then compute missingness per sample
bcftools query -s "$samples" -f '[%GT\t]\n' combined.vcf.gz | \
awk -v OFS="\t" '
{
  for (i=1; i<=NF; i++) if ($i == "./.") missing[i]++
  total++
}
END {
  for (i=1; i<=NF; i++) print i, missing[i]/total
}' > missingness_per_sample.txt

# Combine with sample names for readability:
paste <(bcftools query -l combined.vcf.gz) missingness_per_sample.txt | cut -f1,3 > sample_missingness.txt
```
# ok, 2011_T3 had 99% missing data. 

# I will rerun PoolSNP without this sample to see if we increase the number of SNPs being called.

# first create the new mpileup

5.4 to-mpileup-exc20072011t3.sh
```bash

#!/bin/bash
#SBATCH -D /gxfs_home/geomar/smomw573/work/seasonal_adaptation/log_files/
#SBATCH --mail-type=END
#SBATCH --mail-user=jnascimento@geomar.de
#SBATCH --partition=base
#SBATCH --cpus-per-task=1
#SBATCH --mem=60G
#SBATCH --time=25:00:00
#SBATCH --job-name=to_mpileup20072011
#SBATCH --output=%x.%A.out
#SBATCH --error=%x.%A.err

source ~/miniconda3/bin/activate bwa
cd $WORK/smomw573/seasonal_adaptation/analysis/variants

samtools mpileup -q 15 -Q 0 -d 8000 -R -A -B \
-f $WORK/smomw573/seasonal_adaptation/new_genome/Purple_Eaffinis_genome_assembly.fa \
-b $WORK/smomw573/seasonal_adaptation/analysis/variants/bam-files-new-no2007-2011t3.txt \
-o snps-exc20072011t3.mpileup


````
# and run poolsnp

6.6.poolsnp-no20072011t3.sh

```bash
#!/bin/bash
#SBATCH -D /gxfs_home/geomar/smomw573/work/seasonal_adaptation/log_files/
#SBATCH --mail-type=END
#SBATCH --mail-user=jnascimento@geomar.de
#SBATCH --partition=base
#SBATCH --nodes=1
#SBATCH --tasks-per-node=17
#SBATCH --cpus-per-task=1
#SBATCH --mem=100G
#SBATCH --qos=long
#SBATCH --time=50:00:00
#SBATCH --job-name=PoolSNP-excluding2007
#SBATCH --output=%x.%A.out
#SBATCH --error=%x.%A.err


source ~/miniconda3/bin/activate parallel

bash /gxfs_home/geomar/smomw573/software/PoolSNP/PoolSNP.sh \
mpileup=$WORK/smomw573/seasonal_adaptation/analysis/variants/snps-exc20072011t3.mpileup \
output=$WORK/smomw573/seasonal_adaptation/analysis/PoolSNP/snps-excluding2007-2011t3-new.vcf  \
reference=$WORK/smomw573/seasonal_adaptation/new_genome/Purple_Eaffinis_genome_assembly.fa \
names=EA_2009_T1,EA_2009_T2,EA_2009_T3,EA_2009_T4,EA_2011_T1,EA_2011_T2,EA_2015_T1,EA_2015_T2,EA_2015_T3,EA_2015_T4,EA_2022_T1,EA_2022_T2,EA_2022_T3,EA_2022_T4 \
min-cov=30 \
max-cov=0.95 \
min-count=5 \
min-freq=0.01 \
miss-frac=0.1 \
base-quality=15 \
jobs=14 \
badsites=1 \
allsites=0

```
# 4768624 snps. about double than last time

```bash 
# now let's filter biallelic snps and look at AFs! :D

vcftools --gzvcf snps-excluding2007-2011t3-new.vcf.vcf.gz --min-alleles 2 --max-alleles 2 --recode --out snps-excluding2007-2011t3-new-biallelic

# fixed with python script

# now I need to run greenedalf

cat header-excluding2007and2011-new.txt snps-excluding2007-2011t3-new-biallelic-fixed.recode.vcf > final-excl2007-2011-new.vcf

# tabix file

#tabix the vcf
bgzip final-excl2007-2011-new.vcf
tabix -p vcf final-excl2007-2011-new.vcf.gz

# filtering missing data 


vcftools --gzvcf final-excl2007-2011-new.vcf.gz  --max-missing 1.0  --recode --out final-excl2007-2011-new-no-missing-data.vcf.gz
#kept 138555 

#output genotype depths
vcftools --gzvcf final-excl2007-2011-new-no-missing-data.vcf.gz.recode.vcf --site-mean-depth --out final-excl2007-2011-new-no-missing-data

# Fix vcf for coverage depth
# select SNPs which only have 3x average coverage
#conda vcf env
module load R/4.3.1
```

```r
install.packages("data.table")
install.packages("poolfstat")
library(data.table)
library(poolfstat)
#trial plot
plot(1,1)

#load data 
df <- read.table("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/grenedalf-freq/new-vcf/final-excl2007-2011-new-no-missing-data.ldepth.mean", header = T)
head(df)
hist(df$MEAN_DEPTH)

#check mean and median

mean <- mean(df$MEAN_DEPTH)
69.71
mean3 = mean*3
209.1419

sum(df$MEAN_DEPTH > mean3)
1264
sum(df$MEAN_DEPTH > mean3)/nrow(df)
0.009122731
#not much will be excluded! That is great :D
head(sort(df$MEAN_DEPTH), 10)
36

#poolfstat

#remove all markers with a depth over than three times the mean
dat <- vcf2pooldata(vcf.file="/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/grenedalf-freq/new-vcf/final-excl2007-2011-new-no-missing-data.vcf.gz.recode.vcf",poolsizes=c(rep(100,14)),
                        min.cov.per.pool = 50,max.cov.per.pool=209.1419,min.maf=0.01,nlines.per.readblock=1000000)

# if we use 50/209 there are 726 snps left..

pooldata2genobaypass(dat,writing.dir="/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/grenedalf-freq/vcf-maf0.01-misfrac0.1-maxcov95/")
#ouput is snpdet
```

# something is wrong with 2009 t1

# check missing data per sample

# 16.07.2025

```bash
bcftools query -f '%CHROM\t%POS[\t%AD]\n' final-excl2007-2011-new-no-missing-data.vcf.gz.recode.vcf > allele_depths.txt

vcftools --gzvcf final-excl2007-2011-new.vcf.gz --missing-site --out vcf-with-missing-data
vcftools --gzvcf final-excl2007-2011-new.vcf.gz --missing-ind --out vcf-with-missing-data
#check the bam files
samtools depth -a /gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/aligned_new/mkdups/EA_2009_T2.bam | awk '{sum+=$3} END {print "Average coverage = ", sum/NR}' > EA_2009_T2_coverage.txt
#compare to last analysis
samtools depth -a /gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/aligned/mkdups/EA_2009_T2.bam | awk '{sum+=$3} END {print "Average coverage = ", sum/NR}' > EA_2009_T2_coverage-old.txt
samtools depth -a /gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/aligned/mkdups/EA_2009_T1.bam | awk '{sum+=$3} END {print "Average coverage = ", sum/NR}' > EA_2009_T1_coverage-old.txt

samtools depth -a /gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/aligned_new/mkdups/EA_2009_T2.bam > EA_2009_T2.depth
samtools depth -a /gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/aligned_new/mkdups/EA_2009_T1.bam > EA_2009_T1.depth
module load R/4.3.1
```
```R
# Read in individual missingness file
imiss <- read.table("vcf-with-missing-data.imiss", header=TRUE)

# Check contents
head(imiss)
# Barplot of fraction missing
barplot(imiss$F_MISS, names.arg=imiss$INDV, las=2, col="tomato",
        ylab="Fraction Missing", main="Missing Data per Individual")

lmiss <- read.table("vcf-with-missing-data.lmiss", header=TRUE)

h = hist(lmiss$F_MISS)
h$density = h$counts/sum(h$counts)
plot(h, freq = FALSE)
#check samtools output
# Load the data
cov <- read.table("EA_2009_T1.depth", header=FALSE)
colnames(cov) <- c("CHROM", "POS", "DEPTH")
#Basic coverage histogram (all sites):

hist(cov$DEPTH, breaks=100, col="skyblue", main="Coverage Distribution",
     xlab="Depth per site", ylab="Number of sites",
     xlim=c(0, 2000))
#Coverage per chromosome (optional density plot):
library(ggplot2)

ggplot(cov, aes(x=DEPTH)) +
  geom_density(fill="tomato", alpha=0.5) +
  xlim(0, quantile(cov$DEPTH, 0.99)) + # Avoid long tail
  labs(title="Coverage Density", x="Depth", y="Density")
#Mean coverage per chromosome (optional):

aggregate(DEPTH ~ CHROM, data=cov, FUN=mean)
```
```text 
# ok, so 2009 has loads of missing data
to do list
missing data per individual - figure out which sample is driving. 
depths?
histogram less zoomed? 


ok, given all this, I will remove 2009_t1 from analysis, make a new mpileup file, and call snps again. 
I will also rerun the last poolsnp call and make a new vcf file without a max cov limit
finally, I will remove all old files as I can regenerate them at anypoint and they are just confusing.
5.5 final tompileup
and  6.7.poolsnp-no20072011t3-nomaxcovfilter.sh

# 18.07.25 
call snps again using poolsnp with new mpileup snps-exc20072011t32009t1.mpileup 6.8

# average coverage bam files

for bam in *.bam
do
    avg=$(samtools depth "$bam" | awk '{sum+=$3} END {if (NR>0) print sum/NR; else print 0}')
    echo "$bam average coverage: $avg"
done
```
```bash
#22.07

#I created multiple vcfs, excluding 2007, 2011t3 and 2009t1, using different max-cov filters to check if there is anything weird going on. The decrease in number of snps is linear, 

check Ti/Tv ratio of each of the new vcfs
max-cov 0.9999
max-cov 0.99
max-cov 0.98
max-cov 0.95



vcftools --gzvcf snps-excluding2007-2011t3-2009t1-new-maxcov95.vcf.vcf.gz --min-alleles 2 --max-alleles 2 --recode --out snps-excluding2007-2011t3-2009t1-new-maxcov95_biallelic
vcftools --gzvcf snps-excluding2007-2011t3-2009t1-new-maxcov98.vcf.vcf.gz --min-alleles 2 --max-alleles 2 --recode --out snps-excluding2007-2011t3-2009t1-new-maxcov98_biallelic
vcftools --gzvcf snps-excluding2007-2011t3-2009t1-new-maxcov99.vcf.vcf.gz --min-alleles 2 --max-alleles 2 --recode --out snps-excluding2007-2011t3-2009t1-new-maxcov99_biallelic
vcftools --gzvcf snps-excluding2007-2011t3-2009t1-new-maxcov9999.vcf.vcf.gz --min-alleles 2 --max-alleles 2 --recode --out snps-excluding2007-2011t3-2009t1-new-maxcov9999_biallelic


bcftools stats snps-excluding2007-2011t3-2009t1-new-maxcov95_biallelic.recode.vcf.gz > stats_95.txt
bcftools stats snps-excluding2007-2011t3-2009t1-new-maxcov98_biallelic.recode.vcf.gz > stats_98.txt
bcftools stats snps-excluding2007-2011t3-2009t1-new-maxcov99_biallelic.recode.vcf.gz > stats_99.txt
bcftools stats snps-excluding2007-2011t3-2009t1-new-maxcov9999_biallelic.recode.vcf.gz > stats_9999.txt

check coverage
bcftools query -f '%CHROM\t%POS\t%INFO/DP\n' snps-excluding2007-2011t3-2009t1-new-maxcov9999_biallelic.recode.vcf.gz > depth_raw_9999.tsv
bcftools query -f '%CHROM\t%POS\t%INFO/DP\n' snps-excluding2007-2011t3-2009t1-new-maxcov99_biallelic.recode.vcf.gz > depth_raw_99.tsv
bcftools query -f '%CHROM\t%POS\t%INFO/DP\n' snps-excluding2007-2011t3-2009t1-new-maxcov98_biallelic.recode.vcf.gz > depth_raw_98.tsv
bcftools query -f '%CHROM\t%POS\t%INFO/DP\n' snps-excluding2007-2011t3-2009t1-new-maxcov95_biallelic.recode.vcf.gz > depth_raw_95.tsv
#not much changed, but I will stick with >98 
#fixed files with python script

#fix header

zcat head9999.txt.gz maxcov9999-biallelic-fixed-new.recode.vcf.gz > maxcov9999.vcf.gz
zcat head99.txt.gz maxcov99-biallelic-fixed-new.recode.vcf.gz.gz > maxcov99.vcf.gz
zcat head98.txt.gz maxcov98-biallelic-fixed-new.recode.vcf.gz.gz > maxcov98.vcf.gz
zcat head95.txt.gz maxcov95-biallelic-fixed-new.recode.vcf.gz > maxcov95.vcf.gz

#check missing data

vcftools --gzvcf maxcov9999.vcf.gz --missing-site --out maxcov9999
vcftools --gzvcf maxcov9999.vcf.gz --missing-indv --out maxcov9999
vcftools --gzvcf maxcov99.vcf.gz --missing-indv --out maxcov99
vcftools --gzvcf maxcov98.vcf.gz --missing-indv --out maxcov98
vcftools --gzvcf maxcov95.vcf.gz --missing-indv --out maxcov95


#remove missing data
vcftools --gzvcf maxcov9999.vcf.gz --max-missing 1.0  --recode --out maxcov9999-final.vcf.gz
vcftools --gzvcf maxcov99.vcf.gz --max-missing 1.0  --recode --out maxcov99-final.vcf.gz
vcftools --gzvcf maxcov98.vcf.gz --max-missing 1.0  --recode --out maxcov98-final.vcf.gz
vcftools --gzvcf maxcov95.vcf.gz --max-missing 1.0  --recode --out maxcov95-final.vcf.gz

#ok, for now I will continue the analysis with max-cov 99 

#output genotype depths
vcftools --gzvcf maxcov99-final.vcf.gz.recode.vcf --site-mean-depth --out maxcov99-final

module load R/4.3.1
```
```R
install.packages("data.table")
install.packages("poolfstat")
library(data.table)
library(poolfstat)
#trial plot
plot(1,1)

#load data 
df <- read.table("maxcov99-final.ldepth.mean", header = T)
head(df)
hist(df$MEAN_DEPTH)
summary(df$MEAN_DEPTH)
#not so clear

hist(df$MEAN_DEPTH[df$MEAN_DEPTH <= 500],
     breaks=100,
     col="steelblue",
     main="Mean Depth (≤1000)",
     xlab="Mean Depth",
     ylab="Number of Sites")


#    CHROM                POS             MEAN_DEPTH        VAR_DEPTH
# Length:7575561     Min.   :       6   Min.   :  33.46   Min.   :     8.3
# Class :character   1st Qu.:10602177   1st Qu.:  54.31   1st Qu.:   185.1
# Mode  :character   Median :21468691   Median :  60.69   Median :   265.2
#                    Mean   :22309613   Mean   :  65.78   Mean   :   353.1
#                    3rd Qu.:32842121   3rd Qu.:  73.31   3rd Qu.:   417.0
#                    Max.   :63905152   Max.   :4491.62   Max.   :882536.0


#check mean and median

mean <- mean(df$MEAN_DEPTH)
65.78
mean3 = mean*3
197.3478

sum(df$MEAN_DEPTH > mean_value)
2819400
sum(df$MEAN_DEPTH > mean3)
2412
sum(df$MEAN_DEPTH > mean3)/nrow(df)
0.0003183923

#not much will be excluded! That is great :D
head(sort(df$MEAN_DEPTH), 10)
33

#poolfstat

#remove all markers with a depth over than three times the mean
dat <- vcf2pooldata(vcf.file="maxcov99-final.vcf.gz.recode.vcf",poolsizes=c(rep(100,13)),
                        min.cov.per.pool = 50,max.cov.per.pool=197.3478,min.maf=0.01,nlines.per.readblock=1000000)

# 578435 snps

pooldata2genobaypass(dat,writing.dir="/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/PoolSNP/trial_vcfs/")
#ouput is snpdet

#same for max-cov 98
df <- read.table("maxcov98-final.ldepth.mean", header = T)
head(df)
hist(df$MEAN_DEPTH)
summary(df$MEAN_DEPTH)


hist(df$MEAN_DEPTH[df$MEAN_DEPTH <= 500],
     breaks=100,
     col="steelblue",
     main="Mean Depth (≤1000)",
     xlab="Mean Depth",
     ylab="Number of Sites")

mean <- mean(df$MEAN_DEPTH)
60.98
mean3 = mean*3
182.9401


dat <- vcf2pooldata(vcf.file="maxcov98-final.vcf.gz.recode.vcf",poolsizes=c(rep(100,13)),
                        min.cov.per.pool = 50,max.cov.per.pool=182.9401,min.maf=0.01,nlines.per.readblock=1000000)

# 96522

pooldata2genobaypass(dat,writing.dir="/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/PoolSNP/trial_vcfs/")
#all these snps are also present in max-cov99 out, which has 8k snps more
```

```bash
#correct for depth
awk '{print $1"\t"$2}' snpdet-99 > snp-positions-99.txt

# filter this list with vcf
#screen session
srun --pty --x11 --nodes=1 --cpus-per-task=4 --mem=16000 --time=02:00:00 /bin/bash

vcftools --vcf maxcov99-final.vcf.gz.recode.vcf --positions snp-positions-99.txt --recode --out maxcov99-depth-corrected

#After filtering, kept  578434 snps

#removing MAF errors
bgzip maxcov98-depth-corrected.recode.vcf
tabix -p vcf maxcov98-depth-corrected.recode.vcf.gz

#run grenedalf
module load htslib/1.10.2   
module load bzip2/1.0.8

#format line in vcf file
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">

GRENEDALF=/gxfs_home/geomar/smomw573/software/grenedalf/bin/grenedalf
VCF=/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/PoolSNP/trial_vcfs/maxcov99-depth-corrected.recode.vcf.gz
GENOME=$WORK/smomw573/seasonal_adaptation/new_genome/Purple_Eaffinis_genome_assembly.fa \
OUTPUT=/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/grenedalf-freq/new-genome-maxcov99-no200720092011/

$GRENEDALF frequency --vcf-path $VCF --reference-genome-fasta-file $GENOME --write-sample-alt-freq --allow-file-overwriting --file-suffix maxcov99-sample

#wooohoooo

```
Plot PCA
```R
library(data.table)
library(ggplot2)
df <- read.csv("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/grenedalf-freq/new-genome-maxcov99-no200720092011/frequencymaxcov99-sample.csv", header = T)
head(df)
#df <- noquote(t(frequency_ind_samples[,-(2:4)]))

#colnames(df) <- df[1, ]

# Remove the first row as it's now used for column names
df <- df[,-(2:4)]
#test data
#df2 <- df[0:15,1:14]
#df2
#all okay

pops <- c("EA_2009_T2.FREQ", "EA_2009_T3.FREQ","EA_2009_T4.FREQ", 
        "EA_2011_T1.FREQ", "EA_2011_T2.FREQ", 
        "EA_2015_T1.FREQ", "EA_2015_T2.FREQ", "EA_2015_T3.FREQ", "EA_2015_T4.FREQ", 
        "EA_2022_T1.FREQ", "EA_2022_T2.FREQ", "EA_2022_T3.FREQ", "EA_2022_T4.FREQ") 
varout <- apply(df[, 2:ncol(df)], 1, var)

freqs <- t(df[varout != 0,2:ncol(df)])
nrow(freqs)
ncol(freqs)
#13 

## plot pca

pcaResult <- prcomp(freqs, scale=T)
percentVar <- round(pcaResult$sdev^2/sum(pcaResult$sdev^2)*100, digits=2)

data <- data.frame(
  id = pops,
  Line = substr(pops, 4, 7),
  gen = substr(pops, 9,10),
  PC1 = pcaResult$x[, 1],
  PC2 = pcaResult$x[, 2],
  PC3 = pcaResult$x[, 3]
)

#data$Line <- c(rep("Founding population", 4),
#                rep("Ambient", 4),
#                rep("Acidification", 4),
#                rep("Warming", 4),
#                rep("OWA", 4))

data$Line <- factor(data$Line, levels = c("2009","2011", "2015", "2022"))
data$gen <- factor(data$gen, levels = c("T1","T2", "T3", "T4"))
data$PC2 <- data$PC2*-1

d <- ggplot(data, aes(PC1, PC2, fill=Line, shape=gen)) +
        geom_point(size=4.5) +
        xlab(paste0("PC1: ",percentVar[1],"% variance")) +
        ylab(paste0("PC2: ",percentVar[2],"% variance")) +
        theme_bw() +
       # ylim(-30, 23) + xlim(-50, 65)+
        scale_shape_manual(values=c( 21,22,23,24))+
        scale_color_manual(values=c('black')) +
        #scale_fill_manual(values=c('steelblue1','steelblue','grey45', "darkorchid2", "firebrick3"),
  scale_fill_manual(values=c("#D3DDDC",'#6699CC',"#F2AD00","#00A08A"),
                    labels = c("2009","2011", "2015",
                               "2022"))+
        #theme(legend.position = c(0.83,0.85),
        #    legend.background = element_blank(),
        #legend.box.background = element_rect(colour = "black"),
        theme(legend.title = element_blank()) +
        theme(legend.text=element_text(size=16),
          axis.title.x = element_text(size = 18),
          axis.text.x= element_text(size=18),
          axis.text.y= element_text(size=18),   # Adjust x-axis label size
          axis.title.y = element_text(size = 18))+
        #theme(legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid'))+
#        ggtitle("F1")+
guides(fill=guide_legend(override.aes=list(
        shape=c(21,21,22, 23),
        #fill=c('steelblue1','steelblue','grey45', "darkorchid2", "firebrick3" )),order = 2),
    fill=c("#D3DDDC",'#6699CC',"#F2AD00","#00A08A")),order = 2))

ggsave("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/figures/new-genome-maxcov99-_pca_minorafs_allsamples.pdf",d, w=6, h=4)

#exclude 2015_T3

df <- df[,-c(9)]
pops <- c("EA_2009_T2.FREQ", "EA_2009_T3.FREQ","EA_2009_T4.FREQ", 
        "EA_2011_T1.FREQ", "EA_2011_T2.FREQ", 
        "EA_2015_T1.FREQ", "EA_2015_T2.FREQ",  "EA_2015_T4.FREQ", 
        "EA_2022_T1.FREQ", "EA_2022_T2.FREQ", "EA_2022_T3.FREQ", "EA_2022_T4.FREQ") 
varout <- apply(df[, 2:ncol(df)], 1, var)

freqs <- t(df[varout != 0,2:ncol(df)])
nrow(freqs)
ncol(freqs)
#13 

## plot pca

pcaResult <- prcomp(freqs, scale=T)
percentVar <- round(pcaResult$sdev^2/sum(pcaResult$sdev^2)*100, digits=2)

data <- data.frame(
  id = pops,
  Line = substr(pops, 4, 7),
  gen = substr(pops, 9,10),
  PC1 = pcaResult$x[, 1],
  PC2 = pcaResult$x[, 2],
  PC3 = pcaResult$x[, 3]
)

#data$Line <- c(rep("Founding population", 4),
#                rep("Ambient", 4),
#                rep("Acidification", 4),
#                rep("Warming", 4),
#                rep("OWA", 4))

data$Line <- factor(data$Line, levels = c("2009","2011", "2015", "2022"))
data$gen <- factor(data$gen, levels = c("T1","T2", "T3", "T4"))
data$PC2 <- data$PC2*-1

d <- ggplot(data, aes(PC1, PC2, fill=Line, shape=gen)) +
        geom_point(size=4.5) +
        xlab(paste0("PC1: ",percentVar[1],"% variance")) +
        ylab(paste0("PC2: ",percentVar[2],"% variance")) +
        theme_bw() +
       # ylim(-30, 23) + xlim(-50, 65)+
        scale_shape_manual(values=c( 21,22,23,24))+
        scale_color_manual(values=c('black')) +
        #scale_fill_manual(values=c('steelblue1','steelblue','grey45', "darkorchid2", "firebrick3"),
  scale_fill_manual(values=c("#D3DDDC",'#6699CC',"#F2AD00","#00A08A"),
                    labels = c("2009","2011", "2015",
                               "2022"))+
        #theme(legend.position = c(0.83,0.85),
        #    legend.background = element_blank(),
        #legend.box.background = element_rect(colour = "black"),
        theme(legend.title = element_blank()) +
        theme(legend.text=element_text(size=18),
          axis.title.x = element_text(size = 18),
          axis.text.x= element_text(size=18),
          axis.text.y= element_text(size=18),   # Adjust x-axis label size
          axis.title.y = element_text(size = 18))+
        #theme(legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid'))+
#        ggtitle("F1")+
guides(fill=guide_legend(override.aes=list(
        shape=c(21,21,22, 23),
        #fill=c('steelblue1','steelblue','grey45', "darkorchid2", "firebrick3" )),order = 2),
    fill=c("#D3DDDC",'#6699CC',"#F2AD00","#00A08A")),order = 2))

ggsave("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/figures/new-genome-maxcov99-_pca_minorafs_excl-2015_t3.pdf",d, w=6, h=4)

```

repeat for 98

```bash
#correct for depth
awk '{print $1"\t"$2}' snpdet-98 > snp-positions-98.txt

# filter this list with vcf
#screen session
srun --pty --x11 --nodes=1 --cpus-per-task=4 --mem=16000 --time=02:00:00 /bin/bash

vcftools --vcf maxcov98-final.vcf.gz.recode.vcf --positions snp-positions-98.txt --recode --out maxcov98-depth-corrected

#After filtering, kept  578434 snps

#removing MAF errors
bgzip maxcov99-depth-corrected.recode.vcf
tabix -p vcf maxcov99-depth-corrected.recode.vcf.gz

#run grenedalf
module load htslib/1.10.2   
module load bzip2/1.0.8

#format line in vcf file
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">

GRENEDALF=/gxfs_home/geomar/smomw573/software/grenedalf/bin/grenedalf
VCF=/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/PoolSNP/trial_vcfs/maxcov99-depth-corrected.recode.vcf.gz
GENOME=$WORK/smomw573/seasonal_adaptation/new_genome/Purple_Eaffinis_genome_assembly.fa \
OUTPUT=/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/grenedalf-freq/new-genome-maxcov99-no200720092011/

$GRENEDALF frequency --vcf-path $VCF --reference-genome-fasta-file $GENOME --write-sample-alt-freq --allow-file-overwriting --file-suffix maxcov99-sample

#wooohoooo

```
Plot PCA
```R
library(data.table)
library(ggplot2)
df <- read.csv("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/grenedalf-freq/new-genome-maxcov99-no200720092011/frequencymaxcov99-sample.csv", header = T)
head(df)
#df <- noquote(t(frequency_ind_samples[,-(2:4)]))

#colnames(df) <- df[1, ]

# Remove the first row as it's now used for column names
df <- df[,-(2:4)]
#test data
#df2 <- df[0:15,1:14]
#df2
#all okay

pops <- c("EA_2009_T2.FREQ", "EA_2009_T3.FREQ","EA_2009_T4.FREQ", 
        "EA_2011_T1.FREQ", "EA_2011_T2.FREQ", 
        "EA_2015_T1.FREQ", "EA_2015_T2.FREQ", "EA_2015_T3.FREQ", "EA_2015_T4.FREQ", 
        "EA_2022_T1.FREQ", "EA_2022_T2.FREQ", "EA_2022_T3.FREQ", "EA_2022_T4.FREQ") 
varout <- apply(df[, 2:ncol(df)], 1, var)

freqs <- t(df[varout != 0,2:ncol(df)])
nrow(freqs)
ncol(freqs)
#13 

## plot pca

pcaResult <- prcomp(freqs, scale=T)
percentVar <- round(pcaResult$sdev^2/sum(pcaResult$sdev^2)*100, digits=2)

data <- data.frame(
  id = pops,
  Line = substr(pops, 4, 7),
  gen = substr(pops, 9,10),
  PC1 = pcaResult$x[, 1],
  PC2 = pcaResult$x[, 2],
  PC3 = pcaResult$x[, 3]
)

#data$Line <- c(rep("Founding population", 4),
#                rep("Ambient", 4),
#                rep("Acidification", 4),
#                rep("Warming", 4),
#                rep("OWA", 4))

data$Line <- factor(data$Line, levels = c("2009","2011", "2015", "2022"))
data$gen <- factor(data$gen, levels = c("T1","T2", "T3", "T4"))
data$PC2 <- data$PC2*-1

d <- ggplot(data, aes(PC1, PC2, fill=Line, shape=gen)) +
        geom_point(size=4.5) +
        xlab(paste0("PC1: ",percentVar[1],"% variance")) +
        ylab(paste0("PC2: ",percentVar[2],"% variance")) +
        theme_bw() +
       # ylim(-30, 23) + xlim(-50, 65)+
        scale_shape_manual(values=c( 21,22,23,24))+
        scale_color_manual(values=c('black')) +
        #scale_fill_manual(values=c('steelblue1','steelblue','grey45', "darkorchid2", "firebrick3"),
  scale_fill_manual(values=c("#D3DDDC",'#6699CC',"#F2AD00","#00A08A"),
                    labels = c("2009","2011", "2015",
                               "2022"))+
        #theme(legend.position = c(0.83,0.85),
        #    legend.background = element_blank(),
        #legend.box.background = element_rect(colour = "black"),
        theme(legend.title = element_blank()) +
        theme(legend.text=element_text(size=16),
          axis.title.x = element_text(size = 18),
          axis.text.x= element_text(size=18),
          axis.text.y= element_text(size=18),   # Adjust x-axis label size
          axis.title.y = element_text(size = 18))+
        #theme(legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid'))+
#        ggtitle("F1")+
guides(fill=guide_legend(override.aes=list(
        shape=c(21,21,22, 23),
        #fill=c('steelblue1','steelblue','grey45', "darkorchid2", "firebrick3" )),order = 2),
    fill=c("#D3DDDC",'#6699CC',"#F2AD00","#00A08A")),order = 2))

ggsave("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/figures/new-genome-maxcov99-_pca_minorafs_allsamples.pdf",d, w=6, h=4)

#exclude 2015_T3

df <- df[,-c(9)]
pops <- c("EA_2009_T2.FREQ", "EA_2009_T3.FREQ","EA_2009_T4.FREQ", 
        "EA_2011_T1.FREQ", "EA_2011_T2.FREQ", 
        "EA_2015_T1.FREQ", "EA_2015_T2.FREQ",  "EA_2015_T4.FREQ", 
        "EA_2022_T1.FREQ", "EA_2022_T2.FREQ", "EA_2022_T3.FREQ", "EA_2022_T4.FREQ") 
varout <- apply(df[, 2:ncol(df)], 1, var)

freqs <- t(df[varout != 0,2:ncol(df)])
nrow(freqs)
ncol(freqs)
#13 

## plot pca

pcaResult <- prcomp(freqs, scale=T)
percentVar <- round(pcaResult$sdev^2/sum(pcaResult$sdev^2)*100, digits=2)

data <- data.frame(
  id = pops,
  Line = substr(pops, 4, 7),
  gen = substr(pops, 9,10),
  PC1 = pcaResult$x[, 1],
  PC2 = pcaResult$x[, 2],
  PC3 = pcaResult$x[, 3]
)

#data$Line <- c(rep("Founding population", 4),
#                rep("Ambient", 4),
#                rep("Acidification", 4),
#                rep("Warming", 4),
#                rep("OWA", 4))

data$Line <- factor(data$Line, levels = c("2009","2011", "2015", "2022"))
data$gen <- factor(data$gen, levels = c("T1","T2", "T3", "T4"))
data$PC2 <- data$PC2*-1

d <- ggplot(data, aes(PC1, PC2, fill=Line, shape=gen)) +
        geom_point(size=4.5) +
        xlab(paste0("PC1: ",percentVar[1],"% variance")) +
        ylab(paste0("PC2: ",percentVar[2],"% variance")) +
        theme_bw() +
       # ylim(-30, 23) + xlim(-50, 65)+
        scale_shape_manual(values=c( 21,22,23,24))+
        scale_color_manual(values=c('black')) +
        #scale_fill_manual(values=c('steelblue1','steelblue','grey45', "darkorchid2", "firebrick3"),
  scale_fill_manual(values=c("#D3DDDC",'#6699CC',"#F2AD00","#00A08A"),
                    labels = c("2009","2011", "2015",
                               "2022"))+
        #theme(legend.position = c(0.83,0.85),
        #    legend.background = element_blank(),
        #legend.box.background = element_rect(colour = "black"),
        theme(legend.title = element_blank()) +
        theme(legend.text=element_text(size=18),
          axis.title.x = element_text(size = 18),
          axis.text.x= element_text(size=18),
          axis.text.y= element_text(size=18),   # Adjust x-axis label size
          axis.title.y = element_text(size = 18))+
        #theme(legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid'))+
#        ggtitle("F1")+
guides(fill=guide_legend(override.aes=list(
        shape=c(21,21,22, 23),
        #fill=c('steelblue1','steelblue','grey45', "darkorchid2", "firebrick3" )),order = 2),
    fill=c("#D3DDDC",'#6699CC',"#F2AD00","#00A08A")),order = 2))

ggsave("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/figures/new-genome-maxcov99-_pca_minorafs_excl-2015_t3.pdf",d, w=6, h=4)
```
Okay, so after all this we decided that the coverage (max-cov) filter was not good enough, as it is still doing something funny and removing loci that we might not want to remove.
So we decided to re run the analysis, without using a max-cov filter file. For this, I manually created a max-cov file (in variants), which allows the max cov in all chromossomes to be 10K X. This is really high coverage considering that in the previous analysis, the highest we've seen in the files is 4K X. 
Very good, I reran mpileup and poolsnp using this file as max-cov using the script 6.10
the output was this vcf file snps-excluding2007-2011t3-2009t1-2015t3-new.vcf.vcf.gz with 30961006 variants
We also excluded 2015_T3 from the analysis based on the high percentage of missing data


```bash
#fix vcf


vcftools --gzvcf snps-excluding2007-2011t3-2009t1-2015t3-new.vcf.vcf.gz --min-alleles 2 --max-alleles 2 --recode --out snps-excluding2007-2011t3-2009t1-2015t3-new_biallelic

bcftools stats snps-excluding2007-2011t3-2009t1-2015t3-new_biallelic.recode.vcf.gz > stats_nomaxcov.txt

check coverage
bcftools query -f '%CHROM\t%POS\t%INFO/DP\n' snps-excluding2007-2011t3-2009t1-2015t3-new_biallelic.recode.vcf.gz > depth_raw_nomaxcov.tsv

#remove missing data
vcftools --gzvcf snps-excluding2007-2011t3-2009t1-2015t3-new_biallelic.recode.vcf.gz  --max-missing 1.0  --recode --out final-nomaxcov

cat head-final.txt final-nomaxcov-fixed.vcf > no2007-2009-2011-2015-fixed.vcf


#check missing data
vcftools --gzvcf maxcov9999.vcf.gz --missing-site --out maxcov9999
vcftools --gzvcf maxcov9999.vcf.gz --missing-indv --out maxcov9999

#output genotype depths
vcftools --gzvcf no2007-2009-2011-2015-fixed.vcf --site-mean-depth --out no2007-2009-2011-2015-fixed

```
```R
#same for max-cov 98
df <- read.table("no2007-2009-2011-2015-fixed.ldepth.mean", header = T)
head(df)
hist(df$MEAN_DEPTH)
summary(df$MEAN_DEPTH)


hist(df$MEAN_DEPTH[df$MEAN_DEPTH <= 1000],
     breaks=100,
     col="steelblue",
     main="Mean Depth (≤1000)",
     xlab="Mean Depth",
     ylab="Number of Sites")

mean <- mean(df$MEAN_DEPTH)
68.30948
mean3 = mean*3
204.9284


dat <- vcf2pooldata(vcf.file="no2007-2009-2011-2015-fixed.vcf",poolsizes=c(rep(100,12)),
                        min.cov.per.pool = 50,max.cov.per.pool=204.9284,min.maf=0.01,nlines.per.readblock=1000000)

# 96522

pooldata2genobaypass(dat,writing.dir="/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/PoolSNP/new-excluding2007-2009-2011-2015/")
```
```bash
#correct for depth
awk '{print $1"\t"$2}' snpdet > snp-positions.txt

# filter this list with vcf
#screen session
srun --pty --x11 --nodes=1 --cpus-per-task=4 --mem=16000 --time=02:00:00 /bin/bash

vcftools --vcf no2007-2009-2011-2015-fixed.vcf --positions snp-positions.txt --recode --out no2007-2009-2011-2015-fixed-depth-corrected

#After filtering, kept  1874896  snps
#Perfect!

#run grenedalf
module load hstlib/1.10.2   
module load bzip2/1.0.8

#format line in vcf file
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">

GRENEDALF=/gxfs_home/geomar/smomw573/software/grenedalf/bin/grenedalf
VCF=/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/PoolSNP/new-excluding2007-2009-2011-2015/no2007-2009-2011-2015-fixed-depth-corrected.recode.vcf.gz
GENOME=$WORK/smomw573/seasonal_adaptation/new_genome/Purple_Eaffinis_genome_assembly.fa \
OUTPUT=/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/grenedalf-freq/new-genome-maxcov99-no2007200920112015/

$GRENEDALF frequency --vcf-path $VCF --reference-genome-fasta-file $GENOME --write-sample-alt-freq --allow-file-overwriting --file-suffix altfreq-no2007-2009-2011-2015-depth-corrected

#wooohoooo
frequencyaltfreq-no2007-2009-2011-2015-depth-corrected.csv
```

now I will plot the PCA

```R
library(data.table)
library(ggplot2)
df <- read.csv("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/PoolSNP/new-excluding2007-2009-2011-2015/frequencyaltfreq-no2007-2009-2011-2015-depth-corrected.csv", header = T)
head(df)
#df <- noquote(t(frequency_ind_samples[,-(2:4)]))

#colnames(df) <- df[1, ]

# Remove the first row as it's now used for column names
df <- df[,-(2:4)]
#test data
#df2 <- df[0:15,1:14]
#df2
#all okay

pops <- c("EA_2009_T2.FREQ", "EA_2009_T3.FREQ","EA_2009_T4.FREQ", 
        "EA_2011_T1.FREQ", "EA_2011_T2.FREQ", 
        "EA_2015_T1.FREQ", "EA_2015_T2.FREQ", "EA_2015_T4.FREQ", 
        "EA_2022_T1.FREQ", "EA_2022_T2.FREQ", "EA_2022_T3.FREQ", "EA_2022_T4.FREQ") 
varout <- apply(df[, 2:ncol(df)], 1, var)

freqs <- t(df[varout != 0,2:ncol(df)])
nrow(freqs)
ncol(freqs)
#12 

## plot pca

pcaResult <- prcomp(freqs, scale=T)
percentVar <- round(pcaResult$sdev^2/sum(pcaResult$sdev^2)*100, digits=2)

data <- data.frame(
  id = pops,
  Line = substr(pops, 4, 7),
  gen = substr(pops, 9,10),
  PC1 = pcaResult$x[, 1],
  PC2 = pcaResult$x[, 2],
  PC3 = pcaResult$x[, 3]
)

#data$Line <- c(rep("Founding population", 4),
#                rep("Ambient", 4),
#                rep("Acidification", 4),
#                rep("Warming", 4),
#                rep("OWA", 4))

data$Line <- factor(data$Line, levels = c("2009","2011", "2015", "2022"))
data$gen <- factor(data$gen, levels = c("T1","T2", "T3", "T4"))
data$PC2 <- data$PC2*-1

d <- ggplot(data, aes(PC1, PC2, fill=Line, shape=gen)) +
        geom_point(size=4.5) +
        xlab(paste0("PC1: ",percentVar[1],"% variance")) +
        ylab(paste0("PC2: ",percentVar[2],"% variance")) +
        theme_bw() +
       # ylim(-30, 23) + xlim(-50, 65)+
        scale_shape_manual(values=c( 21,22,23,24))+
        scale_color_manual(values=c('black')) +
        #scale_fill_manual(values=c('steelblue1','steelblue','grey45', "darkorchid2", "firebrick3"),
  scale_fill_manual(values=c("#D3DDDC",'#6699CC',"#F2AD00","#00A08A"),
                    labels = c("2009","2011", "2015",
                               "2022"))+
        #theme(legend.position = c(0.83,0.85),
        #    legend.background = element_blank(),
        #legend.box.background = element_rect(colour = "black"),
        theme(legend.title = element_blank()) +
        theme(legend.text=element_text(size=16),
          axis.title.x = element_text(size = 18),
          axis.text.x= element_text(size=18),
          axis.text.y= element_text(size=18),   # Adjust x-axis label size
          axis.title.y = element_text(size = 18))+
        #theme(legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid'))+
#        ggtitle("F1")+
guides(fill=guide_legend(override.aes=list(
        shape=c(21,21,22, 23),
        #fill=c('steelblue1','steelblue','grey45', "darkorchid2", "firebrick3" )),order = 2),
    fill=c("#D3DDDC",'#6699CC',"#F2AD00","#00A08A")),order = 2))

ggsave("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/figures/new-genome-no2007-2009-2011-2015.pdf",d, w=6, h=4)

#exclude 2015_T3

df <- df[,-c(9)]
pops <- c("EA_2009_T2.FREQ", "EA_2009_T3.FREQ","EA_2009_T4.FREQ", 
        "EA_2011_T1.FREQ", "EA_2011_T2.FREQ", 
        "EA_2015_T1.FREQ", "EA_2015_T2.FREQ",  "EA_2015_T4.FREQ", 
        "EA_2022_T1.FREQ", "EA_2022_T2.FREQ", "EA_2022_T3.FREQ", "EA_2022_T4.FREQ") 
varout <- apply(df[, 2:ncol(df)], 1, var)

freqs <- t(df[varout != 0,2:ncol(df)])
nrow(freqs)
ncol(freqs)
#13 

## plot pca

pcaResult <- prcomp(freqs, scale=T)
percentVar <- round(pcaResult$sdev^2/sum(pcaResult$sdev^2)*100, digits=2)

data <- data.frame(
  id = pops,
  Line = substr(pops, 4, 7),
  gen = substr(pops, 9,10),
  PC1 = pcaResult$x[, 1],
  PC2 = pcaResult$x[, 2],
  PC3 = pcaResult$x[, 3]
)

#data$Line <- c(rep("Founding population", 4),
#                rep("Ambient", 4),
#                rep("Acidification", 4),
#                rep("Warming", 4),
#                rep("OWA", 4))

data$Line <- factor(data$Line, levels = c("2009","2011", "2015", "2022"))
data$gen <- factor(data$gen, levels = c("T1","T2", "T3", "T4"))
data$PC2 <- data$PC2*-1

d <- ggplot(data, aes(PC1, PC2, fill=Line, shape=gen)) +
        geom_point(size=4.5) +
        xlab(paste0("PC1: ",percentVar[1],"% variance")) +
        ylab(paste0("PC2: ",percentVar[2],"% variance")) +
        theme_bw() +
       # ylim(-30, 23) + xlim(-50, 65)+
        scale_shape_manual(values=c( 21,22,23,24))+
        scale_color_manual(values=c('black')) +
        #scale_fill_manual(values=c('steelblue1','steelblue','grey45', "darkorchid2", "firebrick3"),
  scale_fill_manual(values=c("#D3DDDC",'#6699CC',"#F2AD00","#00A08A"),
                    labels = c("2009","2011", "2015",
                               "2022"))+
        #theme(legend.position = c(0.83,0.85),
        #    legend.background = element_blank(),
        #legend.box.background = element_rect(colour = "black"),
        theme(legend.title = element_blank()) +
        theme(legend.text=element_text(size=18),
          axis.title.x = element_text(size = 18),
          axis.text.x= element_text(size=18),
          axis.text.y= element_text(size=18),   # Adjust x-axis label size
          axis.title.y = element_text(size = 18))+
        #theme(legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid'))+
#        ggtitle("F1")+
guides(fill=guide_legend(override.aes=list(
        shape=c(21,21,22, 23),
        #fill=c('steelblue1','steelblue','grey45', "darkorchid2", "firebrick3" )),order = 2),
    fill=c("#D3DDDC",'#6699CC',"#F2AD00","#00A08A")),order = 2))

ggsave("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/figures/new-genome-maxcov99-_pca_minorafs_excl-2015_t3.pdf",d, w=6, h=4)
```

pca was done with scaled data. Need to rerun it for unscaled.
okay, now I am ready to move to the next step.

29.07.2025

- check Fst between populations

```bash
#rename chromossomes so I can run grenedalf
awk '
{
  if (match($0, /^(Chr)([0-9]+)$/, arr)) {
    prefix = arr[1];
    num = arr[2];
    newnum = sprintf("%02d", num);  # 2 digits for Chr
    print $0 "\t" prefix newnum;
  } else if (match($0, /^(Scaffold)([0-9]+)$/, arr)) {
    prefix = arr[1];
    num = arr[2];
    newnum = sprintf("%03d", num);  # 3 digits for Scaffold
    print $0 "\t" prefix newnum;
  } else {
    print $0 "\t" $0;
  }
}' ~/work/seasonal_adaptation/new_genome/chrom-names-new-genome > chr_rename.txt

bcftools annotate --rename-chrs chr_rename.txt -Oz -o renamed.vcf.gz no2007-2009-2011-2015-fixed-depth-corrected.recode.vcf.gz
tabix -p vcf renamed.vcf.gz
bcftools sort -Oz -o final-sorted.vcf.gz renamed.vcf.gz
tabix -p vcf final-sorted.vcf.gz

#run grenedalf
module load htslib/1.10.2    
module load bzip2/1.0.8 

GRENEDALF=/gxfs_home/geomar/smomw573/software/grenedalf/bin/grenedalf
VCF=/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/grenedalf-freq/new-genome-final-vcf/final-sorted.vcf.gz
GENOME=$WORK/smomw573/seasonal_adaptation/new_genome/Purple_Eaffinis_genome_assembly.fa 
OUTPUT=/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/grenedalf-freq/new-genome-final-vcf

#fst_site queue
$GRENEDALF fst \
  --vcf-path $VCF \
  --method unbiased-hudson \
  --allow-file-overwriting \
  --pool-sizes 200 \
  --window-type queue \
  --window-queue-count 1 \
  --omit-na-windows \
  --file-suffix _site \
  --out-dir $OUTPUT/

#fst_site window
$GRENEDALF fst \
  --vcf-path $VCF \
  --method unbiased-hudson \
  --allow-file-overwriting \
  --pool-sizes 200 \
  --window-type sliding \
  --window-sliding-stride 0 \
  --window-sliding-width 1 \
  --omit-na-windows \
  --file-suffix _site_window \
  --out-dir $OUTPUT/
 #ok, these two commands do the exact same thing!

  #fst_window
  $GRENEDALF fst \
  --vcf-path $VCF \
  --method unbiased-hudson \
  --allow-file-overwriting \
  --pool-sizes 200 \
  --window-type sliding \
  --window-sliding-width 5000 \
  --window-sliding-stride 5000 \
  --omit-na-windows \
  --file-suffix _window_5kb \
  --out-dir $OUTPUT/
```

output are two files
fst_window_5kb.csv
fst_site.csv

```R
#plot fst values in R
library(dplyr)
library(ggplot2)
library(tidyr)
library(reshape2)
library(gtools)
library(forcats)
library(corrplot)

fst <- read.csv("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/grenedalf-freq/new-genome-final-vcf/fst_site.csv", header = TRUE, na.strings = "")
#negative values to zeros

fst <- fst %>% mutate_all(function(x) ifelse(x < 0, 0, x))

#means of each colums (are there some higher then others)
# corrplot
mean_fst <- colMeans(fst[, 5:ncol(fst)], na.rm = TRUE)
mean_fst 

# Extract only the Fst columns (assuming first 4 columns are CHROM, POS, etc.)
fst_vals <- fst[, 5:ncol(fst)]

# Calculate column means, ignoring NA values
fst_means <- colMeans(fst_vals, na.rm = TRUE)

# Melt into long format
fst_long <- melt(fst_means)
fst_long$id <- row.names(fst_long)

# Split comparison names into two populations (e.g., "A.B" → "A", "B")
fst_split <- fst_long %>% separate_wider_delim(id, delim = ".", names = c("id1", "id2"))

# Reshape into matrix form (population1 × population2)
fst_matrix <- reshape2::acast(fst_split, id1 ~ id2, value.var = "value")

pdf(file="/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/grenedalf-freq/new-genome-final-vcf/figures/fst_site_final.pdf", h=8, w=8)
corrplot(fst_matrix, method = "color", type = "upper", is.corr = FALSE,
         addCoef.col = "black", number.digits = 4,
         col = colorRampPalette(c("gray95", "firebrick3"))(200),
         tl.col = "black", cl.pos = "n", tl.srt = 45)
dev.off()

####ok saved now another plot

#manhattan plot
# Prepare chromosome ordering
fst$chrom <- factor(fst$chrom, levels = mixedsort(unique(fst$chrom)))

# Compute cumulative base pair position
fst <- fst %>%
  arrange(chrom, start) %>%
  group_by(chrom) %>%
  mutate(chr_len = max(end)) %>%
  ungroup() %>%
  mutate(tot = cumsum(as.numeric(chr_len)) - chr_len) %>%
  mutate(BPcum = start + tot)

# Optional: axis labels centered per chromosome ===
axis_df <- fst %>%
  group_by(chrom) %>%
  summarize(center = (min(BPcum) + max(BPcum)) / 2)

# Melt data to long format for plotting (if multiple comparisons)
fst_long <- fst %>%
  pivot_longer(cols = 5:(ncol(fst)-3), names_to = "comparison", values_to = "Fst")
  # Optional: ensure correct chromosome order
fst_long$chrom <- factor(fst_long$chrom, levels = gtools::mixedsort(unique(fst_long$chrom)))

# Axis label positions per chromosome
axis_df <- fst_long %>%
  group_by(chrom) %>%
  summarize(center = (min(BPcum) + max(BPcum)) / 2)

# Manhattan-style plot with all comparisons
p <- ggplot(fst_filtered, aes(x = BPcum, y = Fst, color = chrom)) +
  geom_point(size = 0.3, alpha = 0.75) +
  facet_wrap(~comparison, ncol = 2) +
  scale_x_continuous(label = axis_df$chrom, breaks = axis_df$center) +
  scale_y_continuous(name = "Fst") +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(angle = 90, size = 5, vjust = 0.5),
    strip.text = element_text(size = 8)
  ) +
  labs(title = "Genome-wide Fst Manhattan Plot", x = "Chromosome", y = "Fst")

# Show plot
print(p)

# Save plot
ggsave("fst_manhattan_all_comparisons.png", plot = p, width = 12, height = 10)

# Plot
# === Manhattan-style FST plot ===
ggplot(fst_long, aes(x = BPcum, y = Fst, color = chrom)) +
  geom_point(alpha = 0.5, size = 0.5) +
  scale_x_continuous(label = axis_df$chrom, breaks = axis_df$center) +
  scale_color_manual(values = rep(c("black", "darkred"), length(unique(fst$chrom)))) +
  facet_wrap(~comparison, ncol = 1) +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, size = 6)
  ) +
  labs(
    x = "Chromosome",
    y = "Fst",
    title = "Manhattan Plot of Per-SNP Fst Values"
  )

ggplot(fst, aes_string(x = "BPcum", y = comp, color = "chrom")) +
    geom_point(size = 0.3, alpha = 0.75) +
    scale_x_continuous(label = axis_df$chrom, breaks = axis_df$center) +
    scale_y_continuous(name = "Fst") +
    theme_bw() +
    theme(
      legend.position = "none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.text.x = element_text(angle = 90, size = 6, vjust = 0.5)
    ) +
    ggtitle(paste("Fst:", comp))
  
# Save each plot as PDF named after the comparison
ggsave(filename = paste0("fst_manhattan_", comp, ".png"), plot = p, width = 10, height = 4)


#plot pair-wise comparisons individually
# Loop over all pairwise comparison columns (from 5th column to last)
comparison_cols <- colnames(fst)[5:(ncol(fst) - 1)]  # excluding cum_len and BPcum if present

for (comp in comparison_cols) {
  p <- ggplot(fst, aes_string(x = "BPcum", y = comp, color = "chrom")) +
    geom_point(size = 0.3, alpha = 0.75) +
    scale_x_continuous(label = axis_df$chrom, breaks = axis_df$center) +
    scale_y_continuous(name = "Fst") +
    theme_bw() +
    theme(
      legend.position = "none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.text.x = element_text(angle = 90, size = 6, vjust = 0.5)
    ) +
    ggtitle(paste("Fst:", comp))
  
  # Save each plot as PDF named after the comparison
  ggsave(filename = paste0("fst_manhattan_", comp, ".png"), plot = p, width = 10, height = 4)
}

#now I will do the comparison scatter plot with only the important comparisons
#for this I will still use the fst_long data used for the manhattan plot! 
# Assign "Line" based on the first timepoint in the comparison
fst_long <- fst_long %>%
  mutate(Line = case_when(
    grepl("^EA_2009", comparison) ~ "2009",
    grepl("^EA_2011", comparison) ~ "2011",
    grepl("^EA_2015", comparison) ~ "2015",
    grepl("^EA_2022", comparison) ~ "2022"
  ))

# Define the desired order of comparisons (x-axis)
comparisons_of_interest <- c(
  "EA_2009_T2.EA_2009_T3",
  "EA_2009_T2.EA_2009_T4",
  "EA_2009_T3.EA_2009_T4",
  "EA_2011_T1.EA_2011_T2",
  "EA_2015_T1.EA_2015_T2",
  "EA_2015_T1.EA_2015_T4",
  "EA_2015_T2.EA_2015_T4",
  "EA_2022_T1.EA_2022_T2",
  "EA_2022_T1.EA_2022_T3",
  "EA_2022_T1.EA_2022_T4",
  "EA_2022_T2.EA_2022_T3",
  "EA_2022_T2.EA_2022_T4",
  "EA_2022_T3.EA_2022_T4"
)

# Ensure comparisons are ordered
fst_long$comparison <- factor(fst_long$comparison, levels = comparisons_of_interest)

# Your PCA colors
line_colors <- c(
  "2009" = "#D3DDDC",
  "2011" = "#6699CC",
  "2015" = "#F2AD00",
  "2022" = "#00A08A"
)
library(stringr)
# Optional: remove missing values or zeros
fst_filtered <- fst_long %>%
  filter(comparison %in% comparisons_of_interest, !is.na(Fst) & Fst > 0)

fst_filtered_final <- fst_subset %>%
  mutate(Line = str_extract(comparison, "20\\d{2}"))

fst_subset <- fst_filtered %>%
  group_by(comparison) %>%
  sample_n(size = min(1000, n()), replace = FALSE) %>%
  ungroup()

line_colors <- c("2009"="#D3DDDC", "2011"="#6699CC", "2015"="#F2AD00", "2022"="#00A08A")

# Extract year from comparison for coloring
fst_subset <- fst_subset %>%
  mutate(Line = str_extract(comparison, "20\\d{2}"))

# Create and assign the plot
d <- ggplot(fst_subset, aes(x = comparison, y = Fst, color = comparison)) +
  geom_jitter(width = 0.3, size = 0.6, alpha = 0.6) +
  scale_color_manual(values = line_colors) +
  theme_bw() +
  labs(title = "Fst single snps (subset)", x = NULL, y = expression(F[ST]), color = "Year") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank()
  )

# Show plot
print(d)

# Plot
d <- ggplot(fst_filtered_final, aes(x = comparison, y = Fst, color = comparison)) +
  geom_jitter(width = 0.3, size = 0.6, alpha = 0.6) +
  scale_color_manual(values = line_colors) +
  theme_bw() +
  labs(title = "Fst single snps", x = NULL, y = expression(F[ST]), color = "Year") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank()
  )

# Display the plot in R
print(d)

# Save the plot to file
ggsave("fst_per_comparison_colored.png", plot = d, width = 10, height = 5)

#old approach
#boxplot to check variation between each pairwise comparison, this way we can see which ones are higher
fst2 <- fst[,-(1:4)]
fst2 <- as.data.frame(fst2)
#remove missing data

#count NAs
sumnas <- sum(is.na(fst2))
#2994703 NAS

# Count NA values in each column
column_na_counts <- colSums(is.na(fst2))

# Print the results
print("Total NA Count:")
print(total_na_count)
print("NA Counts in Each Column:")
print(column_na_counts)

### remove mising data
data_clean <- fst2
data_clean <- na.omit(data_clean)

# Reshape the data into long format
data_long <- tidyr::gather(data_clean, key = "Column", value = "Value")
data_long$Column <- as.factor(data_long$Column)

# Create boxplots using ggplot2
d <- ggplot(data_long, aes(x=Column, y=Value)) + 
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

 ggsave("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/figures/new-genome-boxplotfst.pdf",d, w=7.5, h=4.7)
#8160 rows removed with missing values


#subset for pairwise comparisons of interest

column_names <- colnames(fst2)
print(column_names)
 

fst_subset<- fst2[, c("EA_2009_T1.EA_2009_T2", "EA_2009_T1.EA_2009_T3", "EA_2009_T1.EA_2009_T4", "EA_2009_T2.EA_2009_T3", "EA_2009_T2.EA_2009_T4","EA_2009_T3.EA_2009_T4","EA_2011_T1.EA_2011_T2",  "EA_2015_T1.EA_2015_T2", "EA_2015_T1.EA_2015_T3", 
        "EA_2015_T1.EA_2015_T4", "EA_2015_T2.EA_2015_T3", "EA_2015_T2.EA_2015_T4", "EA_2015_T3.EA_2015_T4","EA_2022_T1.EA_2022_T2",
        "EA_2022_T1.EA_2022_T3", "EA_2022_T1.EA_2022_T4", "EA_2022_T2.EA_2022_T3", "EA_2022_T2.EA_2022_T4", "EA_2022_T3.EA_2022_T4")]


#create new variable year

fst_subset <- as.data.frame(fst_subset)
data_clean <- na.omit(data_clean)
data_long <- tidyr::gather(data_clean, key = "Column", value = "Value") # gather() function from the tidyr package to convert your dataset to long format.
data_long$Column <- as.factor(data_long$Column)

names(data_long)[1] <- "Column"
year<- rep (NA, length(data_long$Column))
year[grep("2009", data_long$Column)] <- "2009"
year[grep("2011", data_long$Column)] <- "2011"
year[grep("2015", data_long$Column)] <- "2015"
year[grep("2022", data_long$Column)] <- "2022"

data_long$year <- year


#plot

d <- ggplot(data_long, aes(x=Column, y=Value)) + 
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

 ggsave("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/figures/boxplot-fst-main-comparisons.pdf",d, w=12, h=4.7)

#removed 1821 rows


#scater plots

# Create a scatter plot grid using facet_grid
d <- ggplot(data_long, aes(x = Column, y = Value)) +
  geom_point() + #~ Column indicates that each facet corresponds to a different column.
  facet_grid(. ~ year, scales = "free") + # The scales = "free" argument allows each facet to have its own y-axis scale.
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Pairwise comparison", y = "Fst")


 ggsave("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/figures/boxplot-fst-scatter-main-comparisons.pdf",d, w=12, h=6)
```

# 31.07.25

# now I want to estimate the genome wide Fst and plot it
# for this, I am using the poolfstat R package

```bash
srun --pty --x11 --nodes=1 --cpus-per-task=4 --mem=32000 --time=06:00:00 /bin/bash
module load R/4.1.1
module load gcc/10.2.0
module load icu4c/67.1

R
```

enter R environment 


```R
library(poolfstat)
#poolfstat for vcf

library(poolfstat)

#import vcf file
ea.readcount50X <- vcf2pooldata(vcf.file="final-sorted.vcf.gz",poolsizes=rep(200,12))
# summary of the resulting pooldata object
ea.readcount50X 
#ha?
selected.snps.idx <- as.numeric(sub("rs","",rownames(ea.readcount50X@snp.info)))
head(selected.snps.idx)

#estimate genome wide Fst over all popuatlions
ea.readcount50X.fst<-computeFST(ea.readcount50X)
ea.readcount50X.fst$FST
#genome wide Fst is -0.003898684, so 0

# Block-Jackknife estimation of FST standard-error and confidence intervals:
ea.readcount50X.fst<-computeFST(ea.readcount50X,nsnp.per.bjack.block = 1000, verbose=TRUE)
ea.readcount50X.fst$FST
#same value  

ea.readcount50X.fst$mean.fst #block-jacknife estimate of s.e.
#-0.003902521
ea.readcount50X.fst$se.fst #s.e. of the genome-wide Fst estimate
#1.855665e-05
ea.readcount50X.fst$mean.fst+c(-1.96,1.96)*ea.readcount50X.fst$se.fst
# -0.02015146 -0.02003757


#Computing multi-locus FST to scan the genome over sliding-windows of SNPs
ea.readcount50X.fst<-computeFST(ea.readcount50X,sliding.window.size=50)
ea.readcount50X.fst<-computeFST(ea.readcount50X,sliding.window.size=100)
ea.readcount50X.fst<-computeFST(ea.readcount50X,sliding.window.size=10)
#we have 42K scaffolds so loosing a lot of data with sliding window size of 50-100 chrom?

#Average (min-max) Window Sizes 17.2 ( 0.1 - 249.7 ) kb
#Average (min-max) Window Sizes 34.8 ( 0.2 - 294 ) kb
#Average (min-max) Window Sizes 3.1 ( 0 - 235.3 ) kb

df <- ea.readcount50X.fst$sliding.windows.fst
plot(df$CumulatedPosition/1e6, df$MultiLocusFst,
     xlab="Cumulated Position (in Mb)",
     ylab="Multi-locus Fst",
     col=as.integer(factor(df$Chr)),  # one color per chromosome
     pch=16)

abline(h=ea.readcount50X.fst$FST, lty=2)

plot(ea.readcount50X.fst$sliding.windows.fst$CumulatedPosition/1e6,
     ea.readcount50X.fst$sliding.windows.fst$MultiLocusFst,
     xlab="Cumulated Position (in Mb)",ylab="Muli-locus Fst",
     col=as.numeric(ea.readcount50X.fst$sliding.windows.fst$Chr),pch=16) 
abline(h=ea.readcount50X.fst$FST,lty=2)

#pairwise FST

ea.pairwisefst<-compute.pairwiseFST(ea.readcount50X,verbose=FALSE)
#heatmap
#Heatmap representing the pairwise-population FST matrix of the 14 populations of the 30XPool-Seq example data set
heatmap(ea.pairwisefst)
#
ea.pairwisefst<-compute.pairwiseFST(ea.pairwise.allelecount,min.maf=0.01,
nsnp.per.bjack.block = 1000,verbose=FALSE)
#it moves pops which are more similar to each other 
```
Generate alt frequency file using grenedalf

```bash
module load htslib/1.10.2  bzip2/1.0.8  

GRENEDALF=/gxfs_home/geomar/smomw573/software/grenedalf/bin/grenedalf
GENOME=/gxfs_work1/geomar/smomw573/seasonal_adaptation/genome/affinis.Atlantic.long_read_draft.Mar22.fasta \
VCF=/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/variants/final_20230810.vcf.gz \
OUTPUT=/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/grenedalf-freq/vcf-maf0.01-misfrac0.1-maxcov95/grenedalfst.log/freq-finalvcf/

$GRENEDALF frequency --vcf-path $VCF --allow-file-overwriting --write-sample-alt-freq --file-suffix sample-frequencies --out-dir $OUTPUT/grenedalftrial.log
```
s
#the fst plot with the new fst dataset fixed for the new genome is on a different file - plot-single-fst.md in my computer
the Fst values are really low between the populations.

 to do list

Reid Brennan
3:41 PM
Jump
just from my notes:

Analysis todo:

From called snps

PCA
Covariance
Mixed models as in drosophila papers
Fst within and across years
From pileup:

pi, waterson’s theta, tajima’s D.
Might need to mask the GC mutations in the vcf. will have to see.