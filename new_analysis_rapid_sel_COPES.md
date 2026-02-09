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
ok I will do a run without the 2007 samples as they are too bad (very low alignment rates)
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
srun --pty --x11 --nodes=1 --cpus-per-task=4 --mem=32000 --time=10:00:00 /bin/bash

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
So we decided to rerun the analysis, without using a max-cov filter file. For this, I manually created a max-cov file (in variants), which allows the max cov in all chromossomes to be 10K X. This is really high coverage considering that in the previous analysis, the highest we've seen in the files is 4K X. 
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

#output genotype depths
vcftools --gzvcf no2007-2009-2011-2015-fixed.vcf --site-mean-depth --out no2007-2009-2011-2015-fixed

```
```R
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
module load htslib/1.10.2    
module load bzip2/1.0.8 

GRENEDALF=/gxfs_home/geomar/smomw573/software/grenedalf/bin/grenedalf
VCF=/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/grenedalf-freq/new-genome-final-vcf/final-sorted.vcf
GENOME=$WORK/smomw573/seasonal_adaptation/new_genome/Purple_Eaffinis_genome_assembly.fa 
OUTPUT=/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/grenedalf-freq/new-genome-final-vcf

$GRENEDALF frequency --vcf-path $VCF  --write-sample-alt-freq --allow-file-overwriting --separator-char tab --file-suffix _alt --out-dir $OUTPUT/
$GRENEDALF frequency --vcf-path $VCF  --write-sample-coverage --allow-file-overwriting --separator-char tab --file-suffix _cov --out-dir $OUTPUT/
$GRENEDALF frequency --vcf-path $VCF --write-total-frequency --allow-file-overwriting --file-suffix _maf --out-dir $OUTPUT
$GRENEDALF frequency --vcf-path $VCF --write-sample-counts --allow-file-overwriting --file-suffix _counts --out-dir $OUTPUT

```
create different files in R - fix data


```R
library(tidyr)

freq <- read.table("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/grenedalf-freq/new-genome-final-vcf/frequency_alt.csv", header=T, sep="\t")

sites <- freq[,c(1,2)]
afmat <- freq[,c(-1:-4)]
afmat <- as.matrix(afmat)
head(afmat)
samps <- read.table("samples.txt", header = TRUE, sep = "\t")
objects <- c("sites", "afmat", "samps")
save(list = "objects", file = "HAFs.Rdata")

data <- freq
#fixing data#

# Remove ".FREQ" from the column names
# Identify the columns containing ".COV" in their names
cov_columns <- grep("EA_\\d+_T\\d+\\.FREQ", names(data))
new_colnames <- colnames(data)
new_colnames[cov_columns] <- gsub("\\.FREQ", "", new_colnames[cov_columns])
colnames(data) <- new_colnames

# Create a data frame with only the "CHROM," "POS," "REF," "ALT," and renamed ".COV" columns
data_FREQ <- data[, c(5:16)]
write.table(data_FREQ, file = "freq_new_vcf.txt", quote = FALSE, sep = "\t")

#cov data
data <- read.table("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/grenedalf-freq/new-genome-final-vcf/frequency_cov.csv", header=T, sep="\t")

# Remove ".COV" from the column names
# Identify the columns containing ".COV" in their names
cov_columns <- grep("EA_\\d+_T\\d+\\.COV", names(data))
new_colnames <- colnames(data)
new_colnames[cov_columns] <- gsub("\\.COV", "", new_colnames[cov_columns])
colnames(data) <- new_colnames

# Create a data frame with only the "CHROM," "POS," "REF," "ALT," and renamed ".COV" columns
data_FREQ <- data[, c(5:16)]

write.table(data_FREQ, file = "cov_new_vcf.txt", quote = FALSE, sep = "\t")
```
get positions
zcat final-sorted.vcf.gz | awk '!/^#/ {print $1"\t"$2}' > chrom_pos.txt

script for GLM
```R
####GLM script trial 1####

cov <- read.table("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/grenedalf-freq/new-genome-final-vcf/cov_new_vcf.txt", header = TRUE, sep = "\t")
freq <- read.table("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/grenedalf-freq/new-genome-final-vcf/freq_new_vcf.txt", header = TRUE, sep = "\t")
pos <- read.table("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/grenedalf-freq/new-genome-final-vcf/chrom-poly-position-seasonal-data.txt", header = FALSE, sep = "\t")
popinfo <- read.table("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/grenedalf-freq/new-genome-final-vcf/popinfo.txt", header = TRUE, sep = "\t")


#subset populations

popnames <- grep("EA_2009_T2|EA_2009_T4|EA_2011_T1|EA_2011_T2|EA_2015_T1|EA_2015_T4|EA_2022_T1|EA_2022_T4", popinfo$pop)

#subset
cov_glm <- cov[,popnames]
freq_glm <- freq[,popnames]
popinfo_glm <- popinfo[popnames,]
popinfo_glm$time <- c("E","L")

#transform into matrix
freq_matrix <- as.matrix(freq_glm)
cov_matrix <- as.matrix(cov_glm)

#fix cov_matrix with Nc

# 50 so 2N = 100
#make new column Nc for the weight
#Nc = (1/N + 1/R) - 1

dp <- (1/100 + 1/cov_matrix)^-1

popinfo_glm$Y <- as.factor(popinfo_glm$Y)

#with year as factor
fileout <- file("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/grenedalf-freq/new-genome-final-vcf/seasonal-glm-output-time-factor.txt", "w")
fileout2 <- file("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/grenedalf-freq/new-genome-final-vcf/seasonal-glm-output-year-factor.txt", "w")
# Loop through each line in freq_matrix
for (i in 1:nrow(freq_matrix)) {
  line <- freq_matrix[i, ]

  # Fit the GLM model for the current line
  out <- summary(glm(line ~ popinfo_glm$time + popinfo_glm$Y, family = binomial, weights = dp[i, ]))

  #out <- summary(glm(freq_matrix[1, ] ~ popinfo_glm$time + popinfo_glm$Y, family = binomial, weights = dp[1, ]))
  # Store the summary
  out2 <- out$coefficient[2,c(1,3,4)]
  out3 <- out$coefficient[3,c(1,3,4)]

  # Concatenate the values into a single line
  output_line2 <- paste(out2,collapse = "\t")
  output_line3 <- paste(out3,collapse = "\t")
  # Write the results for the current line to the file as a single line
  writeLines(output_line2, con = fileout)
  writeLines(output_line3, con = fileout2)
}

#popinfo_glm$Y <- factor(popinfo_glm$Y, ordered = FALSE)

# Close the file
close(fileout)

#time as factor

par(mfrow = c(2, 2))

time <- read.table("seasonal-glm-output-time-factor.txt", header = FALSE)

n <- length(time$V3)
expected_p_values <- seq(0, 1, length.out = n)
# Create a QQ plot
qqplot(-log10(expected_p_values), -log10(time$V3), xlab = "Expected -log10(p)", ylab = "Observed -log10(p)", main = "Time as unordered factor")
# Add a reference line for a perfect uniform distribution
abline(0, 1, col = "red")
#year <- read.table("p-values-year.txt", header = FALSE)
```

the model we ran assumes time as an unordered factor. So each sample stands alone in the analysis. This is what we chose last time.
I will now analyse the qqplot, and then check how many significant SNPs we get out of this model.

```bash 
#file 1
#########
awk -F'\t' '$3 <= 0.05 { count++ } END {print count}' seasonal-glm-output-time-factor.txt 
#6451 (vs 1640 last time)

awk '{print $1"_"$2}' chrom-poly-position-seasonal-data.txt > joint.txt 
awk '{print $3'} seasonal-glm-output-time.txt > p-values-time.txt 
paste chrom-poly-position-seasonal-data.txt joint.txt p-values-time.txt > full-table-ptime.txt

# add header
echo "chrom pos chrompos p" | cat - full-table-ptime.txt > temp && mv temp full-table-ptime.txt
echo "chrom  pos chrompos p" | cat - full-table-pyear.txt > temp && mv temp full-table-pyear.txt

awk -F'\t' '$4 <= 0.05 { count++ } END { print count }' full-table-ptime.txt
awk -F'\t' '$4 <= 0.05 { count++ } END { print count }' full-table-pyear.txt

```bash

#file2

awk '{print $3}' seasonal-glm-output-time-factor-unordered.txt > final-glm-time.txt
awk -F'\t' '$1 <= 0.05 { count++ } END { print count }' final-glm-time.txt
#1684
paste chrom-poly-position-seasonal-data.txt joint.txt final-glm-time.txt> full-table-final-glm-time.txt
# add header
echo "chrohem pos chrompos p" | cat - full-table-final-glm-time.txt > temp && mv temp full-table-final-glm-time.txt

awk -F'\t' '$4 <= 0.05 { print }' full-table-ptime.txt > significant-snps-time.txt

awk 'NR>1' cov_new_vcf.txt > cov_no_head.txt
paste chrom-poly-position-seasonal-data.txt cov_no_head.txt> cov_with_chrom.txt

awk 'NR==FNR{a[$1$2]=$0; next} ($1$2 in a){print a[$1$2], $0}' cov_with_chrom.txt significant-snps-time.txt> cov-sig-snps.txt
```

okay, now we will run baypass 04.08.2025

First I need to convert the files into a baypass object
subset populations of interest
```bash
srun --pty --x11 --nodes=1 --cpus-per-task=4 --mem=50000 --time=07:00:00 /bin/bash
module load bcftools/1.10.2

#subset the vcf for the pops of interest
bcftools view -s EA_2009_T2,EA_2009_T4,EA_2011_T1,EA_2011_T2,EA_2015_T1,EA_2015_T4,EA_2022_T1,EA_2022_T4 /gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/grenedalf-freq/new-genome-final-vcf/final-sorted.vcf.gz -o ~/work/seasonal_adaptation/analysis/new-genome-baypass/new-subset.vcf.gz
```
create genobypass data
```R

library(poolfstat)
dat <- vcf2pooldata(vcf.file="~/work/seasonal_adaptation/analysis/new-genome-baypass/new-subset.vcf.gz",poolsizes=rep(100,8))

pooldata2genobaypass(dat,writing.dir="~/work/seasonal_adaptation/analysis/new-genome-baypass/")
```
now I can run baypass
First, we need to run the core model to generate the omega
 -d0yij should be 1/5 of the smallest pool size
 https://esnielsen.github.io/post/pool-seq-analyses-poolfstat-baypass/

```bash 

#!/bin/bash
#SBATCH -D /gxfs_home/geomar/smomw573/work/seasonal_adaptation/log_files/
#SBATCH --mail-type=END
#SBATCH --mail-user=jnascimento@geomar.de
#SBATCH --partition=base
#SBATCH --nodes=5
#SBATCH --tasks-per-node=2
#SBATCH --cpus-per-task=4
#SBATCH --mem=50G
#SBATCH --time=48:00:00
#SBATCH --job-name=bayenv-stand
#SBATCH --output=baypass-core-model.out
#SBATCH --error=baypass-core-model.err


inputdir=~/work/seasonal_adaptation/analysis/variants
inputdir2=~/work/seasonal_adaptation/analysis/new-genome-baypass
outputdir=~/work/seasonal_adaptation/analysis/new-genome-baypass/core-model
baypassdir=/gxfs_home/geomar/smomw573/software/baypass_public-master/sources

$baypassdir/g_baypass \
-gfile $inputdir2/genobaypass \
-poolsizefile $inputdir/haploid-size-baypass-txt \
-outprefix $outputdir/stand-cov-model -nthreads 10 -seed 3 \
-d0yij 20 \
-npilot 25

```
05.08.2025
Once the core model is done, I will run both, the standard covariate model and the auxiliare covariate model to look at selection.
Both models give intereting outputs, while the first relates more to the strenght and direction of selection in realtion to the covariate (season) How much the allele frequency changes from spring to winter (strength)
and the auxiliate covariate is used to decide if the association is statistically strong enough to be confident.
```bash
#!/bin/bash
#SBATCH -D /gxfs_home/geomar/smomw573/work/seasonal_adaptation/log_files/
#SBATCH --mail-type=END
#SBATCH --mail-user=jnascimento@geomar.de
#SBATCH --partition=base
#SBATCH --nodes=5
#SBATCH --tasks-per-node=2
#SBATCH --cpus-per-task=4
#SBATCH --mem=50G
#SBATCH --time=22:00:00
#SBATCH --array=1-500
#SBATCH --job-name=bpsim-c2
#SBATCH --output=baysim-c2-%a.out
#SBATCH --error=baysim-c2-%a.err


inputdir=~/work/seasonal_adaptation/analysis/variants
omega=~/work/seasonal_adaptation/analysis/baypass
outputdir=~/work/seasonal_adaptation/analysis/new-genome-baypass/c2-model/
baypassdir=~/software/baypass_public-master/sources

~/work/seasonal_adaptation/analysis/baypass/core-model_mat_omega.out

$baypassdir/g_baypass \
-gfile $inputdir/genobaypass \
-efile $inputdir/cov-baypass.txt \
-omegafile $outputdir/core-model_mat_omega.out \
-poolsizefile $inputdir/haploid-size-baypass-txt \
-covmcmc \
-d0yij 20 \
-outprefix $outputdir/stand-cov-model -nthreads 10 -seed 3


# Auxiliary stats model
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


inputdir=~/work/seasonal_adaptation/analysis/variants/genobaypass
outputdir=/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/baypass
baypassdir=/gxfs_home/geomar/smomw573/software/baypass_public-master/sources


$baypassdir/g_baypass \
-gfile $inputdir/genobaypass \
-efile $inputdir/cov-baypass.txt.constrast \
-omegafile $outputdir/core-model_mat_omega.out \
-poolsizefile $inputdir/haploid-size-baypass-txt \
-d0yij 20 \
-outprefix $outputdir/aux-model -nthreads 10 -auxmodel
```
but for now I can already prepare the simulations so I do not waste time.
First I will sum the MAF across all pops
```R
freq <- read.csv("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/grenedalf-freq/new-genome-final-vcf/frequency_maf.csv", header=T)
freq <- read.csv("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/grenedalf-freq/new-genome-final-vcf/header_freq.csv", header=T)

freq$MAF <- ifelse(freq$TOTAL.FREQ > 0.5, 1-freq$TOTAL.FREQ, freq$TOTAL.FREQ)
freq$marker <- 1:nrow(freq)

write.table(freq, file = "maf-new-vcf-corrected.txt", sep="\t", quote=FALSE, row.names=FALSE)
hist(freq$TOTAL.FREQ,breaks=50)
hist(freq$MAF,breaks=50)
#0
sum(freq$MAF < 0.01)
#0
min(freq$MAF)
#0.01
max(freq$MAF)
#0.5
max(freq$TOTAL.FREQ)
#0.99
sum(freq$MAF < 0.5)
#1874822
sum(freq$MAF < 0.05)
#1214086
```
pull counts REF ALT alleles from vcf
pre-sync-file.R
final-sorted.vcf.gz

```bash

cut -d, -f1-4,5,6,9,10,11,12,13,14,15,16,19,20,21,22,27,28 frequency_counts.csv | tr ',' '\t' > selected_frequency_counts.tsv

srun --pty --x11 --nodes=1 --cpus-per-task=4 --mem=20000 --time=07:00:00 /bin/bash
```
```R
#fix columns for ocunts
# Read the tab-delimited file
data <- read.table("selected_frequency_counts.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Columns to keep as is
base_cols <- c("CHROM", "POS", "REF", "ALT")

# List of timepoints you want to keep (change/add as needed)
timepoints <- c("2009_T2", "2009_T4", "2011_T1", "2011_T2", "2015_T1", "2015_T4", "2022_T1", "2022_T4")

# Initialize a new dataframe with the base columns
new_data <- data[, base_cols]

# For each timepoint, combine REF_CNT and ALT_CNT into one column with "REF,ALT" format
for (tp in timepoints) {
  ref_col <- paste0("EA_", tp, ".REF_CNT")
  alt_col <- paste0("EA_", tp, ".ALT_CNT")
  
  # Make sure columns exist in your data
  if (all(c(ref_col, alt_col) %in% colnames(data))) {
    new_data[[tp]] <- paste0(data[[ref_col]], ",", data[[alt_col]])
  } else {
    warning(paste("Columns missing for", tp))
  }
}

# Write the new table to a tab-delimited file without row names or quotes
write.table(new_data, "allele-count-final-new-genome.txt", sep = "\t", quote = FALSE, row.names = FALSE)

```

# to_sync.py

```python
import argparse
from argparse import RawTextHelpFormatter

parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter,description="""description""")
parser.add_argument('--input', required=True, dest='Input', type=str, help="Input file in chrom/pos/ref/alt/allele_counts-comma separated in single columns for each pop format")
parser.add_argument('--output', required=True, dest='Output', type=str, help="Output file in sync format.")

args = parser.parse_args()


InputFile=open(args.Input,"r")
OutputFile=open(args.Output,"w")

base_dict = {'A':0, 'T':1, 'C':2, 'G':3}

###
# define conversion function

def convert_to_sync(ref_allele, var_allele, varscan_input):
    ct_ref=varscan_input.split(",")[0]
    ct_var=varscan_input.split(",")[1]
    sync_line = "0:0:0:0:0:0".split(":")
    sync_line[base_dict.get(ref)]=ct_ref
    sync_line[base_dict.get(var)]=ct_var
    return ":".join(sync_line)


firstline = True

for line in InputFile:

    if firstline:    #skip first line
        firstline = False
        continue

    line=line.split('\n')[0] #drop line break
    cols=line.split('\t')

    ref=cols[2]
    var=cols[3]

    line_out = []

    for pop_input in cols[4:]:
        if len(line_out) == 0:
            line_out = convert_to_sync(ref, var, pop_input)

        else:
            line_out = line_out + "\t" + convert_to_sync(ref, var, pop_input)

    line_out = "\t".join(cols[0:3]) + "\t" + line_out

    OutputFile.write(line_out+'\n')

```

module load python/3.11.5

python3 to.py --input allele-count-final-new-genome.txt --output output_sync.txt


now I can urn the simulations in poolseq

#estimated number of generations in each year
#assuming generational turnover of 2-weeks
#2009: 2 generations
#2011: 1 generations
#2015: 4 generations
#2022: 4 generations


```R
library(poolSeq)
#     allele frequency trajectories. Columns correspond to the time
#     points defined in ‘t’ and rows to individual replicates. If only
#     one replicate should be simulated a numeric vector is returned
#     instead of a matrix.#
#
#     Allele frequencies in the output are guaranteed to be ordered
#     increasing by the number of generations.
#
#Author(s):
#
#     Thomas Taus
#
#Examples:
#
#     # simulate allele frequency trajectories of individual loci using the Wright-Fisher model (diploid individuals)
#     Ne <- 200
#     gen <- seq(0, 100, by=10)
#     alleleFreqs <- wf.traj(p0=rep(0.5, times=500), Ne=Ne, t=gen)
#     # look at a subset of the generated data
#     head(alleleFreqs)
#
#     # plot allele frequency trajectories
#     plot(1, type="n", xlim=c(0, max(gen)), ylim=c(0, 1), main="Neutral Genetic Drift", xlab="Generation", ylab="Allele frequency (%)")
#     for(r in 1:nrow(alleleFreqs)) {
#       lines(gen, alleleFreqs[r,])
#     }
#
#     # simulate allele frequency trajectories including selection
#     alleleFreqs <- wf.traj(p0=rep(0.05, times=500), Ne=Ne, t=gen, s=0.1, h=0.5)
#
#     # plot results
#     plot(1, type="n", xlim=c(0, max(gen)), ylim=c(0, 1), main="Positive Selection", xlab="Generation", ylab="Allele frequency (%)")
#     for(r in 1:nrow(alleleFreqs)) {
#       lines(gen, alleleFreqs[r,])
#     }
#
#     # add the trajectory under selection for a population of infinite size (no random genetic drift)
#     lines(gen, wf.traj(p0=0.05, Ne=NA, t=gen, s=0.1, h=0.5), col="red")

mySync <- read.sync(file="/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/grenedalf-freq/new-genome-final-vcf/output_sync.txt", polarization="reference", gen=c(1,3,1,2,1,4,1,4),repl=c(1,1,2,2,3,3,4,4))

listSync <- list(mySync)
names(listSync) <- c("mySync")

## estimate effective population size:
#reference with explanations: https://pmc.ncbi.nlm.nih.gov/articles/PMC5068858/

af <- as.data.frame(mySync@alleles)
estNe_2009 <- estimateNe(p0=af[,7], pt=af[,9], cov0=af[,8], covt=af[,10], ploidy=2, t=2, poolSize=c(200, 200), method=c("JR.planII", "W.planII", "P.planII", "P.alt.1step.planII"))

estNe_2011 <- estimateNe(p0=af[,11], pt=af[,13], cov0=af[,12], covt=af[,14],ploidy=2, t=1, poolSize=c(200, 200), method=c("JR.planII", "W.planII", "P.planII", "P.alt.1step.planII"))
estNe_2015 <- estimateNe(p0=af[,15], pt=af[,17], cov0=af[,16], covt=af[,18],ploidy=2, t=4, poolSize=c(200, 200), method=c("JR.planII", "W.planII", "P.planII", "P.alt.1step.planII"))

estNe_2022 <- estimateNe(p0=af[,19], pt=af[,21], cov0=af[,20], covt=af[,22], ploidy=2,t=3, poolSize=c(200, 200), method=c("JR.planII", "W.planII", "P.planII", "P.alt.1step.planII"))

estNe_2009
        Nw.planII          Njr.planII Np.alt.1step.planII           Np.planII
        -321.62409        -14715.99298         22401.71524           -98.97588
estNe_2011
        Nw.planII          Njr.planII Np.alt.1step.planII           Np.planII
        -313.04428           332.73255           299.16715           -27.01587
estNe_2015
        Nw.planII          Njr.planII Np.alt.1step.planII           Np.planII
        -1588.33694        150646.28127         15597.16994           -99.96972
estNe_2022
        Nw.planII          Njr.planII Np.alt.1step.planII           Np.planII
        -2225.28315          1096.31407           963.83236           -80.52603


#a bit different than last time... but alt 1 step plan II seems the best? check with Reid

t1Sync <- read.sync(file="/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/grenedalf-freq/new-genome-final-vcf/sync2009.txt", polarization="reference", 
  gen=c(1,3),repl=c(1,2)) #import only the start/end points of interest

t2Sync <- read.sync(file="/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/grenedalf-freq/new-genome-final-vcf/sync2011.txt", 
  gen=c(1,2),repl=c(1,2)) #import only the start/end points of interest

t3Sync <- read.sync(file="/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/grenedalf-freq/new-genome-final-vcf/sync2015.txt", 
  gen=c(1,4),repl=c(1,2)) #import only the start/end points of interest

t4Sync <- read.sync(file="/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/grenedalf-freq/new-genome-final-vcf/sync2022.txt", 
  gen=c(1,4),repl=c(1,2)) #import only the start/end points of interest

```
27.08.2025

So we decided we will use SLiM to calculate drift (or neutral selection) on my populations.
The first step is to run Kwi's code to msprime, so we get the input parameters for SLiM.
For this, I will estimate pi (genetic diversity) of my populations, and with this I can use msprime to find out the representative Ne.

# I will come back to this header once I have read more about msprime and slim and correct the text

that is the script 10.estimate-genome-diversity-final.sh

```
#!/bin/bash
#SBATCH -D /gxfs_home/geomar/smomw573/work/seasonal_adaptation/log_files/
#SBATCH --mail-type=END
#SBATCH --mail-user=jnascimento@geomar.de
#SBATCH --partition=base
#SBATCH --nodes=1
#SBATCH --tasks-per-node=2
#SBATCH --cpus-per-task=1
#SBATCH --mem=50G
#SBATCH --time=3:00:00
#SBATCH --job-name=bayenv-stand
#SBATCH --output=genome-diversity.out
#SBATCH --error=genome-diversity.err
   
module load bzip2/1.0.8 

GRENEDALF=/gxfs_home/geomar/smomw573/software/grenedalf/bin/grenedalf
VCF=/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/grenedalf-freq/new-genome-final-vcf/final-sorted.vcf
MPILEUP=/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/variants/snps-exc20072011t32009t12015t3.mpileup
GENOME=$WORK/smomw573/seasonal_adaptation/new_genome/Purple_Eaffinis_genome_assembly.fa 
OUTPUT=/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/genome-diversity

$GRENEDALF diversity --pileup-path $MPILEUP --filter-sample-min-count 1 --reference-genome-fasta-file $GENOME --allow-file-overwriting --window-type genome --pool-sizes 100 --measure theta-pi --file-prefix whole-genome- --out-dir $OUTPUT

$GRENEDALF diversity --pileup-path $MPILEUP --filter-sample-min-count 1 --reference-genome-fasta-file $GENOME --allow-file-overwriting --window-type sliding  --window-sliding-width 10000 --pool-sizes 100 --measure theta-pi --file-prefix diversity-10kwindow- --out-dir $OUTPUT

$GRENEDALF diversity --pileup-path $MPILEUP --filter-sample-min-count 1 --reference-genome-fasta-file $GENOME --allow-file-overwriting --window-type sliding --window-sliding-width 10000 --window-sliding-stride 10000 --pool-sizes 100 --measure theta-pi --file-prefix diversity-10kwindow- 

$GRENEDALF diversity --pileup-path $MPILEUP --reference-genome-fasta-file $GENOME --allow-file-overwriting --window-type sliding --window-sliding-width 10000 --window-sliding-stride 10000 --pool-sizes 100 --measure theta-pi --filter-sample-min-count 1 --pileup-min-base-qual 15 --filter-sample-min-coverage 30 --filter-sample-max-coverage 205  --file-prefix diversity-10kwindow-filtered-min-count-1 --out-dir $OUTPUT

$GRENEDALF diversity --pileup-path $MPILEUP --reference-genome-fasta-file $GENOME --allow-file-overwriting --window-type sliding --window-sliding-width 10000 --window-sliding-stride 10000 --pool-sizes 100 --measure theta-pi --filter-sample-min-count 2 --pileup-min-base-qual 15 --filter-sample-min-coverage 30 --filter-sample-max-coverage 205  --file-prefix diversity-10kwindow-filtered-min-count-2 --out-dir $OUTPUT


```

16.09 
I spoke with Reid and Kwi, and we decided to analyse the data differently.
First I will estimate pi using vcftools from my vcf file with the actual data.

Then, I will run msprime to produce a new vcf from an original population with no effects of selection. I am doing this so I can estimate the correct effective population size, that I cannot estimate from my data as something is giving me negative values

script 9.run-Burnin-msprime-trial15.sh

In addition, I will estimate pi using msprime 
this will be done using the run-burnin-msprime.sh and .py scripts


```bash
msprime_burnin_n3000_44006078_1.vcf
vcftools --vcf msprime_burnin_n2000_44006078_1.vcf --window-pi 10000 --out window-pi-10kb-msprime_burnin_n2000_44006078_1
vcftools --vcf msprime_burnin_n3000_44006078_1.vcf --window-pi 10000 --out window-pi-10kb-msprime_burnin_n3000_44006078_1
vcftools --vcf final-sorted.vcf --window-pi 10000 --out window-pi-10kb-final-vcf.windowed.pi

#but first I want to compare the diversity across my samples, just to check that everything makes sense, for this I need to separate my files

vcftools --vcf final-sorted.vcf --indv EA_2009_T2 --recode --out EA_2009_T2.vcf
vcftools --vcf final-sorted.vcf --indv EA_2009_T4 --recode --out EA_2009_T4.vcf

vcftools --vcf final-sorted.vcf --indv EA_2011_T1 --recode --out EA_2011_T1.vcf
vcftools --vcf final-sorted.vcf --indv EA_2011_T2 --recode --out EA_2011_T2.vcf

vcftools --vcf final-sorted.vcf --indv EA_2015_T1 --recode --out EA_2015_T1.vcf
vcftools --vcf final-sorted.vcf --indv EA_2015_T4 --recode --out EA_2015_T4.vcf

vcftools --vcf final-sorted.vcf --indv EA_2022_T1 --recode --out EA_2022_T1.vcf
vcftools --vcf final-sorted.vcf --indv EA_2022_T4 --recode --out EA_2022_T4.vcf


for vcf in EA_*.vcf.recode.vcf
do
    prefix=$(basename $vcf .vcf.recode.vcf)
    vcftools --vcf $vcf --window-pi 10000 --out window-pi-10kb-$prefix
done

```
```R
library(ggplot2)
library(dplyr)

# List all .sites.pi files
#files <- list.files(pattern = "*.pi$")

# Read all files and add a column for sample name
year2009.t1 <- read.table("window-pi-10kb-EA_2009_T2.windowed.pi", header = TRUE)
year2009.t2 <- read.table("window-pi-10kb-EA_2009_T4.windowed.pi", header = TRUE)

year2011.t1 <- read.table("window-pi-10kb-EA_2011_T1.windowed.pi", header = TRUE)
year2011.t2 <- read.table("window-pi-10kb-EA_2011_T2.windowed.pi", header = TRUE)

year2015.t1 <- read.table("window-pi-10kb-EA_2015_T1.windowed.pi", header = TRUE)
year2015.t2 <- read.table("window-pi-10kb-EA_2015_T4.windowed.pi", header = TRUE)

year2022.t1 <- read.table("window-pi-10kb-EA_2022_T1.windowed.pi", header = TRUE)
year2022.t2 <- read.table("window-pi-10kb-EA_2022_T4.windowed.pi", header = TRUE)

all <- read.table("window-pi-10kb-final-vcf.windowed.pi", header = TRUE)
msprime.burnin <- read.table("window-pi-10kb-msprime_burnin_n2000_44006078_1.windowed.pi", header = TRUE)
msprime.burnin2 <- read.table("window-pi-10kb-msprime_burnin_n3000_44006078_1.windowed.pi", header = TRUE)
msprime.burnin3 <- read.table("window-pi-10kb-msprime_burnin_n10000_44006078_1.windowed.pi", header = TRUE)
msprime.burnin4 <- read.table("window-pi-10kb-msprime_burnin_n100000_44006078_1.windowed.pi", header = TRUE)
msprime.burnin5 <- read.table("window-pi-10kb-msprime_burnin_n150000_44006078_1.windowed.pi", header = TRUE)
# Suppose your dfs are df1, df2, ..., df10
dfs <- list(year2009.t1, year2009.t2, year2011.t1, year2011.t2, year2015.t1, year2015.t2, year2022.t1, year2022.t2, all, msprime.burnin, msprime.burnin2, msprime.burnin3, msprime.burnin4, msprime.burnin5)

dfs <- list(msprime.burnin, msprime.burnin2, msprime.burnin3, msprime.burnin4, msprime.burnin5)

# Convert CHROM to character in all dfs
dfs <- lapply(dfs, function(x) { x$CHROM <- as.character(x$CHROM); x })

# Add sample name
dfs <- lapply(seq_along(dfs), function(i) { 
  dfs[[i]]$Sample <- paste0("Sample", i)
  dfs[[i]]
})

# Bind together
combined <- bind_rows(dfs)

# combined is your data frame with columns: CHROM, BIN_START, BIN_END, N_VARIANTS, PI, Sample

p <- ggplot(combined, aes(x = Sample, y = PI)) +
  geom_boxplot(fill = "skyblue", outlier.color = "red") +
  theme_bw() +
  labs(title = "Nucleotide Diversity per Sample",
       x = "Sample",
       y = "Nucleotide Diversity (π)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("window-pi-10kb-msprime_burnin_n2000_44006078_1.png", plot = p, width = 8, height = 6, dpi = 300)


#plot diversity from mpileup file using grenedalf whole genome diversity

div1 <- read.csv("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/genome-diversity/diversity-10kwindow-filtered-min-count-1diversity.csv", header = TRUE)

div <- read.csv("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/genome-diversity/diversity-10kwindow-filtered-min-count-2diversity.csv", header = TRUE)

# Absolute theta_pi
abs_cols <- grep("theta_pi_abs", colnames(div1))
theta_abs_mean <- colMeans(div1[, abs_cols], na.rm = TRUE)

# Relative theta_pi
rel_cols <- grep("theta_pi_rel", colnames(div))
theta_rel_mean <- colMeans(div1[, rel_cols], na.rm = TRUE)

theta_abs_mean
theta_rel_mean

#theta_rel_mean for count-1
 snps.exc20072011t32009t12015t3.1.theta_pi_abs
                                      48.16004
 snps.exc20072011t32009t12015t3.2.theta_pi_abs
                                      81.61669
 snps.exc20072011t32009t12015t3.3.theta_pi_abs
                                     113.28190
 snps.exc20072011t32009t12015t3.4.theta_pi_abs
                                      99.85945
 snps.exc20072011t32009t12015t3.5.theta_pi_abs
                                     113.67790
 snps.exc20072011t32009t12015t3.6.theta_pi_abs
                                     117.61982
 snps.exc20072011t32009t12015t3.7.theta_pi_abs
                                     105.47362
 snps.exc20072011t32009t12015t3.8.theta_pi_abs
                                     120.36042
 snps.exc20072011t32009t12015t3.9.theta_pi_abs
                                     117.06856
snps.exc20072011t32009t12015t3.10.theta_pi_abs
                                     112.68734
snps.exc20072011t32009t12015t3.11.theta_pi_abs
                                     118.50635
snps.exc20072011t32009t12015t3.12.theta_pi_abs
                                     119.97110
 snps.exc20072011t32009t12015t3.1.theta_pi_rel
                                    0.02176171
 snps.exc20072011t32009t12015t3.2.theta_pi_rel
                                    0.02211102
 snps.exc20072011t32009t12015t3.3.theta_pi_rel
                                    0.02268525
 snps.exc20072011t32009t12015t3.4.theta_pi_rel
                                    0.02242484
 snps.exc20072011t32009t12015t3.5.theta_pi_rel
                                    0.02257079
 snps.exc20072011t32009t12015t3.6.theta_pi_rel
                                    0.02254696
 snps.exc20072011t32009t12015t3.7.theta_pi_rel
                                    0.02243293
 snps.exc20072011t32009t12015t3.8.theta_pi_rel
                                    0.02265077
 snps.exc20072011t32009t12015t3.9.theta_pi_rel
                                    0.02269361
snps.exc20072011t32009t12015t3.10.theta_pi_rel
                                    0.02258210
snps.exc20072011t32009t12015t3.11.theta_pi_rel
                                    0.02265767
snps.exc20072011t32009t12015t3.12.theta_pi_rel
                                    0.02268685

#the new pi with filtering the mpileup gives a pi of 0.18
 snps.exc20072011t32009t12015t3.1.theta_pi_rel
                                    0.01874406
 snps.exc20072011t32009t12015t3.2.theta_pi_rel
                                    0.01929013
 snps.exc20072011t32009t12015t3.3.theta_pi_rel
                                    0.01998514
 snps.exc20072011t32009t12015t3.4.theta_pi_rel
                                    0.01966919
 snps.exc20072011t32009t12015t3.5.theta_pi_rel
                                    0.01997501
 snps.exc20072011t32009t12015t3.6.theta_pi_rel
                                    0.01993446
 snps.exc20072011t32009t12015t3.7.theta_pi_rel
                                    0.01964826
 snps.exc20072011t32009t12015t3.8.theta_pi_rel
                                    0.01996247
 snps.exc20072011t32009t12015t3.9.theta_pi_rel
                                    0.02005321
snps.exc20072011t32009t12015t3.10.theta_pi_rel
                                    0.02009395
snps.exc20072011t32009t12015t3.11.theta_pi_rel
                                    0.02020233
snps.exc20072011t32009t12015t3.12.theta_pi_rel
                                    0.02019708


#plot diversity
values1 <- as.numeric(theta_rel_mean)
labels <- names(theta_rel_mean)

# Make a barplot
boxplot(values, names.arg = labels,
        las = 2,            # Rotate x-axis labels
        col = "skyblue",
        main = "Theta PI Rel Values",
        ylab = "Theta PI Rel",
        cex.names = 0.6) 

#mean PI for all msprime runs


#this is for all the samples in the slim-simulations/pi-diversity folder
# List all .pi.windowed.pi files
files <- list.files(pattern = ".*pi\\.windowed\\.pi$")

# Compute mean PI for each file
mean_pi <- numeric(length(files))
for(i in seq_along(files)) {
  data <- read.table(files[i], header = TRUE)
  mean_pi[i] <- mean(data$PI, na.rm = TRUE)
}
sample_names <- gsub("\\.pi\\.windowed\\.pi$", "", files)

plot(mean_pi,
     pch = 16,
     col = "blue",
     xaxt = "n",
     xlab = "",
     ylab = "Mean Theta Pi",
     main = "Mean Theta Pi for Each Sample",
     ylim = c(0, max(mean_pi) * 1.1))

# Add tick marks only
axis(1, at = 1:length(mean_pi), labels = FALSE)

# Add 45-degree angled labels
text(x = 1:length(mean_pi), y = par("usr")[3] - 0.001, 
     labels = sample_names, srt = 45, adj = 1, xpd = TRUE)


  ```
> mean_pi
  [1] 4.110437e-05 3.985443e-05 4.125747e-05 4.102248e-05 4.130827e-05
  [6] 4.136886e-05 4.008537e-05 4.234909e-05 3.993301e-05 4.003110e-05
 [11] 4.019036e-05 4.052209e-05 4.089868e-05 3.971625e-05 4.176450e-05
 [16] 4.050267e-05 4.037276e-05 4.116317e-05 4.038356e-05 4.053535e-05
 [21] 3.020412e-03 3.036456e-03 3.024451e-03 3.019392e-03 3.028157e-03
 [26] 3.009479e-03 3.016265e-03 3.025009e-03 3.031782e-03 3.026814e-03
 [31] 3.037835e-03 3.027315e-03 3.002487e-03 3.040671e-03 3.025578e-03
 [36] 3.029020e-03 3.032072e-03 3.024904e-03 3.042467e-03 3.008353e-03
 [41] 3.005264e-03 3.017522e-03 3.017761e-03 3.004536e-03 3.013166e-03
 [46] 2.989575e-03 3.004096e-03 3.000563e-03 3.006003e-03 3.006262e-03
 [51] 3.005566e-03 3.007965e-03 2.995877e-03 3.003328e-03 2.989904e-03
 [56] 3.004867e-03 3.020493e-03 3.003779e-03 3.008372e-03 2.990597e-03
 [61] 1.553114e-05 1.603519e-05 1.607508e-05 1.579066e-05 1.609051e-05
 [66] 1.657744e-05 1.590877e-05 1.568287e-05 1.545750e-05 1.629534e-05
 [71] 1.646302e-05 1.617513e-05 1.618342e-05 1.610755e-05 1.586928e-05
 [76] 1.566955e-05 1.583037e-05 1.602291e-05 1.560014e-05 1.591275e-05
 [81] 7.983775e-05 8.036061e-05 7.894368e-05 8.021408e-05 7.994823e-05
 [86] 8.036035e-05 8.056505e-05 8.098807e-05 7.994363e-05 8.139769e-05
 [91] 8.144618e-05 8.037098e-05 7.976056e-05 8.097686e-05 7.960085e-05
 [96] 7.980144e-05 7.995177e-05 8.009438e-05 7.923251e-05 7.915375e-05
[101] 9.521931e-03 9.713270e-03 9.545206e-03 9.689136e-03 9.793207e-03
[106] 9.684635e-03 9.614764e-03 9.668507e-03 9.583753e-03 9.739784e-03
[111] 9.711075e-03 9.569603e-03 9.672729e-03 9.549037e-03 9.721176e-03
[116] 9.760158e-03 9.667779e-03 9.450604e-03 9.676393e-03 9.455001e-03
[121] 9.676752e-03 9.589606e-03 9.653170e-03
> sample_names
  [1] "msprime_burnin_n10000_44006078_1_pi.windowed.pi"
  [2] "msprime_burnin_n10000_44006078_10_pi.windowed.pi"
  [3] "msprime_burnin_n10000_44006078_11_pi.windowed.pi"
  [4] "msprime_burnin_n10000_44006078_12_pi.windowed.pi"
  [5] "msprime_burnin_n10000_44006078_13_pi.windowed.pi"
  [6] "msprime_burnin_n10000_44006078_14_pi.windowed.pi"
  [7] "msprime_burnin_n10000_44006078_15_pi.windowed.pi"
  [8] "msprime_burnin_n10000_44006078_16_pi.windowed.pi"
  [9] "msprime_burnin_n10000_44006078_17_pi.windowed.pi"
 [10] "msprime_burnin_n10000_44006078_18_pi.windowed.pi"
 [11] "msprime_burnin_n10000_44006078_19_pi.windowed.pi"
 [12] "msprime_burnin_n10000_44006078_2_pi.windowed.pi"
 [13] "msprime_burnin_n10000_44006078_20_pi.windowed.pi"
 [14] "msprime_burnin_n10000_44006078_3_pi.windowed.pi"
 [15] "msprime_burnin_n10000_44006078_4_pi.windowed.pi"
 [16] "msprime_burnin_n10000_44006078_5_pi.windowed.pi"
 [17] "msprime_burnin_n10000_44006078_6_pi.windowed.pi"
 [18] "msprime_burnin_n10000_44006078_7_pi.windowed.pi"
 [19] "msprime_burnin_n10000_44006078_8_pi.windowed.pi"
 [20] "msprime_burnin_n10000_44006078_9_pi.windowed.pi"
 [21] "msprime_burnin_n100000_44006078_1_pi.windowed.pi"
 [22] "msprime_burnin_n100000_44006078_10_pi.windowed.pi"
 [23] "msprime_burnin_n100000_44006078_11_pi.windowed.pi"
 [24] "msprime_burnin_n100000_44006078_12_pi.windowed.pi"
 [25] "msprime_burnin_n100000_44006078_13_pi.windowed.pi"
 [26] "msprime_burnin_n100000_44006078_14_pi.windowed.pi"
 [27] "msprime_burnin_n100000_44006078_15_pi.windowed.pi"
 [28] "msprime_burnin_n100000_44006078_16_pi.windowed.pi"
 [29] "msprime_burnin_n100000_44006078_17_pi.windowed.pi"
 [30] "msprime_burnin_n100000_44006078_18_pi.windowed.pi"
 [31] "msprime_burnin_n100000_44006078_19_pi.windowed.pi"
 [32] "msprime_burnin_n100000_44006078_2_pi.windowed.pi"
 [33] "msprime_burnin_n100000_44006078_20_pi.windowed.pi"
 [34] "msprime_burnin_n100000_44006078_3_pi.windowed.pi"
 [35] "msprime_burnin_n100000_44006078_4_pi.windowed.pi"
 [36] "msprime_burnin_n100000_44006078_5_pi.windowed.pi"
 [37] "msprime_burnin_n100000_44006078_6_pi.windowed.pi"
 [38] "msprime_burnin_n100000_44006078_7_pi.windowed.pi"
 [39] "msprime_burnin_n100000_44006078_8_pi.windowed.pi"
 [40] "msprime_burnin_n100000_44006078_9_pi.windowed.pi"
 [41] "msprime_burnin_n120000_44006078_1_pi.windowed.pi"
 [42] "msprime_burnin_n120000_44006078_10_pi.windowed.pi"
 [43] "msprime_burnin_n120000_44006078_11_pi.windowed.pi"
 [44] "msprime_burnin_n120000_44006078_12_pi.windowed.pi"
 [45] "msprime_burnin_n120000_44006078_13_pi.windowed.pi"
 [46] "msprime_burnin_n120000_44006078_14_pi.windowed.pi"
 [47] "msprime_burnin_n120000_44006078_15_pi.windowed.pi"
 [48] "msprime_burnin_n120000_44006078_16_pi.windowed.pi"
 [49] "msprime_burnin_n120000_44006078_17_pi.windowed.pi"
 [50] "msprime_burnin_n120000_44006078_18_pi.windowed.pi"
 [51] "msprime_burnin_n120000_44006078_19_pi.windowed.pi"
 [52] "msprime_burnin_n120000_44006078_2_pi.windowed.pi"
 [53] "msprime_burnin_n120000_44006078_20_pi.windowed.pi"
 [54] "msprime_burnin_n120000_44006078_3_pi.windowed.pi"
 [55] "msprime_burnin_n120000_44006078_4_pi.windowed.pi"
 [56] "msprime_burnin_n120000_44006078_5_pi.windowed.pi"
 [57] "msprime_burnin_n120000_44006078_6_pi.windowed.pi"
 [58] "msprime_burnin_n120000_44006078_7_pi.windowed.pi"
 [59] "msprime_burnin_n120000_44006078_8_pi.windowed.pi"
 [60] "msprime_burnin_n120000_44006078_9_pi.windowed.pi"
 [61] "msprime_burnin_n2000_44006078_1_pi.windowed.pi"
 [62] "msprime_burnin_n2000_44006078_10_pi.windowed.pi"
 [63] "msprime_burnin_n2000_44006078_11_pi.windowed.pi"
 [64] "msprime_burnin_n2000_44006078_12_pi.windowed.pi"
 [65] "msprime_burnin_n2000_44006078_13_pi.windowed.pi"
 [66] "msprime_burnin_n2000_44006078_14_pi.windowed.pi"
 [67] "msprime_burnin_n2000_44006078_15_pi.windowed.pi"
 [68] "msprime_burnin_n2000_44006078_16_pi.windowed.pi"
 [69] "msprime_burnin_n2000_44006078_17_pi.windowed.pi"
 [70] "msprime_burnin_n2000_44006078_18_pi.windowed.pi"
 [71] "msprime_burnin_n2000_44006078_19_pi.windowed.pi"
 [72] "msprime_burnin_n2000_44006078_2_pi.windowed.pi"
 [73] "msprime_burnin_n2000_44006078_20_pi.windowed.pi"
 [74] "msprime_burnin_n2000_44006078_3_pi.windowed.pi"
 [75] "msprime_burnin_n2000_44006078_4_pi.windowed.pi"
 [76] "msprime_burnin_n2000_44006078_5_pi.windowed.pi"
 [77] "msprime_burnin_n2000_44006078_6_pi.windowed.pi"
 [78] "msprime_burnin_n2000_44006078_7_pi.windowed.pi"
 [79] "msprime_burnin_n2000_44006078_8_pi.windowed.pi"
 [80] "msprime_burnin_n2000_44006078_9_pi.windowed.pi"
 [81] "msprime_burnin_n20000_44006078_1_pi.windowed.pi"
 [82] "msprime_burnin_n20000_44006078_10_pi.windowed.pi"
 [83] "msprime_burnin_n20000_44006078_11_pi.windowed.pi"
 [84] "msprime_burnin_n20000_44006078_12_pi.windowed.pi"
 [85] "msprime_burnin_n20000_44006078_13_pi.windowed.pi"
 [86] "msprime_burnin_n20000_44006078_14_pi.windowed.pi"
 [87] "msprime_burnin_n20000_44006078_15_pi.windowed.pi"
 [88] "msprime_burnin_n20000_44006078_16_pi.windowed.pi"
 [89] "msprime_burnin_n20000_44006078_17_pi.windowed.pi"
 [90] "msprime_burnin_n20000_44006078_18_pi.windowed.pi"
 [91] "msprime_burnin_n20000_44006078_19_pi.windowed.pi"
 [92] "msprime_burnin_n20000_44006078_2_pi.windowed.pi"
 [93] "msprime_burnin_n20000_44006078_20_pi.windowed.pi"
 [94] "msprime_burnin_n20000_44006078_3_pi.windowed.pi"
 [95] "msprime_burnin_n20000_44006078_4_pi.windowed.pi"
 [96] "msprime_burnin_n20000_44006078_5_pi.windowed.pi"
 [97] "msprime_burnin_n20000_44006078_6_pi.windowed.pi"
 [98] "msprime_burnin_n20000_44006078_7_pi.windowed.pi"
 [99] "msprime_burnin_n20000_44006078_8_pi.windowed.pi"
[100] "msprime_burnin_n20000_44006078_9_pi.windowed.pi"
[101] "msprime_burnin_n200000_44006078_10_pi.windowed.pi"
[102] "msprime_burnin_n200000_44006078_10_run2_pi.windowed.pi"
[103] "msprime_burnin_n200000_44006078_11_pi.windowed.pi"
[104] "msprime_burnin_n200000_44006078_11_run2_pi.windowed.pi"
[105] "msprime_burnin_n200000_44006078_12_pi.windowed.pi"
[106] "msprime_burnin_n200000_44006078_12_run2_pi.windowed.pi"
[107] "msprime_burnin_n200000_44006078_14_pi.windowed.pi"
[108] "msprime_burnin_n200000_44006078_14_run2_pi.windowed.pi"
[109] "msprime_burnin_n200000_44006078_15_pi.windowed.pi"
[110] "msprime_burnin_n200000_44006078_18_pi.windowed.pi"
[111] "msprime_burnin_n200000_44006078_18_run2_pi.windowed.pi"
[112] "msprime_burnin_n200000_44006078_19_pi.windowed.pi"
[113] "msprime_burnin_n200000_44006078_19_run2_pi.windowed.pi"
[114] "msprime_burnin_n200000_44006078_20_pi.windowed.pi"
[115] "msprime_burnin_n200000_44006078_20_run2_pi.windowed.pi"
[116] "msprime_burnin_n200000_44006078_5_pi.windowed.pi"
[117] "msprime_burnin_n200000_44006078_5_run2_pi.windowed.pi"
[118] "msprime_burnin_n200000_44006078_6_pi.windowed.pi"
[119] "msprime_burnin_n200000_44006078_6_run2_pi.windowed.pi"
[120] "msprime_burnin_n200000_44006078_7_pi.windowed.pi"
[121] "msprime_burnin_n200000_44006078_7_run2_pi.windowed.pi"
[122] "msprime_burnin_n200000_44006078_9_pi.windowed.pi"
[123] "msprime_burnin_n200000_44006078_9_run2_pi.windowed.pi"


Date 29.07.2025

Ok, summing up what I have done so far

I am playing around with msprime to understand how much can we lower the chrom lenght without messing up with the pi measurement. the average genome pi calculated for our mpileup files was of 0.018 to 0.02, whilst for the msprime burn runs we got an average pi of 0.010. whilst this is a bit lower, it is still within the same unit. I spoke with Kwi and she said that this should work fine, and she has done the same for her population.
I also tested chromlenghts of 3Mb, 5Mb and 44Mb (the average chrom length for our species), and it did not really make a difference. Need to chat with Reid why, but it is somehow because we are estimating diversity in certain regions and as it is a representation of the genome, kind of divided in smaller chuncks, it should not make a difference. 
I check and the diversity levels I get with Ne 200000 and chrom length of either 5Mb or 3Mb is the same. For this reason I decided to run the simulations with 3Mb chrom length to save memory. 
I used the mutation rates and recombinations rates from these two papers recombination https://www.nature.com/articles/s41467-022-31622-8 Highly structured populations of copepods at risk to deep (Arbizu 2024)
I calculated the average number of generations per year (4.4) and estimated as good as I could how many generations were in between the years. I ran slim with an output of 100 individuals, as this is what we sampled from the populations for poolseq.

Ok, I ran the simulations and it is super fast, so I decided to have a bit longer chrom length (5Mb). (details are in the script 11.slim-neutral.sh)

script 11.slim-neutral.sh and 01_neutral_run_slim_copes.sh

I will now create pooled vcf files from the ones generated from the simulations 

files from 20th of october in slim-simulations/simulations directory
```bash
in_dir="/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/slim-simulations-new-genome/simulations"
out_dir="${in_dir}/vcfs/compressed"


for f in "${in_dir}"/*.vcf; do
    fname=$(basename "$f")                # extract file name only (no path)
    out="${out_dir}/${fname}.gz"          # where to save compressed file
    echo "Compressing and indexing $fname ..."
    bgzip -c "$f" > "$out"
    tabix -p vcf "$out"
done

#into vcf_compressed directory 
```

Questions to Reid:

1- Is it correct to assume that the msprime coalescent population will be inputed in slim at 2009_t1, to see how this changed through time? YES. 
#need to produce a vcf with the 100 individuals without simulating. DONE
2- Do we assume no reproduction has occured? YES. so it is a late event
3- can he recommend a not so hard paper for the lecture?

what is evolutionary rescue, why do we care? 
bumpus
experimental evolution - lensky

20.10.2025
I moved all the files into /simulations/vcfs

Okay, now that I have finished the mpileup burnin runs, and have ran slim (script 11.slim-neutral.sh), I am ready to simulate everything with baypass. However, I first need to do a few adjustments.
Change everything to a baypass friendly format (.sync format) and then adjust the script with the actual baypass simulation, adjusting the coverage part for my samples and removing the random sampling from a population, as this has already been done by slim. I ran 100 simulations per time point in slim, with 100 individuals from the 200000 that were estimated to reach the same levels of diversity.

okay. I will search in my old scripts how I transformed the vcf file into a sync file.

First, I transformed the vcfs into pooled vcfs as in the code i pool-vcfs.py in the folder. I also had to fix the header of the vcf files with the script fix-vcfs.sh that is in each of the folders. There was a rpboelm with the AD field.
Then, with grenedalf, I transformed the vcfs in sync files


OK- some files had duplicate variants from slim. I removed them. Ask reid if thats okay?

#check the duplicates
bcftools query -f '%CHROM\t%POS\n' neutral_76_sampled_2022_late_sorted.vcf.gz | sort | uniq -d

```bash
#!/bin/bash

# Create a folder for compressed VCFs
mkdir -p ../vcfs_compressed

# Loop over all uncompressed VCFs
for VCF in *.vcf; do
    echo "Compressing $VCF ..."
    bgzip -c "$VCF" > "/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/slim-simulations-new-genome/simulations/vcfs/vcfs_compressed/${VCF}.gz"
done

for f in *.vcf.gz; do
    tabix -p vcf "$f"
done
#run pool-vcfs.py
#moved to pooled

#sort vcfs
#!/bin/bash
for f in /gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/slim-simulations-new-genome/simulations/vcfs/pooled/*.vcf; do
    echo "Compressing $f ..."
    bgzip -f "$f"
done
for f in /gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/slim-simulations-new-genome/simulations/vcfs/pooled/*.vcf.gz; do
    echo "Indexing $f ..."
    tabix -f -p vcf "$f"
done

in_dir="/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/slim-simulations-new-genome/simulations/vcfs/vcfs_compressed"
out_dir="/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/slim-simulations-new-genome/simulations/vcfs/sorted"

for VCF in "${in_dir}"/*.vcf.gz; do
    echo "Processing $VCF ..."
    BASE=$(basename "$VCF" .vcf.gz)
    OUT="${out_dir}/${BASE}_sorted.vcf.gz"

    # Sort and compress
    bcftools sort -Oz -o "$OUT" "$VCF"

    # Index the sorted file
    tabix -p vcf "$OUT"

    echo "✓ Done: $OUT"
done

# Loop over all sorted VCFs
#!/bin/bash

input="/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/slim-simulations-new-genome/simulations/vcfs/sorted"
output="/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/slim-simulations-new-genome/simulations/vcfs/no_dups_2"

mkdir -p "$output"

for VCF in "$input"/*sorted.vcf.gz; do
    echo "Removing duplicates from $VCF ..."

    # Clean up the weird suffix (.vcf.gz_sorted.vcf.gz → just _sorted)
    BASE=$(basename "$VCF" .sorted.vcf.gz)

    # Remove duplicate sites
    bcftools norm -d all -Oz \
        -o "$output/${BASE}_nodup.vcf.gz" "$VCF"

    # Index the cleaned file
    tabix -p vcf "$output/${BASE}_nodup.vcf.gz"

    echo "Finished $VCF → $output/${BASE}_nodup.vcf.gz"
done
```
I will run grenedalf in the files inside this folder no_dups

```bash
srun --pty --x11 --nodes=1 --cpus-per-task=4 --mem=32000 --time=10:00:00 /bin/bash
module load htslib/1.10.2    
module load bzip2/1.0.8 
#format line in vcf file
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed"> with pool-vcsf.py script (wrong name)

GRENEDALF=/gxfs_home/geomar/smomw573/software/grenedalf/bin/grenedalf
VCF=/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/slim-simulations-new-genome/simulations/vcfs/pooled 
GENOME=$WORK/smomw573/seasonal_adaptation/new_genome/Purple_Eaffinis_genome_assembly.fa \
OUTPUT=/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/slim-simulations-new-genome/simulations/vcfs/sync_2

$GRENEDALF sync --vcf-path $VCF --allow-file-overwriting  --out-dir ../sync_2
#loop
for VCF in /gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/slim-simulations-new-genome/simulations/vcfs/pooled/*.vcf; do
    BASE=$(basename "$VCF" .vcf)
    echo "Processing $VCF ..."

    $GRENEDALF sync \
        --vcf-path "$VCF" \
        --allow-file-overwriting \
        --out-dir "$OUTPUT" \
        --file-prefix "$BASE"

    echo "Finished $VCF → $OUTPUT/${BASE}.sync"
done

for VCF in /gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/slim-simulations-new-genome/simulations/vcfs/pooled/*.vcf; do
    BASE=$(basename "$VCF" .vcf)

    echo "Processing $VCF ..."

    $GRENEDALF sync \
        --vcf-path "$VCF" \
        --allow-file-overwriting \
        --out-dir "$OUTPUT" \
        --file-prefix "$BASE"

    echo
done


```

23.10.25

Now, as discussed with reid, i will combine all timepoints from one simulation run (eg_sim 78, 2009,2011,2015,2022), but in a way that the missing lines (because we excluded the snps) will have an NA. I will then substitute this NA for 0:0:0:0:0 as the other lines. This way all files will have the same number. I can then run the simulation on this combined files (100)

check how many columns dont match in a single simulation
```bash 
awk '
{
    file=FILENAME
    pos[$2][file]=1
    files[file]=1
}
END {
    total=0
    mismatch=0
    nfiles=0
    for (f in files) nfiles++
    for (p in pos) {
        count=0
        for (f in files) if (pos[p][f]) count++
        total++
        if (count < nfiles) mismatch++
    }
    print "Total unique positions:", total
    print "Positions missing in at least one file:", mismatch
    print "Positions present in all files:", total - mismatch
}' neutral_100_sampled_*_sorted_pooled_nodupsync.sync
```
now I want to merge all simulations (/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/slim-simulations-new-genome/simulations/vcfs/sync directory)

```bash
for i in {101..500}; do
    echo "Merging simulation $i..."

    awk -v i="$i" '
    FNR==NR { base[$1"\t"$2"\t"$3]=$4; next }  # 2009_early
    FNR!=NR {
        file=FILENAME
        val[file][$1"\t"$2"\t"$3]=$4
        allkeys[$1"\t"$2"\t"$3]=1
    }
    END {
        files[1]="neutral_"i"_sampled_2009_late_sorted_pooled.vcf_nodup_pooledsync.sync"
        files[2]="neutral_"i"_sampled_2011_early_sorted_pooled.vcf_nodup_pooledsync.sync"
        files[3]="neutral_"i"_sampled_2011_late_sorted_pooled.vcf_nodup_pooledsync.sync"
        files[4]="neutral_"i"_sampled_2015_early_sorted_pooled.vcf_nodup_pooledsync.sync"
        files[5]="neutral_"i"_sampled_2015_late_sorted_pooled.vcf_nodup_pooledsync.sync"
        files[6]="neutral_"i"_sampled_2022_early_sorted_pooled.vcf_nodup_pooledsync.sync"
        files[7]="neutral_"i"_sampled_2022_late_sorted_pooled.vcf_nodup_pooledsync.sync"

        for (key in allkeys) {
            split(key, arr, "\t")
            printf "%s\t%s\t%s\t%s", arr[1], arr[2], arr[3], (key in base)?base[key]:"NA"
            for (j=1; j<=7; j++) {
                f=files[j]
                printf "\t%s", (key in val[f])?val[f][key]:"NA"
            }
            print ""
        }
    }' \
    neutral_${i}_sampled_2009_early_sorted_pooled.vcf_nodup_pooledsync.sync \
    neutral_${i}_sampled_2009_late_sorted_pooled.vcf_nodup_pooledsync.sync \
    neutral_${i}_sampled_2011_early_sorted_pooled.vcf_nodup_pooledsync.sync \
    neutral_${i}_sampled_2011_late_sorted_pooled.vcf_nodup_pooledsync.sync \
    neutral_${i}_sampled_2015_early_sorted_pooled.vcf_nodup_pooledsync.sync \
    neutral_${i}_sampled_2015_late_sorted_pooled.vcf_nodup_pooledsync.sync \
    neutral_${i}_sampled_2022_early_sorted_pooled.vcf_nodup_pooledsync.sync \
    neutral_${i}_sampled_2022_late_sorted_pooled.vcf_nodup_pooledsync.sync \
    > /gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/slim-simulations-new-genome/simulations/vcfs/final_syncs/neutral_${i}_joined.sync

done

````

#this bit has to be skipped so I can run Reids script.
replace the NA's for 0:0:0:0:0:0 

```bash 
for f in /gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/slim-simulations-new-genome/simulations/vcfs/final_syncs/neutral_*_joined.sync; do
    sed -i 's/\bNA\b/0:0:0:0:0:0/g' "$f"
done
```

28.10.25

okay, now I will run baypass in all my 100 files.
I need to run the null (standart) to get the first file. 
But first I need to figure out how to run it and how to do the simulations.
Right now, we just want to simulate the different coverages across the time points.

29.10.2025
Meeting Reid 
rerun the analysis with vcftools to merge vcfs. Set missing data for homozygous (0/0) call.
verify in the sync files if the first position is the reference? make sure it is consistent.
until the sync step our coverages are all 200x. We need to simulate a binomial distribution of our coverages (mean coverages for all samples).
leave NA's in the final sync as NA's - /final_sync directory so they can be replaced as fixed during the simulation.

06.11.2025
Reid helped fixing the script (binomial sampling) to fix for the coverage issue in our samples. The script is in the folder, but also in the email. This was done as the NAs are actually invariant (not SNPs) in some of the samples in which they are snps.
```R
# modifed 2025-10-30 by RSB. 

library(poolSeq)

# List all neutral_*.sync files
sync_files_list <- list.files("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/slim-simulations-new-genome/simulations/vcfs/final_syncs",
                         pattern="neutral_.*_joined.sync$", full.names=TRUE)

for(file_in in sync_files_list){

  # Pick the file for this array task
  sync_file <- file_in
  
  cat("Processing file:", sync_file, "\n")
  # get sync name:
  sync_name <- basename(sync_file)
  rep_name <- gsub(".sync", "", sync_name)
  
  subSync <- read.sync(file=sync_file, 
                     polarization="reference",
                     gen=c(1,2,1,2,1,2,1,2), repl=c(1,1,2,2,3,3,4,4))


alleledf <- subSync@alleles

# there are some that are NA. these are actually invariant. 
# we polarized by reference allele when reading in. So set these NA to 1

# loop through each freq column and replace NAs
for(colidx in grep("freq", names(alleledf))) {
  col_name <- names(alleledf)[colidx]
  alleledf[is.na(alleledf[[colidx]]), (col_name) := 1]
}

# syntax of af:
# af(sync, chr, pos, repl, gen)

# get all allele freqs
af_all <- af(subSync, , , c(1,2,3,4), c(1,2))

#head(af_all)

# drop maf 
# @ jenny - what did you use?

# get mean allele freq
meanaf <- rowMeans(af_all)
# get minor allele freq
meanmaf <- ifelse(meanaf > 0.5, 1-meanaf, meanaf)

#sum(meanmaf > 0.01)
# about %75 of loci
# keep only those with minor allele freq > threshold here.
filtaf <- af_all[meanmaf > 0.01,]
nrow(filtaf)


#########################################
###now simulate the data for each year###
#########################################

cov_in <- 68 


  #binomial sampling, using the SLiM AF simulations as input. 
  # Those already have subsampled 100 individuals from the 200000Ne population.
  # basically just adding noise from sequencing

simulate_allele_freq <- function(af_in, cov_in){
  n_sites <- nrow(af_in)  
  baypass_out <- data.frame(matrix(ncol = 0, nrow = n_sites))  
  for(i in 1:ncol(af_in)){
    # sample ref allele
    sampled_ref <- rbinom(n = n_sites, 
                          size = cov_in, 
                          prob = af_in[,i])
    
    # get alt allele
    sampled_alt <- cov_in-sampled_ref
    
    # return alleles
    baypass_out[, paste0("pop", i, "_ref")] <- sampled_ref
    baypass_out[, paste0("pop", i, "_alt")] <- sampled_alt
  }
  return(baypass_out)

}
  
  
out <- simulate_allele_freq(af_in = filtaf, cov_in = 68)  

#write.table(out, file=paste('/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/baypass-new-genome/slim-neutral-baypass/',rep_name, '.binomAF.baypass',sep=""), 
#                            sep = '\t',
#                            col.names = FALSE, row.names = FALSE, quote=FALSE)

write.table(out, file=paste('/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/slim-simulations-new-genome/simulations/vcfs/baypass_output/',rep_name, '.binomAF.baypass',sep=""), 
            sep = '\t',
            col.names = FALSE, row.names = FALSE, quote=FALSE)

}
```


24.11.2025
Okay, now I managed to finish the 500 simulations with slim. I showed everything to Kwi and she said it is alright. I will now start with the baypass runs.
It is a bit messy with SLim cause I first did simulations 1-100 and then 101-500. So I had to create a few folders in the simulations directory (no_dups and no_dups_2 for example). When rerunning the code above, I need to fix this for each directory, but those are minor issues. 
The final 500 simulatins are in the sym_final directory

#25.11.2025

Okay, now I have 500 simulations, already corrected for sampling error using Reids script.
The next step is to run baypass in those simulations.
For this, I will need to run several baypass models. For today, I will read the baypass manual and aftewards run the first model (standard)

15.12.2025
Everything has been a bit slow because I am on holidays, but today I will run the baypass on the 500 files
first, I will run the core model.

```bash

#!/bin/bash
#SBATCH -D /gxfs_home/geomar/smomw573/work/seasonal_adaptation/log_files/
#SBATCH --mail-type=END
#SBATCH --mail-user=jnascimento@geomar.de
#SBATCH --partition=base
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=36G
#SBATCH --time=5:00:00
#SBATCH --job-name=baypass-core
#SBATCH --output=baypass_slim_%A_%a.out
#SBATCH --error=baypass_slim_%A_%a.err
#SBATCH --array=1-500%20
# Directories
inputdir=/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/baypass-new-genome/
inputdir2=/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/slim-simulations-new-genome/simulations/vcfs/baypass_output
outputdir=/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/baypass-new-genome/core-model-slim
baypassdir=/gxfs_home/geomar/smomw573/software/baypass_public-master/sources


# Pick the .baypass file for this array task
INPUT_FILE=$(ls $inputdir2/neutral_*_joined.binomAF.baypass | sort -V  | sed -n "${SLURM_ARRAY_TASK_ID}p")
BASENAME=$(basename $INPUT_FILE .baypass)

# Run BAYPASS core model
$baypassdir/g_baypass \
  -gfile $INPUT_FILE \
  -poolsizefile ${inputdir}/poolsize \
  -d0yij 40 \
  -npilot 25 \
  -outprefix ${outputdir}/${BASENAME}_core-model \
  -nthreads 4 \
  -seed ${SLURM_ARRAY_TASK_ID}
  ```


I also need to rerun baypass on the new vcf 
First, I need to get a ref and alternate count for the samples as an input for baypass

```bash
#run grenedalf
module load htslib/1.10.2    
module load bzip2/1.0.8 

GRENEDALF=/gxfs_home/geomar/smomw573/software/grenedalf/bin/grenedalf
VCF=/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/grenedalf-freq/new-genome-final-vcf/final-sorted.vcf
GENOME=$WORK/smomw573/seasonal_adaptation/new_genome/Purple_Eaffinis_genome_assembly.fa 
OUTPUT=/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/grenedalf-freq/new-genome-final-vcf

#fst_site queue
$GRENEDALF frequency \
  --vcf-path $VCF \
  --write-sample-counts \
  --allow-file-overwriting \
  --file-suffix _frequency_ref_alt_table\
  --out-dir $OUTPUT/
```

fix ouptut

awk 'NR>1 { $1=""; print }' frequency_ref_alt_table.csv | sed 's/^ //' | tr ',' ' ' > fix-genotype-baypass-original and became fixed-genotypes-baypass-original

now I run baypass 3 times:

```bash
#!/bin/bash
#SBATCH -D /gxfs_home/geomar/smomw573/work/seasonal_adaptation/log_files/
#SBATCH --mail-type=END
#SBATCH --mail-user=jnascimento@geomar.de
#SBATCH --qos=long
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=50G
#SBATCH --time=60:00:00
#SBATCH --job-name=baypass-core
#SBATCH --output=baypass_original_%A_%a.out
#SBATCH --error=baypass_original_%A_%a.err
#SBATCH --array=1-3
# Directories

# Directories
inputdir=/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/baypass-new-genome
outputdir=/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/baypass-new-genome/core-model-original
baypassdir=/gxfs_home/geomar/smomw573/software/baypass_public-master/sources

# Create a unique basename for output
BASENAME="original_population_run_${SLURM_ARRAY_TASK_ID}"

# Run BAYPASS core model
$baypassdir/g_baypass \
  -gfile ${inputdir}/fixed-genotypes-baypass-original \
  -poolsizefile ${inputdir}/poolsize \
  -d0yij 40 \
  -npilot 20 \
  -outprefix ${outputdir}/${BASENAME}_core-model \
  -nthreads 4 \
  -seed ${SLURM_ARRAY_TASK_ID}
```

07.01.2026

Core model of baypass is ready!
I will now run the follow up model using the matrix for both the real data and the simulated data!

c2 model for 
```bash
#!/bin/bash
#SBATCH -D /gxfs_home/geomar/smomw573/work/seasonal_adaptation/log_files/
#SBATCH --mail-type=END
#SBATCH --mail-user=jnascimento@geomar.de
#SBATCH --qos=long
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=50G
#SBATCH --time=60:00:00
#SBATCH --job-name=baypass-core
#SBATCH --output=baypass_c2original_%A_%a.out
#SBATCH --error=baypass_c2original_%A_%a.err
#SBATCH --array=1-3
# Directories

# Directories
inputdir=/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/baypass-new-genome
inputdir2=/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/baypass-new-genome/core-model-original/
outputdir=/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/baypass-new-genome/c2-model-original
baypassdir=/gxfs_home/geomar/smomw573/software/baypass_public-master/sources

# Create a unique basename for output
BASENAME="original_population_run_${SLURM_ARRAY_TASK_ID}"

# Run BAYPASS c2 model
$baypassdir/g_baypass \
  -gfile ${inputdir}/fixed-genotypes-baypass-original \
  -efile $inputdir/cov-baypass.txt \
  -omegafile $inputdir2/original_population_run_1_core-model_mat_omega.out \
  -poolsizefile ${inputdir}/poolsize \
  -contrastfile $inputdir/cov-baypass.txt
  -d0yij 40 \
  -outprefix ${outputdir}/${BASENAME}_c2-model \
  -nthreads 4 \
  -seed ${SLURM_ARRAY_TASK_ID}
```
now I need to do the same for the simulated data!

```bash
#!/bin/bash
#SBATCH -D /gxfs_home/geomar/smomw573/work/seasonal_adaptation/log_files/
#SBATCH --mail-type=END
#SBATCH --mail-user=jnascimento@geomar.de
#SBATCH --partition=base
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=36G
#SBATCH --time=5:00:00
#SBATCH --job-name=baypass-core
#SBATCH --output=baypass_slim_%A_%a.out
#SBATCH --error=baypass_slim_%A_%a.err
#SBATCH --array=1-500%20
# Directories
inputdir=/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/baypass-new-genome/
inputdir2=/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/slim-simulations-new-genome/simulations/vcfs/baypass_output
inputdir3=/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/baypass-new-genome/core-model-slim
outputdir=/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/baypass-new-genome/c2-model-slim
baypassdir=/gxfs_home/geomar/smomw573/software/baypass_public-master/sources


# Pick the .baypass file for this array task
INPUT_FILE=$(ls $inputdir2/neutral_*_joined.binomAF.baypass | sort -V | sed -n "${SLURM_ARRAY_TASK_ID}p")
BASENAME=$(basename $INPUT_FILE .baypass)
INPUT_FILE2=$(ls $inputdir3/neutral_*_joined.binomAF_core-model_mat_omega.out | sort -V  | sed -n "${SLURM_ARRAY_TASK_ID}p")

# Run BAYPASS c2 model
$baypassdir/g_baypass \
  -gfile $INPUT_FILE \
  -efile $inputdir/cov-baypass.txt \
  -omegafile $INPUT_FILE2 \
  -poolsizefile ${inputdir}/poolsize \
  -contrastfile $inputdir/cov-baypass.txt
  -d0yij 40 \
  -outprefix ${outputdir}/${BASENAME}_c2-model \
  -nthreads 4 \
  -seed ${SLURM_ARRAY_TASK_ID}

```
ls /gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/slim-simulations-new-genome/simulations/vcfs/baypass_output/neutral_*_joined.binomAF.baypass | sort -V 


ls /gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/baypass-new-genome/core-model-slim/neutral_*_joined.binomAF_core-model_mat_omega.out | sort -V 
cd 


16.01.2026
Today I will first check if the Omega matrices are similar from the core model in the original data

```R
source("/gxfs_home/geomar/smomw573/software/baypass_public-master/utils/baypass_utils.R")
omega1 <- as.matrix(read.table("original_population_run_1_core-model_mat_omega.out")) 
omega2 <- as.matrix(read.table("original_population_run_2_core-model_mat_omega.out")) 
omega3 <- as.matrix(read.table("original_population_run_3_core-model_mat_omega.out")) 

pop.names=c("2009_start",
        "2009_end",
        "2011_start",
        "2011_end",
        "2015_start",
        "2015_end",
        "2022_start",
        "2022_end"
        )
plot.omega(omega=omega1, names=pop.names)
#they are basically the same

require(corrplot)
31
cor.mat=cov2cor(omega)
corrplot(cor.mat,method="color",mar=c(2,1,2,2)+0.1,
main=expression("Correlation map based on"~hat(Omega)))
 
#
#> fmd.dist(omega1,omega2)
#[1] 0.04381255
#> fmd.dist(omega1,omega2,omega3)
#Error in fmd.dist(omega1, omega2, omega3) : unused argument (omega3)
#> fmd.dist(omega1,omega3)
#[1] 0.03421107
#> fmd.dist(omega3,omega2)
#[1] 0.04151861

```

All good.

Now I need to check the P values for each simulation.
For the C2 model, we look at the c2 values.
how to check if results converge?

```R

source("/gxfs_home/geomar/smomw573/software/baypass_public-master/utils/baypass_utils.R")

c2_1 <- read.table("original_population_run_1_c2-model_summary_contrast.out", h=TRUE)
c2_2 <- read.table("original_population_run_2_c2-model_summary_contrast.out", h=TRUE)
c2_3 <- read.table("original_population_run_3_c2-model_summary_contrast.out", h=TRUE)

c2v_1 <- c2_1$C2
c2v_2 <- c2_2$C2
c2v_3 <- c2_3$C2

c2_median <- apply(cbind(c2v_1, c2v_2, c2v_3), 1, median)

hist(10**(-1*lsa.ecotype.C2$log10.1.pval.),freq=F,breaks=50)
hist(10**(-1*lsa.ecotype.C2$log10.1.pval.),freq=F,breaks=50)
hist(10**(-1*lsa.ecotype.C2$log10.1.pval.),freq=F,breaks=50)

cor(c2v_1,c2v_2)
cor(c2v_1,c2v_3)
cor(c2v_3,c2v_2)

#they dont converge well, but not sure what I can do (70%)
M_C2_mat <- cbind(c2_1$M_C2, c2_2$M_C2, c2_3$M_C2)
# Median C2 across chains
M_C2_median <- apply(M_C2_mat, 1, median)

# Standard deviation across chains (optional, for reference)
M_C2_sd <- apply(M_C2_mat, 1, sd)

# Combine M_C2 and log10 p-values
M_C2_mat <- cbind(c2_1$M_C2, c2_2$M_C2, c2_3$M_C2)
logp_mat <- cbind(c2_1$log10.1.pval., c2_2$log10.1.pval., c2_3$log10.1.pval.)

# Function to get median and corresponding p-value
get_median_with_p <- function(c2_row, logp_row) {
  # Which chain gives the median
  idx <- which(c2_row == median(c2_row))[1]  # first if tie
  return(c(median = c2_row[idx], log10p = logp_row[idx]))
}

# Apply row-wise
res <- t(apply(cbind(M_C2_mat, logp_mat), 1, function(x) {
  c2_row <- x[1:3]
  logp_row <- x[4:6]
  get_median_with_p(c2_row, logp_row)
}))

# Build final table
c2_median_table <- data.frame(
  CONTRAST = c2_1$CONTRAST,
  MRK      = c2_1$MRK,
  M_C2_med = res[, "median"],
  log10.1.pval_med = res[, "log10p"]
)

head(c2_median_table)

# Plot median C2 across SNPs
plot(c2_median_table$M_C2_med, type="h",
     main="Median C2 per SNP",
     xlab="SNP index",
     ylab="Median C2",
     col="black")

# Convert log10 p-values to actual p-values
p_values <- 10^(-c2_median_table$log10.1.pval_med)

# Plot using points
plot(p_values, 
     pch=19,            # solid circles
     cex=0.5,           # smaller points
     col="black",
     xlab="SNP index",
     ylab="p-value",
     main="Median p-values per SNP")
     
# Optional: add threshold line, e.g., 0.001
abline(h=0.001, col="red", lty=2)

# Convert log10 p-value to actual p-value
c2_median_table$p_value <- 10^(-c2_median_table$log10.1.pval_med)

# Filter for p < 0.05
sig_snps <- subset(c2_median_table, p_value < 0.05)

# C2 points plot
plot(sig_snps$M_C2_med, 
     pch=19, 
     cex=0.7, 
     col="steelblue",
     xlab="Significant SNP index",
     ylab="Median C2",
     main="Median C2 for significant SNPs (p < 0.05)")

plot(-log10(sig_snps$p_value), 
     pch=19, 
     cex=0.7, 
     col="steelblue",
     xlab="Significant SNP index",
     ylab="-log10(p-value)",
     main="-log10 p-values for significant SNPs (p < 0.05)")
abline(h=-log10(0.05), col="red", lty=2)  # threshold line

# Filter significant SNPs (already done)
sig_snps <- subset(c2_median_table, p_value < 0.05)

# Histogram of p-values
hist(sig_snps$p_value,
     breaks=20,                  # number of bins
     col="skyblue",
     border="white",
     main="Histogram of p-values for significant SNPs",
     xlab="p-value",
     ylab="Frequency")


# check how many snps are significant

c2_median_table$p_value <- 10^(-c2_median_table$log10.1.pval_med)
# Total number of SNPs
total_snps <- nrow(c2_median_table)

# Number of significant SNPs
num_sig_snps <- nrow(subset(c2_median_table, p_value < 0.05))

# Percentage of significant SNPs
percent_sig <- (num_sig_snps / total_snps) * 100

cat("Total SNPs:", total_snps, "\n")
cat("Significant SNPs (p < 0.05):", num_sig_snps, "\n")
cat("Percentage of significant SNPs:", round(percent_sig, 2), "%\n")

#1.3 is significant
```

19.01.2026

So I reran the c2 model using the original data and trying out a more powerful repetition to see if the models converge a bit better

I will check this now

```r

source("/gxfs_home/geomar/smomw573/software/baypass_public-master/utils/baypass_utils.R") 
c2_1 <- read.table("original_population_run_1_c2-model_summary_contrast.out", h=TRUE)
c2_2 <- read.table("original_population_run_2_c2-model_summary_contrast.out", h=TRUE)
c2_3 <- read.table("original_population_run_3_c2-model_summary_contrast.out", h=TRUE)

c2v_1 <- c2_1$C2
c2v_2 <- c2_2$C2
c2v_3 <- c2_3$C2

c2_median <- apply(cbind(c2v_1, c2v_2, c2v_3), 1, median)

hist(10**(-1*lsa.ecotype.C2$log10.1.pval.),freq=F,breaks=50)
hist(10**(-1*lsa.ecotype.C2$log10.1.pval.),freq=F,breaks=50)
hist(10**(-1*lsa.ecotype.C2$log10.1.pval.),freq=F,breaks=50)

cor(c2v_1,c2v_2)
cor(c2v_1,c2v_3)
cor(c2v_3,c2v_2)
```
okay, this will have to be redone later this week as I am rerunning the models - needed more time to finish the runs

I will now check the SNPs in the baypass from the simulations

```R
# Path to the directory containing the files
simulations <- "/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/baypass-new-genome/c2-model-slim"

# List all neutral C2 files
slim <- list.files(
  path = simulations,
  pattern = "neutral_.*_c2-model_summary_contrast.out",
  full.names = TRUE
)

length(slim)
get_percent_significant <- function(file, alpha = 0.05) {
  
  dat <- read.table(file, header = TRUE)
  
  # extract log10(1/p) values
  logp <- as.numeric(dat$log10.1.pval.)
  
  # convert to p-values
  pvals <- 10^(-logp)
  
  # percentage significant
  mean(pvals < alpha) * 100
}

percent_sig_slim <- sapply(slim, get_percent_significant)

length(percent_sig_slim)


hist(percent_sig_slim,
     breaks = 30,
     col = "lightgray",
     border = "white",
     main = "Neutral simulations (SLiM)\n% significant SNPs (p < 0.05)",
     xlab = "Percentage of significant SNPs")

extract_significant_snps <- function(file, alpha = 0.05) {
  
  dat <- read.table(file, header = TRUE)
  
  # Convert log10(1/p) to p-values
  pvals <- 10^(-dat$log10.1.pval.)
  
  # Identify significant SNPs
  sig_idx <- which(pvals < alpha)
  
  if (length(sig_idx) == 0) return(NULL)
  
  # Build output table
  out <- dat[sig_idx, ]
  out$p_value <- pvals[sig_idx]
  out$simulation <- basename(file)
  
  return(out)
}
sig_snps_list <- lapply(slim, extract_significant_snps)


sig_snps_all <- do.call(rbind, sig_snps_list)


dim(sig_snps_all)
head(sig_snps_all)

write.table(
  sig_snps_all,
  file = "slim_all_significant_snps_p05.txt",
  row.names = FALSE,
  quote = FALSE
)
sig_snps_all$log10p <- -log10(sig_snps_all$p_value)

# numeric simulation index for plotting
sig_snps_all$sim_id <- as.numeric(
  factor(sig_snps_all$simulation)
)
plot(
  sig_snps_all$sim_id,
  sig_snps_all$log10p,
  pch = 16,
  cex = 0.6,
  col = "black",
  xlab = "Simulation",
  ylab = expression(-log[10](p)),
  main = "Significant SNPs (p < 0.05)\nAcross neutral simulations"
)

# significance threshold
abline(h = -log10(0.05), col = "red", lty = 2)
```


okay, so I saw that the overall significance is similar in both. 

I will now see the distribution of values, or how extreme those are, and then compare to the real data baypass.

```r
sig_snps_sorted <- sig_snps_all[order(sig_snps_all$p_value), ]
n <- nrow(sig_snps_sorted)

quantiles <- c(0.5, 0.25, 0.10, 0.05, 0.03, 0.01)
names(quantiles) <- c("top50", "top25", "top10", "top5", "top3", "top1")

# List to store results
extreme_snps_list <- list()

for (q in names(quantiles)) {
  n_snps <- ceiling(n * quantiles[q])
  extreme_snps_list[[q]] <- sig_snps_sorted[1:n_snps, ]
}

for (q in names(extreme_snps_list)) {
  write.table(
    extreme_snps_list[[q]],
    file = paste0("sig_snps_", q, ".txt"),
    row.names = FALSE,
    quote = FALSE
  )
}

sapply(extreme_snps_list, nrow)



# Combine quantiles into one data frame for plotting
plot_data <- do.call(rbind, lapply(names(extreme_snps_list), function(q) {
  df <- extreme_snps_list[[q]]
  df$quantile <- q
  df$log10p_plot <- df$log10p
  df$log10p_plot[!is.finite(df$log10p_plot)] <- max_finite + 1
  return(df)
}))


# Replace Inf with a large number for plotting
max_finite <- max(sig_snps_sorted$log10p[is.finite(sig_snps_sorted$log10p)])
sig_snps_sorted$log10p_plot <- sig_snps_sorted$log10p
sig_snps_sorted$log10p_plot[!is.finite(sig_snps_sorted$log10p_plot)] <- max_finite + 1

# Simple dot plot
library(ggplot2)

ggplot(plot_data, aes(x = quantile, y = log10p_plot)) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 1) +
  geom_hline(yintercept = -log10(0.05), color = "red", linetype = "dashed", size = 1) +
  labs(
    x = "Quantile of extreme SNPs",
    y = expression(-log[10](p)),
    title = "Extremeness of significant SNPs across quantiles"
  ) +
  theme_minimal()


ggplot(plot_data, aes(x = quantile, y = log10p_plot)) +
  geom_boxplot(fill = "lightblue", outlier.size = 1) +
  geom_hline(yintercept = -log10(0.05), color = "red", linetype = "dashed", size = 1) +
  labs(
    x = "Quantile of extreme SNPs",
    y = expression(-log[10](p)),
    title = "Distribution of extremeness across SNP quantiles"
  ) +
  theme_minimal()

  > sapply(extreme_snps_list, nrow)
 top50  top25  top10   top5   top3   top1
171719  85860  34344  17172  10304   3435

# Get top 1% SNPs
top1_snps <- extreme_snps_list$top1

# Replace Inf values for plotting
max_finite <- max(sig_snps_sorted$log10p[is.finite(sig_snps_sorted$log10p)])
top1_snps$log10p_plot <- top1_snps$log10p
top1_snps$log10p_plot[!is.finite(top1_snps$log10p_plot)] <- max_finite + 1

hist(
  top1_snps$log10p_plot,
  breaks = 30,
  col = "lightcoral",
  border = "white",
  main = "Histogram of -log10(p) for top 1% SNPs",
  xlab = expression(-log[10](p)),
  ylab = "Frequency"
)

# Significance threshold line
abline(v = -log10(0.05), col = "red", lty = 2, lwd = 2)

# Top 5% SNPs
top5_snps <- extreme_snps_list$top5

# Replace Inf values for plotting
max_finite <- max(sig_snps_sorted$log10p[is.finite(sig_snps_sorted$log10p)])
top5_snps$log10p_plot <- top5_snps$log10p
top5_snps$log10p_plot[!is.finite(top5_snps$log10p_plot)] <- max_finite + 1

hist(
  top5_snps$log10p_plot,
  breaks = 30,
  col = "lightgreen",
  border = "white",
  main = "Histogram of -log10(p) for top 5% SNPs",
  xlab = expression(-log[10](p)),
  ylab = "Frequency"
)
boxplot(
  top5_snps$log10p_plot,
  main = "Boxplot of -log10(p) for top 5% SNPs",
  ylab = expression(-log[10](p)),
  col = "lightgreen",
  border = "darkgreen"
)


```
I also produced these files
[4001] "sig_snps_top1.txt"
[4002] "sig_snps_top10.txt"
[4003] "sig_snps_top25.txt"
[4004] "sig_snps_top3.txt"
[4005] "sig_snps_top5.txt"
[4006] "sig_snps_top50.txt"
[4007] "slim_all_significant_snps_p05.txt"


Yes! now we got to the most exciting bit! lets find those snpppppsssssssssssssss!

I will start with the simulations

First I will try with one of the simulations

```R
library(tidyr)
libryr(ggplot2)
library(dyplr)

df <- read.table("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/baypass-new-genome/c2-model-slim/neutral_499_joined.binomAF_c2-model_summary_contrast.out", h=TRUE)
#first add an extra table converting the log.1.pval. into normal pvalues.

df$pvals <- 10^(-df$log10.1.pval.)

df_sig <- df[df$pvals < 0.05, ]
nrow(df_sig)


#load corresponding gfile
gfile <- read.table("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/slim-simulations-new-genome/simulations/vcfs/baypass_output/neutral_499_joined.binomAF.baypass", header=F)

# Sanity check
stopifnot(ncol(gfile) %% 2 == 0)

# Indices of REF and ALT columns
ref_idx <- seq(1, ncol(gfile), by = 2)
alt_idx <- seq(2, ncol(gfile), by = 2)

colnames(gfile) <- c(
  "2009_start_REF", "2009_start_ALT",
  "2009_end_REF",   "2009_end_ALT",
  "2011_start_REF",       "2011_start_ALT",
  "2011_end_REF",       "2011_end_ALT",
  "2015_start_REF",       "2015_start_ALT",
  "2015_end_REF",       "2015_end_ALT",
  "2022_start_REF",       "2022_start_ALT",
  "2022_end_REF",       "2022_end_ALT"
)

# Calculate ALT allele frequencies
alt_freq <- gfile[, alt_idx] / (gfile[, ref_idx] + gfile[, alt_idx])

# Name columns (remove _REF/_ALT)
colnames(alt_freq) <- sub("_ALT$", "", colnames(gfile)[alt_idx])

# Convert to data frame
alt_freq <- as.data.frame(alt_freq)

# Inspect
head(alt_freq)
alt_freq$marker <- 1:nrow(alt_freq)

# Move it to the first column (optional)
alt_freq <- alt_freq[, c("marker", setdiff(names(alt_freq), "marker"))]

# Check
head(alt_freq)


#save file
outdir <- "/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/baypass-new-genome/ALT-freq-slim"

write.table(
  alt_freq,
  file = file.path(outdir, "alt_freq_slim_499.txt"),
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)


#pull the significant snsps from the AF table

# Assume df_sig$MRK contains significant SNP indices
sig_markers <- df_sig$MRK

# Subset alt_freq
alt_sig <- alt_freq[alt_freq$marker %in% sig_markers, ]

#OKAY, now is the time
# Select only the frequency columns
freq_cols <- names(alt_sig)[-1]  

# remove markers
start_cols <- freq_cols[grep("_start", freq_cols)]
end_cols   <- freq_cols[grep("_end", freq_cols)]
# Logical vectors
is_positive <- apply(alt_sig, 1, function(x) {
  all(x[end_cols] > x[start_cols])
})

is_negative <- apply(alt_sig, 1, function(x) {
  all(x[end_cols] < x[start_cols])
})

# Subset datasets
fluctuating_positive <- alt_sig[is_positive, ]
fluctuating_negative <- alt_sig[is_negative, ]


colnames(fluctuating_positive) <- c("marker", "2009.start", "2009.end",
                   "2011.start", "2011.end",
                   "2015.start", "2015.end",
                   "2022.start", "2022.end")

colnames(fluctuating_negative) <- c("marker", "2009.start", "2009.end",
                   "2011.start", "2011.end",
                   "2015.start", "2015.end",
                   "2022.start", "2022.end")

norm_positive <- fluctuating_positive %>%
  mutate(
    '2009.end' = `2009.end` - `2009.start`,
    '2011.end' = `2011.end` - `2011.start`,
    '2015.end' = `2015.end` - `2015.start`,
    '2022.end' = `2022.end` - `2022.start`
  )

norm_positive <- norm_positive %>%
  mutate(
    # for start columns, change relative to 2009
    `2009.start` = 0,
    `2011.start` = 0,
    `2015.start` = 0,
    `2022.start` = 0
  )

norm_long <- norm_positive %>%
  pivot_longer(
    cols = -marker,
    names_to = "timepoint",
    values_to = "af"
  ) %>%
  mutate(timepoint = factor(timepoint,
                            levels = c("2009.start","2009.end",
                                       "2011.start","2011.end",
                                       "2015.start","2015.end",
                                       "2022.start","2022.end")))

ggplot(norm_long, aes(x = timepoint, y = af, group = marker, color = marker)) +
  geom_line(alpha = 0.4) +
  geom_point(size = 1, position = position_jitter(width = 0.1, height = 0)) +
  theme_minimal() +
  labs(
    title = "Normalized AF Changes (Starts = 0)",
    x = "Timepoint",
    y = "AF"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(color = "none") +
  scale_y_continuous(limits = c(0, 1))

  #same for negative

colnames(fluctuating_negative) <- c("marker", "2009.start", "2009.end",
                   "2011.start", "2011.end",
                   "2015.start", "2015.end",
                   "2022.start", "2022.end")


norm_negative <- fluctuating_negative %>%
  mutate(
    '2009.end' = (`2009.end` - `2009.start`) * -1,
    '2011.end' = (`2011.end` - `2011.start`) * -1,
    '2015.end' = (`2015.end` - `2015.start`) * -1,
    '2022.end' = (`2022.end` - `2022.start`) * -1
  )

norm_negative <- norm_negative %>%
  mutate(
    # for start columns, change relative to 2009
    `2009.start` = 0,
    `2011.start` = 0,
    `2015.start` = 0,
    `2022.start` = 0
  )

norm_long <- norm_negative %>%
  pivot_longer(
    cols = -marker,
    names_to = "timepoint",
    values_to = "af"
  ) %>%
  mutate(timepoint = factor(timepoint,
                            levels = c("2009.start","2009.end",
                                       "2011.start","2011.end",
                                       "2015.start","2015.end",
                                       "2022.start","2022.end")))

ggplot(norm_long, aes(x = timepoint, y = af, group = marker, color = marker)) +
  geom_line(alpha = 0.4) +
  geom_point(size = 1, position = position_jitter(width = 0.1, height = 0)) +
  theme_minimal() +
  labs(
    title = "Normalized AF Changes (Starts = 0)",
    x = "Timepoint",
    y = "AF"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(color = "none") +
  scale_y_continuous(limits = c(0, 1))


#finally I will combine norm_positive and norm_negative and make the stats
norm_all <- bind_rows(norm_positive, norm_negative)
norm_all <- norm_all %>%
  mutate(avg_change = (`2009.end` + `2011.end` + `2015.end` + `2022.end`) / 4)

head(norm_all)

outdir_stats <- "/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/baypass-new-genome/ALT-freq-slim"
write.table(
  norm_all,
  file = file.path(
    outdir_stats,
    paste0("norm_all_slim_", run_id, ".txt")
  ),
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)

median_value <- median(norm_all$avg_change, na.rm = TRUE)

outdir_stats <- "/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/baypass-new-genome/ALT-freq-slim"

write.table(
  median_table,
  file = file.path(outdir_stats, "median_avg_change_slim.txt"),
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)


p_box <- ggplot(norm_all, aes(y = avg_change)) +
  geom_boxplot(fill = "steelblue", alpha = 0.6) +
  theme_minimal() +
  labs(
    title = "Distribution of Average Changes",
    y = "Average Change"
  ) +
  geom_hline(yintercept = median_value, color = "red", linetype = "dashed") +
  annotate("text", x = 1.2, y = median_value, label = paste("Median =", round(median_value, 4)), color = "red", hjust = 0)
ggplot(norm_all, aes(x = avg_change)) +
  geom_histogram(binwidth = 0.02, fill = "steelblue", color = "black", alpha = 0.7) +
  theme_minimal() +
  labs(
    title = "Histogram of Average Changes",
    x = "Average Change",
    y = "Count"
  ) +
  geom_vline(xintercept = median(norm_all$avg_change, na.rm = TRUE),
             color = "red", linetype = "dashed") +
  annotate("text", x = median(norm_all$avg_change, na.rm = TRUE) + 0.02, 
           y = 5, label = paste("Median =", round(median(norm_all$avg_change, 4))),
           color = "red")

ggsave(
  filename = "average_change_slim_499_histogram.pdf",
  plot = p_hist,
  width = 7,
  height = 5
)

#now a plot of all changes

  # Keep only the end columns
end_long <- norm_all %>%
  select(marker, `2009.end`, `2011.end`, `2015.end`, `2022.end`) %>%
  pivot_longer(
    cols = -marker,
    names_to = "year",
    values_to = "end_value"
  )

d <-   ggplot(end_long, aes(x = end_value)) +
  geom_histogram(binwidth = 0.02, fill = "tomato", color = "black", alpha = 0.7) +
  theme_minimal() +
  labs(
    title = "Histogram of All End Values",
    x = "End Value",
    y = "Count"
  ) +
  geom_vline(xintercept = median(end_long$end_value, na.rm = TRUE),
             color = "blue", linetype = "dashed") +
  annotate("text", x = median(end_long$end_value, na.rm = TRUE) + 0.02, 
           y = 5, label = paste("Median =", round(median(end_long$end_value, 4))),
           color = "blue")

ggsave(
  filename = "all_changes_slim_499_histogram.pdf",
  plot = d,
  width = 7,
  height = 5
)

median_value_all <- median(end_long$end_value, na.rm = TRUE)

```

try to automize it

```R
library(tidyr)
library(ggplot2)
library(dplyr)

# Directories
indir <- "/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/baypass-new-genome/c2-model-slim"
gfile_dir <- "/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/slim-simulations-new-genome/simulations/vcfs/baypass_output"
outdir_freq <- "/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/baypass-new-genome/ALT-freq-slim"

# Initialize a table to store medians for all simulations
median_table <- data.frame(simulation = integer(), median_avg_change = numeric())

for (sim in 1:500) {

  message("Processing simulation ", sim, "...")

  # File paths
  df_file <- file.path(indir, paste0("neutral_", sim, "_joined.binomAF_c2-model_summary_contrast.out"))
  gfile_file <- file.path(gfile_dir, paste0("neutral_", sim, "_joined.binomAF.baypass"))

  # Read data
  df <- read.table(df_file, header = TRUE)
  df$pvals <- 10^(-df$log10.1.pval.)
  df_sig <- df[df$pvals < 0.05, ]

  if (nrow(df_sig) == 0) {
    message("No significant SNPs in simulation ", sim, ", skipping...")
    next
  }

  gfile <- read.table(gfile_file, header = FALSE)
  stopifnot(ncol(gfile) %% 2 == 0)

  # Indices of REF and ALT columns
  ref_idx <- seq(1, ncol(gfile), by = 2)
  alt_idx <- seq(2, ncol(gfile), by = 2)

  # Set column names
  colnames(gfile) <- c(
    "2009_start_REF", "2009_start_ALT",
    "2009_end_REF", "2009_end_ALT",
    "2011_start_REF", "2011_start_ALT",
    "2011_end_REF", "2011_end_ALT",
    "2015_start_REF", "2015_start_ALT",
    "2015_end_REF", "2015_end_ALT",
    "2022_start_REF", "2022_start_ALT",
    "2022_end_REF", "2022_end_ALT"
  )

  # Calculate ALT allele frequencies
  alt_freq <- gfile[, alt_idx] / (gfile[, ref_idx] + gfile[, alt_idx])
  colnames(alt_freq) <- sub("_ALT$", "", colnames(gfile)[alt_idx])
  alt_freq <- as.data.frame(alt_freq)

  # Add marker column
  alt_freq$marker <- 1:nrow(alt_freq)
  alt_freq <- alt_freq[, c("marker", setdiff(names(alt_freq), "marker"))]

  # Save ALT frequencies
  write.table(
    alt_freq,
    file = file.path(outdir_freq, paste0("alt_freq_slim_", sim, ".txt")),
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE
  )

  # Subset significant SNPs
  sig_markers <- df_sig$MRK
  alt_sig <- alt_freq[alt_freq$marker %in% sig_markers, ]

  # Identify start/end columns
  freq_cols <- names(alt_sig)[-1]
  start_cols <- freq_cols[grep("_start", freq_cols)]
  end_cols <- freq_cols[grep("_end", freq_cols)]

  # Logical vectors for positive/negative changes
  is_positive <- apply(alt_sig, 1, function(x) all(x[end_cols] > x[start_cols]))
  is_negative <- apply(alt_sig, 1, function(x) all(x[end_cols] < x[start_cols]))

  # Subset
  fluctuating_positive <- alt_sig[is_positive, ]
  fluctuating_negative <- alt_sig[is_negative, ]

  # Rename columns
  colnames(fluctuating_positive) <- c("marker", "2009.start", "2009.end",
                                      "2011.start", "2011.end",
                                      "2015.start", "2015.end",
                                      "2022.start", "2022.end")
  colnames(fluctuating_negative) <- colnames(fluctuating_positive)

  # Normalize positive
  norm_positive <- fluctuating_positive %>%
    mutate(
      '2009.end' = `2009.end` - `2009.start`,
      '2011.end' = `2011.end` - `2011.start`,
      '2015.end' = `2015.end` - `2015.start`,
      '2022.end' = `2022.end` - `2022.start`,
      `2009.start` = 0, `2011.start` = 0, `2015.start` = 0, `2022.start` = 0
    )

  # Normalize negative
  norm_negative <- fluctuating_negative %>%
    mutate(
      '2009.end' = (`2009.end` - `2009.start`) * -1,
      '2011.end' = (`2011.end` - `2011.start`) * -1,
      '2015.end' = (`2015.end` - `2015.start`) * -1,
      '2022.end' = (`2022.end` - `2022.start`) * -1,
      `2009.start` = 0, `2011.start` = 0, `2015.start` = 0, `2022.start` = 0
    )

  # Combine
  norm_all <- bind_rows(norm_positive, norm_negative) %>%
    mutate(avg_change = (`2009.end` + `2011.end` + `2015.end` + `2022.end`) / 4)

  # Save normalized table
  write.table(
    norm_all,
    file = file.path(outdir_freq, paste0("norm_all_slim_", sim, ".txt")),
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE
  )

  # Compute median
  median_sim <- median(norm_all$avg_change, na.rm = TRUE)
  median_table <- rbind(median_table, data.frame(simulation = sim, median_avg_change = median_sim))

  # Boxplot
  p_box <- ggplot(norm_all, aes(y = avg_change)) +
    geom_boxplot(fill = "steelblue", alpha = 0.6) +
    theme_minimal() +
    labs(title = paste0("Distribution of Average Changes - Simulation ", sim),
         y = "Average Change") +
    geom_hline(yintercept = median_sim, color = "red", linetype = "dashed") +
    annotate("text", x = 1.2, y = median_sim, label = paste("Median =", round(median_sim, 4)), color = "red", hjust = 0)

  ggsave(
    filename = file.path(outdir_freq, paste0("average_change_slim_", sim, "_boxplot.pdf")),
    plot = p_box,
    width = 7, height = 5
  )

  # Histogram
  p_hist <- ggplot(norm_all, aes(x = avg_change)) +
    geom_histogram(binwidth = 0.02, fill = "steelblue", color = "black", alpha = 0.7) +
    theme_minimal() +
    labs(title = paste0("Histogram of Average Changes - Simulation ", sim),
         x = "Average Change", y = "Count") +
    geom_vline(xintercept = median_sim, color = "red", linetype = "dashed") +
    annotate("text", x = median_sim + 0.02, y = 5, label = paste("Median =", round(median_sim, 4)), color = "red")

  ggsave(
    filename = file.path(outdir_freq, paste0("average_change_slim_", sim, "_histogram.pdf")),
    plot = p_hist,
    width = 7, height = 5
  )
}

# Save median table for all simulations
write.table(
  median_table,
  file = file.path(outdir_freq, "median_avg_change_all_sims.txt"),
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)
```
and final plot
```R
all_end_values <- data.frame()

for (sim in 1:500) {
  norm_file <- paste0("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/baypass-new-genome/ALT-freq-slim/norm_all_slim_", sim, ".txt")
  if (!file.exists(norm_file)) next
  
  norm_all <- read.table(norm_file, header = TRUE, sep = "\t", check.names = TRUE) # check.names=TRUE is default
  
  # Get all end columns (works with X2009.end, X2011.end, etc.)
  end_cols <- grep("end$", colnames(norm_all), value = TRUE)
  
  end_long <- norm_all %>%
    select(marker, all_of(end_cols)) %>%
    pivot_longer(
      cols = -marker,
      names_to = "year",
      values_to = "end_value"
    ) %>%
    mutate(simulation = sim)
  
  all_end_values <- bind_rows(all_end_values, end_long)
}

# Plot histogram
p_all_end <- ggplot(all_end_values, aes(x = end_value)) +
  geom_histogram(binwidth = 0.02, fill = "tomato", color = "black", alpha = 0.7) +
  theme_minimal() +
  labs(
    title = "Histogram of All End Values Across Simulations",
    x = "End Value",
    y = "Count"
  ) +
  geom_vline(xintercept = median(all_end_values$end_value, na.rm = TRUE),
             color = "blue", linetype = "dashed") +
  annotate(
    "text",
    x = median(all_end_values$end_value, na.rm = TRUE) + 0.02,
    y = 50,
    label = paste("Global Median =", round(median(all_end_values$end_value), 4)),
    color = "blue"
  )

# Save PDF
ggsave(
  filename = "/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/baypass-new-genome/ALT-freq-slim/all_simulations_end_values_histogram.pdf",
  plot = p_all_end,
  width = 7,
  height = 5
)
```

03.02.2026
Now I have to do the same for the real data!

```R
library(tidyr)
library(ggplot2)
library(dplyr)

#first I need to pick the output that has the median values.

df1 <- read.table("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/baypass-new-genome/c2-model-original/original_population_run_1_c2-model_summary_contrast.out", h=TRUE)

df2 <- read.table("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/baypass-new-genome/c2-model-original/original_population_run_2_c2-model_summary_contrast.out", h=TRUE)

df3 <- read.table("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/baypass-new-genome/c2-model-original/original_population_run_3_c2-model_summary_contrast.out", h=TRUE)

combined_df <- df1 %>%
  select(CONTRAST, MRK,
         C2_std_df1 = C2_std,
         log10_p_df1 = log10.1.pval.) %>%
  left_join(
    df2 %>%
      select(CONTRAST, MRK,
             C2_df2 = M_C2,
             log10_p_df2 = log10.1.pval.),
    by = c("CONTRAST", "MRK")
  ) %>%
  left_join(
    df3 %>%
      select(CONTRAST, MRK,
             C2_df3 = M_C2,
             log10_p_df3 = log10.1.pval.),
    by = c("CONTRAST", "MRK")
  )

#median

combined_df <- combined_df %>%
  mutate(
    median_C2 = apply(
      select(., C2_std_df1, C2_df2, C2_df3),
      1,
      median,
      na.rm = TRUE
    ),
    median_log10_p = apply(
      select(., log10_p_df1, log10_p_df2, log10_p_df3),
      1,
      median,
      na.rm = TRUE
    )
  )
#first add an extra table converting the log.1.pval. into normal pvalues.

combined_df <- combined_df %>%
  mutate(
    median_pvalue = 10^(-median_log10_p)
  )

significant_df <- combined_df %>%
  filter(median_pvalue < 0.05)

#> nrow(significant_df)
# 92360

# Plot
p <- ggplot(significant_df, aes(x = median_C2)) +
  geom_histogram(bins = 50, fill = "steelblue", color = "black") +
  theme_minimal() +
  xlab("Median C2") +
  ylab("Count") +
  ggtitle("Distribution of Median C2 for Significant SNPs")

# Show the plot
print(p)

# Save as PNG
ggsave("histogram_median_C2.png", plot = p, width = 6, height = 4, dpi = 300)

p <- ggplot(significant_df, aes(x = median_pvalue)) +
  geom_histogram(bins = 50, fill = "steelblue", color = "black") +
  theme_minimal() +
  xlab("Median pvalue") +
  ylab("Count") +
  ggtitle("Distribution of Median pvalue for Significant SNPs")

# Show the plot
print(p)

# Save as PNG
ggsave("histogram_median_pvalue.png", plot = p, width = 6, height = 4, dpi = 300)


#okay, now I want to filter for the ones which are fluctuating


freq_df <- read.table("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/grenedalf-freq/new-genome-final-vcf/frequency_alt.csv", header=TRUE)

# Add MRK column
freq_df$MRK <- 1:nrow(freq_df)

# Check
head(freq_df)

# Keep only SNPs that are in significant_df
freq_signif <- freq_df %>%
  filter(MRK %in% significant_df$MRK)

# Quick check
nrow(freq_signif)
head(freq_signif)

library(dplyr)

# Keep only the columns you want
freq_signif_filtered <- freq_signif %>%
  select(
    EA_2009_T2.FREQ, EA_2009_T4.FREQ,
    EA_2011_T1.FREQ, EA_2011_T2.FREQ,
    EA_2015_T1.FREQ, EA_2015_T4.FREQ,
    EA_2022_T1.FREQ, EA_2022_T4.FREQ, MRK
  )

freq_signif_filtered <- freq_signif_filtered %>%
  select(MRK, everything())

# Extract marker IDs
sig_markers <- freq_signif_filtered$MRK

alt_sig <- freq_signif_filtered %>%
  filter(MRK %in% sig_markers)

alt_sig_diff <- alt_sig %>%
  mutate(
    diff_2009 = as.numeric(EA_2009_T4.FREQ) - as.numeric(EA_2009_T2.FREQ),
    diff_2011 = as.numeric(EA_2011_T2.FREQ) - as.numeric(EA_2011_T1.FREQ),
    diff_2015 = as.numeric(EA_2015_T4.FREQ) - as.numeric(EA_2015_T1.FREQ),
    diff_2022 = as.numeric(EA_2022_T4.FREQ) - as.numeric(EA_2022_T1.FREQ)
  )

# Positive SNPs: all differences > 0
is_positive <- apply(alt_sig_diff[, c("diff_2009","diff_2011","diff_2015","diff_2022")], 1, function(x) all(x > 0))

# Negative SNPs: all differences < 0
is_negative <- apply(alt_sig_diff[, c("diff_2009","diff_2011","diff_2015","diff_2022")], 1, function(x) all(x < 0))

fluctuating_positive <- alt_sig_diff[is_positive, ]
fluctuating_negative <- alt_sig_diff[is_negative, ]

# For positive SNPs
colnames(fluctuating_positive) <- c("MRK", "2009.start","2011.start","2015.start","2022.start",
                                    "2009.end","2011.end","2015.end","2022.end",
                                    "diff_2009","diff_2011","diff_2015","diff_2022")[1:ncol(fluctuating_positive)]

colnames(fluctuating_negative) <- c("MRK", 
                                    "2009.start","2011.start","2015.start","2022.start",
                                    "2009.end","2011.end","2015.end","2022.end",
                                    "diff_2009","diff_2011","diff_2015","diff_2022")[1:ncol(fluctuating_negative)]

norm_positive <- fluctuating_positive %>%
  mutate(
    # Start = 0
    `2009.start` = 0, `2011.start` = 0, `2015.start` = 0, `2022.start` = 0,
    # End = difference (already diff_* columns)
    `2009.end` = diff_2009,
    `2011.end` = diff_2011,
    `2015.end` = diff_2015,
    `2022.end` = diff_2022
  ) %>%
  # Keep only columns needed for plotting
  select(MRK,
         `2009.start`, `2009.end`,
         `2011.start`, `2011.end`,
         `2015.start`, `2015.end`,
         `2022.start`, `2022.end`)

# -------------------------
# 2️⃣ Pivot to long format
# -------------------------
library(tidyr)
norm_long_positive <- norm_positive %>%
  pivot_longer(
    cols = -MRK,
    names_to = "timepoint",
    values_to = "af"
  ) %>%
  mutate(
    # Factor levels ensure correct plotting order
    timepoint = factor(timepoint,
                       levels = c("2009.start","2009.end",
                                  "2011.start","2011.end",
                                  "2015.start","2015.end",
                                  "2022.start","2022.end"))
  )

# Check for NA
any(is.na(norm_long_positive$timepoint))  # should be FALSE

# -------------------------
# 3️⃣ Plot positive trajectories
# -------------------------
library(ggplot2)
ggplot(norm_long_positive, aes(x = timepoint, y = af, group = MRK, color = as.factor(MRK))) +
  geom_line(alpha = 0.3) +
#  geom_point(size = 1, position = position_jitter(width = 0.1, height = 0)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "Positive SNP Trajectories (Normalized, Start = 0)",
    x = "Timepoint",
    y = "AF Change"
  ) +
  guides(color = "none")



norm_negative <- fluctuating_negative %>%
  mutate(
    # Start = 0
    `2009.start` = 0, `2011.start` = 0, `2015.start` = 0, `2022.start` = 0,
    # End = -diff to flip to positive
    `2009.end` = -diff_2009,
    `2011.end` = -diff_2011,
    `2015.end` = -diff_2015,
    `2022.end` = -diff_2022
  ) %>%
  # Keep only columns needed for plotting
  select(MRK,
         `2009.start`, `2009.end`,
         `2011.start`, `2011.end`,
         `2015.start`, `2015.end`,
         `2022.start`, `2022.end`)

# -------------------------
# 2️⃣ Pivot to long format
# -------------------------
norm_long_negative <- norm_negative %>%
  pivot_longer(
    cols = -MRK,
    names_to = "timepoint",
    values_to = "af"
  ) %>%
  mutate(
    timepoint = factor(timepoint,
                       levels = c("2009.start","2009.end",
                                  "2011.start","2011.end",
                                  "2015.start","2015.end",
                                  "2022.start","2022.end"))
  )

# Check for NA
any(is.na(norm_long_negative$timepoint))  # should be FALSE


norm_long_all <- bind_rows(norm_long_positive, norm_long_negative)

# -------------------------
# 4️⃣ Plot all trajectories together
# -------------------------
ggplot(norm_long_all, aes(x = timepoint, y = af, group = MRK)) +
  geom_line(alpha = 0.3, color = "black") +
 # geom_point(size = 1, position = position_jitter(width = 0.1, height = 0), color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "Normalized SNP Trajectories (Positive and Negative)",
    x = "Timepoint",
    y = "AF Change"
  )

# okay, now lets look at the p-values of those snps

# Assuming 'significant_df' has columns: MRK, median_pvalue
pvals_positive <- significant_df %>%
  filter(MRK %in% norm_positive$MRK) %>%
  select(MRK, median_pvalue)

pvals_negative <- significant_df %>%
  filter(MRK %in% norm_negative$MRK) %>%
  select(MRK, median_pvalue)

# Quick summary of p-values
summary(pvals_positive$median_pvalue)
summary(pvals_negative$median_pvalue)

pvals_combined <- bind_rows(
  pvals_positive %>% mutate(type = "positive"),
  pvals_negative %>% mutate(type = "negative")
)

#plot p values fluctuating snps
ggplot(
  pvals_combined %>% mutate(logp = -log10(median_pvalue)),
  aes(x = logp)
) +
  geom_histogram(
    bins = 50,
    fill = "steelblue",
    color = "black",
    alpha = 0.7
  ) +
  theme_minimal() +
  labs(
    title = "Distribution of –log10(Median p-values)\n(Positive + Negative SNPs)",
    x = expression(-log[10](median~p)),
    y = "Count"
  )



ggplot(pvals_combined, aes(x = median_pvalue)) +
  geom_histogram(
    bins = 50,
    fill = "steelblue",
    color = "black",
    alpha = 0.7
  ) +
  theme_minimal() +
  labs(
    title = "Distribution of Median p-values\n(Positive + Negative SNPs)",
    x = "Median p-value",
    y = "Count"
  )


#median
median_all <- median(pvals_combined$median_pvalue, na.rm = TRUE)
median_all

0.02170757

mean_all <- median(pvals_combined$me_pvalue, na.rm = TRUE)
median_all

#top 1%

# Quantile thresholds
thr_5  <- quantile(pvals_combined$median_pvalue, 0.05, na.rm = TRUE)
thr_1  <- quantile(pvals_combined$median_pvalue, 0.01, na.rm = TRUE)

thr_5
thr_1
#         5%
#0.001754513
#          1%
#0.0003396386


top5 <- pvals_combined %>%
  filter(median_pvalue <= thr_5)

top1 <- pvals_combined %>%
  filter(median_pvalue <= thr_1)

nrow(top5)
nrow(top1)

#[1] 2596
#[1] 520

summary_pvals <- tibble(
  group = c("all", "top_5_percent", "top_1_percent"),
  mean_pvalue = c(
    mean(pvals_combined$median_pvalue, na.rm = TRUE),
    mean(top5$median_pvalue, na.rm = TRUE),
    mean(top1$median_pvalue, na.rm = TRUE)
  ),
  median_pvalue = c(
    median(pvals_combined$median_pvalue, na.rm = TRUE),
    median(top5$median_pvalue, na.rm = TRUE),
    median(top1$median_pvalue, na.rm = TRUE)
  )
)

summary_pvals
# A tibble: 3 × 3
#  group         mean_pvalue median_pvalue
#  <chr>               <dbl>         <dbl>
#1 all              0.0227        0.0217
#2 top_5_percent    0.000848      0.000834
#3 top_1_percent    0.000163      0.000165


summary_logp <- tibble(
  group = c("all", "top_5_percent", "top_1_percent"),
  mean_logp = c(
    mean(-log10(pvals_combined$median_pvalue), na.rm = TRUE),
    mean(-log10(top5$median_pvalue), na.rm = TRUE),
    mean(-log10(top1$median_pvalue), na.rm = TRUE)
  ),
  median_logp = c(
    median(-log10(pvals_combined$median_pvalue), na.rm = TRUE),
    median(-log10(top5$median_pvalue), na.rm = TRUE),
    median(-log10(top1$median_pvalue), na.rm = TRUE)
  )
)

summary_logp


#plot top 1% and 5%

ggplot(top5, aes(x = median_pvalue)) +
  geom_histogram(bins = 40, fill = "steelblue", color = "black") +
  theme_minimal() +
  labs(
    title = "Top 5% Median p-values",
    x = "Median p-value",
    y = "Count"
  )

ggplot(top1, aes(x = median_pvalue)) +
  geom_histogram(bins = 40, fill = "firebrick", color = "black") +
  theme_minimal() +
  labs(
    title = "Top 1% Median p-values",
    x = "Median p-value",
    y = "Count"
  )

###okay, now I will estimate the average changes in AF 

norm_all <- bind_rows(
  norm_positive,
  norm_negative
)

nrow(norm_all)
head(norm_all)
norm_all <- norm_all %>%
  mutate(
    mean_change = rowMeans(
      select(., `2009.end`, `2011.end`, `2015.end`, `2022.end`),
      na.rm = TRUE
    ),
    median_change = apply(
      select(., `2009.end`, `2011.end`, `2015.end`, `2022.end`),
      1,
      median,
      na.rm = TRUE
    )
  )

#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
#0.005785 0.020448 0.029022 0.036321 0.045809 0.327886
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
#0.002492 0.018525 0.026974 0.034065 0.043074 0.343678

#plot changes in AF of all fluctuating SNPs

ggplot(norm_all, aes(x = mean_change)) +
  geom_histogram(bins = 50, fill = "grey70", color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_minimal() +
  labs(
    title = "Distribution of mean allele-frequency change",
    x = "Mean AF change across years",
    y = "Number of SNPs"
  )

#now 1% and 5%
#add p value to norm_all
norm_all_p <- norm_all %>%
  left_join(
    significant_df %>% select(MRK, median_pvalue),
    by = "MRK"
  )

top5_p <- norm_all_p %>%
  filter(median_pvalue <= thr_5)

top1_p <- norm_all_p %>%
  filter(median_pvalue <= thr_1)

ggplot(top5_p, aes(x = mean_change)) +
  geom_histogram(bins = 40, fill = "steelblue", color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_minimal() +
  labs(
    title = "Mean AF change – SNPs with top 5% p-values",
    x = "Mean AF change",
    y = "Number of SNPs"
  )


ggplot(top1_p, aes(x = mean_change)) +
  geom_histogram(bins = 30, fill = "firebrick", color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_minimal() +
  labs(
    title = "Mean AF change – SNPs with top 1% p-values",
    x = "Mean AF change",
    y = "Number of SNPs"
  )


#Save pvalues fluctuating snps
write.table(norm_all, "p_values_fluctuating_snps.tsv", sep = "\t", row.names = FALSE, quote = FALSE)


#plot only those with a change higher than 10%
# Filter SNP trajectories with changes >= 0.1 at least in one timepoint
norm_long_all_filtered <- norm_long_all %>%
  group_by(MRK) %>%
  filter(max(af, na.rm = TRUE) >= 0.1) %>%
  ungroup()

ggplot(norm_long_all_filtered, aes(x = timepoint, y = af, group = MRK)) +
  geom_line(alpha = 0.3, color="black") +
 # geom_point(size = 1, position = position_jitter(width = 0.1, height = 0)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "Normalized SNP Trajectories (Changes >= 0.1)",
    x = "Timepoint",
    y = "AF Change",
    color = "Trajectory"
  )

# Count number of SNPs with at least one AF change >= 0.1
markers_above_threshold <- norm_long_all %>%
  group_by(MRK) %>%
  filter(max(af, na.rm = TRUE) >= 0.1) %>%
  summarise()  # keep only MRK column

n_markers <- nrow(markers_above_threshold)
n_markers

### okay, now only those which fluctuated more than 10%

fluctuating_positive <- alt_sig_diff %>%
  filter(diff_2009 > 0.1 & diff_2011 > 0.1 & diff_2015 > 0.1 & diff_2022 > 0.1)

fluctuating_negative <- fluctuating_negative %>%
  filter(abs(diff_2009) >= 0.1 & abs(diff_2011) >= 0.1 &
         abs(diff_2015) >= 0.1 & abs(diff_2022) >= 0.1)

# -------------------------
# 4️⃣ Normalize and flip negative SNPs
# -------------------------
norm_positive <- fluctuating_positive %>%
  mutate(
    `2009.start` = 0, `2009.end` = diff_2009,
    `2011.start` = 0, `2011.end` = diff_2011,
    `2015.start` = 0, `2015.end` = diff_2015,
    `2022.start` = 0, `2022.end` = diff_2022
  ) %>%
  select(MRK, `2009.start`,`2009.end`, `2011.start`,`2011.end`, `2015.start`,`2015.end`, `2022.start`,`2022.end`) 

norm_negative <- fluctuating_negative %>%
  mutate(
    `2009.start` = 0, `2009.end` = -diff_2009,
    `2011.start` = 0, `2011.end` = -diff_2011,
    `2015.start` = 0, `2015.end` = -diff_2015,
    `2022.start` = 0, `2022.end` = -diff_2022
  ) %>%
  select(MRK, `2009.start`,`2009.end`, `2011.start`,`2011.end`, `2015.start`,`2015.end`, `2022.start`,`2022.end`)

# -------------------------
# 5️⃣ Pivot long and combine
# -------------------------
norm_long_positive <- norm_positive %>%
  pivot_longer(-MRK, names_to = "timepoint", values_to = "af") %>%
  mutate(timepoint = factor(timepoint,
                            levels = c("2009.start","2009.end",
                                       "2011.start","2011.end",
                                       "2015.start","2015.end",
                                       "2022.start","2022.end")),
         trajectory = "Positive")

norm_long_negative <- norm_negative %>%
  pivot_longer(-MRK, names_to = "timepoint", values_to = "af") %>%
  mutate(timepoint = factor(timepoint,
                            levels = c("2009.start","2009.end",
                                       "2011.start","2011.end",
                                       "2015.start","2015.end",
                                       "2022.start","2022.end")),
         trajectory = "Negative")

# Combine
norm_long_all <- bind_rows(norm_long_positive, norm_long_negative)

ggplot(norm_long_all, aes(x = timepoint, y = af, group = MRK, color = trajectory)) +
  geom_line(alpha = 0.3) +
#  geom_point(size = 1, position = position_jitter(width = 0.1, height = 0)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Normalized SNP Trajectories (|Change| ≥ 0.1)",
       x = "Timepoint", y = "AF Change", color = "Trajectory")

# Assuming 'significant_df' has columns: MRK, median_pvalue
pvals_positive <- significant_df %>%
  filter(MRK %in% fluctuating_positive$MRK) %>%
  select(MRK, median_pvalue)

pvals_negative <- significant_df %>%
  filter(MRK %in% fluctuating_negative$MRK) %>%
  select(MRK, median_pvalue)

# Quick summary of p-values
summary(pvals_positive$median_pvalue)
summary(pvals_negative$median_pvalue)

#     Min.   1st Qu.    Median      Mean   3rd Qu.      Max.
#0.0001378 0.0003416 0.0005240 0.0027816 0.0007494 0.0143985
#     Min.   1st Qu.    Median      Mean   3rd Qu.      Max.
#1.000e-08 2.031e-04 8.363e-04 6.703e-03 4.101e-03 4.945e-02

# Combine positive and negative p-values
pvals_all <- bind_rows(
  pvals_positive %>% mutate(trajectory = "Positive"),
  pvals_negative %>% mutate(trajectory = "Negative")
)

# Histogram of p-values
ggplot(pvals_all, aes(x = median_pvalue)) +
  geom_histogram(binwidth = 0.005, color = "black", alpha = 0.7) +
  theme_minimal() +
  labs(
    title = "Histogram of Median p-values for Extreme SNPs (|ΔAF| ≥ 0.1)",
    x = "Median p-value",
    y = "Count",
    fill = "Trajectory"
  ) +
  scale_fill_manual(values = c("Positive" = "steelblue", "Negative" = "tomato"))

```
okay I had a meeting with Reid this week, and need to re-do a few things.
First correct q-values, then get a table with mean p values for 0.05 of the significant snps and 0.01 and 0.001 and then fluctuating snps. and then mean AF in each of those frequencies and then do the simulations... WOHOOO

lets do it in R

```R
install.packages("remotes")
remotes::install_github("jdstorey/qvalue")
library(qvalue)
library(dplyr)

df1 <- read.table("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/baypass-new-genome/c2-model-original/original_population_run_1_c2-model_summary_contrast.out", h=TRUE)

df2 <- read.table("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/baypass-new-genome/c2-model-original/original_population_run_2_c2-model_summary_contrast.out", h=TRUE)

df3 <- read.table("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/baypass-new-genome/c2-model-original/original_population_run_3_c2-model_summary_contrast.out", h=TRUE)

pval1 <- df1$log10.1.pval
#range(df1$log10.1.pval)
#plot(df1$log10.1.pval)

pval2 <- df2$log10.1.pval
pval3 <- df3$log10.1.pval


pval1.2 <- 10^(-pval1)
pval2.2 <- 10^(-pval2)
pval3.2 <- 10^(-pval3)

qobj1 <- qvalue(p = pval1.2)
qobj2 <- qvalue(p = pval2.2)
qobj3 <- qvalue(p = pval3.2)

pdf("q-values.pdf", width=8, height=6)
par(mfrow = c(3,2))
plot(qobj1$qvalues)
plot(qobj2$qvalues)
plot(qobj3$qvalues)



#now calculate the median and then backtransform to -log10
qvalues1 <- qobj1$qvalues
qvalues2 <- qobj2$qvalues
qvalues3 <- qobj3$qvalues


combinedq <- data.frame(q1=qvalues1, q2=qvalues2, q3=qvalues3)

combinedq <- cbind(marker = 1:nrow(combinedq), combinedq)
head(combinedq)


#check for convergence
q_columns <- grep("^q", names(combinedq), value = TRUE)

# Create an empty matrix to store correlation values
cor_matrix <- matrix(NA, nrow = length(q_columns), ncol = length(q_columns),
                    dimnames = list(q_columns, q_columns))

# Loop through the pairs of columns
for (i in seq_along(q_columns)) {
 for (j in seq_along(q_columns)) {
   cor_matrix[i, j] <- cor(x = combinedq[[q_columns[i]]], 
                           y = combinedq[[q_columns[j]]], 
                           method = "pearson", 
                           use = "complete.obs")
 }
}

# Print the correlation matrix
print(cor_matrix)
#0.68 for all

# median
head(combinedq)
combinedq$median_q <- apply(combinedq[, c("q1", "q2", "q3")], 1, median)

combinedq$log_transformed <- -log10(combinedq$median_q)

#save table
write.csv(
  combinedq,
  file = "/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/baypass-new-genome/c2-original-output-files/combinedq.csv",
  row.names = FALSE
)

#okay now I have the correct p values. lets redo the calculation of significant snps

#okay, now I want to filter for the ones which are fluctuating


freq_df <- read.table("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/grenedalf-freq/new-genome-final-vcf/frequency_alt.csv", header=TRUE)

# Add MRK column
freq_df$MRK <- 1:nrow(freq_df)

# Check
head(freq_df)

#select significant markers
sum(qvalues2 < 0.05)

sig_combinedq <- subset(combinedq, median_q <= 0.05)

# Keep only SNPs that are in significant_df
freq_signif <- freq_df %>%
  filter(MRK %in% sig_combinedq$marker)

#plot the trajectory of the n_snps
# Quick check
nrow(freq_signif)
head(freq_signif)

library(dplyr)

# Keep only the columns you want
freq_signif_filtered <- freq_signif %>%
  select(
    EA_2009_T2.FREQ, EA_2009_T4.FREQ,
    EA_2011_T1.FREQ, EA_2011_T2.FREQ,
    EA_2015_T1.FREQ, EA_2015_T4.FREQ,
    EA_2022_T1.FREQ, EA_2022_T4.FREQ, MRK
  )

freq_signif_filtered <- freq_signif_filtered %>%
  select(MRK, everything())

# Extract marker IDs
sig_markers <- freq_signif_filtered$MRK

#alt_sig <- freq_signif_filtered %>%
#  filter(MRK %in% sig_markers)

alt_sig_diff <- freq_signif_filtered %>%
  mutate(
    diff_2009 = as.numeric(EA_2009_T4.FREQ) - as.numeric(EA_2009_T2.FREQ),
    diff_2011 = as.numeric(EA_2011_T2.FREQ) - as.numeric(EA_2011_T1.FREQ),
    diff_2015 = as.numeric(EA_2015_T4.FREQ) - as.numeric(EA_2015_T1.FREQ),
    diff_2022 = as.numeric(EA_2022_T4.FREQ) - as.numeric(EA_2022_T1.FREQ)
  )

#all significant snps trajectory

all_sig <- alt_sig_diff

colnames(all_sig) <- c("MRK", "2009.start","2011.start","2015.start","2022.start",
                                    "2009.end","2011.end","2015.end","2022.end",
                                    "diff_2009","diff_2011","diff_2015","diff_2022")[1:ncol(all_sig)]

norm_all_sig <- all_sig %>%
  mutate(
    # Start = 0
    `2009.start` = 0, `2011.start` = 2011.start - 2009.start , `2015.start` = 2015.start - 2009.start, `2022.start` = 2022.start - 2009.start,
    # End = difference (already diff_* columns)
    `2009.end` = diff_2009,
    `2011.end` = diff_2011,
    `2015.end` = diff_2015,
    `2022.end` = diff_2022
  ) %>%
  # Keep only columns needed for plotting
  select(MRK,
         `2009.start`, `2009.end`,
         `2011.start`, `2011.end`,
         `2015.start`, `2015.end`,
         `2022.start`, `2022.end`)

library(tidyr)
norm_long_all <- norm_all_sig %>%
  pivot_longer(
    cols = -MRK,
    names_to = "timepoint",
    values_to = "af"
  ) %>%
  mutate(
    # Factor levels ensure correct plotting order
    timepoint = factor(timepoint,
                       levels = c("2009.start","2009.end",
                                  "2011.start","2011.end",
                                  "2015.start","2015.end",
                                  "2022.start","2022.end"))
  )

# Check for NA
any(is.na(norm_long_all$timepoint))  # should be FALSE

# -------------------------
# 3️⃣ Plot positive trajectories
# -------------------------
library(ggplot2)
ggplot(norm_long_all, aes(x = timepoint, y = af, group = MRK, color = as.factor(MRK))) +
  geom_line(alpha = 0.3) +
#  geom_point(size = 1, position = position_jitter(width = 0.1, height = 0)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "Positive SNP Trajectories (Normalized, Start = 0)",
    x = "Timepoint",
    y = "AF Change"
  ) +
  guides(color = "none")

#look at fluctuating ones...

positive_snps <- norm_all_sig[is_positive, ]
library(tidyr)
library(ggplot2)

norm_long_positive <- positive_snps %>%
  pivot_longer(
    cols = -MRK,
    names_to = "timepoint",
    values_to = "af"
  ) %>%
  mutate(
    timepoint = factor(timepoint,
                       levels = c("2009.start","2009.end",
                                  "2011.start","2011.end",
                                  "2015.start","2015.end",
                                  "2022.start","2022.end"))
  )

#no positive ones!

# Subset negative
negative_snps <- norm_all_sig[is_negative, ]
#117 snps

norm_long_negative <- negative_snps %>%
  pivot_longer(
    cols = -MRK,
    names_to = "timepoint",
    values_to = "af"
  ) %>%
  mutate(
    timepoint = factor(timepoint,
                       levels = c("2009.start","2009.end",
                                  "2011.start","2011.end",
                                  "2015.start","2015.end",
                                  "2022.start","2022.end"))
  )

ggplot(norm_long_negative, aes(x = timepoint, y = af, group = MRK, color = as.factor(MRK))) +
  geom_line(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Negative SNP Trajectories", x = "Timepoint", y = "AF Change") +
  guides(color = "none")

# okay, now lets look at the p-values of those snps

# Assuming 'significant_df' has columns: MRK, median_pvalue

pvals_positive <- sig_combinedq %>%
  filter(marker %in% positive_snps$MRK) %>%
  select(marker, median_q)

pvals_negative <- sig_combinedq %>%
  filter(marker %in% negative_snps$MRK) %>%
  select(marker, median_q)

# Quick summary of p-values
summary(pvals_positive$median_q)
summary(pvals_negative$median_q)

pvals_combined <- bind_rows(
  pvals_positive %>% mutate(type = "positive"),
  pvals_negative %>% mutate(type = "negative")
)


#median
median_all <- median(pvals_combined$median_q, na.rm = TRUE)
median_all

#0.0006504364

# q < 0.05
q_005 <- pvals_combined %>% filter(median_q < 0.05)

# q < 0.01
q_001 <- pvals_combined %>% filter(median_q < 0.01)

# q < 0.001
q_0001 <- pvals_combined %>% filter(median_q < 0.001)

# Set your output folder
outdir <- "/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/baypass-new-genome/c2-original-output-files/"

# Save as CSV
write.csv(q_005, file = file.path(outdir, "negative_snps_q005.csv"), row.names = FALSE)
write.csv(q_001, file = file.path(outdir, "negative_snps_q001.csv"), row.names = FALSE)
write.csv(q_0001, file = file.path(outdir, "negative_snps_q0001.csv"), row.names = FALSE)

ggplot(pvals_combined, aes(x = median_q, fill = type)) +
  geom_histogram(alpha = 0.6, position = "identity", bins = 30) +
  theme_minimal() +
  labs(
    title = "Distribution of Median Q values",
    x = "Median Q",
    y = "Count"
  ) +
  scale_fill_manual(values = c("positive" = "green", "negative" = "red")) +
  theme(legend.title = element_blank())

#filter for top 0.01 p values
q_001 <- pvals_combined %>%
  filter(median_q < 0.001)

# Quick check
nrow(q_001)
table(q_001$type)
#90 snps

q_0001 <- pvals_combined %>%
  filter(median_q < 0.001)

# Quick check
nrow(q_0001)
table(q_0001$type)
#64 snps

# For q < 0.01
q_001_changes <- q_001 %>%
  left_join(alt_sig_diff %>% select(MRK, diff_2009, diff_2011, diff_2015, diff_2022),
            by = c("marker" = "MRK")) %>%
  rowwise() %>%
  mutate(mean_change = mean(c(diff_2009, diff_2011, diff_2015, diff_2022), na.rm = TRUE)) %>%
  ungroup()

# For q < 0.001
q_0001_changes <- q_0001 %>%
  left_join(alt_sig_diff %>% select(MRK, diff_2009, diff_2011, diff_2015, diff_2022),
            by = c("marker" = "MRK")) %>%
  rowwise() %>%
  mutate(mean_change = mean(c(diff_2009, diff_2011, diff_2015, diff_2022), na.rm = TRUE)) %>%
  ungroup()

# q < 0.01
ggplot(q_001_changes, aes(x = mean_change)) +
  geom_histogram(binwidth = 0.05, fill = "red", color = "black") +
  theme_minimal() +
  labs(title = "Histogram of Mean AF Changes (q < 0.01)",
       x = "Mean AF Change per SNP", y = "Count")

# q < 0.001
ggplot(q_0001_changes, aes(x = mean_change)) +
  geom_histogram(binwidth = 0.05, fill = "blue", color = "black") +
  theme_minimal() +
  labs(title = "Histogram of Mean AF Changes (q < 0.001)",
       x = "Mean AF Change per SNP", y = "Count")

#for all snps

negative_snps <- negative_snps %>%
  mutate(
    diff_2009 = `2009.end` - `2009.start`,
    diff_2011 = `2011.end` - `2011.start`,
    diff_2015 = `2015.end` - `2015.start`,
    diff_2022 = `2022.end` - `2022.start`
  )

negative_snps_changes <- negative_snps %>%
  rowwise() %>%
  mutate(mean_change = mean(c(diff_2009, diff_2011, diff_2015, diff_2022), na.rm = TRUE)) %>%
  ungroup()

ggplot(negative_snps_changes, aes(x = mean_change)) +
  geom_histogram(binwidth = 0.05, fill = "blue", color = "black") +
  theme_minimal() +
  labs(title = "Histogram of Mean AF Changes (All Negative SNPs)",
       x = "Mean AF Change per SNP", y = "Count")

#same for median

# For all negative SNPs
negative_snps <- negative_snps %>%
  mutate(
    diff_2009 = `2009.end` - `2009.start`,
    diff_2011 = `2011.end` - `2011.start`,
    diff_2015 = `2015.end` - `2015.start`,
    diff_2022 = `2022.end` - `2022.start`
  )

# For q < 0.01
q_001 <- q_001 %>%
  left_join(negative_snps %>% select(MRK, diff_2009, diff_2011, diff_2015, diff_2022),
            by = c("marker" = "MRK"))

# For q < 0.001
q_0001 <- q_0001 %>%
  left_join(negative_snps %>% select(MRK, diff_2009, diff_2011, diff_2015, diff_2022),
            by = c("marker" = "MRK"))
# All negative SNPs
negative_snps_median <- negative_snps %>%
  rowwise() %>%
  mutate(median_change = median(c(diff_2009, diff_2011, diff_2015, diff_2022), na.rm = TRUE)) %>%
  ungroup()

# q < 0.01
q_001_median <- q_001 %>%
  rowwise() %>%
  mutate(median_change = median(c(diff_2009, diff_2011, diff_2015, diff_2022), na.rm = TRUE)) %>%
  ungroup()

# q < 0.001
q_0001_median <- q_0001 %>%
  rowwise() %>%
  mutate(median_change = median(c(diff_2009, diff_2011, diff_2015, diff_2022), na.rm = TRUE)) %>%
  ungroup()

# All negative SNPs
ggplot(negative_snps_median, aes(x = median_change)) +
  geom_histogram(binwidth = 0.05, fill = "blue", color = "black") +
  theme_minimal() +
  labs(title = "Median AF Changes (All Negative SNPs)", x = "Median AF Change per SNP", y = "Count")

ggplot(q_001_median, aes(x = median_change)) +
  geom_histogram(binwidth = 0.05, fill = "red", color = "black") +
  theme_minimal() +
  labs(title = "Median AF Changes (q < 0.01)", x = "Median AF Change per SNP", y = "Count")
ggplot(q_0001_median, aes(x = median_change)) +
  geom_histogram(binwidth = 0.05, fill = "green", color = "black") +
  theme_minimal() +
  labs(title = "Median AF Changes (q < 0.001)", x = "Median AF Change per SNP", y = "Count")

  ```

#1 q-value correction



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


ASK REID: should we consider year as a factor or not in GLM


convert vcf from indivudal to pooled.
calculate pop gen stats from it (check what Fst) calculate AF from each site, use line 117 from reids scripts AF trajectory.


























```bash
#!/bin/bash
#SBATCH -D /gxfs_home/geomar/smomw573/work/seasonal_adaptation/log_files/
#SBATCH --mail-type=END
#SBATCH --mail-user=jnascimento@geomar.de
#SBATCH --partition=base
#SBATCH --nodes=3
#SBATCH --tasks-per-node=3
#SBATCH --cpus-per-task=1
#SBATCH --mem=50G
#SBATCH --time=05:00:00
#SBATCH --array=1-$(wc -l < vcf_list.txt)    # Auto set to number of files
#SBATCH --job-name=vcftools_div
#SBATCH --output=vcftools_div_%A_%a.out
#SBATCH --error=vcftools_div_%A_%a.err


# Activate your environment with vcftools
source ~/miniconda3/bin/activate parallel

INPUT=/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/slim-simulations-new-genome/vcfs
VCF_LIST=/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/slim-simulations-new-genome/vcfs/vcf_list.txt
OUTDIR=/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/slim-simulations-new-genome/pi_results

# Get the VCF file for this array task
VCF=$(sed -n "${SLURM_ARRAY_TASK_ID}p" vcf_list.txt)

# Get filename without extension
BASENAME=$(basename "$VCF" .vcf)


# Run VCFtools for windowed nucleotide diversity
vcftools --vcf "$INPUT/$VCF" \
  --window-pi 10000 \
  --out "$OUTDIR/${BASENAME}_pi"
  ```
vcftools --vcf msprime_burnin_n200000_44006078_20.vcf \
  --window-pi 10000 \
  --out /gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/slim-simulations-new-genome/pi_results
15:18 msprime_burnin_n200000_44006078_20.vcf
