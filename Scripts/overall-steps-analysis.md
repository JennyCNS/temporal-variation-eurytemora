# create local path to R library

```bash
module load R/4.1.1
module load gcc/10.2.0
cd ~
mkdir R_libs
export R_LIBS=$HOME/R_libs:$R_LIBS
```

# check quality of trimmed reads
```
#!/bin/bash
#SBATCH -D /gxfs_work1/geomar/smomw573/seasonal_adaptation/raw_data/trimmed/
#SBATCH --mail-type=END
#SBATCH --mail-user=jnascimento@geomar.de
#SBATCH --partition=cluster
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=47:00:00
#SBATCH --job-name=trim
#SBATCH --output=seqkit.out
#SBATCH --error=seqkit.err


seqkit stats /gxfs_work1/geomar/smomw573/seasonal_adaptation/raw_data/trimmed/*.gz
```
# call SNPs using PoolSNP

```bash
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
output=/gxfs_work1/geomar/smomw573/seasonal_adaptation/analysis/PoolSNP/final-maf0.01-mincov50-mincount10-allsites0-allsamples_b.vcf  \
reference=/gxfs_work1/geomar/smomw573/seasonal_adaptation/genome/Eaffinis.Baltic.PseudoRef.Mar22.fasta \
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

#trial 1 wrong sample name input
#trial 2 did not have scripts directory in the Poolsnp directory
#trial 3 ran together with nohup so missed header
#trial 4 fixed all of the above but ran with all sites = 1
#trial 5 = trial 4 but no sites = 0
#trial 6 = trial 5 but different working direcotry (softare) and all sites = 1
#trial 7 = min-count = 10

```
#I tried this with multiple parameters but the final one were these
min-cov=50 \
max-cov=0.95 \
min-count=10 \
min-freq=0.01 \
miss-frac=0.1 \
base-quality 15 \
jobs=17 \
badsites=1 \
allsites=0

# Now we checked the quality of the vcf file using all samples 

# explore data statistics with plink
### this is wrong so please IGNORE THIS OUPUT
```bash script 
conda install -c bioconda plink

# Output directory where processed files will be stored
output=/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/PoolSNP/stats

plink --file final-maf0.01-mincov50-mincount10-allsites0-allsamples_b.vcf --allow-extra-chr --missing --hardy --freq --out /$output/final-allsamples

```
# plot ind missing data in R

```R
###Poolseq data 
###plink output data statistics
##July 2023


library(ggplot2)
library(tidyverse)
library(data.table)


###################
#all samples
#MAF 0.01
#MIN COV 50
#MIN COUNT 10
###################

#MAF
frq_data <- fread("final-allsamples.frq")

#in ascending order
ggplot(frq_data, aes(x = MAF)) +
  geom_density(color="dodgerblue", fill="lightblue")+
  theme_bw()
ggsave("allsamplesMAF.pdf", width = 7, height = 7, units =c("cm"), dpi =3000)


#simple plot
h = hist(frq_data$MAF)
h$density = h$counts/sum(h$counts)
plot(h, freq = FALSE)

### lmiss
lmiss <- fread("final-allsamples.lmiss")
a <- ggplot(lmiss, aes(F_MISS)) + geom_density(colour = "dodgerblue1", fill = "lightblue", alpha = 0.3)
a + theme_light()+
  labs(x="frequency LMiss",
       y="density") +
  theme_classic() +
  theme(axis.text.x=element_text(hjust = 1)) +
  theme(text = element_text(size = 12)) 


h = hist(lmiss$F_MISS)
h$density = h$counts/sum(h$counts)
plot(h, freq = FALSE)

#individual missing

imiss <- fread("final-allsamples.imiss")
str(imiss)
a <- ggplot(imiss, aes(F_MISS)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()

ggplot(imiss, aes(x=FID, y=F_MISS)) + geom_bar(stat="identity", fill = "dodgerblue1", colour = "black", alpha = 0.3) +
  labs(x="Sample",
       y="Frequency missing data") +
  theme_classic() +
  theme(axis.text.x=element_text(angle = 45, hjust = 1)) +
  theme(text = element_text(size = 8)) 
#density

ggplot(imiss, aes(F_MISS)) + geom_density( fill = "dodgerblue1", colour = "black", alpha = 0.3) +
  labs(x="missing data (individual)",
       y="Density") +
  theme_classic() +
  theme(axis.text.x=element_text(angle = 45, hjust = 1)) +
  theme(text = element_text(size = 8)) 

ggsave("freq_missing_data_allvariants_maf0.01-mincount20-mincov50.pdf", width = 14, height = 6, units =c("cm"), dpi =3000)


#Individual missing

h = hist(lmiss$F_MISS)
h$density = h$counts/sum(h$counts)
plot(h, freq = FALSE)

###################
#excluding bad samples
#MAF 0.01
#MIN COV 50
#MIN COUNT 10
###################

#MAF
frq_data <- fread("10mincount_filtered.frq")

#in ascending order
ggplot(frq_data, aes(x = MAF)) +
  geom_density(color="dodgerblue", fill="lightblue")+
  theme_bw()
ggsave("allsamplesMAF.pdf", width = 7, height = 7, units =c("cm"), dpi =3000)


#simple plot
h = hist(frq_data$MAF)
h$density = h$counts/sum(h$counts)
plot(h, freq = FALSE)

### lmiss
lmiss <- fread("10mincount_filtered.lmiss")
a <- ggplot(lmiss, aes(F_MISS)) + geom_density(colour = "dodgerblue1", fill = "lightblue", alpha = 0.3)
a + theme_light()+
  labs(x="frequency LMiss",
       y="density") +
  theme_classic() +
  theme(axis.text.x=element_text(hjust = 1)) +
  theme(text = element_text(size = 12)) 


h = hist(lmiss$F_MISS)
h$density = h$counts/sum(h$counts)
plot(h, freq = FALSE)

#individual missing

imiss <- fread("10mincount_filtered.imiss")
str(imiss)
a <- ggplot(imiss, aes(F_MISS)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()

ggplot(imiss, aes(x=FID, y=F_MISS)) + geom_bar(stat="identity", fill = "dodgerblue1", colour = "black", alpha = 0.3) +
  labs(x="Sample",
       y="Frequency missing data") +
  theme_classic() +
  theme(axis.text.x=element_text(angle = 45, hjust = 1)) +
  theme(text = element_text(size = 8)) 
#density

ggplot(imiss, aes(F_MISS)) + geom_density( fill = "dodgerblue1", colour = "black", alpha = 0.3) +
  labs(x="missing data (individual)",
       y="Density") +
  theme_classic() +
  theme(axis.text.x=element_text(angle = 45, hjust = 1)) +
  theme(text = element_text(size = 8)) 

ggsave("freq_missing_data_allvariants_maf0.01-mincount20-mincov50.pdf", width = 14, height = 6, units =c("cm"), dpi =3000)


#Individual missing

h = hist(lmiss$F_MISS)
h$density = h$counts/sum(h$counts)
plot(h, freq = FALSE)


###################
#excluding bad samples
#MAF 0.01
#MIN COV 50
#MIN COUNT 5
###################

#MAF
frq_data <- fread("10mincount_filtered.frq")

#in ascending order
ggplot(frq_data, aes(x = MAF)) +
  geom_density(color="dodgerblue", fill="lightblue")+
  theme_bw()
ggsave("allsamplesMAF.pdf", width = 7, height = 7, units =c("cm"), dpi =3000)


#simple plot
h = hist(frq_data$MAF)
h$density = h$counts/sum(h$counts)
plot(h, freq = FALSE)

### lmiss
lmiss <- fread("5mincount_filtered.lmiss")
a <- ggplot(lmiss, aes(F_MISS)) + geom_density(colour = "dodgerblue1", fill = "lightblue", alpha = 0.3)
a + theme_light()+
  labs(x="frequency LMiss",
       y="density") +
  theme_classic() +
  theme(axis.text.x=element_text(hjust = 1)) +
  theme(text = element_text(size = 12)) 


h = hist(lmiss$F_MISS)
h$density = h$counts/sum(h$counts)
plot(h, freq = FALSE)

#individual missing

imiss <- fread("5mincount_filtered.imiss")
str(imiss)
a <- ggplot(imiss, aes(F_MISS)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()

ggplot(imiss, aes(x=FID, y=F_MISS)) + geom_bar(stat="identity", fill = "dodgerblue1", colour = "black", alpha = 0.3) +
  labs(x="Sample",
       y="Frequency missing data") +
  theme_classic() +
  theme(axis.text.x=element_text(angle = 45, hjust = 1)) +
  theme(text = element_text(size = 8)) 
#density

ggplot(imiss, aes(F_MISS)) + geom_density( fill = "dodgerblue1", colour = "black", alpha = 0.3) +
  labs(x="missing data (individual)",
       y="Density") +
  theme_classic() +
  theme(axis.text.x=element_text(angle = 45, hjust = 1)) +
  theme(text = element_text(size = 8)) 

ggsave("freq_missing_data_allvariants_maf0.01-mincount20-mincov50.pdf", width = 14, height = 6, units =c("cm"), dpi =3000)


#Individual missing

h = hist(lmiss$F_MISS)
h$density = h$counts/sum(h$counts)
plot(h, freq = FALSE)

```

# From the output we decided to remove samples 2007 T1 and T2 as well as 2011 T3.

# make new mpileup file

```bash
#!/bin/bash
#SBATCH -D .
#SBATCH --mail-type=END
#SBATCH --mail-user=jnascimento@geomar.de
#SBATCH --partition=cluster
#SBATCH --cpus-per-task=1
#SBATCH --mem=60G
#SBATCH --time=18:00:00
#SBATCH --job-name=to_mpileup
#SBATCH --output=tompileup.out
#SBATCH --error=tompileup.err

module load samtools/1.10

cd $WORK/seasonal_adaptation/analysis/variants

samtools mpileup -q 15 -Q 0 -d 8000 -R -A -B \
        -f $WORK/seasonal_adaptation/genome/Eaffinis.Baltic.PseudoRef.Mar22.fasta \
        -b $WORK/seasonal_adaptation/analysis/variants/bamfiles-excluding2007-2011.txt \
        -o all.mpileup-excluding20072011
```
```bash

# rerun PoolSNP with same parameters but trying min-count10 and 5
# kept the min count 10
bash /gxfs_home/geomar/smomw573/software/PoolSNP/PoolSNP.sh \
mpileup=/gxfs_work1/geomar/smomw573/seasonal_adaptation/analysis/PoolSNP/all.mpileup-excluding20072011  \
output=/gxfs_work1/geomar/smomw573/seasonal_adaptation/analysis/PoolSNP/out-maf0.01-mincov50-mincount10-allsites-excluding20072011 \
reference=/gxfs_work1/geomar/smomw573/seasonal_adaptation/genome/Eaffinis.Baltic.PseudoRef.Mar22.fasta \
names=EA_2009_T1,EA_2009_T2,EA_2009_T3,EA_2009_T4,EA_2011_T1,EA_2011_T2,EA_2015_T1,EA_2015_T2,EA_2015_T3,EA_2015_T4,EA_2022_T1,EA_2022_T2,EA_2022_T3,EA_2022_T4 \
min-cov=50 \
max-cov=0.95 \
min-count=10 \
min-freq=0.01 \
miss-frac=0.1 \
base-quality 15 \
jobs=17 \
badsites=1 \
allsites=1
```
# Check number of SNPs in each vcf output file

```bash
cd /gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/PoolSNP

#count snps
#all samples included, maf 0.01, min cov 50, min count 10
zgrep -v "^##" final-maf0.01-mincov50-mincount10-allsites0-allsamples_b.vcf.gz | wc -l
# SNPs 21779

#removing the 3 bad samples, maf 0.01, min cov 50, min count 10
zgrep -v "^#" final-maf0.01-mincov50-mincount10-allsites0-excluding20072011.vcf.gz | wc -l
#SNPs 2424258

#removing the 3 bad samples, maf 0.01, min cov 50, min count 5
zgrep -v "^#" final-maf0.01-mincov50-mincount5-allsites0-excluding20072011.vcf.gz | wc -l
#SNPs 2424280

```
# ran script fix_vcf.py to fix it for grenedalf

# keep biallelic snps

```bash
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
VCF=/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/PoolSNP/finalfile.vcf.gz
OUTPUT=/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/PoolSNP


vcftools --gzvcf $VCF --min-alleles 2 --max-alleles 2 --recode --out $OUT/finalfile_output_biallelic

```

# fix file so it can be read by grenedalf using python script

```bash
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

#add in final vcf header 
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
# file needs to be indexed for rest of the analysis
```
# tabix file

```bash
bgzip finalfile.vcf
tabix -p vcf finalfile.vcf

```

# calculate missing data per indivudal sample

```bash
module load vcftools/0.1.14
vcftools --vcf final-maf0.01-mincov50-mincount10-allsites0-allsamples_b.vcf.gz --missing-indv
 
#plot outputs in R
#R scripts are saved in computer
#check snps in the output files

module load bcftools/1.10.2
bcftools view -v snps input.vcf.gz -o output_snps.vcf.gzls
```

# filtering missing data 

```bash 
module load vcftools/0.1.14

vcftools --gzvcf finalfile.vcf.gz  --max-missing 1.0  --recode --out finalfile.nomissingdata

#output genotype depths
vcftools --vcf finalfile.nomissingdata.recode.vcf --site-mean-depth --out finalfile.nomissingdata

module load R/4.1.1
```
# Fix vcf for coverage depth
# select SNPs which only have 3x average coverage

```r
install.packages("data.table")
install.packages("poolfstat")
library(data.table)
library(poolfstat)
#trial plot
plot(1,1)

#load data 
df <- read.table("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/grenedalf-freq/vcf-maf0.01-misfrac0.1-maxcov95/finalfile.nomissingdata.ldepth.mean", header = T)

head(df)
hist(df$MEAN_DEPTH)

#check mean and median

mean <- mean(df$MEAN_DEPTH)
497.8522
mean3 = mean*3
1493.557

sum(df$MEAN_DEPTH > mean3)
81366
sum(df$MEAN_DEPTH > mean3)/nrow(df)
0.05871821

#poolfstat
#did not do this step with individual datanames
inNames <- read.table("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/grenedalf-freq/vcf-maf0.01-misfrac0.1-maxcov95/finalfile.nomissingdata.ldepth.mean", header=T, sep="\t", nrow=1)
inNames <- colnames(inNames)[grep("FREQ",colnames(inNames))]
inNames <- inNames[1:length(inNames)-1]
################################

#remove all markers with a depth over than three times the mean
dat <- vcf2pooldata(vcf.file="/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/grenedalf-freq/vcf-maf0.01-misfrac0.1-maxcov95/finalfile.nomissingdata.recode.vcf",poolsizes=c(rep(100,14)),
			min.cov.per.pool = 50,max.cov.per.pool=1493.5570,min.maf=0.01,nlines.per.readblock=1000000)

pooldata2genobaypass(dat,writing.dir="/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/grenedalf-freq/vcf-maf0.01-misfrac0.1-maxcov95/")
#ouput is snpdet
```

# correct for depth 
# first make file with chrom + position of SNPs that have the proper depth 
```bash
awk '{print $1"\t"$2}' snpdet > snp-positions.txt
```
# filter this list with vcf

```bash
srun --pty --x11 --nodes=1 --cpus-per-task=4 --mem=16000 --time=02:00:00 /bin/bash

vcftools --vcf finalfile.nomissingdata.recode.vcf --positions snp-positions.txt --recode --out depth-corrected

#After filtering, kept 1258282 out of a possible 1385703 Sites


#removing MAF errors
bgzip depth-corrected.recode.vcf
tabix -p vcf depth-corrected.recode.vcf.gz


#run grenedalf
module load htslib/1.10.2  bzip2/1.0.8  

GRENEDALF=/gxfs_home/geomar/smomw573/software/grenedalf/bin/grenedalf
VCF=/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/grenedalf-freq/vcf-maf0.01-misfrac0.1-maxcov95/depth-corrected.recode.vcf.gz
GENOME=/gxfs_work1/geomar/smomw573/seasonal_adaptation/genome/Eaffinis.Baltic.PseudoRef.Mar22.fasta \
OUTPUT=/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/grenedalf-freq/vcf-maf0.01-misfrac0.1-maxcov95/

$GRENEDALF frequency --vcf-path $VCF  --reference-genome-fasta-file $GENOME --write-total-frequency --allow-file-overwriting --file-suffix depth-frequencies > $OUTPUT/grenedalftrial.log



```


# open file in R
```r
freq <- read.csv("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/grenedalf-freq/vcf-maf0.01-misfrac0.1-maxcov95/frequencydepth-frequencies.csv", header=T)

freq$MAF <- ifelse(freq$TOTAL.FREQ > 0.5, 1-freq$TOTAL.FREQ, freq$TOTAL.FREQ)

hist(freq$TOTAL.FREQ,breaks=50)
hist(freq$MAF,breaks=50)

sum(freq$TOTAL.FREQ < 0.01)
sum(freq$TOTAL.FREQ > 0.99)
sum(freq$MAF < 0.01)
sum(freq$MAF > 0.99)

min(freq$MAF)
max(freq$TOTAL.FREQ)
sum(freq$TOTAL.FREQ < 0.01)
#0
sum(freq$TOTAL.FREQ > 0.99)
#0
sum(freq$MAF < 0.01)
#0
sum(freq$MAF > 0.99)
#0
min(freq$MAF)
#0.01
max(freq$MAF)
#0.5
max(freq$TOTAL.FREQ)
#0.99
sum(freq$MAF < 0.5)
#1258257
sum(freq$MAF < 0.05)
#848949

```
# provide Reid access to my folders
```bash
setfacl -m u:smomw504:rx /gxfs_work1/geomar/smomw573
```

# generate individual frequencies grenedalf

```bash
module load htslib/1.10.2  bzip2/1.0.8  

GRENEDALF=/gxfs_home/geomar/smomw573/software/grenedalf/bin/grenedalf
VCF=/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/grenedalf-freq/vcf-maf0.01-misfrac0.1-maxcov95/depth-corrected.recode.vcf.gz
GENOME=/gxfs_work1/geomar/smomw573/seasonal_adaptation/genome/Eaffinis.Baltic.PseudoRef.Mar22.fasta \
OUTPUT=/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/grenedalf-freq/vcf-maf0.01-misfrac0.1-maxcov95/

$GRENEDALF frequency --vcf-path $VCF  --reference-genome-fasta-file $GENOME --write-sample-alt-freq --file-suffix sample-frequencies > $OUTPUT/grenedalftrial.log
```

# plot pca in R

```bash 
module load R/4.1.1
module load gcc/10.2.0
```

```R
library(data.table)
library(ggplot2)
df<- read.csv("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/grenedalf-freq/vcf-maf0.01-misfrac0.1-maxcov95/frequencysample-frequencies.csv", header = T)
head(df)
#df <- noquote(t(frequency_ind_samples[,-(2:4)]))

#colnames(df) <- df[1, ]

# Remove the first row as it's now used for column names
df <- df[,-(2:4)]
#test data
#df2 <- df[0:15,1:14]
#df2
#all okay

pops <- c("EA_2009_T1.FREQ", "EA_2009_T2.FREQ", "EA_2009_T3.FREQ",
          "EA_2009_T4.FREQ", "EA_2011_T1.FREQ", "EA_2011_T2.FREQ", 
          "EA_2015_T1.FREQ", "EA_2015_T2.FREQ", "EA_2015_T3.FREQ", 
          "EA_2015_T4.FREQ", "EA_2022_T1.FREQ", "EA_2022_T2.FREQ", 
          "EA_2022_T3.FREQ", "EA_2022_T4.FREQ") 
varout <- apply(df[, 2:ncol(df)], 1, var)

freqs <- t(df[varout != 0,2:ncol(df)])


nrow(freqs)
ncol(freqs)
#14 

##
## plot pca
##

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
          axis.title.x = element_text(size = 28),
          axis.text.x= element_text(size=26),
          axis.text.y= element_text(size=26),   # Adjust x-axis label size
          axis.title.y = element_text(size = 28))+
        #theme(legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid'))+
#        ggtitle("F1")+
guides(fill=guide_legend(override.aes=list(
        shape=c(21,21,22, 23),
        #fill=c('steelblue1','steelblue','grey45', "darkorchid2", "firebrick3" )),order = 2),
    fill=c("#D3DDDC",'#6699CC',"#F2AD00","#00A08A")),order = 2))

ggsave("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/figures/pca_afs_allsamples.pdf",d, w=6, h=4)


#exclude samples 2009 T1 and 2015 T3

df<- read.csv("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/grenedalf-freq/vcf-maf0.01-misfrac0.1-maxcov95/frequencysample-frequencies.csv", header = T)
df <- df[,-c(2:5, 13)]



pops <- c("EA_2009_T2.FREQ", "EA_2009_T3.FREQ",
          "EA_2009_T4.FREQ", "EA_2011_T1.FREQ", 
          "EA_2011_T2.FREQ", "EA_2015_T1.FREQ", 
          "EA_2015_T2.FREQ", "EA_2015_T4.FREQ", 
          "EA_2022_T1.FREQ", "EA_2022_T2.FREQ", 
          "EA_2022_T3.FREQ", "EA_2022_T4.FREQ") 
varout <- apply(df[, 2:ncol(df)], 1, var)
freqs <- t(df[varout != 0,2:ncol(df)])


nrow(freqs)
ncol(freqs)


##
## plot pca
##

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
        geom_point(size=16) +
        xlab(paste0("PC1: ",percentVar[1],"% variance")) +
        ylab(paste0("PC2: ",percentVar[2],"% variance")) +
        theme_bw() +
       # ylim(-30, 23) + xlim(-50, 65)+
        scale_shape_manual(values=c( 21,22,23,24))+
        scale_color_manual(values=c('black')) +
        #scale_fill_manual(values=c('steelblue1','steelblue','grey45', "darkorchid2", "firebrick3"),
  scale_fill_manual(values=c("#D3DDDC",'#6699CC',"#F2AD00","#00A08A"),
                    labels = c("2009","2011", "2015","2022"))+
        #theme(legend.position = c(0.83,0.85),
        #    legend.background = element_blank(),
        #legend.box.background = element_rect(colour = "black"),
         theme(legend.title = element_blank()) +
         theme(legend.text=element_text(size=18),
          axis.title.x = element_text(size = 28),
          axis.text.x= element_text(size=28),
          axis.text.y= element_text(size=28),   # Adjust x-axis label size
          axis.title.y = element_text(size = 28))+
       # theme(legend.text=element_text(size=8))+
        #theme(legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid'))+
#        ggtitle("F1")+
guides(fill=guide_legend(override.aes=list(
        shape=c(21,21, 21,21),
        #fill=c('steelblue1','steelblue','grey45', "darkorchid2", "firebrick3" )),order = 2),
    fill=c("#D3DDDC",'#6699CC',"#F2AD00","#00A08A")),order = 2))

ggsave("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/figures/pca_afs_excluding20092015.pdf",d, w=10, h=7.8)
ggsave("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/figures/pca_afs_excluding20092015.jpg",d, dpi=300, w=10, h=7.8)

######################################
#repeat analysis with non-scaled data#
######################################

df1<- read.csv("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/grenedalf-freq/vcf-maf0.01-misfrac0.1-maxcov95/frequencysample-frequencies.csv", header = T)

#df <- noquote(t(frequency_ind_samples[,-(2:4)]))

#colnames(df) <- df[1, ]

# Remove the first row as it's now used for column names
df <- df1[,-(2:3)]
df <- df[,-(2)]
#test data
#df2
#all okay

pops <- c("EA_2009_T1.FREQ", "EA_2009_T2.FREQ", "EA_2009_T3.FREQ",
          "EA_2009_T4.FREQ", "EA_2011_T1.FREQ", "EA_2011_T2.FREQ", 
          "EA_2015_T1.FREQ", "EA_2015_T2.FREQ", "EA_2015_T3.FREQ", 
          "EA_2015_T4.FREQ", "EA_2022_T1.FREQ", "EA_2022_T2.FREQ", 
          "EA_2022_T3.FREQ", "EA_2022_T4.FREQ") 
varout <- apply(df[, 2:ncol(df)], 1, var)

freqs <- t(df[varout != 0,2:ncol(df)])


nrow(freqs)
ncol(freqs)
#14 

##
## plot pca
##

pcaResult <- prcomp(freqs, scale=F)
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
        scale_shape_manual(values=c( 21,21,22,23))+
        scale_color_manual(values=c('black')) +
        #scale_fill_manual(values=c('steelblue1','steelblue','grey45', "darkorchid2", "firebrick3"),
  scale_fill_manual(values=c("#D3DDDC",'#6699CC',"#F2AD00","#00A08A"),
                    labels = c("2009","2011", "2015",
                               "2022"))+
        #theme(legend.position = c(0.83,0.85),
        #    legend.background = element_blank(),
        #legend.box.background = element_rect(colour = "black"),
        theme(legend.title = element_blank()) +
       # theme(legend.text=element_text(size=8))+
        #theme(legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid'))+
#        ggtitle("F1")+
guides(fill=guide_legend(override.aes=list(
        shape=c(21,21,22, 23),
        #fill=c('steelblue1','steelblue','grey45', "darkorchid2", "firebrick3" )),order = 2),
    fill=c("#D3DDDC",'#6699CC',"#F2AD00","#00A08A")),order = 2))

ggsave("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/figures/pca_afs-nonscaled_allsamples.pdf",d, w=5.5, h=3.7)


#exclude samples 2009 T1 and 2015 T3

df<- read.csv("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/grenedalf-freq/vcf-maf0.01-misfrac0.1-maxcov95/frequencysample-frequencies.csv", header = T)
df <- df[,-c(2:5, 13)]

pops <- c("EA_2009_T2.FREQ", "EA_2009_T3.FREQ",
          "EA_2009_T4.FREQ", "EA_2011_T1.FREQ", 
          "EA_2011_T2.FREQ", "EA_2015_T1.FREQ", 
          "EA_2015_T2.FREQ", "EA_2015_T4.FREQ", 
          "EA_2022_T1.FREQ", "EA_2022_T2.FREQ", 
          "EA_2022_T3.FREQ", "EA_2022_T4.FREQ") 
varout <- apply(df[, 2:ncol(df)], 1, var)
freqs <- t(df[varout != 0,2:ncol(df)])


nrow(freqs)
ncol(freqs)


##
## plot pca
##

pcaResult <- prcomp(freqs, scale=F)
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
                    labels = c("2009","2011", "2015","2022"))+
        #theme(legend.position = c(0.83,0.85),
        #    legend.background = element_blank(),
        #legend.box.background = element_rect(colour = "black"),
        theme(legend.title = element_blank()) +
       # theme(legend.text=element_text(size=8))+
        #theme(legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid'))+
#        ggtitle("F1")+
guides(fill=guide_legend(override.aes=list(
        shape=c(21,21, 21,21),
        #fill=c('steelblue1','steelblue','grey45', "darkorchid2", "firebrick3" )),order = 2),
    fill=c("#D3DDDC",'#6699CC',"#F2AD00","#00A08A")),order = 2))

ggsave("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/figures/pca_afs-nonscaled_excluding20092015.pdf",d, w=5.5, h=3.7)
```

# calculate Fst between different populations
 as discussed with Reid, we should see if there are any chances in Fst during the course of the year - so comparing early vs late populations could give us an idea of this.
 the comparisons will be:
 EA_2009_T1 vs EA_2009_T4
 EA_2011_T1 vs EA_2011_T3
 EA_2015_T1 vs EA_2015_T4
 EA_2022_T1 vs EA_2015_T4

# Trial Fst comparisons with grenedalf
```bash
#run grenedalf
module load htslib/1.10.2  bzip2/1.0.8  

GRENEDALF=/gxfs_home/geomar/smomw573/software/grenedalf/bin/grenedalf
VCF=/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/grenedalf-freq/vcf-maf0.01-misfrac0.1-maxcov95/depth-corrected.recode.vcf.gz
GENOME=/gxfs_work1/geomar/smomw573/seasonal_adaptation/genome/Eaffinis.Baltic.PseudoRef.Mar22.fasta \
OUTPUT=/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/grenedalf-freq/vcf-maf0.01-misfrac0.1-maxcov95/

$GRENEDALF fst --vcf-path $VCF --reference-genome-fasta-file $GENOME --method unbiased-hudson --pool-sizes 50 --window-type sliding --window-sliding-width 1000 --file-suffix all-parwise  --out-dir $OUTPUT/grenedalfst.log

$GRENEDALF fst --vcf-path $VCF --reference-genome-fasta-file $GENOME --method unbiased-hudson --pool-sizes 50 --window-type sliding --window-sliding-width 5000 --file-suffix all-parwise-5000  --out-dir $OUTPUT/grenedalfst.log

$GRENEDALF fst --vcf-path $VCF --reference-genome-fasta-file $GENOME --method unbiased-hudson --pool-sizes 50 --window-type sliding --window-sliding-width 10000 --file-suffix all-parwise-10000  --out-dir $OUTPUT/grenedalfst.log

#single snp by snp FST calculation
#genomic fst

$GRENEDALF fst --vcf-path $VCF --reference-genome-fasta-file $GENOME --method unbiased-hudson --pool-sizes 50 --window-type single --file-suffix single-snps2 --out-dir $OUTPUT/grenedalfst.log

#sort vcf
#we now have a new vcf file that has the SNPs alligned to the newer version of the genome
#!!!!!!!!!!!!!!!!!!!
############
#rerun some analyisis

##not needed cause Reid had already sorted it
module load bcftools/1.10.2

GENOME=/gxfs_work1/geomar/smomw573/seasonal_adaptation/genome/affinis.Atlantic.long_read_draft.Mar22.fasta \
VCF=/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/variants/final_20230810.vcf.gz
OUTPUT=/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/variants/final_20230810_sorted.vcf.gz


bcftools sort -Ov -o $OUTPUT -T tmp_dir -m 1000M $VCF
##not needed cause Reid had already sorted it

GRENEDALF=/gxfs_home/geomar/smomw573/software/grenedalf/bin/grenedalf
VCF=/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/variants/final_20230810_sorted.vcf.gz \
GENOME=/gxfs_work1/geomar/smomw573/seasonal_adaptation/genome/affinis.Atlantic.long_read_draft.Mar22.fasta \
OUTPUT=/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/grenedalf-freq/vcf-maf0.01-misfrac0.1-maxcov95/

$GRENEDALF fst --vcf-path $VCF --method unbiased-hudson --allow-file-overwriting --pool-sizes 100 --window-type queue --window-queue-count 1 --file-suffix single-snps-newvcf-queue --out-dir $OUTPUT/grenedalfst.log

#rerun analyisis with transformed genome

```

transform all negative values to zero

Fix csv file and plot pairwise comparisons in R.

```bash
module load R/4.1.1
R
```
```R
library(dplyr)
library(ggplot2)
library(tidyr)
library(reshape2)

fst <- read.csv("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/grenedalf-freq/vcf-maf0.01-misfrac0.1-maxcov95/grenedalfst.log/fst-all-pairwise.csv", header = TRUE, na.strings = "")
#negative values to zeros

fst <- fst %>% mutate_all(function(x) ifelse(x < 0, 0, x))

#means of each colums (are the some higher then others)

mean_fst <- colMeans(fst[, 5:ncol(fst)], na.rm = TRUE)

#scatterplot number of snps per window

mean_fst 
> means
EA_2009_T1.EA_2009_T2 EA_2009_T1.EA_2009_T3 EA_2009_T1.EA_2009_T4
         6.522382e-04          9.119229e-04          8.192032e-04
EA_2009_T1.EA_2011_T1 EA_2009_T1.EA_2011_T2 EA_2009_T1.EA_2015_T1
         9.224172e-04          1.136448e-03          1.360173e-03
EA_2009_T1.EA_2015_T2 EA_2009_T1.EA_2015_T3 EA_2009_T1.EA_2015_T4
         1.306930e-03          7.707226e-04          1.405240e-03
EA_2009_T1.EA_2022_T1 EA_2009_T1.EA_2022_T2 EA_2009_T1.EA_2022_T3
         1.432560e-03          1.388207e-03          1.539748e-03
EA_2009_T1.EA_2022_T4 EA_2009_T2.EA_2009_T3 EA_2009_T2.EA_2009_T4
         1.408387e-03          2.212941e-04          1.727231e-04
EA_2009_T2.EA_2011_T1 EA_2009_T2.EA_2011_T2 EA_2009_T2.EA_2015_T1
         3.530734e-04          3.538277e-04          5.111806e-04
EA_2009_T2.EA_2015_T2 EA_2009_T2.EA_2015_T3 EA_2009_T2.EA_2015_T4
         5.006985e-04          4.776478e-04          5.687172e-04
EA_2009_T2.EA_2022_T1 EA_2009_T2.EA_2022_T2 EA_2009_T2.EA_2022_T3
         5.954697e-04          5.563285e-04          6.646338e-04
EA_2009_T2.EA_2022_T4 EA_2009_T3.EA_2009_T4 EA_2009_T3.EA_2011_T1
         5.654983e-04          9.451288e-05          3.842776e-04
EA_2009_T3.EA_2011_T2 EA_2009_T3.EA_2015_T1 EA_2009_T3.EA_2015_T2
         2.044344e-04          2.918076e-04          3.178607e-04
EA_2009_T3.EA_2015_T3 EA_2009_T3.EA_2015_T4 EA_2009_T3.EA_2022_T1
         5.442187e-04          3.234878e-04          4.064266e-04
EA_2009_T3.EA_2022_T2 EA_2009_T3.EA_2022_T3 EA_2009_T3.EA_2022_T4
         4.129459e-04          4.696007e-04          3.566869e-04
EA_2009_T4.EA_2011_T1 EA_2009_T4.EA_2011_T2 EA_2009_T4.EA_2015_T1
         3.001030e-04          1.674710e-04          2.809781e-04
EA_2009_T4.EA_2015_T2 EA_2009_T4.EA_2015_T3 EA_2009_T4.EA_2015_T4
         2.926337e-04          4.850598e-04          3.177162e-04
EA_2009_T4.EA_2022_T1 EA_2009_T4.EA_2022_T2 EA_2009_T4.EA_2022_T3
         3.578791e-04          3.527122e-04          4.222247e-04
EA_2009_T4.EA_2022_T4 EA_2011_T1.EA_2011_T2 EA_2011_T1.EA_2015_T1
         3.444598e-04          3.314538e-04          3.071385e-04
EA_2011_T1.EA_2015_T2 EA_2011_T1.EA_2015_T3 EA_2011_T1.EA_2015_T4
         2.882376e-04          3.725812e-04          3.199376e-04
EA_2011_T1.EA_2022_T1 EA_2011_T1.EA_2022_T2 EA_2011_T1.EA_2022_T3
         3.075624e-04          2.813889e-04          3.398815e-04
EA_2011_T1.EA_2022_T4 EA_2011_T2.EA_2015_T1 EA_2011_T2.EA_2015_T2
         3.092691e-04          1.928528e-04          1.731642e-04
EA_2011_T2.EA_2015_T3 EA_2011_T2.EA_2015_T4 EA_2011_T2.EA_2022_T1
         6.036474e-04          2.024127e-04          2.407018e-04
EA_2011_T2.EA_2022_T2 EA_2011_T2.EA_2022_T3 EA_2011_T2.EA_2022_T4
         2.340229e-04          2.840422e-04          2.042804e-04
EA_2015_T1.EA_2015_T2 EA_2015_T1.EA_2015_T3 EA_2015_T1.EA_2015_T4
         7.216167e-05          5.353443e-04          4.000097e-05
EA_2015_T1.EA_2022_T1 EA_2015_T1.EA_2022_T2 EA_2015_T1.EA_2022_T3
         1.371440e-04          1.315774e-04          1.163237e-04
EA_2015_T1.EA_2022_T4 EA_2015_T2.EA_2015_T3 EA_2015_T2.EA_2015_T4
         7.212490e-05          5.151115e-04          8.000057e-05
EA_2015_T2.EA_2022_T1 EA_2015_T2.EA_2022_T2 EA_2015_T2.EA_2022_T3
         1.242289e-04          1.416661e-04          1.333631e-04
EA_2015_T2.EA_2022_T4 EA_2015_T3.EA_2015_T4 EA_2015_T3.EA_2022_T1
         1.043142e-04          6.021291e-04          6.011546e-04
EA_2015_T3.EA_2022_T2 EA_2015_T3.EA_2022_T3 EA_2015_T3.EA_2022_T4
         5.817207e-04          6.458708e-04          6.006643e-04
EA_2015_T4.EA_2022_T1 EA_2015_T4.EA_2022_T2 EA_2015_T4.EA_2022_T3
         1.287545e-04          1.227300e-04          1.031789e-04
EA_2015_T4.EA_2022_T4 EA_2022_T1.EA_2022_T2 EA_2022_T1.EA_2022_T3
         7.271113e-05          1.041804e-04          1.297219e-04
EA_2022_T1.EA_2022_T4 EA_2022_T2.EA_2022_T3 EA_2022_T2.EA_2022_T4
         1.114017e-04          1.310401e-04          1.054126e-04
EA_2022_T3.EA_2022_T4
         8.650827e-05

unique_values <- as.list(unique(fst$end))
unique_values
#[[1]] 1000
#[[2]]3000
#[[3]]2000
#[[4]] 4000
#[[5]] 5000
#[[6]] 6000
#[[7]] 7000
#[[8]] 8000
#[[9]] 9000
#[[10]] 10000
#[[11]] 11000
#[[12]] 12000
#[[13]] 13000



#boxplot to check variation between each pairwise comparison, this way we can see which ones are higher
fst2 <- fst[,-(1:4)]
fst2 <- as.data.frame(fst2)
#remove missing data

#count NAs
sumnas <- sum(is.na(fst2))
#it has 8160 NAs

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

 ggsave("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/figures/boxplotfst.pdf",d, w=7.5, h=4.7)
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


#trial heatmap


mean_fst_subset <- colMeans(fst_subset, na.rm = TRUE)
mean_fst_subset

EA_2009_T1.EA_2009_T2 EA_2009_T1.EA_2009_T3 EA_2009_T1.EA_2009_T4
         6.522382e-04          9.119229e-04          8.192032e-04
EA_2009_T2.EA_2009_T3 EA_2009_T2.EA_2009_T4 EA_2009_T3.EA_2009_T4
         2.212941e-04          1.727231e-04          9.451288e-05
EA_2011_T1.EA_2011_T2 EA_2015_T1.EA_2015_T2 EA_2015_T1.EA_2015_T3
         3.314538e-04          7.216167e-05          5.353443e-04
EA_2015_T1.EA_2015_T4 EA_2015_T2.EA_2015_T3 EA_2015_T2.EA_2015_T4
         4.000097e-05          5.151115e-04          8.000057e-05
EA_2015_T3.EA_2015_T4 EA_2022_T1.EA_2022_T2 EA_2022_T1.EA_2022_T3
         6.021291e-04          1.041804e-04          1.297219e-04
EA_2022_T1.EA_2022_T4 EA_2022_T2.EA_2022_T3 EA_2022_T2.EA_2022_T4
         1.114017e-04          1.310401e-04          1.054126e-04
EA_2022_T3.EA_2022_T4
         8.650827e-05

heatmap_matrix <- mean_fst_subset %>% rownames_to_column(var = "Pairwise") %>% pivot_longer(cols = -Pairwise, names_to = "Variable", values_to = "Mean")

#fix data cause something went wrong

mean_fst_matrix <- matrix(
  c(
    6.522382e-04, 9.119229e-04, 8.192032e-04,
    2.212941e-04, 1.727231e-04, 9.451288e-05,
    3.314538e-04, 7.216167e-05, 5.353443e-04,
    4.000097e-05, 5.151115e-04, 8.000057e-05,
    6.021291e-04, 1.041804e-04, 1.297219e-04,
    1.114017e-04, 1.310401e-04, 1.054126e-04,
    8.650827e-05),
 nrow = 19,
 ncol=19,   # Number of pairwise comparisons (rows)
 byrow = TRUE,
 dimnames = list(
    c(
      "EA_2009_T1.EA_2009_T2", "EA_2009_T1.EA_2009_T3", "EA_2009_T1.EA_2009_T4",
      "EA_2009_T2.EA_2009_T3", "EA_2009_T2.EA_2009_T4", "EA_2009_T3.EA_2009_T4",
      "EA_2011_T1.EA_2011_T2", "EA_2015_T1.EA_2015_T2", "EA_2015_T1.EA_2015_T3",
      "EA_2015_T1.EA_2015_T4", "EA_2015_T2.EA_2015_T3", "EA_2015_T2.EA_2015_T4",
      "EA_2015_T3.EA_2015_T4", "EA_2022_T1.EA_2022_T2", "EA_2022_T1.EA_2022_T3",
      "EA_2022_T1.EA_2022_T4", "EA_2022_T2.EA_2022_T3", "EA_2022_T2.EA_2022_T4",
      "EA_2022_T3.EA_2022_T4"),
    c(
      "EA_2009_T1.EA_2009_T2", "EA_2009_T1.EA_2009_T3", "EA_2009_T1.EA_2009_T4",
      "EA_2009_T2.EA_2009_T3", "EA_2009_T2.EA_2009_T4", "EA_2009_T3.EA_2009_T4",
      "EA_2011_T1.EA_2011_T2", "EA_2015_T1.EA_2015_T2", "EA_2015_T1.EA_2015_T3",
      "EA_2015_T1.EA_2015_T4", "EA_2015_T2.EA_2015_T3", "EA_2015_T2.EA_2015_T4",
      "EA_2015_T3.EA_2015_T4", "EA_2022_T1.EA_2022_T2", "EA_2022_T1.EA_2022_T3",
      "EA_2022_T1.EA_2022_T4", "EA_2022_T2.EA_2022_T3", "EA_2022_T2.EA_2022_T4",
      "EA_2022_T3.EA_2022_T4")
  	)
  )

# Plot a heatmap using the heatmap() function

# Get the row and column names from dimnames
row_names <- rownames(mean_fst_matrix)
col_names <- colnames(mean_fst_matrix)

# Calculate the number of rows and columns in the matrix
n_rows <- nrow(mean_fst_matrix)
n_cols <- ncol(mean_fst_matrix)

#list diagonal values
 diagonal_values <- diag(mean_fst_matrix)

# Create a data frame with the diagonal values
diagonal_data <- data.frame(
  Pairwise = rownames(mean_fst_matrix),
  Value = diagonal_values
)

# Plot the heatmap of diagonal values
d <- ggplot(diagonal_data, aes(x = Pairwise, y = Pairwise, fill = Value)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "blue") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    axis.ticks = element_blank()
  ) +
  labs(x = NULL, y = NULL)

print(heatmap_plot)
ggsave("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/figures/heatmap-fst-window1000-maincomparisons.pdf",d, w=10, h=8)




#plot full fst
packages <- c("ComplexHeatmap", "grid", "Rcpp", "stats", "base", "plotrix", "wesanderson", "textshape")
install.packages(packages)
#install complexheatmap for older version
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("ComplexHeatmap")
library("ComplexHeatmap")
library("grid")
library("Rcpp")
library("stats")
library("base")
library("plotrix")
library("wesanderson")
library("textshape")
install.packages("pheatmap")
library(pheatmap)

fst_cleaned <- na.omit(fst)
missing_values <- any(is.na(fst))
if (missing_values) {
  cat("Missing values detected in the data.\n")
} else {
  cat("No missing values found in the data.\n")
}


# Remove row names and column names for heatmap
data_for_heatmap <- fst_cleaned[, 5:ncol(fst)]  # Exclude the first column (chrom)

fst_cleaned$countdown <- seq(nrow(fst_cleaned), 1)
fst_cleaned$full_name <- paste0(fst_cleaned$chrom, fst_cleaned$countdown)
rownames(data_for_heatmap) <- fst_cleaned$full_name


# Convert data to a matrix
mat <- as.matrix(data_for_heatmap)

# Create a heatmap using pheatmap
d <- pheatmap(mat,
         scale = "none",  # You can customize scaling if needed
         col = colorRampPalette(c("white", "blue"))(100),
         fontsize_row = 8,
         fontsize_col = 8,
      
         annotation_legend = TRUE, annotation_names_row = FALSE, annotation_names_col = TRUE,
         drop_levels = FALSE, show_rownames = F, show_colnames = T, 
         main = "Heatmap of Data")

ggsave("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/figures/trial-heatmap-fst-window1000.pdf",d, w=10, h=8)

########################snps per bin####################
#####comparison between 1000 and 5000 bin sizes#########
########################################################


fst <- read.csv("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/grenedalf-freq/vcf-maf0.01-misfrac0.1-maxcov95/grenedalfst.log/fst-all-pairwise.csv", header = TRUE, na.strings = "")
#negative values to zeros

fst <- fst %>% mutate_all(function(x) ifelse(x < 0, 0, x))
fst_cleaned <- na.omit(fst)

sum_smaller_than_10 <- sum(fst_cleaned$snps < 2)

# Print the result
cat("Number of entries smaller than 10:", sum_smaller_than_10, "\n")
#9596
#1349 smaller than 2

#create a new colum for each bin
#fst_cleaned$bins <- paste(fst_cleaned$start, fst_cleaned$end, sep = "-")
fst_cleaned$bins <- paste(fst_cleaned$chrom, fst_cleaned$start, fst_cleaned$end,sep = "-")
fst_cleaned$bins<- factor(fst_cleaned$bins)

#number of bins
num_different_values <- length(unique(fst_cleaned$bins))
num_different_values
#42728 chroms
#2.5% of the bins only have 1 SNPs, perhaps using 5K windows makes more sense 

#plot scatter plot per bin
pdf(file="/gxfs_home/geomar/smomw573/work/seasonal_adaptation/figures/trial-snps-per.bin-1000.pdf", height=6,width=10)
plot(fst_cleaned$bins, fst_cleaned$snps, main = "Scatter Plot", xlab = "Bins", ylab = "N SNPs", ylim = c(0, 800), pch = 16, col = "blue",  xaxt = "n")
dev.off()


###########################
#for 5000K window size #
###########################

fst <- read.csv("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/grenedalf-freq/vcf-maf0.01-misfrac0.1-maxcov95/grenedalfst.log/fstall-parwise-5000.csv", header = TRUE, na.strings = "")
#negative values to zeros

fst <- fst %>% mutate_all(function(x) ifelse(x < 0, 0, x))

#means of each colums (are the some higher then others)

mean_fst <- colMeans(fst[, 5:ncol(fst)], na.rm = TRUE)

#scatterplot number of snps per window

mean_fst 


#boxplot to check variation between each pairwise comparison, this way we can see which ones are higher

fst2 <- as.data.frame(fst2)
#remove missing data

#count NAs
sumnas <- sum(is.na(fst2))
#2411 less missing data than 1K windows

### remove mising data
data_clean<- na.omit(fst2)

# Reshape the data into long format
data_long <- tidyr::gather(data_clean, key = "Column", value = "Value")
data_long$Column <- as.factor(data_long$Column)

# Create boxplots using ggplot2
d <- ggplot(data_long, aes(x=Column, y=Value)) + 
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

 ggsave("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/figures/boxplotfst-5Kwindows.pdf",d, w=7.5, h=4.7)

#subset for pairwise comparisons of interest

column_names <- colnames(fst2)
print(column_names)
 

fst_subset<- data_clean[, c("EA_2009_T1.EA_2009_T2", "EA_2009_T1.EA_2009_T3", "EA_2009_T1.EA_2009_T4", "EA_2009_T2.EA_2009_T3", "EA_2009_T2.EA_2009_T4","EA_2009_T3.EA_2009_T4","EA_2011_T1.EA_2011_T2",  "EA_2015_T1.EA_2015_T2", "EA_2015_T1.EA_2015_T3", 
	"EA_2015_T1.EA_2015_T4", "EA_2015_T2.EA_2015_T3", "EA_2015_T2.EA_2015_T4", "EA_2015_T3.EA_2015_T4","EA_2022_T1.EA_2022_T2",
	"EA_2022_T1.EA_2022_T3", "EA_2022_T1.EA_2022_T4", "EA_2022_T2.EA_2022_T3", "EA_2022_T2.EA_2022_T4", "EA_2022_T3.EA_2022_T4")]


#create new variable year

fst_subset <- as.data.frame(fst_subset)
#data_clean <- fst_subset
data_long <- tidyr::gather(fst_subset, key = "Column", value = "Value") # gather() function from the tidyr package to convert your dataset to long format.
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

 ggsave("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/figures/boxplot-fst-main-comparisons-5Kwindows.pdf",d, w=12, h=4.7)

#scater plots
d <- ggplot(data_long, aes(x = Column, y = Value)) +
  geom_point() + #~ Column indicates that each facet corresponds to a different column.
  facet_grid(. ~ year, scales = "free") + # The scales = "free" argument allows each facet to have its own y-axis scale.
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Pairwise comparison", y = "Fst")


 ggsave("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/figures/boxplot-fst-scatter-main-comparisons-5kwindows.pdf",d, w=12, h=6)


# Create a scatter plot grid using facet_grid


#trial heatmap


mean_fst_subset <- colMeans(fst_subset, na.rm = TRUE)
mean_fst_subset

heatmap_matrix <- mean_fst_subset %>% rownames_to_column(var = "Pairwise") %>% pivot_longer(cols = -Pairwise, names_to = "Variable", values_to = "Mean")

#fix data cause something went wrong

mean_fst_matrix <- matrix(
  c(
     6.251438e-04, 8.726515e-04, 7.699255e-04,
     2.069383e-04, 1.547305e-04, 8.344243e-05,
     3.222143e-04, 6.337224e-05, 5.170126e-04,
     3.289386e-05, 4.931738e-04, 7.707632e-05, 
     5.848110e-04, 9.421378e-05, 1.233745e-04, 
     1.060758e-04, 1.230012e-04, 9.586955e-05, 8.112256e-05),
 nrow = 19,
 ncol=19,   # Number of pairwise comparisons (rows)
 byrow = TRUE,
 dimnames = list(
    c(
      "EA_2009_T1.EA_2009_T2", "EA_2009_T1.EA_2009_T3", "EA_2009_T1.EA_2009_T4",
      "EA_2009_T2.EA_2009_T3", "EA_2009_T2.EA_2009_T4", "EA_2009_T3.EA_2009_T4",
      "EA_2011_T1.EA_2011_T2", "EA_2015_T1.EA_2015_T2", "EA_2015_T1.EA_2015_T3",
      "EA_2015_T1.EA_2015_T4", "EA_2015_T2.EA_2015_T3", "EA_2015_T2.EA_2015_T4",
      "EA_2015_T3.EA_2015_T4", "EA_2022_T1.EA_2022_T2", "EA_2022_T1.EA_2022_T3",
      "EA_2022_T1.EA_2022_T4", "EA_2022_T2.EA_2022_T3", "EA_2022_T2.EA_2022_T4",
      "EA_2022_T3.EA_2022_T4"),
    c(
      "EA_2009_T1.EA_2009_T2", "EA_2009_T1.EA_2009_T3", "EA_2009_T1.EA_2009_T4",
      "EA_2009_T2.EA_2009_T3", "EA_2009_T2.EA_2009_T4", "EA_2009_T3.EA_2009_T4",
      "EA_2011_T1.EA_2011_T2", "EA_2015_T1.EA_2015_T2", "EA_2015_T1.EA_2015_T3",
      "EA_2015_T1.EA_2015_T4", "EA_2015_T2.EA_2015_T3", "EA_2015_T2.EA_2015_T4",
      "EA_2015_T3.EA_2015_T4", "EA_2022_T1.EA_2022_T2", "EA_2022_T1.EA_2022_T3",
      "EA_2022_T1.EA_2022_T4", "EA_2022_T2.EA_2022_T3", "EA_2022_T2.EA_2022_T4",
      "EA_2022_T3.EA_2022_T4")
  	)
  )

# Plot a heatmap using the heatmap() function

# Get the row and column names from dimnames
row_names <- rownames(mean_fst_matrix)
col_names <- colnames(mean_fst_matrix)

# Calculate the number of rows and columns in the matrix
n_rows <- nrow(mean_fst_matrix)
n_cols <- ncol(mean_fst_matrix)

#list diagonal values
 diagonal_values <- diag(mean_fst_matrix)

# Create a data frame with the diagonal values
diagonal_data <- data.frame(
  Pairwise = rownames(mean_fst_matrix),
  Value = diagonal_values
)

# Plot the heatmap of diagonal values
d <- ggplot(diagonal_data, aes(x = Pairwise, y = Pairwise, fill = Value)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "blue") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    axis.ticks = element_blank()
  ) +
  labs(x = NULL, y = NULL)

print(heatmap_plot)
ggsave("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/figures/heatmap-fst-window5000-maincomparisons.pdf",d, w=10, h=8)


#bin calculations

sum_smaller_than_10 <- sum(data_clean$snps < 10)
sum_smaller_than_10
# Print the result
cat("Number of entries smaller than 10:", sum_smaller_than_10, "\n")
#0
data_clean<- na.omit(fst)

#create a new colum for each bin
#fst_cleaned$bins <- paste(fst_cleaned$start, fst_cleaned$end, sep = "-")
data_clean$bins <- paste(data_clean$chrom, data_clean$start, data_clean$end,sep = "-")
data_clean$bins<- factor(data_clean$bins)

#number of bins
num_different_values <- length(unique(data_clean$bins))
num_different_values
#41161
#0 % of the bins has less than 10 snps per bin

#plot scatter plot per bin
pdf(file="/gxfs_home/geomar/smomw573/work/seasonal_adaptation/figures/trial-snps-per.bin-5000.pdf", height=6,width=10)
plot(data_clean$bins, data_clean$snps, main = "Scatter Plot", xlab = "Bins", ylab = "N SNPs", pch = 16, col = "blue",  xaxt = "n")
dev.off()
```

# Manhattan plot of single snp by snp Fst calculation

``` bash
srun --pty --x11 --nodes=1 --cpus-per-task=4 --mem=32000 --time=06:00:00 /bin/bash
module load R/4.1.1
module load gcc/10.2.0
R
```
```R
library(pcadapt)
library(data.table)
library(qqman)
library(ggplot2)
library(dplyr)
library(tidyr)



#input data must be a matrix with samples as row names and fst values as colums
freq <- read.csv("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/grenedalf-freq/vcf-maf0.01-misfrac0.1-maxcov95/grenedalfst.log/fstsingle-snps.csv", header = TRUE)

# Count rows with NaN values
sum(apply(freq, 1, function(row) any(is.nan(row))))
freq_clean <- na.omit(freq)
save(freq_clean, file = "/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/grenedalf-freq/vcf-maf0.01-misfrac0.1-maxcov95/grenedalfst.log/freq_clean.RData")

#########################################################
######trial boxplot for relevant pairwise comparisons####
#########################################################

#Processed 41599 chromosomes with 1258282 (non-filtered) positions in 6872734 windows.
#Total filter summary (after applying all sample filters):
#Passed:               1258282

freq_clean<- na.omit(freq_clean)
freq_clean2  <- freq_clean  %>% mutate_all(function(x) ifelse(x < 0, 0, x))



fst_subset<- freq_clean2[, c("EA_2009_T1.EA_2009_T2", "EA_2009_T1.EA_2009_T3", "EA_2009_T1.EA_2009_T4", "EA_2009_T2.EA_2009_T3", "EA_2009_T2.EA_2009_T4","EA_2009_T3.EA_2009_T4","EA_2011_T1.EA_2011_T2",  "EA_2015_T1.EA_2015_T2", "EA_2015_T1.EA_2015_T3", 
  "EA_2015_T1.EA_2015_T4", "EA_2015_T2.EA_2015_T3", "EA_2015_T2.EA_2015_T4", "EA_2015_T3.EA_2015_T4","EA_2022_T1.EA_2022_T2",
  "EA_2022_T1.EA_2022_T3", "EA_2022_T1.EA_2022_T4", "EA_2022_T2.EA_2022_T3", "EA_2022_T2.EA_2022_T4", "EA_2022_T3.EA_2022_T4")]
fst2 <- fst_subset[,-(1:4)]
fst2  <- fst2  %>% mutate_all(function(x) ifelse(x < 0, 0, x))
##############################################continuation of analyisis


df_long <- freq_clean %>%
  pivot_longer(cols = 5:95, names_to = "Chrom", values_to = "Value")
df_long$combined <- paste(df_long$chrom, df_long$start)



#ad a numeric column for chrom
data <- df_long %>%
  mutate(number = match(chrom, unique(chrom)))
save("/gxfs_work1/geomar/smomw573/seasonal_adaptation/analysis/grenedalf-freq/vcf-maf0.01-misfrac0.1-maxcov95/grenedalfst.log/formated_table_FSTSINGLE.RData")

#simple plot
x <- seq_along(df_long$chrom)

pdf(file="/gxfs_home/geomar/smomw573/work/seasonal_adaptation/figures/singlesnpfst.pdf", height=6,width=10)
plot(data$number, data$Value, xlab = "Index", ylab = "Fst", main = "Single SNP Fst", xaxt = "n")
dev.off()

#doesnt work with FST data?

library(qqman)
load("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/grenedalf-freq/vcf-maf0.01-misfrac0.1-maxcov95/grenedalfst.log/longFST.RData")
#plot manhattan plot long data
pdf(file="/gxfs_home/geomar/smomw573/work/seasonal_adaptation/figures/manhattan-singlesnps.pdf", height=6,width=10)
manhattan(data, chr="number", bp="start", snp="combined", p="Value" )
dev.off()

###

load("/gxfs_work1/geomar/smomw573/seasonal_adaptation/analysis/grenedalf-freq/vcf-maf0.01-misfrac0.1-maxcov95/grenedalfst.log/formated_table_FSTSINGLE.RData")
d <- ggplot(data, aes(x=number, y=Value)) +
  geom_point(size=2, shape=23)+
    theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

ggsave("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/figures/trialscatterfst.pdf",d, w=7.5, h=4.7)
####
#calculate number of SNPs per chrom



# Define the file path for the new CSV file (same folder)
output_file <- "/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/grenedalf-freq/vcf-maf0.01-misfrac0.1-maxcov95/grenedalfst.log/fstsingle-snps-cleaned.csv"

# Save the cleaned data frame as a new CSV file
write.csv(freq_clean, file = output_file, row.names = FALSE)

cat("Cleaned data frame has been saved as:", output_file, "\n")

#why dont we have an FST for all snps?
nrow(freq_clean)
# FST values 1098611 from 1258282 SNPs in original vcf file, why?


freq_clean <- freq_clean[, !(names(freq_clean) %in% c("start", "end", "snps"))]
freq_clean  <- freq_clean  %>% mutate_all(function(x) ifelse(x < 0, 0, x))
save(freq_clean, file = "nonanozeros-allsnps.rds")
df <- freq_clean[-1,-1]

count_zero_columns <- sum(rowSums(df) == 0)
print(count_zero_columns)
# 555450

trans<- t(freq_clean)
mtrans <- as.matrix(trans)
mtrans<-matrix(mtrans, ncol= ncol(trans), dimnames =NULL)


#create pcadapt file
output_file <- "/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/grenedalf-freq/vcf-maf0.01-misfrac0.1-maxcov95/grenedalfst.log/pcadapt-input.csv"

# Save the cleaned data frame as a new CSV file
write.table(mtrans, file = output_file, row.names = FALSE, col.names = FALSE, quote=FALSE)

#heatmap FST snp per snp

library("ComplexHeatmap")
library("grid")
library("Rcpp")
library("stats")
library("base")
library("plotrix")
library("wesanderson")
library("textshape")
library(pheatmap)

load ("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/grenedalf-freq/vcf-maf0.01-misfrac0.1-maxcov95/grenedalfst.log/freq_clean.RData")


freq_clean  <- freq_clean  %>% mutate_all(function(x) ifelse(x < 0, 0, x))
check <- freq_clean[c(1:100),c(1:95)]

# Remove row names and column names for heatmap
data_for_heatmap <- freq_clean[, 5:ncol(fst)]  # Exclude the first column (chrom)

freq_clean$countdown <- seq(nrow(freq_clean), 1)
freq_clean$full_name <- paste0(freq_clean$chrom, freq_clean$countdown)
rownames(data_for_heatmap) <- freq_clean$full_name

# Convert data to a matrix
mat <- as.matrix(data_for_heatmap)
subset_mat <- mat[,c(1,2,3,64,65,66,86,87,88,89,90,91)]

#freq_clean_sub <- data_for_heatmap[,c(1,2,3,64,65,66,86,87,88,89,90,91)]
#plot mean of single snps across all samples

mean_fst <- rowMeans(data_for_heatmap)
colnames(mean_fst) <- c("chrom", "mean_fst")
pdf(file="/gxfs_home/geomar/smomw573/work/seasonal_adaptation/figures/singlesnpfst_mean.pdf", height=6,width=10)
plot(mean_fst, xlab = "Index", ylab = "Fst", main = "Mean Single SNP Fst")
dev.off()

#plot histogram

#transform data columns into rows
load ("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/grenedalf-freq/vcf-maf0.01-misfrac0.1-maxcov95/grenedalfst.log/freq_clean.RData")
freq_clean  <- freq_clean  %>% mutate_all(function(x) ifelse(x < 0, 0, x))
slot1 <- freq_clean[, 5:ncol(fst)]  # Exclude the first column (chrom)


#frombefore


new_row_names <- paste0(row.names(slot1), "_", colnames(slot1)[-1])
# Assign the new row names to the dataframe
row.names(data) <- new_row_names

# Reshape the data
reshaped_data <- slot1 %>%
  tibble::rownames_to_column(var = "row_name") %>%
  pivot_longer(-row_name, names_to = "Column", values_to = "Value")
#plot histogram
pdf(file="/gxfs_home/geomar/smomw573/work/seasonal_adaptation/figures/histogram_singlesnp_fst.pdf", height=6,width=10)
hist_data <- hist(reshaped_data$Value, xlab="breaks", ylab= "Frequency", breaks=30)
# Calculate midpoints of each bin
midpoints <- (hist_data$breaks[-1] + hist_data$breaks[-length(hist_data$breaks)]) / 2
# Get counts in each bin
counts <- hist_data$counts
text(hist_data$mids, hist_data$counts, labels = hist_data$counts, pos = 3)
dev.off()

sum(reshaped_data$Value > 0.95)


#freq table
slot1 <- freq_clean[, c(1, (5:ncol(freq_clean)))]
transformed_data <- slot1 %>%
  pivot_longer(cols = -chrom, names_to = "Column", values_to = "Value") %>%
  select(chrom, Column, Value)

freq <- as.data.frame(table(transformed_data$chrom))
sum(transformed_data$Value == 0)
#93805542
sum(transformed_data$Value > 0)
#2872226


#plot SNPs per chrom
pdf(file="/gxfs_home/geomar/smomw573/work/seasonal_adaptation/figures/nsnpsperchrom.pdf", height=6,width=10)
plot(freq$Freq, pch = 16, col = "black", xlab = "Index", ylab = "n SNPs", main = "SNPs per chrom")
dev.off()
mean(freq$Freq)
#2361.277
```

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
ea.readcount30X <- vcf2pooldata(vcf.file="depth-corrected.recode.vcf.gz",poolsizes=rep(50,14))
# summary of the resulting pooldata object
ea.readcount30X

#ha?
selected.snps.idx <- as.numeric(sub("rs","",rownames(ea.readcount30X@snp.info)))
head(selected.snps.idx)

#estimate genome wide Fst across all the popuatlions
ea.readcount30X.fst<-computeFST(ea.readcount30X)
ea.readcount30X.fst$FST
#genome wide Fst is -0.01960076, so 0

# Block-Jackknife estimation of FST standard-error and confidence intervals:
ea.readcount30X.fst<-computeFST(ea.readcount30X,nsnp.per.bjack.block = 1000, verbose=FALSE)
ea.readcount30X.fst$FST
#same value  -0.01960076

ea.readcount30X.fst$mean.fst #block-jacknife estimate of s.e.
#-0.02009452
ea.readcount30X.fst$se.fst #s.e. of the genome-wide Fst estimate
#2.90547e-05
ea.readcount30X.fst$mean.fst+c(-1.96,1.96)*ea.readcount30X.fst$se.fst
# -0.02015146 -0.02003757


#Computing multi-locus FST to scan the genome over sliding-windows of SNPs

ea.readcount30X.fst<-computeFST(ea.readcount30X,sliding.window.size=50)
ea.readcount30X.fst<-computeFST(ea.readcount30X,sliding.window.size=100)
ea.readcount30X.fst<-computeFST(ea.readcount30X,sliding.window.size=10)
#we have 42K scaffolds so loosing a lot of data with sliding window size of 50-100 chrom?


#6722 chromosomes scanned (with more than 50 SNPs)
#Average (min-max) Window Sizes 0.3 ( 0.1 - 4.9 ) kb


#931 chromosomes scanned (with more than 100 SNPs)
#Average (min-max) Window Sizes 0.8 ( 0.2 - 7.4 ) kb


#31251 chromosomes scanned (with more than 10 SNPs)
#Average (min-max) Window Sizes 0 ( 0 - 6.9 ) kb

#I will just play around with 100 as I dont really understand
ea.readcount30X.fst<-computeFST(ea.readcount30X,sliding.window.size=100)

plot(ea.readcount30X.fst$sliding.windows.fst$CumulatedPosition/1e6,
     ea.readcount30X.fst$sliding.windows.fst$MultiLocusFst,
     xlab="Cumulated Position (in Mb)",ylab="Muli-locus Fst")
     #col=as.numeric(ea.readcount30X.fst$sliding.windows.fst$Chr),pch=16) Doesnt work as we dont have chromossome numbers
abline(h=ea.readcount30X.fst$FST,lty=2)

head(ea.readcount30X.fst$sliding.windows.fst$CumulatedPosition/1e6)
head(ea.readcount30X.fst$sliding.windows.fst$MultiLocusFst)
head(ea.readcount30X.fst$sliding.windows.fst$Chr)

#Manhattan plot of the multi-locus FST computed over sliding-windows of 50 SNPs on the PoolSeq example data. The dashed line indicates the estimated overall genome-wide FST . The 20 simulated
#chromosomes are represented by alternate colors


#pairwise FST

ea.pairwisefst<-compute.pairwiseFST(ea.readcount30X,verbose=FALSE)
#heatmap
#Heatmap representing the pairwise-population FST matrix of the 14 populations of the 30XPool-Seq example data set
heatmap(ea.pairwisefst)
#it moves pops which are more similar to each other 

#Block-Jackknife estimation of FST standard-error and visualisation of confidence intervals
ea.pairwisefst@PairwiseFSTmatrix
plot(ea.pairwisefst)

#the fst plot with the new fst dataset fixed for the new genome is on a different file - plot-single-fst.md in my computer
```

## Trial Fst Boxplot and T test

```bash
srun --pty --x11 --nodes=1 --cpus-per-task=4 --mem=32000 --time=06:00:00 /bin/bash
module load R/4.1.1
module load gcc/10.2.0
module load icu4c/67.1

R
```

```R
library(ggplot2)
library(dplyr)
library(ggpubr)

fst <- read.csv("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/grenedalf-freq/vcf-maf0.01-misfrac0.1-maxcov95/grenedalfst.log/fstsingle-snps-newvcf-queue.csv", header = TRUE, na.strings = "")
#negative values to zeros

fst <- fst %>% mutate_all(function(x) ifelse(x < 0, 0, x))


fst2009a <- fst[,c("chrom", "EA_2009_T1.EA_2009_T2", "EA_2009_T1.EA_2009_T3")]
fst2009b <- fst[,c("chrom", "EA_2009_T1.EA_2009_T2", "EA_2009_T1.EA_2009_T4")]
fst2009c <- fst[,c("chrom", "EA_2009_T1.EA_2009_T3", "EA_2009_T1.EA_2009_T4")]
fst2009all <- fst[,c("chrom", "EA_2009_T1.EA_2009_T2","EA_2009_T1.EA_2009_T3", "EA_2009_T1.EA_2009_T4")]

fst2011 <- fst[,c("chrom", "EA_2011_T1.EA_2011_T2")]

fst2015a <- fst[,c("chrom", "EA_2015_T1.EA_2015_T2", "EA_2015_T1.EA_2015_T3")]
fst2015b <- fst[,c("chrom", "EA_2015_T1.EA_2015_T2", "EA_2015_T1.EA_2015_T4")]
fst2015c <- fst[,c("chrom", "EA_2015_T1.EA_2015_T3", "EA_2015_T1.EA_2015_T4")]
fst2015all <- fst[,c("chrom", "EA_2015_T1.EA_2015_T2", "EA_2015_T1.EA_2015_T3", "EA_2015_T1.EA_2015_T4")]

fst2022a <- fst[,c("chrom", "EA_2022_T1.EA_2022_T2","EA_2022_T1.EA_2022_T3")]
fst2022b <- fst[,c("chrom", "EA_2022_T1.EA_2022_T2", "EA_2022_T1.EA_2022_T4")]
fst2022c <- fst[,c("chrom", "EA_2022_T1.EA_2022_T3", "EA_2022_T1.EA_2022_T4")]
fst2022all <- fst[,c("chrom", "EA_2022_T1.EA_2022_T2", "EA_2022_T1.EA_2022_T3", "EA_2022_T1.EA_2022_T4")]

#prep data for t-test
#2009a 
t2009a <- tidyr::gather(fst2009a, "pairwise", "fst", 2:3)
t2009$pairwise <- as.factor(t2009$pairwise)
t2009a <- na.omit(t2009a)

#2009b
t2009b <- tidyr::gather(fst2009b, "pairwise", "fst", 2:3)
t2009b$pairwise <- as.factor(t2009b$pairwise)
t2009b <- na.omit(t2009b)

#2009c
t2009c <- tidyr::gather(fst2009c, "pairwise", "fst", 2:3)
t2009c$pairwise <- as.factor(t2009c$pairwise)
t2009c <- na.omit(t2009c)

t2009all <- tidyr::gather(fst2009all, "pairwise", "fst", 2:4)

#t2011 is okay and cant to the t-test

t2011all <- tidyr::gather(fst2011, "pairwise", "fst", 2)

#2015
t2015a <- tidyr::gather(fst2015a, "pairwise", "fst", 2:3)
t2015a$pairwise <- as.factor(t2015a$pairwise)
t2015a <- na.omit(t2015a)

t2015b <- tidyr::gather(fst2015b, "pairwise", "fst", 2:3)
t2015b$pairwise <- as.factor(t2015b$pairwise)
t2015b <- na.omit(t2015b)

t2015c <- tidyr::gather(fst2015c, "pairwise", "fst", 2:3)
t2015c$pairwise <- as.factor(t2015c$pairwise)
t2015c <- na.omit(t2015c)

t2015all <- tidyr::gather(fst2015all, "pairwise", "fst", 2:4)

#2022
t2022a <- tidyr::gather(fst2022a, "pairwise", "fst", 2:3)
t2022a$pairwise <- as.factor(t2022a$pairwise)
t2022a <- na.omit(t2022a)

t2022b <- tidyr::gather(fst2022b, "pairwise", "fst", 2:3)
t2022b$pairwise <- as.factor(t2022b$pairwise)
t2022b <- na.omit(t2022b)

t2022c <- tidyr::gather(fst2022c, "pairwise", "fst", 2:3)
t2022c$pairwise <- as.factor(t2022c$pairwise)
t2022c <- na.omit(t2022c)

t2022all <- tidyr::gather(fst2022all, "pairwise", "fst", 2:4)

#check if all pairwise comparisons are there
unique(new_dataset$pairwise)

#t-test for each of the pairwise comparisons and years

#2009a
t.test(t2009a$fst ~ t2009a$pairwise)

#        Welch Two Sample t-test

#data:  t2009a$fst by t2009a$pairwise
#t = 2.9031, df = 2360868, p-value = 0.003695
#alternative hypothesis: true difference in means between group EA_2009_T1.EA_2009_T2 and group EA_2009_T1.EA_2009_T3 is not equal to 0
#95 percent confidence interval:
# 1.089348e-05 5.617110e-05
#sample estimates:
#mean in group EA_2009_T1.EA_2009_T2 mean in group EA_2009_T1.EA_2009_T3
#                        0.002211764                         0.002178232

#2009b
t.test(t2009b$fst ~ t2009b$pairwise)

#        Welch Two Sample t-test

#data:  t2009b$fst by t2009b$pairwise
#t = 19.21, df = 2370704, p-value < 2.2e-16
#alternative hypothesis: true difference in means between group EA_2009_T1.EA_2009_T2 and group EA_2009_T1.EA_2009_T4 is not equal to 0
#95 percent confidence interval:
# 0.0001926584 0.0002364383
#sample estimates:
#mean in group EA_2009_T1.EA_2009_T2 mean in group EA_2009_T1.EA_2009_T4
#                        0.002211764                         0.001997216


t.test(t2009c$fst ~ t2009c$pairwise)

#2015a
#2015a
t.test(t2015a$fst ~ t2015a$pairwise)

#        Welch Two Sample t-test

#data:  t2015a$fst by t2015a$pairwise
#t = -103.22, df = 1692658, p-value < 2.2e-16
#alternative hypothesis: true difference in means between group EA_2015_T1.EA_2015_T2 and group EA_2015_T1.EA_2015_T3 is not equal to 0
#95 percent confidence interval:
# -0.0007057102 -0.0006794093
#sample estimates:
#mean in group EA_2015_T1.EA_2015_T2 mean in group EA_2015_T1.EA_2015_T3
#                       0.0005338422                        0.0012264019


#2015b
t.test(t2015b$fst ~ t2015b$pairwise)

#        Welch Two Sample t-test

#data:  t2015b$fst by t2015b$pairwise
#t = 39.981, df = 2287976, p-value < 2.2e-16
#alternative hypothesis: true difference in means between group EA_2015_T1.EA_2015_T2 and group EA_2015_T1.EA_2015_T4 is not equal to 0
#95 percent confidence interval:
# 0.0001376141 0.0001518021
#sample estimates:
#mean in group EA_2015_T1.EA_2015_T2 mean in group EA_2015_T1.EA_2015_T4
#                       0.0005338422                        0.0003891341

t.test(t2015c$fst ~ t2015c$pairwise)

#2022a
t.test(t2022a$fst ~ t2022a$pairwise)

#        Welch Two Sample t-test


#t = 10.641, df = 2394592, p-value < 2.2e-16
#alternative hypothesis: true difference in means between group EA_2022_T1.EA_2022_T2 and group EA_2022_T1.EA_2022_T3 is not equal to 0
#95 percent confidence interval:
# 3.986161e-05 5.786058e-05
#sample estimates:
#mean in group EA_2022_T1.EA_2022_T2 mean in group EA_2022_T1.EA_2022_T3
#                       0.0006455680                        0.0005967069


#2022b
t.test(t2022b$fst ~ t2022b$pairwise)

#        Welch Two Sample t-test

#data:  t2022b$fst by t2022b$pairwise
#t = 27.349, df = 2367138, p-value < 2.2e-16
#alternative hypothesis: true difference in means between group EA_2022_T1.EA_2022_T2 and group EA_2022_T1.EA_2022_T4 is not equal to 0
#95 percent confidence interval:
# 0.0001121237 0.0001294348
#sample estimates:
#mean in group EA_2022_T1.EA_2022_T2 mean in group EA_2022_T1.EA_2022_T4
#                       0.0006507826                        0.0005300033

t.test(t2022c$fst ~ t2022c$pairwise)


#fst_subset2<- fst2[, c( "EA_2009_T2.EA_2009_T3", "EA_2009_T2.EA_2009_T4","EA_2009_T3.EA_2009_T4",  
#                        " EA_2015_T2.EA_2015_T3", "EA_2015_T2.EA_2015_T4", "EA_2015_T3.EA_2015_T4", 
#                        "EA_2022_T2.EA_2022_T3", "EA_2022_T2.EA_2022_T4", "EA_2022_T3.EA_2022_T4")]

#boxplots

a <- ggplot(t2009all, aes(x=pairwise, y=fst)) + 
  geom_boxplot(notch=TRUE, outlier.colour="#D3DDDC", outlier.shape=8,
                outlier.size=2) + 
  coord_flip() +
    theme_bw()
#ggsave("fst2009comparisons.pdf", width = 7, height = 7, units =c("cm"), dpi =3000)

b <- ggplot(t2015all, aes(x=pairwise, y=fst)) + 
  geom_boxplot(notch=TRUE, outlier.colour="#F2AD00", outlier.shape=8,
                outlier.size=2) + 
  coord_flip() +
    theme_bw()
#ggsave("fst2015comparisons.pdf", width = 7, height = 7, units =c("cm"), dpi =3000)

c <- ggplot(t2022all, aes(x=pairwise, y=fst)) + 
  geom_boxplot(notch=TRUE, outlier.colour="#00A08A", outlier.shape=8,outlier.size=2) + 
  coord_flip() +
    theme_bw()
#ggsave("fst2022comparisons.pdf", width = 7, height = 7, units =c("cm"), dpi =3000)

d<- ggplot(t2011all, aes(x=pairwise, y=fst)) + 
  geom_boxplot(notch=TRUE, outlier.colour="#6699CC", outlier.shape=8, outlier.size=2) + 
  coord_flip() +
    theme_bw()


#merge 4 plots
figure <- ggarrange(a, d, c, b,
                    labels = c("2009", "2011", "2015", "2022"),
                    ncol = 1, nrow = 4)
figure
ggsave("fst-allyearpairwisecomparisons.jpeg", width = 16, height = 20, units =c("cm"))

```

# tips 
##install local CRAN miror repository

wget
local_cran_mirror <- "file:///gxfs_home/geomar/smomw573/R/CRAN_mirrors.csv"
options(repos = structure(c(CRAN = local_cran_mirror)))
# Load the CRAN_mirrors.csv file
file_path <- "/gxfs_home/geomar/smomw573/R/CRAN_mirrors.csv"
mirrors_data <- read.csv(file_path)

# Extract the URL column from the data
mirror_urls <- mirrors_data$URL

# Set the repos option using the extracted URLs
options(repos = structure(c(CRAN = mirror_urls)))

# Now you can install packages from the specified CRAN mirrors
install.packages("dplyr")
#something went wrong so I can just use this function
install.packages("tidyr",        # Using repos argument
                 repos = "https://cran.uni-muenster.de/")

#####


######

# Generate new Freq estimation using the final vcf file produced by Reid

# generate individual frequencies grenedalf

```bash
module load htslib/1.10.2  bzip2/1.0.8  

GRENEDALF=/gxfs_home/geomar/smomw573/software/grenedalf/bin/grenedalf
GENOME=/gxfs_work1/geomar/smomw573/seasonal_adaptation/genome/affinis.Atlantic.long_read_draft.Mar22.fasta \
VCF=/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/variants/final_20230810.vcf.gz \
OUTPUT=/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/grenedalf-freq/vcf-maf0.01-misfrac0.1-maxcov95/grenedalfst.log/freq-finalvcf/

$GRENEDALF frequency --vcf-path $VCF --allow-file-overwriting --write-sample-alt-freq --file-suffix sample-frequencies --out-dir $OUTPUT/grenedalftrial.log

#need to fix the new vcf
#this did not work so there is something going on on the vcf but we dont care

grep "^##" final_20230810.vcf > final_20230810.vcf.header.txt
head final_20230810.vcf.header.txt
grep -v "^##" final_20230810.vcf > final_20230810.vcf.body.txt
head final_20230810.vcf.body.txt
awk '!/^#/ { print NF; exit }' final_20230810.vcf.body.txt
awk 'BEGIN{OFS="\t"} {temp=$4; $4=$5; $5=temp; print}' final_20230810.vcf.body.txt > final_20230810.vcf.body.sorted.txt
head final_20230810.vcf.header.txt
nano final_20230810.vcf.body.sorted.txt
#change alt ref names
cat final_20230810.vcf.header.txt final_20230810.vcf.body.sorted.txt > final_20230810.fixed.vcf
head -n 20 final_20230810.fixed.vcf
```
# create environment with the three datasets
module load R/4.1.1
module load gcc/10.2.0
```R

freq <- read.csv("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/grenedalf-freq/vcf-maf0.01-misfrac0.1-maxcov95/grenedalfst.log/freq-finalvcf/frequencysample-frequencies.csv", header=T)

sites <- freq[,c(1,2)]
afmat <- freq[,c(-1:-4)]
afmat <- as.matrix(afmat)
head(afmat)
samps <- read.table("samples.txt", header = TRUE, sep = "\t")
objects <- c("sites", "afmat", "samps")
save(list = "objects", file = "HAFs.Rdata")
```
```bash
# everything has been transfered to  ~/work/seasonal_adaptation/analysis/GLM
# I will now run the GLM from there

#extract depth information from mpileup file
module load bcftools/1.10.2
bcftools depth -a -o depth.txt /gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/variants/all.mpileup-excluding20072011


Rscript calc_GLM.R /gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/GLM/HAFs.Rdata --effectiveCov 200 -o /gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/output_GLM


#errors to fix
either --readDepth or --effectiveCov must be supplied.
***EXITING***

loading HAFs
Error in `$<-.data.frame`(`*tmp*`, cage, value = integer(0)) :
  replacement has 0 rows, data has 14
Calls: $<- -> $<-.data.frame
Execution halted
(base) [smomw573@nes


```

# re-calculate allele freq for the new vcf file

```bash 

library(tidyr)

module load htslib/1.10.2  bzip2/1.0.8  

GRENEDALF=/gxfs_home/geomar/smomw573/software/grenedalf/bin/grenedalf
VCF=/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/variants/final_20230810.vcf.gz
GENOME=/gxfs_work1/geomar/smomw573/seasonal_adaptation/genome/Eaffinis.Baltic.PseudoRef.Mar22.fasta \
OUTPUT=/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/GLM

$GRENEDALF frequency --vcf-path $VCF  --write-sample-alt-freq  --separator-char tab --file-suffix alt-new-vcf --out-dir $OUTPUT
$GRENEDALF frequency --vcf-path $VCF  --write-sample-coverage --separator-char tab --file-suffix cov-new-vcf --out-dir $OUTPUT

#output freq-new-vcf.csv - not a csv

#create file in the rigth format in R
srun --pty --x11 --nodes=1 --cpus-per-task=4 --mem=20000 --time=07:00:00 /bin/bash

module load R/4.1.1
module load gcc/10.2.0
R
```

# GLM


# run Rscript

```R
library(tidyr)
data <- read.table("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/GLM/frequencyalt-new-vcf.txt", header = TRUE, sep = "\t")


#fixing data#

# Remove ".FREQ" from the column names
# Identify the columns containing ".COV" in their names
cov_columns <- grep("EA_\\d+_T\\d+\\.FREQ", names(data))
new_colnames <- colnames(data)
new_colnames[cov_columns] <- gsub("\\.FREQ", "", new_colnames[cov_columns])
colnames(data) <- new_colnames

# Create a data frame with only the "CHROM," "POS," "REF," "ALT," and renamed ".COV" columns
data_FREQ <- data[, c(5:18)]
write.table(data_FREQ, file = "freq_new_vcf.txt", quote = FALSE, sep = "\t")

#cov data
data <- read.table("frequencycov-new-vcf.txt", header = TRUE, sep = "\t")

# Remove ".COV" from the column names
# Identify the columns containing ".COV" in their names
cov_columns <- grep("EA_\\d+_T\\d+\\.COV", names(data))
new_colnames <- colnames(data)
new_colnames[cov_columns] <- gsub("\\.COV", "", new_colnames[cov_columns])
colnames(data) <- new_colnames

# Create a data frame with only the "CHROM," "POS," "REF," "ALT," and renamed ".COV" columns
data_FREQ <- data[, c(5:18)]

write.table(data_FREQ, file = "cov_new_vcf.txt", quote = FALSE, sep = "\t")
```
# Bash GLM

```bash 
sbatch run-glm.sh 

#!/bin/bash
#SBATCH -D .
#SBATCH --mail-type=END
#SBATCH --mail-user=jnascimento@geomar.de
#SBATCH --partition=cluster
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=20G
#SBATCH --time=48:00:00
#SBATCH --job-name=glm
#SBATCH --output=glm.out
#SBATCH --error=glm.err

module load R/4.1.1
module load gcc/10.2.0

Rscript glm-r.R

```

```R
####GLM#####

####GLM script trial 1####

cov <- read.table("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/GLM/cov_new_vcf.txt", header = TRUE, sep = "\t")
freq <- read.table("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/GLM/freq_new_vcf.txt", header = TRUE, sep = "\t")
pos <- read.table("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/GLM/chrom-poly-position-seasonal-data.txt", header = FALSE, sep = "\t")
popinfo <- read.table("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/GLM/popinfo.txt", header = TRUE, sep = "\t")


#subset populations

popnames <- grep("EA_2009_T1|EA_2009_T2|EA_2011_T1|EA_2011_T2|EA_2015_T1|EA_2015_T4|EA_2022_T1|EA_2022_T4", popinfo$pop)

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
#in case we want unordered factors
#popinfo_glm$Y <- factor(popinfo_glm$Y, ordered = FALSE)
#GLM loop#
# Open the file for writing (creates a new file if it doesn't exist)
#fileout <- file("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/GLM/seasonal-glm-output.txt", "w")

# Open the file for writing (creates a new file if it doesn't exist)
fileout <- file("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/GLM/seasonal-glm-output-time-factor-unordered.txt", "w")
fileout2 <- file("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/GLM/seasonal-glm-output-year-factor-unordered.txt", "w")
# Loop through each line in freq_matrix
for (i in 1:nrow(freq_matrix)) {
  # Fit the GLM model for the current line
  out <- summary(glm(freq_matrix[i, ] ~ popinfo_glm$time + popinfo_glm$Y, family = binomial, weights = dp[i, ]))
  #glm_model <- summary(glm(freq_matrix[1, ] ~ popinfo_popnames$time, family = binomial, weights = dp[1, ]))
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


#if I would like to see the global effect of year that would be the code
#glm <- glm(line ~ popinfo_glm$time + popinfo_glm$Y, family = binomial, weights = dp[i, ]) 
#drop1(glm, test = "F")

# Close the file
close(fileout)

#with year as factor
fileout <- file("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/GLM/seasonal-glm-output-time-factor.txt", "w")
fileout2 <- file("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/GLM/seasonal-glm-output-year-factor.txt", "w")
# Loop through each line in freq_matrix
for (i in 1:nrow(freq_matrix)) {
  line <- freq_matrix[i, ]
  
  # Fit the GLM model for the current line
  out <- summary(glm(line ~ popinfo_glm$time + popinfo_glm$Y, family = binomial, weights = dp[i, ]))
  #glm_model <- summary(glm(freq_matrix[1, ] ~ popinfo_popnames$time, family = binomial, weights = dp[1, ]))
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

#with year as unordered factor

#with year as factor
fileout <- file("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/GLM/seasonal-glm-output-time-factor-un.txt", "w")
fileout2 <- file("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/GLM/seasonal-glm-output-year-factor-un.txt", "w")
# Loop through each line in freq_matrix
for (i in 1:nrow(freq_matrix)) {
  line <- freq_matrix[i, ]
  
  # Fit the GLM model for the current line
  out <- summary(glm(line ~ popinfo_glm$time + popinfo_glm$Y, family = binomial, weights = dp[i, ]))
  #glm_model <- summary(glm(freq_matrix[1, ] ~ popinfo_popnames$time, family = binomial, weights = dp[1, ]))
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

# Close the file
close(fileout)

```
# fix output files to get qqplots

```b
#file 1
#########
awk '{print $3}' seasonal-glm-output.txt | sort > seasonal-glm-output-*.txt
awk -F'\t' '$1 <= 0.05 { count++ } seasonal-glm-output-p-sorted.txt
#1641

awk '{print $1"_"$2}' chrom-poly-position-seasonal-data.txt > joint.txt 
awk '{print $3'} seasonal-glm-output-time.txt > p-values-time.txt 
paste chrom-poly-position-seasonal-data.txt joint.txt p-values-time.txt > full-table-ptime.txt

# add header
echo "chrom pos chrompos p" | cat - full-table-ptime.txt > temp && mv temp full-table-ptime.txt
echo "chrom  pos chrompos p" | cat - full-table-pyear.txt > temp && mv temp full-table-pyear.txt

awk -F'\t' '$4 <= 0.05 { count++ } END { print count }' full-table-ptime.txt
#1684
awk -F'\t' '$4 <= 0.05 { count++ } END { print count }' full-table-pyear.txt
#3730

awk -F'\t' '$4 <= 0.05 { print }' full-table-ptime.txt > significant-snps-time.txt
awk -F'\t' '$4 <= 0.05 { print }' full-table-pyear.txt > significant-snps-year.txt

#####


awk '{print $3'} seasonal-glm-output-time-factor.txt > p-values-time.txt 
paste chrom-poly-position-seasonal-data.txt joint.txt p-values-time.txt > full-table-ptime-factor.txt

awk '{print $3'} seasonal-glm-output-year-factor.txt > p-values-season.txt 
paste chrom-poly-position-seasonal-data.txt joint.txt p-values-year.txt > full-table-pyear-factor.txt

# add header
echo "chrom pos chrompos p" | cat - full-table-ptime-factor.txt > temp && mv temp full-table-ptime-factor.txt
echo "chrom  pos chrompos p" | cat - full-table-pyear-factor.txt > temp && mv temp full-table-pyear-factor.txt

awk -F'\t' '$4 <= 0.05 { count++ } END { print count }' full-table-ptime-factor.txt
#1684
awk -F'\t' '$4 <= 0.05 { count++ } END { print count }' full-table-pyear-factor.txt
#3730


#########
#file2

awk '{print $3}' seasonal-glm-output-time-factor-unordered-noweights.txt > time-noweights.txt 
awk -F'\t' '$1 <= 0.05 { count++ } END { print count }' time-noweights.txt
#1
paste chrom-poly-position-seasonal-data.txt joint.txt time-noweights.txt> full-table-time-noweights.txt
# add header
echo "chrom pos chrompos p" | cat - full-table-time-noweights.txt > temp && mv temp full-table-time-noweights.txt
awk -F'\t' '$4 <= 0.05 { count++ } END { print count }' full-table-time-noweights.txt
#1

awk -F'\t' '$4 <= 0.05 { print }' full-table-ptime.txt > significant-snps-time.txt
awk -F'\t' '$4 <= 0.05 { print }' full-table-pyear.txt > significant-snps-year.txt


#########
#file3


awk '{print $3}' seasonal-glm-output-time-factor-unordered.txt > time-unordered.txt 
awk -F'\t' '$1 <= 0.05 { count++ } END { print count }' time-unordered.txt
#1648
paste chrom-poly-position-seasonal-data.txt joint.txt time-unordered.txt > full-table-time-unordered.txt

# add header
echo "chrom pos chrompos p" | cat - full-table-time-unordered.txt > temp && mv temp full-table-time-unordered.txt

awk -F'\t' '$4 <= 0.05 { print }' full-table-time-unordered.txt> significant-full-table-time-unordered.txt

#########
#file4


awk '{print $3}' seasonal-glm-output-time-uncorrected-coverage.txt > time-covnonc.txt 
awk -F'\t' '$1 <= 0.05 { count++ } END { print count }' time-covnonc.txt 
#67292
paste chrom-poly-position-seasonal-data.txt joint.txt time-covnonc.txt > full-table-time-covnonc.txt




```
module load R/4.1.1
module load gcc/10.2.0
```R
#qqplots

#time as factor

par(mfrow = c(2, 2))

time <- read.table("p-values-time.txt", header = FALSE)

n <- length(time$V1)
expected_p_values <- seq(0, 1, length.out = n)
# Create a QQ plot
qqplot(-log10(expected_p_values), -log10(time$V1), xlab = "Expected -log10(p)", ylab = "Observed -log10(p)", main = "ordered factor (OF)")
# Add a reference line for a perfect uniform distribution
abline(0, 1, col = "red")
#year <- read.table("p-values-year.txt", header = FALSE)


#file2
#time unordered factor
time2 <- read.table("time-unordered.txt", header = FALSE)
n <- length(time2$V1)
expected_p_values <- seq(0, 1, length.out = n)
# Create a QQ plot
qqplot(-log10(expected_p_values), -log10(time2$V1), xlab = "Expected -log10(p)", ylab = "Observed -log10(p)", main = "unordered factor (UF)")
# Add a reference line for a perfect uniform distribution
abline(0, 1, col = "red")


#file3
#time no weights unordered factor
time3 <- read.table("time-noweights.txt", header = FALSE)
n <- length(time3$V1)
expected_p_values <- seq(0, 1, length.out = n)
# Create a QQ plot
qqplot(-log10(expected_p_values), -log10(time3$V1), xlab = "Expected -log10(p)", ylab = "Observed -log10(p)", main = "UF no weights")
# Add a reference line for a perfect uniform distribution
abline(0, 1, col = "red")

#file4
#time no weights unordered factor with uncorrected cov
time4 <- read.table("time-covnonc.txt", header = FALSE)
n <- length(time4$V1)
expected_p_values <- seq(0, 1, length.out = n)
# Create a QQ plot
qqplot(-log10(expected_p_values), -log10(time4$V1), xlab = "Expected -log10(p)", ylab = "Observed -log10(p)", main = "UF uncorrected cov")
# Add a reference line for a perfect uniform distribution
abline(0, 1, col = "red")
```
## 05-10-2023 -> based on the qqplot output, we decided to check what the variables that have a significant p value are doing.
## first I will check the changes in AF of these snps in each population, and second I will check if the ones who have a significant p value actually have a high coverage (mean cov from the matrix or pull it from grenedalf)

# 13-10-2023 -> okay, no we decided we want to have year as an unordered factor, and I ran the glm model with this, the output was the file seasonal-glm-output-time-factor-unordered.txt with 1202736 SNPs (what is what we expected) and the same file with p-values for years seasonal-glm-output-year-factor-unordered.txt
# now I will create the final file so we can check which SNPs are doing what

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
# now I will plot in R the the output
module load R/4.1.1
module load gcc/10.2.0
```R
library(tidyr)
library(ggplot2)
library(matrixStats)

df<- read.table(file="cov-sig-snps-txt", header= F)
df <- df[-3]
df$chrompos <- paste(df[,V1], df[,V2])
df$chrompos <- paste(df$V1, df$V2)
colnames(df) <- c("chrom","pos",
                  "2009t1","2009t2","2009t3","2009t4",
                  "2011t1","2011t2",
                  "2015t1","2015t2","2015t3","2015t4",
                  "2022t1","2022t2","2022t3","2022t4",
                  "chrompos")

df <- df[, -c(1, 2)]
df_long <- df %>%
  gather(key = "sample", value = "coverage", -chrompos)

ggplot(df_long, aes(x = chrompos, y = coverage, color = sample)) +
  geom_point(size = 3) +
  labs(title = "coverage", x = "chromosome", y = "coverage") +
  scale_x_discrete(labels = NULL) +
  theme_minimal()
ggsave("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/figures/coverage-individual-snps.pdf", width = 8, height = 6, units = "in")

#how many SNPs have the average coverage (500x)

sum(df_long$coverage <= 500)
 #21528 
 sum(df_long$coverage >= 500)                                                  
 #2055


df$mean_cov <- rowMeans(df[, -15])
df$se <- rowSds(as.matrix(df[, -15])) / sqrt(ncol(df) - 1)

ggplot(df, aes(x = chrompos, y = mean_cov)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_cov - se, ymax = mean_cov + se), width = 0.2) +
  labs(title = "coverage", x = "chromosome", y = "coverage") +
  scale_x_discrete(labels = NULL) +
  theme_minimal()

  #not yet

df <- separate(df, chrompos, into = c("chrom", "pos"), sep = " ", remove=FALSE)
library(dplyr)
distinct_count <- df %>% 
  count(chrom)

                      chrom   n
#1    Scz6wRH_101;HRSCAF=195   2
#2      Scz6wRH_12;HRSCAF=69   1
#3  Scz6wRH_1478;HRSCAF=1610   5
#4  Scz6wRH_1684;HRSCAF=1838   1
#5  Scz6wRH_1685;HRSCAF=1840 457
#6  Scz6wRH_1688;HRSCAF=1872   2
#7  Scz6wRH_1691;HRSCAF=1912   1
#8  Scz6wRH_1693;HRSCAF=1917 299
#9    Scz6wRH_185;HRSCAF=284   1
#10     Scz6wRH_19;HRSCAF=90   3
#11   Scz6wRH_234;HRSCAF=334   1
#12     Scz6wRH_23;HRSCAF=98 436
#13   Scz6wRH_267;HRSCAF=368   1
#14   Scz6wRH_277;HRSCAF=379   2
#15      Scz6wRH_2;HRSCAF=24 447
#16   Scz6wRH_306;HRSCAF=408   1
#17   Scz6wRH_345;HRSCAF=449   2
#18   Scz6wRH_357;HRSCAF=465   1
#19    Scz6wRH_37;HRSCAF=122   2
#20   Scz6wRH_407;HRSCAF=518   2
#21   Scz6wRH_430;HRSCAF=541   2
#22    Scz6wRH_76;HRSCAF=166   3
#23    Scz6wRH_85;HRSCAF=179   7
#24      Scz6wRH_8;HRSCAF=60   3
#25    Scz6wRH_98;HRSCAF=192   2

#plot SNPs by chromossome
ggplot(distinct_count, aes(x = reorder(chrom, -n), y = n)) +
  geom_bar(stat = "identity") +
  labs(title = "SNP counts", x = "Chromosome", y = "Count") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/figures/snps-by-choromosome-glm.pdf", width = 10, height = 6, units = "in")

df$nchrom <- as.integer(factor(df$chrom))
ggplot(df, aes(x = factor(chrom), y = mean_cov, color=chrom)) +
  geom_boxplot(width = 0.5, show.legend = FALSE) +  # Adjust width as needed
  labs(title = "Mean Coverage", x = "Chromosome", y = "Coverage") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/figures/mean-coverage-snps-chrom.pdf", width = 8, height = 6, units = "in")


#now I will find out which SNPs were deviating the norm in the QQplot, and see how they behave

#file2
#time unordered factor
pvalues <- read.table("full-table-final-glm-time.txt", header = TRUE)
n <- length(pvalues$p)
expected_p_values <- seq(0, 1, length.out = n)
# Create a QQ plot
qqplot(-log10(expected_p_values), -log10(pvalues$p), xlab = "Expected -log10(p)", ylab = "Observed -log10(p)", main = "unordered factor (UF)")

abline(0, 1, col = "red")

# Create the QQ plot
qqplot(-log10(expected_p_values), -log10(pvalues$p), xlab = "Expected -log10(p)", ylab = "Observed -log10(p)", main = "unordered factor (UF)")

# Identify SNPs with expected > 4.5 and observed > 5
selected_snps <- pvalues$chrompos[(-log10(pvalues$p) > 5)]
text(-log10(expected_p_values[pvalues$chrompos %in% selected_snps]), -log10(pvalues$p[pvalues$chrompos %in% selected_snps]), labels = selected_snps, pos = 1, cex = 0.7)
# Assuming you have logged values
logged_expected_p_values <- -log10(expected_p_values)
logged_p_values <- -log10(pvalues$p)

# Create a data frame with logged values
logged_data <- data.frame(
  Logged_Expected = logged_expected_p_values,
  Logged_Observed = logged_p_values
)

#okay now I want to plot cov in relation with pvalues
#########

covsig <- read.table("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/GLM/cov-sig-snps.txt", header = FALSE)
covsig$chrompos <- paste(covsig$V1, covsig$V2, sep="_")
covsig <- c[, -c(2, 3)]
pvalues <- read.table("/gxfs_work1/geomar/smomw573/seasonal_adaptation/analysis/GLM/significant-snps-time.txt", header = TRUE)
pvalues <- pvalues[, -c(1,2)]
colnames(pvalues) <- c("chrompos", "p")
merged_data <- merge(pvalues, covsig, by = "chrompos")
merged_data <- merged_data[,-c(3,4,5)]
colnames(merged_data) <- c("chrompos", "p", "2009t1", "2009t2", "2009t3", "2009t4", "2011t1","2011t2","2015t1","2015t2","2015t3","2015t4","2022t1","2022t2","2022t3","2022t4")
merged_data$mean_cov <- rowMeans(merged_data[, 3:16], na.rm = TRUE)

merged_data$log <- -log10(merged_data$p)


ggplot(merged_data, aes(x = log, y = mean_cov)) +
  geom_point(size = 3, show.legend = FALSE) +  # Adjust size as needed
  labs(title = "Mean Coverage", x = "p value (-log10)", y = "Mean Coverage") +
  geom_smooth(method='lm')
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/figures/mean-cov-p-value.pdf", width = 8, height = 6, units = "in")

model <- lm(mean_cov ~ log, data = merged_data)
summary(model)
```

# now I will check what the allele frequencies are doing!
# 16.10.2023
# will do everything in R

```R
library(tidyr)
library(ggplot2)
library(dplyr)

chrom <- read.table(file="/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/GLM/chroms.txt", header= T)
freq <- read.table(file="/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/GLM/freq_new_vcf.txt", header= T)
sig <- read.table(file="/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/GLM/significant-snps-time.txt", header= F)

af <- freq[, c("EA_2009_T1","EA_2009_T4" ,"EA_2011_T1", "EA_2011_T2","EA_2015_T1",
 "EA_2015_T4", "EA_2022_T1", "EA_2022_T4")]
af2 <- cbind(chrom, af)
af2 <- af2 %>%
  inner_join(sig, by = c("chrom" = "V3")) %>%
  select(all_of(names(af2)))

transaf<- gather(af2, key = "Sample", value = "af", -chrom)


# new columns
transaf <- transaf %>%
  mutate(time_point = ifelse(grepl("T1", Sample), "start", "end"))
transaf <- transaf %>%
  mutate(year = as.integer(sub("EA_(\\d{4})_.*", "\\1", Sample)))

transaf <- transaf %>%
  group_by(chrom, year) %>%
  mutate(change_af = af - lag(af, default = first(af)))
#this graph gives no information at all.

ggplot(data = transaf, aes(x = interaction(year, time_point, lex.order = TRUE), y = af, color = year)) +
  geom_point(position = position_dodge(width = 0.2), size = 3) +
  geom_line(aes(group = chrom), position = position_dodge(width = 0.2), linewidth = 1) 
  labs(x = "Year & Time Point", y = "AF", color = "Year") +
  ggtitle("Dot Plot of AF with Connecting Lines")

# first I will check per year what is going on 
#lets try with 2009


af2009 <- af2 %>%
  mutate(af_difference = EA_2009_T1 - EA_2009_T4) %>%
  select(chrom, af_difference) %>%
  arrange(af_difference)

ggplot(data = af2009, aes(x = af_difference, y = 1)) +
  geom_point(size = 1) +
  labs(x = "af_difference", y = "Count") +
  ggtitle("Scatter Plot of af_difference")

summary_counts <- af2009 %>%
  summarize(positive_count = sum(af_difference > 0),
            negative_count = sum(af_difference < 0))
summary_counts$positive_count
#1019
summary_counts$negative_count
#656

positive2009<- af2009 %>%
  filter(af_difference > 0)

negative2009<- af2009 %>%
  filter(af_difference < 0)


#so here I go back to af2 and fetch the chroms that were decreasing from start to end in 2009
filtered_positive_2009 <- af2 %>%
  filter(chrom %in% positive2009$chrom)

#trial crazy plot

filtered_af2 <- filtered_af2 %>%
  arrange(EA_2009_T4)


#first ill check if the snps that are increasing in one year behave the same in others (visually, tomorrow ill do it manually with the overlap)
#reorder chroms,
filtered_positive_2009$chrom <- factor(filtered_positive_2009 $chrom, levels = filtered_positive_2009 $chrom[order(filtered_positive_2009 $EA_2009_T4)])
#again....
# new columns

transafp2009<- gather(filtered_positive_2009, key = "Sample", value = "af", -chrom)
transafp2009 <- transafp2009 %>%
  mutate(time_point = ifelse(grepl("T1", Sample), "start", "end"))
transafp2009 <- transafp2009 %>%
  mutate(year = as.integer(sub("EA_(\\d{4})_.*", "\\1", Sample)))


#this graph gives no information at all.

ggplot(data = transafp2009, aes(x = interaction(year, time_point, lex.order = TRUE), y = af, color = year)) +
  geom_point(position = position_dodge(width = 0.2), size = 1) +
  geom_line(aes(group = chrom), position = position_dodge(width = 0.2), linewidth = 1) 
  labs(x = "Year & Time Point", y = "AF", color = "Year") +
  ggtitle("Dot Plot of AF with Connecting Lines")


# Create a plot with the specified order, light gray lines, and thin points
transafp2009$int <- interaction(transafp2009$year, transafp2009$time_point)
# Check levels
levels(transafp2009$int)
levels <- c("2009.start", "2009.end", "2011.start", "2011.end", "2015.start", "2015.end", "2022.start", "2022.end")

# Reorder the levels of the factor
transafp2009$int <- factor(transafp2009$int, levels = levels)

# Reverse the order of "time_point" factor levels
#transafp2009$time_point <- factor(transafp2009$time_point, levels = (time_point_order))
transafp2009$year <- factor(transafp2009$year, levels = (year_order))
# Create a plot
ggplot(data = transafp2009, aes(x = int, y = af, color = year)) +
  geom_line(aes(group = chrom), color = "lightgray", size = 0.5) +
  labs(x = "Year & Time Point", y = "AF", color = "Year") +
  ggtitle("Dot Plot of AF with Connecting Lines") +
  theme_bw()
 # scale_x_discrete(limits = int, labels = int, sep = " "))


#same for negative
filtered_positive_2009$chrom <- factor(filtered_positive_2009 $chrom, levels = filtered_positive_2009 $chrom[order(filtered_positive_2009 $EA_2009_T4)])
#again....
# new columns
#so here I go back to af2 and fetch the chroms that were decreasing from start to end in 2009
filtered_negative_2009 <- af2 %>%
  filter(chrom %in% negative2009$chrom)

filtered_negative_2009 <- filtered_negative_2009 %>%
  arrange(EA_2009_T4)


#first ill check if the snps that are increasing in one year behave the same in others (visually, tomorrow ill do it manually with the overlap)
#reorder chroms,
transafpneg2009<- gather(filtered_negative_2009, key = "Sample", value = "af", -chrom)
transafpneg2009 <- transafpneg2009 %>%
  mutate(time_point = ifelse(grepl("T1", Sample), "start", "end"))
transafpneg2009 <- transafpneg2009 %>%
  mutate(year = as.integer(sub("EA_(\\d{4})_.*", "\\1", Sample)))

# Create a plot with the specified order, light gray lines, and thin points
transafpneg2009$int <- interaction(transafpneg2009$year, transafpneg2009$time_point)
# Check levels
levels(transafpneg2009$int)
levels <- c("2009.start", "2009.end", "2011.start", "2011.end", "2015.start", "2015.end", "2022.start", "2022.end")
# Reorder the levels of the factor
transafpneg2009$int <- factor(transafpneg2009$int, levels = levels)

# Create a plot
ggplot(data = transafpneg2009, aes(x = int, y = af, color = year)) +
  geom_line(aes(group = chrom), color = "lightgray", size = 0.5) +
  labs(x = "Year & Time Point", y = "AF", color = "Year") +
  ggtitle("Dot Plot of AF with Connecting Lines") +
  theme_bw()
```

# okay, now I will check the top 10% and 1% of the snps in terms of pvalues, what are they doing
```R

library(tidyr)
library(ggplot2)
library(dplyr)

chrom <- read.table(file="/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/GLM/chroms.txt", header= T)
freq <- read.table(file="/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/GLM/freq_new_vcf.txt", header= T)
sig <- read.table(file="/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/GLM/significant-snps-time.txt", header= F)

af <- freq[, c("EA_2009_T1","EA_2009_T4" ,"EA_2011_T1", "EA_2011_T2","EA_2015_T1",
 "EA_2015_T4", "EA_2022_T1", "EA_2022_T4")]
af2 <- cbind(chrom, af)
af2 <- af2 %>%
  inner_join(sig, by = c("chrom" = "V3")) %>%
  select(names(af2), V4)

#what are the average p values
mean(af2$V4)
# 0.02527572
quantile(af2$V4, 0.9)
#       90%
#0.04639943
threshold <- quantile(af2$V4, 0.01)
lower_10_percent_data <- af2[af2$V4 <= threshold, ]

# new columns

low10 <- lower_10_percent_data[, -which(names(lower_10_percent_data) == "V4")]
low10 <- gather(low10, key = "Sample", value = "af", -chrom)
low10 <- low10 %>%
  mutate(time_point = ifelse(grepl("T1", Sample), "start", "end"))
low10 <- low10 %>%
  mutate(year = as.integer(sub("EA_(\\d{4})_.*", "\\1", Sample)))
# Create a plot with the specified order, light gray lines, and thin points
low10$int <- interaction(low10$year, low10$time_point)
# Check levels
levels(low10$int)
levels <- c("2009.start", "2009.end", "2011.start", "2011.end", "2015.start", "2015.end", "2022.start", "2022.end")
# Reorder the levels of the factor
low10$int <- factor(low10$int, levels = levels)


ggplot(data = low10, aes(x = int, y = af, color = year)) +
  geom_point(position = position_dodge(width = 0.2), linewidth = 1) +
  geom_line(aes(group = chrom),  linewidth = 0.5) 
  labs(x = "years", y = "AF", color = "Year") +
  ggtitle("top 10% p values") +
  theme_bw()

ggplot(data = low10, aes(x = int, y = af, color = year)) +
  geom_line(aes(group = chrom), color = "lightgray", linewidth = 0.5) +
  labs(x = "years", y = "AF", color = "Year") +
  ggtitle("top 10% p values") +
  theme_bw()
```
# 18-10-2023

we went bak to the PCA data - the original one had samples 2009 T1 and 2015 T3 grouped very far apart... 
for this reason, I will rerun the analysis using 2009 t2 instead of t1 (to see if there are any changes..)
wondering if I should use t2 for all other years as well... but dont think this one week should do that much although it did a lot for 2011...

arghhh
I'll try both..

# rerun the GLM analysis, first using 2009 t2 as input, second using all years apart from 2011 starting at t2.

module load R/4.1.1
module load gcc/10.2.0
```R
cov <- read.table("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/GLM/cov_new_vcf.txt", header = TRUE, sep = "\t")
freq <- read.table("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/GLM/freq_new_vcf.txt", header = TRUE, sep = "\t")
pos <- read.table("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/GLM/chrom-poly-position-seasonal-data.txt", header = FALSE, sep = "\t")
popinfo <- read.table("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/GLM/popinfo.txt", header = TRUE, sep = "\t")
#subset populations
popnames <- grep("EA_2009_T1|EA_2009_T4|EA_2011_T1|EA_2011_T2|EA_2015_T1|EA_2015_T4|EA_2022_T1|EA_2022_T4", popinfo$pop)

#subset 
cov_glm <- cov[,popnames]
freq_glm <- freq[,popnames]
popinfo_glm <- popinfo[popnames,]
popinfo_glm$time <- c("E","L")

#transform into matrix
freq_matrix <- as.matrix(freq_glm)
cov_matrix <- as.matrix(cov_glm)
dp <- (1/100 + 1/cov_matrix)^-1

popinfo_glm$Y <- as.factor(popinfo_glm$Y)
popinfo_glm$Y <- factor(popinfo_glm$Y, ordered = FALSE)

#with year as factor
fileout <- file("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/GLM/glm-181023-time.txt", "w")
fileout2 <- file("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/GLM/glm-181023-year.txt", "w")

# Loop through each line in freq_matrix
for (i in 1:nrow(freq_matrix)) {
  # Fit the GLM model for the current line
  out <- summary(glm(freq_matrix[i,] ~ popinfo_glm$time + popinfo_glm$Y, family = binomial, weights = dp[i, ]))
  #glm_model <- summary(glm(freq_matrix[1, ] ~ popinfo_popnames$time, family = binomial, weights = dp[1, ]))
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
# Close the file
close(fileout)
close(fileout2)


#same GLM but for all t1 t2 comparisons
#######################################
cov <- read.table("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/GLM/cov_new_vcf.txt", header = TRUE, sep = "\t")
freq <- read.table("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/GLM/freq_new_vcf.txt", header = TRUE, sep = "\t")
pos <- read.table("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/GLM/chrom-poly-position-seasonal-data.txt", header = FALSE, sep = "\t")
popinfo <- read.table("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/GLM/popinfo.txt", header = TRUE, sep = "\t")
popnames <- grep("EA_2009_T1|EA_2009_T2|EA_2011_T1|EA_2011_T2|EA_2015_T1|EA_2015_T2|EA_2022_T1|EA_2022_T2", popinfo$pop)

#subset 
cov_glm <- cov[,popnames]
freq_glm <- freq[,popnames]
popinfo_glm <- popinfo[popnames,]
popinfo_glm$time <- c("E","L")

#transform into matrix
freq_matrix <- as.matrix(freq_glm)
cov_matrix <- as.matrix(cov_glm)
dp <- (1/100 + 1/cov_matrix)^-1

popinfo_glm$Y <- as.factor(popinfo_glm$Y)
popinfo_glm$Y <- factor(popinfo_glm$Y, ordered = FALSE)

#with year as factor
fileout <- file("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/GLM/glm-181023-t2-time.txt", "w")
fileout2 <- file("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/GLM/glm-181023-t2-year.txt", "w")

# Loop through each line in freq_matrix
for (i in 1:nrow(freq_matrix)) {
  # Fit the GLM model for the current line
  out <- summary(glm(freq_matrix[i,] ~ popinfo_glm$time + popinfo_glm$Y, family = binomial, weights = dp[i, ]))
  #glm_model <- summary(glm(freq_matrix[1, ] ~ popinfo_popnames$time, family = binomial, weights = dp[1, ]))
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
close(fileout)
close(fileout2)


#same GLM but for all 2009t2 instead of t1 comparisons
#######################################
cov <- read.table("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/GLM/cov_new_vcf.txt", header = TRUE, sep = "\t")
freq <- read.table("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/GLM/freq_new_vcf.txt", header = TRUE, sep = "\t")
pos <- read.table("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/GLM/chrom-poly-position-seasonal-data.txt", header = FALSE, sep = "\t")
popinfo <- read.table("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/GLM/popinfo.txt", header = TRUE, sep = "\t")
popnames <- grep("EA_2009_T2|EA_2009_T4|EA_2011_T1|EA_2011_T2|EA_2015_T1|EA_2015_T4|EA_2022_T1|EA_2022_T4", popinfo$pop)

#subset 
cov_glm <- cov[,popnames]
freq_glm <- freq[,popnames]
popinfo_glm <- popinfo[popnames,]
popinfo_glm$time <- c("E","L")

#transform into matrix
freq_matrix <- as.matrix(freq_glm)
cov_matrix <- as.matrix(cov_glm)
dp <- (1/100 + 1/cov_matrix)^-1

popinfo_glm$Y <- as.factor(popinfo_glm$Y)
popinfo_glm$Y <- factor(popinfo_glm$Y, ordered = FALSE)

#with year as factor
fileout <- file("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/GLM/glm-181023-2009t2-time.txt", "w")
fileout2 <- file("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/GLM/glm-181023-2009t2-year.txt", "w")

# Loop through each line in freq_matrix
for (i in 1:nrow(freq_matrix)) {
  # Fit the GLM model for the current line
  out <- summary(glm(freq_matrix[i,] ~ popinfo_glm$time + popinfo_glm$Y, family = binomial, weights = dp[i, ]))
  #glm_model <- summary(glm(freq_matrix[1, ] ~ popinfo_popnames$time, family = binomial, weights = dp[1, ]))
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
close(fileout)
close(fileout2)
```

two output files that I need to look at
glm-181023-t2-time.txt
glm-181023-time.txt

```bash
awk '{print $3}' glm-181023-time.txt > glm-181023-time-p.txt
awk -F'\t' '$1 <= 0.05 { count++ } END { print count }' glm-181023-time-p.txt
#1924
paste chrom-poly-position-seasonal-data.txt joint.txt glm-181023-time-p.txt > full-table-glm-181023-time-p.txt
# add header
echo "chrom pos chrompos p" | cat - full-table-glm-181023-time-p.txt > temp && mv temp full-table-glm-181023-time-p.txt

#nano and fix header

awk -F'\t' '$4 <= 0.05 { print }' full-table-glm-181023-time-p.txt > significant-glm-181023-time-p.txt

awk 'NR==FNR{a[$1$2]=$0; next} ($1$2 in a){print a[$1$2]}' cov_with_chrom.txt significant-glm-181023-time-p.txt > cov-sig-snps-glm-181023-time-p.txt
```
# now I will plot in R the the output
module load R/4.1.1
module load gcc/10.2.0
```R
library(tidyr)
library(ggplot2)
library(matrixStats)

df<- read.table(file="/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/GLM/cov-sig-snps-glm-181023-time-p.txt", header= F)
df <- df[-3]
df$chrompos <- paste(df$V1, df$V2)
colnames(df) <- c("chrom","pos",
                  "2009t1","2009t2","2009t3","2009t4",
                  "2011t1","2011t2",
                  "2015t1","2015t2","2015t3","2015t4",
                  "2022t1","2022t2","2022t3","2022t4",
                  "chrompos")

df <- df[, -c(2)]
df_long <- df %>%
  gather(key = "sample", value = "coverage", -chrompos, -chrom)

ggplot(df_long, aes(x = chrom, y = coverage, color = sample)) +
  geom_point(position = position_dodge(width = 0.2),size = 3) +
  labs(title = "coverage", x = "chromosome", y = "coverage") +
  #scale_x_discrete(labels = df_long$chrom) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
ggsave("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/figures/coverage-individual-snps-fixed-glm.pdf", width = 8, height = 6, units = "in")

#how many SNPs have the average coverage (500x)
sum(df_long$coverage <= 500)
 # 24720 
 sum(df_long$coverage >= 500)                                                  
  #2222

df$mean_cov <- rowMeans(df[, -c(1,16)])
df$se <- rowSds(as.matrix(df[, -c(1,16)])) / sqrt(ncol(df) - 2)

ggplot(df, aes(x = chrompos, y = mean_cov)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_cov - se, ymax = mean_cov + se), width = 0.2) +
  labs(title = "coverage", x = "chromosome", y = "coverage") +
  scale_x_discrete(labels = NULL) +
  theme_minimal()

  #not yet

df <- separate(df, chrompos, into = c("chrom", "pos"), sep = " ", remove=FALSE)
library(dplyr)
distinct_count <- df %>% 
  count(chrom)

#                      chrom   n
#1    Scz6wRH_101;HRSCAF=195   1
#2    Scz6wRH_108;HRSCAF=204   1
#3  Scz6wRH_1478;HRSCAF=1610   8
#4  Scz6wRH_1684;HRSCAF=1838   1
#5  Scz6wRH_1685;HRSCAF=1840 516
#6  Scz6wRH_1688;HRSCAF=1872   2
#7  Scz6wRH_1693;HRSCAF=1917 354
#8    Scz6wRH_177;HRSCAF=276   1
#9    Scz6wRH_185;HRSCAF=284   1
#10     Scz6wRH_19;HRSCAF=90   4
#11     Scz6wRH_23;HRSCAF=98 483
#12   Scz6wRH_267;HRSCAF=368   1
#13   Scz6wRH_276;HRSCAF=378   1
#14      Scz6wRH_2;HRSCAF=24 526
#15   Scz6wRH_306;HRSCAF=408   1
#16   Scz6wRH_345;HRSCAF=449   2
#17   Scz6wRH_360;HRSCAF=468   1
#18    Scz6wRH_37;HRSCAF=122   3
#19   Scz6wRH_430;HRSCAF=541   2
#20   Scz6wRH_508;HRSCAF=620   1
#21   Scz6wRH_682;HRSCAF=800   1
#22    Scz6wRH_76;HRSCAF=166   1
#23    Scz6wRH_85;HRSCAF=179  10
#24    Scz6wRH_98;HRSCAF=192   1
#25      Scz6wRH_9;HRSCAF=61   1

#plot SNPs by chromossome
ggplot(distinct_count, aes(x = reorder(chrom, -n), y = n)) +
  geom_bar(stat = "identity") +
  labs(title = "SNP counts", x = "Chromosome", y = "Count") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/figures/snps-by-choromosome-glm-fixed.pdf", width = 10, height = 6, units = "in")

df$nchrom <- as.integer(factor(df$chrom))
ggplot(df, aes(x = factor(chrom), y = mean_cov, color=chrom)) +
  geom_boxplot(width = 0.5, show.legend = FALSE) +  # Adjust width as needed
  labs(title = "Mean Coverage", x = "Chromosome", y = "Coverage") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/figures/mean-coverage-snps-chrom-fixed.pdf", width = 8, height = 6, units = "in")


#now I will find out which SNPs were deviating the norm in the QQplot, and see how they behave

#file2
#time unordered factor
pvalues <- read.table("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/GLM/full-table-glm-181023-time-p.txt", header = TRUE)
n <- length(pvalues$p)
expected_p_values <- seq(0, 1, length.out = n)
# Create a QQ plot
qqplot(-log10(expected_p_values), -log10(pvalues$p), xlab = "Expected -log10(p)", ylab = "Observed -log10(p)", main = "unordered factor (UF)")

abline(0, 1, col = "red")

#okay now I want to plot cov in relation with pvalues
#########

covsig <- read.table(file="/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/GLM/cov-sig-snps-glm-181023-time-p.txt", header= F)
covsig$chrompos <- paste(covsig$V1, covsig$V2, sep="_")
covsig <- covsig[, -c(1,2)]
pvalues <- read.table("/gxfs_work1/geomar/smomw573/seasonal_adaptation/analysis/GLM/significant-glm-181023-time-p.txt", header = FALSE)
pvalues <- pvalues[, -c(1,2)]
colnames(pvalues) <- c("chrompos", "p")
merged_data <- merge(pvalues, covsig, by = "chrompos")
merged_data <- merged_data[,-c(3)]
colnames(merged_data) <- c("chrompos", "p", "2009t1", "2009t2", "2009t3", "2009t4", "2011t1","2011t2","2015t1","2015t2","2015t3","2015t4","2022t1","2022t2","2022t3","2022t4")
merged_data$mean_cov <- rowMeans(merged_data[, 3:16], na.rm = TRUE)

merged_data$log <- -log10(merged_data$p)


ggplot(merged_data, aes(x = log, y = mean_cov)) +
  geom_point(size = 3, show.legend = FALSE) +  # Adjust size as needed
  labs(title = "Mean Coverage", x = "p value (-log10)", y = "Mean Coverage") +
  geom_smooth(method='lm')
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/figures/mean-cov-p-value-fixed.pdf", width = 8, height = 6, units = "in")

model <- lm(mean_cov ~ log, data = merged_data)
summary(model)


# now I will check what the allele frequencies are doing!
# 16.10.2023
# will do everything in R

chrom <- read.table(file="/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/GLM/chroms.txt", header= T)
freq <- read.table(file="/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/GLM/freq_new_vcf.txt", header= T)
sig <- read.table(file="/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/GLM/cov-sig-snps-glm-181023-time-p.txt", header= F) 

sig$chrompos <- paste(sig$V1, sig$V2, sep="_")
af <- freq[, c("EA_2009_T1","EA_2009_T4" ,"EA_2011_T1", "EA_2011_T2","EA_2015_T1",
 "EA_2015_T4", "EA_2022_T1", "EA_2022_T4")]
af2 <- cbind(chrom, af)
af2 <- af2 %>%
  inner_join(sig, by = c("chrom" = "chrompos")) %>%
  select(all_of(names(af2)))

transaf<- gather(af2, key = "Sample", value = "af", -chrom)


# new columns
transaf <- transaf %>%
  mutate(time_point = ifelse(grepl("T1", Sample), "start", "end"))
transaf <- transaf %>%
  mutate(year = as.integer(sub("EA_(\\d{4})_.*", "\\1", Sample)))

transaf <- transaf %>%
  group_by(chrom, year) %>%
  mutate(change_af = af - lag(af, default = first(af)))
#this graph gives no information at all.

ggplot(data = transaf, aes(x = interaction(year, time_point, lex.order = TRUE), y = af, color = year)) +
  geom_point(position = position_dodge(width = 0.2), size = 3) +
  geom_line(aes(group = chrom), position = position_dodge(width = 0.2), linewidth = 1) +
  labs(x = "Year & Time Point", y = "AF", color = "Year") +
  ggtitle("changes af t1-t4")

# first I will check per year what is going on 
#lets try with 2009


af2009 <- af2 %>%
  mutate(af_difference = EA_2009_T1 - EA_2009_T4) %>%
  select(chrom, af_difference) %>%
  arrange(af_difference)

ggplot(data = af2009, aes(x = af_difference, y = 1)) +
  geom_point(size = 1) +
  labs(x = "af_difference", y = "Count") +
  ggtitle("Scatter Plot of af_difference")

summary_counts <- af2009 %>%
  summarize(positive_count = sum(af_difference > 0),
            negative_count = sum(af_difference < 0))
summary_counts$positive_count
#1247
summary_counts$negative_count
#671

positive2009<- af2009 %>%
  filter(af_difference > 0)

negative2009<- af2009 %>%
  filter(af_difference < 0)


#so here I go back to af2 and fetch the chroms that were decreasing from start to end in 2009
filtered_positive_2009 <- af2 %>%
  filter(chrom %in% positive2009$chrom)

#trial crazy plot

filtered_af2 <- filtered_positive_2009 %>%
  arrange(EA_2009_T4)


#first ill check if the snps that are increasing in one year behave the same in others (visually, tomorrow ill do it manually with the overlap)
#reorder chroms,
filtered_positive_2009$chrom <- factor(filtered_positive_2009 $chrom, levels = filtered_positive_2009 $chrom[order(filtered_positive_2009 $EA_2009_T4)])
#again....
# new columns

transafp2009<- gather(filtered_positive_2009, key = "Sample", value = "af", -chrom)
transafp2009 <- transafp2009 %>%
  mutate(time_point = ifelse(grepl("T1", Sample), "start", "end"))
transafp2009 <- transafp2009 %>%
  mutate(year = as.integer(sub("EA_(\\d{4})_.*", "\\1", Sample)))


#this graph gives no information at all.

ggplot(data = transafp2009, aes(x = interaction(year, time_point, lex.order = TRUE), y = af, color = year)) +
  geom_point(position = position_dodge(width = 0.2), size = 1) +
  geom_line(aes(group = chrom), position = position_dodge(width = 0.2), linewidth = 1) 
  labs(x = "Year & Time Point", y = "AF", color = "Year") +
  ggtitle("Dot Plot of AF with Connecting Lines")


# Create a plot with the specified order, light gray lines, and thin points
transafp2009$int <- interaction(transafp2009$year, transafp2009$time_point)
# Check levels
levels(transafp2009$int)
levels <- c("2009.start", "2009.end", "2011.start", "2011.end", "2015.start", "2015.end", "2022.start", "2022.end")

# Reorder the levels of the factor
transafp2009$int <- factor(transafp2009$int, levels = levels)

ggplot(data = transafp2009, aes(x = int, y = af, color = year)) +
  geom_line(aes(group = chrom), color = "lightgray", linewidth = 0.5) +
  labs(x = "Year & Time Point", y = "AF", color = "Year") +
  ggtitle("") +
  theme_bw()
 # scale_x_discrete(limits = int, labels = int, sep = " "))


#same for negative
filtered_positive_2009$chrom <- factor(filtered_positive_2009 $chrom, levels = filtered_positive_2009 $chrom[order(filtered_positive_2009 $EA_2009_T4)])
#again....
# new columns
#so here I go back to af2 and fetch the chroms that were decreasing from start to end in 2009
filtered_negative_2009 <- af2 %>%
  filter(chrom %in% negative2009$chrom)
filtered_negative_2009 <- filtered_negative_2009[,-c(10)]
filtered_negative_2009 <- filtered_negative_2009 %>%
  arrange(EA_2009_T4)


#first ill check if the snps that are increasing in one year behave the same in others (visually, tomorrow ill do it manually with the overlap)
#reorder chroms,
transafpneg2009<- gather(filtered_negative_2009, key = "Sample", value = "af", -chrom)
transafpneg2009 <- transafpneg2009 %>%
  mutate(time_point = ifelse(grepl("T1", Sample), "start", "end"))
transafpneg2009 <- transafpneg2009 %>%
  mutate(year = as.integer(sub("EA_(\\d{4})_.*", "\\1", Sample)))

# Create a plot with the specified order, light gray lines, and thin points
transafpneg2009$int <- interaction(transafpneg2009$year, transafpneg2009$time_point)
# Check levels
levels(transafpneg2009$int)
levels <- c("2009.start", "2009.end", "2011.start", "2011.end", "2015.start", "2015.end", "2022.start", "2022.end")
# Reorder the levels of the factor
transafpneg2009$int <- factor(transafpneg2009$int, levels = levels)

# Create a plot
ggplot(data = transafpneg2009, aes(x = int, y = af, color = year)) +
  geom_line(aes(group = chrom), color = "lightgray", size = 0.5) +
  labs(x = "Year & Time Point", y = "AF", color = "Year") +
  ggtitle("") +
  theme_bw()


# okay, now I will check the top 10% and 1% of the snps in terms of pvalues, what are they doing


sig <- read.table(file="/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/GLM/significant-glm-181023-time-p.txt", header= F)

af <- freq[, c("EA_2009_T1","EA_2009_T4" ,"EA_2011_T1", "EA_2011_T2","EA_2015_T1",
 "EA_2015_T4", "EA_2022_T1", "EA_2022_T4")]
af2 <- cbind(chrom, af)
af2 <- af2 %>%
  inner_join(sig, by = c("chrom" = "V3")) %>%
  select(names(af2), V4)

#whath are the average p values
mean(af2$V4)
#  0.02283937
quantile(af2$V4, 0.9)
#       90%
# 0.04594068
threshold <- quantile(af2$V4, 0.01)
lower_10_percent_data <- af2[af2$V4 <= threshold, ]

# new columns
low10 <- lower_10_percent_data[, -which(names(lower_10_percent_data) == "V4")]
low10 <- gather(low10, key = "Sample", value = "af", -chrom)
low10 <- low10 %>%
  mutate(time_point = ifelse(grepl("T1", Sample), "start", "end"))
low10 <- low10 %>%
  mutate(year = as.integer(sub("EA_(\\d{4})_.*", "\\1", Sample)))
# Create a plot with the specified order, light gray lines, and thin points
low10$int <- interaction(low10$year, low10$time_point)
# Check levels
levels(low10$int)
levels <- c("2009.start", "2009.end", "2011.start", "2011.end", "2015.start", "2015.end", "2022.start", "2022.end")
# Reorder the levels of the factor
low10$int <- factor(low10$int, levels = levels)

ggplot(data = low10, aes(x = int, y = af, color = chrom)) +
  geom_line(aes(group = chrom), color = "lightgray", linewidth = 0.5) +
  labs(x = "years", y = "AF", color = "Year") +
  ggtitle("top 10% p values") +
  theme_bw()
```
now for second file (t2 glm output)

```bash
awk '{print $3}' glm-181023-t2-time.txt > glm-181023-t2-time-p.txt
awk -F'\t' '$1 <= 0.05 { count++ } END { print count }' glm-181023-t2-time-p.txt
#1830
paste chrom-poly-position-seasonal-data.txt joint.txt glm-181023-t2-time-p.txt > full-table-glm-181023-t2-time-p.txt
# add header
echo "chrom pos chrompos p" | cat - full-table-glm-181023-t2-time-p.txt > temp && mv temp full-table-glm-181023-t2-time-p.txt

#nano and fix header

awk -F'\t' '$4 <= 0.05 { print }' full-table-glm-181023-t2-time-p.txt > significant-glm-181023-t2-time-p.txt

awk 'NR==FNR{a[$1$2]=$0; next} ($1$2 in a){print a[$1$2]}' cov_with_chrom.txt significant-glm-181023-t2-time-p.txt > cov-sig-snps-glm-181023-t2-time-p.txt
```
# now I will plot in R the the output
#that is the final set for the GLM - do not get confused
module load R/4.1.1
module load gcc/10.2.0
```R
library(tidyr)
library(ggplot2)
library(matrixStats)

df<- read.table(file="/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/GLM/cov-sig-snps-glm-181023-t2-time-p.txt", header= F)
df <- df[-3]
df$chrompos <- paste(df$V1, df$V2)
colnames(df) <- c("chrom","pos",
                  "2009t1","2009t2","2009t3","2009t4",
                  "2011t1","2011t2",
                  "2015t1","2015t2","2015t3","2015t4",
                  "2022t1","2022t2","2022t3","2022t4",
                  "chrompos")

df <- df[, -c(2)]
df_long <- df %>%
  gather(key = "sample", value = "coverage", -chrompos, -chrom)

ggplot(df_long, aes(x = chrom, y = coverage, color = sample)) +
  geom_point(position = position_dodge(width = 0.2),size = 3) +
  labs(title = "coverage", x = "chromosome", y = "coverage") +
  #scale_x_discrete(labels = df_long$chrom) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
ggsave("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/figures/coverage-individual-snps-fixed-glmt2.pdf", width = 8, height = 6, units = "in")

#how many SNPs have the average coverage (500x)
sum(df_long$coverage <= 500)
 # 23586 
 sum(df_long$coverage >= 500)
#2041

df$mean_cov <- rowMeans(df[, -c(1,16)])
df$se <- rowSds(as.matrix(df[, -c(1,16)])) / sqrt(ncol(df) - 2)

ggplot(df, aes(x = chrompos, y = mean_cov)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_cov - se, ymax = mean_cov + se), width = 0.2) +
  labs(title = "coverage", x = "chromosome", y = "coverage") +
  scale_x_discrete(labels = NULL) +
  theme_minimal()
#not yet

df <- separate(df, chrompos, into = c("chrom", "pos"), sep = " ", remove=FALSE)
library(dplyr)
distinct_count <- df %>% 
  count(chrom)
#                      chrom   n
#1  Scz6wRH_1092;HRSCAF=1221   1
#2  Scz6wRH_1132;HRSCAF=1261   1
#3      Scz6wRH_12;HRSCAF=69   2
#4  Scz6wRH_1478;HRSCAF=1610   7
#5  Scz6wRH_1657;HRSCAF=1792   1
#6  Scz6wRH_1684;HRSCAF=1838   2
#7  Scz6wRH_1685;HRSCAF=1840 525
#8  Scz6wRH_1688;HRSCAF=1872   2
#9  Scz6wRH_1693;HRSCAF=1917 347
#10 Scz6wRH_1699;HRSCAF=1948   1
#11 Scz6wRH_1702;HRSCAF=1974   1
#12     Scz6wRH_19;HRSCAF=90   1
#13   Scz6wRH_201;HRSCAF=300   1
#14     Scz6wRH_20;HRSCAF=91   4
#15     Scz6wRH_23;HRSCAF=98 457
#16   Scz6wRH_267;HRSCAF=368   2
#17   Scz6wRH_277;HRSCAF=379   1
#18      Scz6wRH_2;HRSCAF=24 451
#19   Scz6wRH_306;HRSCAF=408   1
#20   Scz6wRH_345;HRSCAF=449   2
#21    Scz6wRH_38;HRSCAF=123   2
#22   Scz6wRH_407;HRSCAF=518   1
#23   Scz6wRH_430;HRSCAF=541   2
#24    Scz6wRH_43;HRSCAF=130   2
#25   Scz6wRH_528;HRSCAF=640   1
#26   Scz6wRH_667;HRSCAF=783   1
#27   Scz6wRH_691;HRSCAF=811   1
#28    Scz6wRH_76;HRSCAF=166   2
#29    Scz6wRH_85;HRSCAF=179   7

#these are the longest choromosomes
#1703:Scz6wRH_1693;HRSCAF=1917   100968129
#1704:Scz6wRH_23;HRSCAF=98       128296010
#1705:Scz6wRH_2;HRSCAF=24        140955377
#1706:Scz6wRH_1685;HRSCAF=1840   147976548


#plot SNPs by chromossome
ggplot(distinct_count, aes(x = reorder(chrom, -n), y = n)) +
  geom_bar(stat = "identity") +
  labs(title = "SNP counts", x = "Chromosome", y = "Count") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/figures/snps-by-choromosome-glm-fixed.pdf", width = 10, height = 6, units = "in")

df$nchrom <- as.integer(factor(df$chrom))
ggplot(df, aes(x = factor(chrom), y = mean_cov, color=chrom)) +
  geom_boxplot(width = 0.5, show.legend = FALSE) +  # Adjust width as needed
  labs(title = "Mean Coverage", x = "Chromosome", y = "Coverage") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/figures/mean-coverage-snps-chrom-fixed-t2.pdf", width = 8, height = 6, units = "in")


#now I will find out which SNPs were deviating the norm in the QQplot, and see how they behave

#time unordered factor
pvalues <- read.table("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/GLM/full-table-glm-181023-t2-time-p.txt", header = TRUE)
n <- length(pvalues$p)
expected_p_values <- seq(0, 1, length.out = n)
# Create a QQ plot
qqplot(-log10(expected_p_values), -log10(pvalues$p), xlab = "Expected -log10(p)", ylab = "Observed -log10(p)", main = "unordered factor (UF)")

abline(0, 1, col = "red")

#okay now I want to plot cov in relation with pvalues
#########

covsig <- read.table(file="/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/GLM/cov-sig-snps-glm-181023-t2-time-p.txt", header= F)
covsig$chrompos <- paste(covsig$V1, covsig$V2, sep="_")
covsig <- covsig[, -c(1,2,3)]
pvalues <- read.table("/gxfs_work1/geomar/smomw573/seasonal_adaptation/analysis/GLM/significant-glm-181023-t2-time-p.txt", header = FALSE)
pvalues <- pvalues[, -c(1,2)]
colnames(pvalues) <- c("chrompos", "p")
merged_data <- merge(pvalues, covsig, by = "chrompos")

colnames(merged_data) <- c("chrompos", "p", "2009t1", "2009t2", "2009t3", "2009t4", "2011t1","2011t2","2015t1","2015t2","2015t3","2015t4","2022t1","2022t2","2022t3","2022t4")
merged_data$mean_cov <- rowMeans(merged_data[, 3:16], na.rm = TRUE)

merged_data$log <- -log10(merged_data$p)


ggplot(merged_data, aes(x = log, y = mean_cov)) +
  geom_point(size = 3, show.legend = FALSE) +  # Adjust size as needed
  labs(title = "Mean Coverage", x = "p value (-log10)", y = "Mean Coverage") +
  geom_smooth(method='lm')
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/figures/mean-cov-p-value-fixed-t2.pdf", width = 8, height = 6, units = "in")

model <- lm(mean_cov ~ log, data = merged_data)
summary(model)


# now I will check what the allele frequencies are doing!

chrom <- read.table(file="/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/GLM/chroms.txt", header= T)
freq <- read.table(file="/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/GLM/freq_new_vcf.txt", header= T)
sig <- read.table(file="/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/GLM/cov-sig-snps-glm-181023-t2-time-p.txt", header= F) 

sig$chrompos <- paste(sig$V1, sig$V2, sep="_")
af <- freq[, c("EA_2009_T1","EA_2009_T2" ,"EA_2011_T1", "EA_2011_T2","EA_2015_T1",
 "EA_2015_T2", "EA_2022_T1", "EA_2022_T2")]
af2 <- cbind(chrom, af)
af2 <- af2 %>%
  inner_join(sig, by = c("chrom" = "chrompos")) %>%
  select(all_of(names(af2)))

transaf<- gather(af2, key = "Sample", value = "af", -chrom)


# new columns
transaf <- transaf %>%
  mutate(time_point = ifelse(grepl("T1", Sample), "start", "end"))
transaf <- transaf %>%
  mutate(year = as.integer(sub("EA_(\\d{4})_.*", "\\1", Sample)))

transaf <- transaf %>%
  group_by(chrom, year) %>%
  mutate(change_af = af - lag(af, default = first(af)))

transaf$int <- interaction(transaf$year, transaf$time_point)
# Check levels
levels(transafp2009$int)
levels <- c("2009.start", "2009.end", "2011.start", "2011.end", "2015.start", "2015.end", "2022.start", "2022.end")
# Reorder the levels of the factor
transaf$int <- factor(transaf$int, levels = levels)
transaf2 <- transaf %>%
  arrange(af)


#this graph gives no information at all.
ggplot(data = transaf2, aes(x = int, y = af, color = year)) +
  geom_point(position = position_dodge(width = 0.2), size = 3) +
  geom_line(aes(group = chrom), position = position_dodge(width = 0.2), linewidth = 1) +
  labs(x = "Year & Time Point", y = "AF", color = "Year") +
  ggtitle("changes af t1-t2")

# first I will check per year what is going on 
#lets try with 2009
af2009 <- af2 %>%
  mutate(af_difference = EA_2009_T1 - EA_2009_T2) %>%
  select(chrom, af_difference) %>%
  arrange(af_difference)

ggplot(data = af2009, aes(x = af_difference, y = 1)) +
  geom_point(size = 1) +
  labs(x = "af_difference", y = "Count") +
  ggtitle("Scatter Plot of af_difference")

summary_counts <- af2009 %>%
  summarize(positive_count = sum(af_difference > 0),
            negative_count = sum(af_difference < 0))
summary_counts$positive_count
#1105
summary_counts$negative_count
#716

positive2009<- af2009 %>%
  filter(af_difference > 0)

negative2009<- af2009 %>%
  filter(af_difference < 0)


#so here I go back to af2 and fetch the chroms that were increasing from start to end in 2009
filtered_positive_2009 <- af2 %>%
  filter(chrom %in% positive2009$chrom)

#first ill check if the snps that are increasing in one year behave the same in others (visually, tomorrow ill do it manually with the overlap)
#reorder chroms,
filtered_positive_2009$chrom <- factor(filtered_positive_2009 $chrom, levels = filtered_positive_2009 $chrom[order(filtered_positive_2009 $EA_2009_T2)])
#again....
# new columns

transafp2009<- gather(filtered_positive_2009, key = "Sample", value = "af", -chrom)
transafp2009 <- transafp2009 %>%
  mutate(time_point = ifelse(grepl("T1", Sample), "start", "end"))
transafp2009 <- transafp2009 %>%
  mutate(year = as.integer(sub("EA_(\\d{4})_.*", "\\1", Sample)))


#this graph gives no information at all.

ggplot(data = transafp2009, aes(x = interaction(year, time_point, lex.order = TRUE), y = af, color = year)) +
  geom_point(size = 1) +
  geom_line(aes(group = chrom), linewidth = 1) 
  labs(x = "Year & Time Point", y = "AF", color = "Year") +
  ggtitle("")


# Create a plot with the specified order, light gray lines, and thin points
transafp2009$int <- interaction(transafp2009$year, transafp2009$time_point)
# Check levels
levels(transafp2009$int)
levels <- c("2009.start", "2009.end", "2011.start", "2011.end", "2015.start", "2015.end", "2022.start", "2022.end")

# Reorder the levels of the factor
transafp2009$int <- factor(transafp2009$int, levels = levels)

ggplot(data = transafp2009, aes(x = int, y = af, color = year)) +
  geom_line(aes(group = chrom), color = "lightgray", linewidth = 0.5) +
  labs(x = "Year / Time Point", y = "AF", color = "Year") +
  ggtitle("") +
  theme_bw()
 # scale_x_discrete(limits = int, labels = int, sep = " "))


#same for negative
filtered_positive_2009$chrom <- factor(filtered_positive_2009 $chrom, levels = filtered_positive_2009 $chrom[order(filtered_positive_2009 $EA_2009_T4)])
#again....
# new columns
#so here I go back to af2 and fetch the chroms that were decreasing from start to end in 2009
filtered_negative_2009 <- af2 %>%
  filter(chrom %in% negative2009$chrom)
filtered_negative_2009 <- filtered_negative_2009[,-c(10)]
#filtered_negative_2009 <- filtered_negative_2009 %>%
#  arrange(EA_2009_T4)


#first ill check if the snps that are increasing in one year behave the same in others (visually, tomorrow ill do it manually with the overlap)
#reorder chroms,
transafpneg2009<- gather(filtered_negative_2009, key = "Sample", value = "af", -chrom)
transafpneg2009 <- transafpneg2009 %>%
  mutate(time_point = ifelse(grepl("T1", Sample), "start", "end"))
transafpneg2009 <- transafpneg2009 %>%
  mutate(year = as.integer(sub("EA_(\\d{4})_.*", "\\1", Sample)))

# Create a plot with the specified order, light gray lines, and thin points
transafpneg2009$int <- interaction(transafpneg2009$year, transafpneg2009$time_point)
# Check levels
levels(transafpneg2009$int)
levels <- c("2009.start", "2009.end", "2011.start", "2011.end", "2015.start", "2015.end", "2022.start", "2022.end")
# Reorder the levels of the factor
transafpneg2009$int <- factor(transafpneg2009$int, levels = levels)

# Create a plot
ggplot(data = transafpneg2009, aes(x = int, y = af, color = year)) +
  geom_line(aes(group = chrom), color = "lightgray", size = 0.5) +
  labs(x = "Year & Time Point", y = "AF", color = "Year") +
  ggtitle("") +
  theme_bw()

# okay, now I will check the top 10% and 1% of the snps in terms of pvalues, what are they doing

sig <- read.table(file="/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/GLM/significant-glm-181023-t2-time-p.txt", header= F) 

af <- freq[, c("EA_2009_T1","EA_2009_T2" ,"EA_2011_T1", "EA_2011_T2","EA_2015_T1",
 "EA_2015_T2", "EA_2022_T1", "EA_2022_T2")]
af2 <- cbind(chrom, af)
af2 <- af2 %>%
  inner_join(sig, by = c("chrom" = "V3")) %>%
  select(names(af2), V4)

#whath are the average p values
mean(af2$V4)
#0.02710126
quantile(af2$V4, 0.9)
#       90%
#0.04696166
threshold <- quantile(af2$V4, 0.1)
lower_10_percent_data <- af2[af2$V4 <= threshold, ]

# new columns
low10 <- lower_10_percent_data[, -which(names(lower_10_percent_data) == "V4")]
low10 <- gather(low10, key = "Sample", value = "af", -chrom)
low10 <- low10 %>%
  mutate(time_point = ifelse(grepl("T1", Sample), "start", "end"))
low10 <- low10 %>%
  mutate(year = as.integer(sub("EA_(\\d{4})_.*", "\\1", Sample)))
# Create a plot with the specified order, light gray lines, and thin points
low10$int <- interaction(low10$year, low10$time_point)
# Check levels
levels(low10$int)
levels <- c("2009.start", "2009.end", "2011.start", "2011.end", "2015.start", "2015.end", "2022.start", "2022.end")
# Reorder the levels of the factor
low10$int <- factor(low10$int, levels = levels)

ggplot(data = low10, aes(x = int, y = af, color = chrom)) +
  geom_line(aes(group = chrom), color = "lightgray", linewidth = 0.5) +
  labs(x = "years", y = "AF", color = "Year") +
  ggtitle("top 10% p values") +
  theme_bw()
```
# now for the third file using 

#THAT IS THE FINAL FILE
```bash
awk '{print $3}' glm-181023-2009t2-time.txt > glm-181023-2009t2-time-p.txt
awk -F'\t' '$1 <= 0.05 { count++ } END { print count }' glm-181023-2009t2-time-p.txt
#1099
paste chrom-poly-position-seasonal-data.txt joint.txt glm-181023-2009t2-time-p.txt > full-table-glm-181023-2009t2-time-p.txt
# add header
echo "chrom pos chrompos p" | cat - full-table-glm-181023-2009t2-time-p.txt > temp && mv temp full-table-glm-181023-2009t2-time-p.txt

#nano and fix header

awk -F'\t' '$4 <= 0.05 { print }' full-table-glm-181023-2009t2-time-p.txt > significant-glm-181023-2009t2-time-p.txt

awk 'NR==FNR{a[$1$2]=$0; next} ($1$2 in a){print a[$1$2]}' cov_with_chrom.txt significant-glm-181023-2009t2-time-p.txt > cov-sig-snps-glm-181023-2009t2-time-p.txt
```
# now I will plot in R the the output
module load R/4.1.1
module load gcc/10.2.0
```R
library(tidyr)
library(ggplot2)
library(matrixStats)

df<- read.table(file="/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/GLM/cov-sig-snps-glm-181023-2009t2-time-p.txt", header= F)
df <- df[-3]
df$chrompos <- paste(df$V1, df$V2)
colnames(df) <- c("chrom","pos",
                  "2009t1","2009t2","2009t3","2009t4",
                  "2011t1","2011t2",
                  "2015t1","2015t2","2015t3","2015t4",
                  "2022t1","2022t2","2022t3","2022t4",
                  "chrompos")

df <- df[, -c(2)]
df_long <- df %>%
  gather(key = "sample", value = "coverage", -chrompos, -chrom)

ggplot(df_long, aes(x = chrom, y = coverage, color = sample)) +
  geom_point(position = position_dodge(width = 0.2),size = 3) +
  labs(title = "coverage", x = "chromosome", y = "coverage") +
  #scale_x_discrete(labels = df_long$chrom) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
ggsave("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/figures/coverage-individual-snps-fixed-glm-2009t2.pdf", width = 8, height = 6, units = "in")

#how many SNPs have the average coverage (500x)
sum(df_long$coverage <= 500)
 # 14382 
 sum(df_long$coverage >= 500)
#1005

df$mean_cov <- rowMeans(df[, -c(1,16)])
df$se <- rowSds(as.matrix(df[, -c(1,16)])) / sqrt(ncol(df) - 2)

ggplot(df, aes(x = chrompos, y = mean_cov)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_cov - se, ymax = mean_cov + se), width = 0.2) +
  labs(title = "coverage", x = "chromosome", y = "coverage") +
  scale_x_discrete(labels = NULL) +
  theme_minimal()

df <- separate(df, chrompos, into = c("chrom", "pos"), sep = " ", remove=FALSE)
library(dplyr)
distinct_count <- df %>% 
  count(chrom)

#plot SNPs by chromossome
ggplot(distinct_count, aes(x = reorder(chrom, -n), y = n)) +
  geom_bar(stat = "identity") +
  labs(title = "SNP counts", x = "Chromosome", y = "Count") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
#ggsave("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/figures/snps-by-choromosome-glm-fixed.pdf", width = 10, height = 6, units = "in")

df$nchrom <- as.integer(factor(df$chrom))
ggplot(df, aes(x = factor(chrom), y = mean_cov, color=chrom)) +
  geom_boxplot(width = 0.5, show.legend = FALSE) +  # Adjust width as needed
  labs(title = "Mean Coverage", x = "Chromosome", y = "Coverage") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/figures/mean-coverage-snps-chrom-fixed-t2.pdf", width = 8, height = 6, units = "in")


#now I will find out which SNPs were deviating the norm in the QQplot, and see how they behave

#time unordered factor
pvalues <- read.table("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/GLM/full-table-glm-181023-2009t2-time-p.txt", header = TRUE)
n <- length(pvalues$p)
expected_p_values <- seq(0, 1, length.out = n)
# Create a QQ plot
qqplot(-log10(expected_p_values), -log10(pvalues$p), xlab = "Expected -log10(p)", ylab = "Observed -log10(p)", main = "unordered factor (UF)")
abline(0, 1, col = "red")

#okay now I want to plot cov in relation with pvalues
#########

covsig <- read.table(file="/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/GLM/cov-sig-snps-glm-181023-2009t2-time-p.txt", header= F)
covsig$chrompos <- paste(covsig$V1, covsig$V2, sep="_")
covsig <- covsig[, -c(1,2,3)]
pvalues <- read.table("/gxfs_work1/geomar/smomw573/seasonal_adaptation/analysis/GLM/significant-glm-181023-2009t2-time-p.txt", header = FALSE)
pvalues <- pvalues[, -c(1,2)]
colnames(pvalues) <- c("chrompos", "p")
merged_data <- merge(pvalues, covsig, by = "chrompos")

colnames(merged_data) <- c("chrompos", "p", "2009t1", "2009t2", "2009t3", "2009t4", "2011t1","2011t2","2015t1","2015t2","2015t3","2015t4","2022t1","2022t2","2022t3","2022t4")
merged_data$mean_cov <- rowMeans(merged_data[, 3:16], na.rm = TRUE)
merged_data$log <- -log10(merged_data$p)

ggplot(merged_data, aes(x = log, y = mean_cov)) +
  geom_point(size = 3, show.legend = FALSE) +  # Adjust size as needed
  labs(title = "Mean Coverage", x = "p value (-log10)", y = "Mean Coverage") +
  geom_smooth(method='lm')
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/figures/mean-cov-p-value-fixed-2009t2.pdf", width = 8, height = 6, units = "in")

model <- lm(mean_cov ~ log, data = merged_data)
summary(model)

###that is the final and correct one, the others before were not right!
# now I will check what the allele frequencies are doing!

chrom <- read.table(file="/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/GLM/chroms.txt", header= T)
freq <- read.table(file="/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/GLM/freq_new_vcf.txt", header= T)
sig <- read.table(file="/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/GLM/cov-sig-snps-glm-181023-2009t2-time-p.txt", header= F) 

sig$chrompos <- paste(sig$V1, sig$V2, sep="_")
af <- freq[, c("EA_2009_T2","EA_2009_T4" ,"EA_2011_T1", "EA_2011_T2","EA_2015_T1",
 "EA_2015_T4", "EA_2022_T1", "EA_2022_T4")]
af2 <- cbind(chrom, af)
af2 <- af2 %>%
  inner_join(sig, by = c("chrom" = "chrompos")) %>%
  select(all_of(names(af2)))

transaf<- gather(af2, key = "Sample", value = "af", -chrom)

# new columns
transaf <- transaf %>%
  mutate(time_point = ifelse(grepl("T1|2009_T2", Sample), "start", "end"))
transaf <- transaf %>%
  mutate(year = as.integer(sub("EA_(\\d{4})_.*", "\\1", Sample)))

transaf <- transaf %>%
  group_by(chrom, year) %>%
  mutate(change_af = af - lag(af, default = first(af)))

transaf$int <- interaction(transaf$year, transaf$time_point)
# Check levels
levels(transafp2009$int)
levels <- c("2009.start", "2009.end", "2011.start", "2011.end", "2015.start", "2015.end", "2022.start", "2022.end")
# Reorder the levels of the factor
transaf$int <- factor(transaf$int, levels = levels)
transaf2 <- transaf %>%
  arrange(af)


#this graph gives no information at all.
ggplot(data = transaf2, aes(x = int, y = af, color = year)) +
  geom_point(position = position_dodge(width = 0.2), size = 3) +
  geom_line(aes(group = chrom), position = position_dodge(width = 0.2), linewidth = 1) +
  labs(x = "Year & Time Point", y = "AF", color = "Year") +
  ggtitle("changes af ")

# first I will check per year what is going on 
#lets try with 2009
af2009 <- af2 %>%
  mutate(af_difference = EA_2009_T2 - EA_2009_T4) %>%
  select(chrom, af_difference) %>%
  arrange(af_difference)

ggplot(data = af2009, aes(x = af_difference, y = 1)) +
  geom_point(size = 1) +
  labs(x = "af_difference", y = "Count") +
  ggtitle("Scatter Plot of af_difference")

summary_counts <- af2009 %>%
  summarize(positive_count = sum(af_difference > 0),
            negative_count = sum(af_difference < 0))
summary_counts$positive_count
#652
summary_counts$negative_count
#443
positive2009<- af2009 %>%
  filter(af_difference > 0)
negative2009<- af2009 %>%
  filter(af_difference < 0)

#so here I go back to af2 and fetch the chroms that were increasing from start to end in 2009
filtered_positive_2009 <- af2 %>%
  filter(chrom %in% positive2009$chrom)

#first ill check if the snps that are increasing in one year behave the same in others (visually, tomorrow ill do it manually with the overlap)
#reorder chroms,
filtered_positive_2009$chrom <- factor(filtered_positive_2009$chrom, levels = filtered_positive_2009$chrom[order(filtered_positive_2009$EA_2009_T4)])
#again....
# new columns

transafp2009<- gather(filtered_positive_2009, key = "Sample", value = "af", -chrom)
transafp2009 <- transafp2009 %>%
  mutate(time_point = ifelse(grepl("T1|2009_T2", Sample), "start", "end"))
transafp2009 <- transafp2009 %>%
  mutate(year = as.integer(sub("EA_(\\d{4})_.*", "\\1", Sample)))


#this graph gives no information at all.

ggplot(data = transafp2009, aes(x = interaction(year, time_point, lex.order = TRUE), y = af, color = year)) +
  geom_point(size = 1) +
  geom_line(aes(group = chrom), linewidth = 1) 
  labs(x = "Year & Time Point", y = "AF", color = "Year") +
  ggtitle("")

# Create a plot with the specified order, light gray lines, and thin points
transafp2009$int <- interaction(transafp2009$year, transafp2009$time_point)
# Check levels
levels(transafp2009$int)
levels <- c("2009.start", "2009.end", "2011.start", "2011.end", "2015.start", "2015.end", "2022.start", "2022.end")

# Reorder the levels of the factor
transafp2009$int <- factor(transafp2009$int, levels = levels)

ggplot(data = transafp2009, aes(x = int, y = af, color = year)) +
  geom_line(aes(group = chrom), color = "lightgray", linewidth = 0.5) +
  labs(x = "Year / Time Point", y = "AF", color = "Year") +
  ggtitle("") +
  theme_bw()
 # scale_x_discrete(limits = int, labels = int, sep = " "))

#same for negative

#so here I go back to af2 and fetch the chroms that were decreasing from start to end in 2009
filtered_negative_2009 <- af2 %>%
  filter(chrom %in% negative2009$chrom)
filtered_negative_2009 <- filtered_negative_2009[,-c(10)]
#filtered_negative_2009 <- filtered_negative_2009 %>%
#  arrange(EA_2009_T4)


#first ill check if the snps that are increasing in one year behave the same in others (visually, tomorrow ill do it manually with the overlap)
#reorder chroms,
transafpneg2009<- gather(filtered_negative_2009, key = "Sample", value = "af", -chrom)
transafpneg2009 <- transafpneg2009 %>%
  mutate(time_point = ifelse(grepl("T1|2009_T2", Sample), "start", "end"))
transafpneg2009 <- transafpneg2009 %>%
  mutate(year = as.integer(sub("EA_(\\d{4})_.*", "\\1", Sample)))

# Create a plot with the specified order, light gray lines, and thin points
transafpneg2009$int <- interaction(transafpneg2009$year, transafpneg2009$time_point)
# Check levels
levels(transafpneg2009$int)
levels <- c("2009.start", "2009.end", "2011.start", "2011.end", "2015.start", "2015.end", "2022.start", "2022.end")
# Reorder the levels of the factor
transafpneg2009$int <- factor(transafpneg2009$int, levels = levels)

# Create a plot
ggplot(data = transafpneg2009, aes(x = int, y = af, color = year)) +
  geom_line(aes(group = chrom), color = "lightgray", size = 0.5) +
  labs(x = "Year & Time Point", y = "AF", color = "Year") +
  ggtitle("") +
  theme_bw()

# okay, now I will check the top 10% and 1% of the snps in terms of pvalues, what are they doing

sig <- read.table(file="/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/GLM/significant-glm-181023-2009t2-time-p.txt", header= F) 

af <- freq[, c("EA_2009_T2","EA_2009_T4" ,"EA_2011_T1", "EA_2011_T2","EA_2015_T1",
 "EA_2015_T4", "EA_2022_T1", "EA_2022_T4")]
af2 <- cbind(chrom, af)
af2 <- af2 %>%
  inner_join(sig, by = c("chrom" = "V3")) %>%
  select(names(af2), V4)

#whath are the average p values
mean(af2$V4)
#0.02710126
quantile(af2$V4, 0.9)
#       90%
#0.04696166
threshold <- quantile(af2$V4, 0.1)
lower_10_percent_data <- af2[af2$V4 <= threshold, ]

# new columns
low10 <- lower_10_percent_data[, -which(names(lower_10_percent_data) == "V4")]
low10 <- gather(low10, key = "Sample", value = "af", -chrom)
low10 <- low10 %>%
  mutate(time_point = ifelse(grepl("T1|2009_T2", Sample), "start", "end"))
low10 <- low10 %>%
  mutate(year = as.integer(sub("EA_(\\d{4})_.*", "\\1", Sample)))
# Create a plot with the specified order, light gray lines, and thin points
low10$int <- interaction(low10$year, low10$time_point)
# Check levels
levels(low10$int)
levels <- c("2009.start", "2009.end", "2011.start", "2011.end", "2015.start", "2015.end", "2022.start", "2022.end")
# Reorder the levels of the factor
low10$int <- factor(low10$int, levels = levels)

ggplot(data = low10, aes(x = int, y = af, color = chrom)) +
  geom_line(aes(group = chrom), color = "lightgray", linewidth = 0.5) +
  labs(x = "years", y = "AF", color = "Year") +
  ggtitle("top 10% p values") +
  theme_bw()
```
# play around with allele freqs to figure out what is going on

first I will normalize the data and plot everything for T2 2009, and then every start end

will do it in R

```bash
module load R/4.1.1
module load gcc/10.2.0
```
```R
library(tidyr)
library(ggplot2)
library(matrixStats)

freq <- read.table(file="/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/GLM/freq_new_vcf.txt", header= T)
chrom <- read.table(file="/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/GLM/chroms.txt", header= T)
sig <- read.table(file="/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/GLM/cov-sig-snps-glm-181023-2009t2-time-p.txt", header= F) 

sig$chrompos <- paste(sig$V1, sig$V2, sep="_")

af <- freq[, c("EA_2009_T2","EA_2009_T4" ,"EA_2011_T1", "EA_2011_T2","EA_2015_T1", "EA_2015_T4", "EA_2022_T1", "EA_2022_T4")]
af2 <- cbind(chrom, af)

colnames(af2) <- c("chrom", "2009.start", "2009.end",
  "2011.start", "2011.end", 
  "2015.start", "2015.end",
  "2022.start","2022.end")

#normalize everything for 2009
normalized_af <- af2 %>%
  mutate(across(starts_with("EA_"), ~ . - 2009.start))

#select significant snps

normaf2 <- normalized_af %>%
  inner_join(sig, by = c("chrom" = "chrompos")) %>%
  select(all_of(names(normalized_af))) #do not include columns in sig

colnames(normaf2) <- c("chrom", "2009.start", "2009.end",
  "2011.start", "2011.end", 
  "2015.start", "2015.end",
  "2022.start","2022.end")

transnormaf<- gather(normaf2, key = "year", value = "af", -chrom)
#now 
levels <- c("2009.start", "2009.end", "2011.start", "2011.end", "2015.start", "2015.end", "2022.start", "2022.end")
# Reorder the levels of the factor
transnormaf$year <- factor(transnormaf$year, levels = levels)

d <- ggplot(data = transnormaf, aes(x = year, y = af, color = year)) +
  geom_line(aes(group = chrom), color = "lightgray", linewidth = 0.5) +
  labs(x = "Year / Time Point", y = "AF", color = "Year") +
  ggtitle("") +
  theme_bw()

ggsave("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/figures/af-normalized2009.pdf",d, w=6, h=6)

#doesnt help
ggplot(transnormaf, aes(x = year, y = af, fill = year)) +
  geom_boxplot() +
  stat_summary(fun = "mean", geom = "point", shape = 8,
               size = 2, color = "white")

# Calculate mean and SE by year
sumdata <- transnormaf %>%
  group_by(year) %>%
  summarize(mean_af = mean(af), se = sd(af) / sqrt(n()))

#percentage positive SNPs (increasing AF) per year

percpos <- colSums(result[, -1] > 0) / nrow(result) * 100

# Print the percentage of positive values for each column
print(percpos)

# Create the dot plot
  ggplot(sumdata, aes(x = year, y = mean_af)) +
  geom_dotplot(binaxis = "y", binwidth = 0.05, dotsize = 0.1, stackdir = "center") +
  geom_errorbar(aes(ymin = mean_af - se, ymax = mean_af + se), width = 0.2) +
  labs(x = "Year", y = "Mean AF") +
  ggtitle("Average change in AF") +
  ylim(-0.05, 0.05)

#estimate percentage of positive values (or SNPs that went up from 2009 start)

result <- normaf2 %>%
  mutate(across(starts_with("2009.") | starts_with("2011.") | starts_with("2015.") | starts_with("2022."), ~ . > 0))

# Print the result
print(result)

#okay, now from this I want to keep the values which are increasing 

result2 <- result %>%
  filter_at(vars(ends_with(".end")), all_vars(. == TRUE))
#366 snps

#create the dataset with these snps
merged <- merge(result2, normaf2, by = "chrom", suffixes = c(".result2", ".normaf2"))
# Find the column names ending with "normaf2"
normaf2_columns <- grep("normaf2", names(merged), value = TRUE)

# Select the columns with those names
positivesnps <- merged %>%
  select(chrom, all_of(normaf2_columns))

positivesnpst <- gather(positivesnps, key = "year", value = "af", -chrom)

levels <- c("2009.start.normaf2", "2009.end.normaf2", "2011.start.normaf2", "2011.end.normaf2", "2015.start.normaf2", "2015.end.normaf2", "2022.start.normaf2", "2022.end.normaf2")
# Reorder the levels of the factor
positivesnpst$year <- factor(positivesnpst$year, levels = levels)

ggplot(data = positivesnpst, aes(x = year, y = af, color = year)) +
  geom_line(aes(group = chrom), color = "lightgray", linewidth = 0.5) +
  labs(x = "Year / Time Point", y = "AF", color = "Year") +
  ggtitle("") +
  theme_bw()

#negative snps
neg <- result %>%
  filter_at(vars(ends_with(".end")), all_vars(. == FALSE))
mergedneg <- merge(neg, normaf2, by = "chrom", suffixes = c(".result2", ".normaf2"))
negsnps <- mergedneg %>%
  select(chrom, all_of(normaf2_columns))
negsnpst <- gather(negsnps, key = "year", value = "af", -chrom)
negsnpst$year <- factor(negsnpst$year, levels = levels)

ggplot(data = negsnpst, aes(x = year, y = af, color = year)) +
  geom_line(aes(group = chrom), color = "lightgray", linewidth = 0.5) +
  labs(x = "Year / Time Point", y = "AF", color = "Year") +
  ggtitle("") +
  theme_bw()


###now I'll normalize things by end start every year

normafindyear <- af2 %>%
  mutate(
    '2009.end.normalized' = `2009.end` - `2009.start`,
    '2011.end.normalized' = `2011.end` - `2011.start`,
    '2015.end.normalized' = `2015.end` - `2015.start`,
    '2022.end.normalized' = `2022.end` - `2022.start`
  )

normafindyear <- normafindyear %>%
  mutate(
    '2009.start' = 0,
    '2011.start' = 0,
    '2015.start' = 0,
    '2022.start' = 0
  )
# Select the columns you want to keep
selected_columns <- normafindyear %>%
  select(chrom, , ends_with(".start"), ends_with(".normalized"))

# Create a new dataset with the selected columns
normafindyear <- selected_columns

# Rename the columns if needed
colnames(normafindyear) <- gsub(".normalized", "", colnames(normafindyear))

# View the new dataset
head(normafindyear)

normaf2 <- normafindyear %>%
  inner_join(sig, by = c("chrom" = "chrompos")) %>%
  select(all_of(names(normafindyear))) #do not include columns in sig

normindyear<- gather(normaf2, key = "year", value = "af", -chrom)
#now 
levels <- c("2009.start", "2009.end", "2011.start", "2011.end", "2015.start", "2015.end", "2022.start", "2022.end")
# Reorder the levels of the factor
normindyear$year <- factor(normindyear$year, levels = levels)

ggplot(data = normindyear, aes(x = year, y = af, color = year)) +
  geom_line(aes(group = chrom), color = "lightgray", linewidth = 0.5)+
  labs(x = "Year / Time Point", y = "AF", color = "Year") +
  ggtitle("") +
  theme_bw()

#now the pos and negatives###
#############################

#estimate percentage of positive values (or SNPs that went up from 2009 start)

result <- normaf2 %>%
  mutate(across(starts_with("2009.") | starts_with("2011.") | starts_with("2015.") | starts_with("2022."), ~ . > 0))

# Print the result
print(result)

#okay, now from this I want to keep the values which are increasing 

result2 <- result %>%
  filter_at(vars(ends_with(".end")), all_vars(. == TRUE))
#335 snps

#create the dataset with these snps
merged <- merge(result2, normaf2, by = "chrom", suffixes = c(".result2", ".normaf2"))
# Find the column names ending with "normaf2"
normaf2_columns <- grep("normaf2", names(merged), value = TRUE)

# Select the columns with those names
positivesnps <- merged %>%
  select(chrom, all_of(normaf2_columns))

positivesnpst <- gather(positivesnps, key = "year", value = "af", -chrom)
colnames(positivesnpst) <- gsub(".normalized", "", colnames(positivesnpst))

levels <- c("2009.start", "2009.end", "2011.start", "2011.end", "2015.start", "2015.end", "2022.start", "2022.end")
# Reorder the levels of the factor
positivesnpst$year <- factor(positivesnpst$year, levels = levels)

d <- ggplot(data = positivesnpst, aes(x = year, y = af, color = year)) +
  geom_line(aes(group = chrom), color = "lightgray", linewidth = 0.5) +
  labs(x = "Year / Time Point", y = "AF", color = "Year") +
  ggtitle("") +
  theme_bw()

ggsave("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/figures/af-normalizedperyear-pos.pdf",d, w=7, h=5)

#so now I have the list of the significant snps which constantly increased from start to end,
#I will pull the original AF from these snps and plot it

finalpos <- af2 %>%
  inner_join(positivesnps, by = c("chrom" = "chrom")) %>%
  select(all_of(names(af2)))
finalpost <- gather(finalpos, key = "year", value = "af", -chrom)
levels <- c("2009.start", "2009.end", "2011.start", "2011.end", "2015.start", "2015.end", "2022.start", "2022.end")
# Reorder the levels of the factor
finalpost$year <- factor(finalpost$year, levels = levels)
ggplot(data = finalpost, aes(x = year, y = af, color = year)) +
  geom_line(aes(group = chrom), color = "lightgray", linewidth = 0.5) +
  labs(x = "Year / Time Point", y = "AF", color = "Year") +
  ggtitle("") +
  theme_bw()


#negative snps
neg <- result %>%
  filter_at(vars(ends_with(".end")), all_vars(. == FALSE))
#474 snps

mergedneg <- merge(neg, normaf2, by = "chrom", suffixes = c(".result2", ".normaf2"))
normaf2_columns <- grep("normaf2", names(mergedneg), value = TRUE)
negsnps <- mergedneg %>%
  select(chrom, all_of(normaf2_columns))
colnames(negsnps) <- gsub(".normaf2", "", colnames(negsnps))
negsnpst <- gather(negsnps, key = "year", value = "af", -chrom)
negsnpst$year <- factor(negsnpst$year, levels = levels)

ggplot(data = negsnpst, aes(x = year, y = af, color = year)) +
  geom_line(aes(group = chrom), color = "lightgray", linewidth = 0.5) +
  labs(x = "Year / Time Point", y = "AF", color = "Year") +
  ggtitle("") +
  theme_bw()

#so now I have the list of the significant snps which constantly decreased from start to end,
#I will pull the original AF from these snps and plot it
mergedneg <- merge(neg, normaf2, by = "chrom", suffixes = c(".result2", ".normaf2"))
normaf2_columns <- grep("normaf2", names(mergedneg), value = TRUE)
negsnps <- mergedneg %>%
  select(chrom, all_of(normaf2_columns))
colnames(negsnps) <- gsub(".normaf2", "", colnames(negsnps))
finalneg <- af2 %>%
  inner_join(negsnps, by = c("chrom" = "chrom")) %>%
  select(all_of(names(af2)))
finalnegt <- gather(finalneg, key = "year", value = "af", -chrom)
levels <- c("2009.start", "2009.end", "2011.start", "2011.end", "2015.start", "2015.end", "2022.start", "2022.end")
# Reorder the levels of the factor
finalnegt$year <- factor(finalnegt$year, levels = levels)
d <- ggplot(data = finalnegt, aes(x = year, y = af, color = year)) +
  geom_line(aes(group = chrom), color = "lightgray", linewidth = 0.5) +
  labs(x = "Year / Time Point", y = "AF", color = "Year") +
  ggtitle("") +
  theme_bw()

ggsave("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/figures/af-normalizedperyear-neg.pdf",d, w=7, h=5)

#now I am going to export the datasets with the list of snps

write.table(positivesnps, "/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/GLM/positive-sig-seasonal-snps-list.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(negsnps, "/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/GLM/neg-sig-seasonal-snps-list.txt", sep = "\t", quote = FALSE, row.names = FALSE)


```
# 26-10-2023
# Running bypass and byenv
```bash
#in the end I dowloaded PGDspider with wget
#went to the software directory and ran 
srun --pty --x11 --nodes=1 --cpus-per-task=4 --mem=50000 --time=07:00:00 /bin/bash
module load bcftools/1.10.2

#subset the vcf for the pops of interest
bcftools view -s EA_2009_T2,EA_2009_T4,EA_2011_T1,EA_2011_T2,EA_2015_T1,EA_2015_T4,EA_2022_T1,EA_2022_T4 /gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/variants/final_20230810.vcf -o /gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/variants/final_20230810_subset.vcf.gz

cd /gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/

#baypass
git clone https://gitlab.com/YDorant/Toolbox


#ok lets try it with a subsample data
bcftools view -Ov -o /gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/variants/final_20230810_subset_corrected2.vcf -m2 -M2 -v snps -n 100000 /gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/variants/final_20230810_subset.vcf.gz
head -n 50000 final_20230810_subset.vcf > first_50000_lines.vcf
grep -v "^#" first_50000_lines.vcf |wc -l
499743

#what is going on in the example data?
python Toolbox/reshaper_baypass.py ~/canada.vcf.gz ~/popmap_canada.txt /gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/variants/canada.baypass

#okay

python Toolbox/reshaper_baypass.py /gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/variants/ /gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/variants/popmap2.txt /gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/variants/final_20230810_subset_trimmed2.baypass
#works so I will run a job (below)
#but first try bayenv dataprep
#PDGSpider
#if needed conda activate java
#create spid file on GUI
java -Xmx1024m -Xms512m -jar PGDSpider2.jar
#run command line
java -Xmx1024m -Xms512m -jar PGDSpider2-cli.jar -inputfile /gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/variants/final_20230810_subset.vcf -inputformat VCF -outputfile /gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/variants/final_20230810_trimmed.bayenv -outputformat BAYENV -spid /gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/variants/template_VCF_BAYENV.spid

#okay, it works with the trimmed data
#so I will also run a job
#ran the job - does not work! something to do with vcf format but Im in a hurry and will try poolfstat (below scripts)

module load R/4.3.1
module load gcc/10.2.0
```
run it in a job
convert-baypass.sh

#!/bin/bash
#SBATCH -D /gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/
#SBATCH --mail-type=END
#SBATCH --mail-user=jnascimento@geomar.de
#SBATCH --partition=cluster
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=5:00:00
#SBATCH --job-name=trim
#SBATCH --output=baypass_convert.out
#SBATCH --error=baypass_convert.err


python Toolbox/reshaper_baypass.py /gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/variants/final_20230810_subset.vcf /gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/variants/popmap2.txt /gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/variants/final_20230810_subset.baypass



###############################

#job 2
convert-bayenv.sh

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

conda activate java

java -Xmx1024m -Xms1024m -jar PGDSpider2-cli.jar -inputfile /gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/variants/final_20230810_subset.vcf -inputformat VCF -outputfile /gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/variants/final_20230810_trimmed.bayenv -outputformat BAYENV -spid /gxfs_home/geomar/smomw573/software/PGDSpider_2.1.1.5/spid-bayenv.spid


```bash 

#check if data is okay with python
with open('/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/variants/popmap2.txt', 'r') as file:
    for line_number, line in enumerate(file, start=1):
        # Remove leading and trailing whitespace and split by comma
        parts = line.strip().split(',')
        
        if len(parts) == 2:
            ind, pop = parts
            print(f"Line {line_number}: ind = {ind}, pop = {pop}")
        else:
            print(f"Error in line {line_number}: Invalid format - {line}")
```

```
```R
library(vcfR)
library(adegenet)
library(hierfstat)

#transform vcf into bayescan input data
vcf <- read.vcfR("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/variants/final_20230810_subset.vcf.gz", verbose=FALSE)
pop_map <- read.table("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/variants/popmap.txt", header=TRUE, stringsAsFactors = TRUE)
show(vcf)
head(vcf)
genind <- vcfR2genind(vcf)

genind@pop <- pop_map$STRATA
save(genind, file="final_vcf_subset.RData")

pop(genind)
#no pops
pop_names <- indNames(genind)
pop2 <- rep (NA, length(pop_names))

pop2[grep("2009_T2", pop_names)] <- "2009.start"
pop2[grep("9_T4", pop_names)] <- "2009.end"
pop2[grep("11_T1", pop_names)] <- "2011.start"
pop2[grep("11_T2", pop_names)] <- "2011.end"
pop2[grep("5_T1", pop_names)] <- "2015.start"
pop2[grep("5_T4", pop_names)] <- "2015.end"
pop2[grep("2_T1", pop_names)] <- "2022.start"
pop2[grep("2_T4", pop_names)] <- "2022.end"

pop2

pop(genind) 
#nothing in the pop field at the moment for this genind...
# try attaching this vector of popdata to the genind object...
#?strata
strata(genind) <- data.frame(pop2)
strata(genind)
setPop(genind) <- ~pop2
genind

```
# nothing worked so lets try poolfstat

```R

library(poolfstat)
dat <- vcf2pooldata(vcf.file="/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/variants/final_20230810_subset.vcf",poolsizes=rep(50,8))

pooldata2genobaypass(dat,writing.dir="/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/variants/")

```
# that did work
# 30-10-2023

#install baypass following instructions https://forgemia.inra.fr/mathieu.gautier/baypass_public
#produced gfile using poolfstats, efile and poolsize file manually

#trial run
```bash
/gxfs_home/geomar/smomw573/software/baypass_public-master/sources/g_baypass \
-gfile /gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/variants/headgenobaypass \
-efile /gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/variants/cov-baypass.txt \
-poolsizefile /gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/variants/haploid-size-baypass-txt \
-outprefix /gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/baypass/trial1 -nthreads 2

#okay works now I'll run a job with the whole file 
```
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


# dec 16th

# core model

#!/bin/bash
#SBATCH -D /gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/
#SBATCH --mail-type=END
#SBATCH --mail-user=jnascimento@geomar.de
#SBATCH --partition=base
#SBATCH --nodes=5
#SBATCH --tasks-per-node=2
#SBATCH --cpus-per-task=4
#SBATCH --mem=50G
#SBATCH --time=30:00:00
#SBATCH --job-name=bayenv-conv
#SBATCH --output=baypass-core.out
#SBATCH --error=baypass-core.err

inputdir=/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/variants
outputdir=/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/baypass
baypassdir=/gxfs_home/geomar/smomw573/software/baypass_public-master/sources

$baypassdir/g_baypass \
-gfile $inputdir/genobaypass \
-poolsizefile $inputdir/haploid-size-baypass-txt \
-d0yij 20 \
-npilot 100 \
-outprefix $outputdir/core-model -nthreads 10


# Jan 08
# standard covariate model

#!/bin/bash
#SBATCH -D /gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/
#SBATCH --mail-type=END
#SBATCH --mail-user=jnascimento@geomar.de
#SBATCH --partition=base
#SBATCH --nodes=5
#SBATCH --tasks-per-node=2
#SBATCH --cpus-per-task=4
#SBATCH --mem=50G
#SBATCH --time=48:00:00
#SBATCH --job-name=bayenv-stand
#SBATCH --output=baypass-stand.out
#SBATCH --error=baypass-stand.err


inputdir=/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/variants
outputdir=/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/baypass
baypassdir=/gxfs_home/geomar/smomw573/software/baypass_public-master/sources

$baypassdir/g_baypass \
-gfile $inputdir/genobaypass \
-efile $inputdir/cov-baypass.txt \
-omegafile $outputdir/core-model_mat_omega.out \
-poolsizefile $inputdir/haploid-size-baypass-txt \
-covmcmc \
-d0yij 20 \
-outprefix $outputdir/stand-cov-model -nthreads 10 -seed 3


# auxiliary model
#!/bin/bash
#SBATCH -D /gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/
#SBATCH --mail-type=END
#SBATCH --mail-user=jnascimento@geomar.de
#SBATCH --partition=base
#SBATCH --nodes=5
#SBATCH --tasks-per-node=2
#SBATCH --cpus-per-task=4
#SBATCH --mem=50G
#SBATCH --time=48:00:00
#SBATCH --job-name=bayenv-aux
#SBATCH --output=baypass-aux.out
#SBATCH --error=baypass-aux.err


inputdir=/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/variants
outputdir=/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/baypass
baypassdir=/gxfs_home/geomar/smomw573/software/baypass_public-master/sources

$baypassdir/g_baypass \
-gfile $inputdir/genobaypass \
-efile $inputdir/cov-baypass.txt.constrast \
-omegafile $outputdir/core-model_mat_omega.out \
-poolsizefile $inputdir/haploid-size-baypass-txt \
-d0yij 20 \
-outprefix $outputdir/aux-model -nthreads 10 -auxmodel



# analysis of core model in R
# Jan 09
# following baypass manual

 module load gcc12-env
 module load R/4.3.1

```R
library(corrplot)
library(ape)
library(WriteXLS)
library(RColorBrewer)
#source the baypass R functions (check PATH)
source("baypass_utils.R")
#upload the estimated Omega matrix
omega=as.matrix(read.table("core-model_mat_omega.out"))
pop.names=c("EA_2009_T1","EA_2009_T2","EA_2011_T1","EA_2011_T2","EA_2015_T1","EA_2015_T2","EA_2022_T1","EA_2022_T2")
dimnames(omega)=list(pop.names,pop.names)
# Visualization of the matrix

# Using SVD decomposition
pdf("omega_plot_core-model.pdf")
plot.omega(omega=omega,pop.names=pop.names)
dev.off()
# as a correlation plot
require(corrplot)
# Define your custom color palette (for example, using RColorBrewer palette)
COL2(diverging = c("RdBu", "BrBG", "PiYG", "PRGn", "PuOr", "RdYlBu"), n = 200)

# Calculate the correlation matrix from omega
cor.mat = cov2cor(omega)

# Plot the correlation matrix using corrplot with custom colors
pdf("hm_omega_plot_core-model.pdf")
corrplot(cor.mat, method = "color", mar = c(2, 1, 2, 2) + 0.1, col.lim=c(0,1),
         main = expression("Correlation map based on" ~ hat(Omega)),
         col = COL2('RdBu', 20))
dev.off()
# as a heatmap and hierarchical clustering tree (using the average agglomeration method)
hclust.ave <- function(x) hclust(x, method="average")
pdf("hm2_core-model.pdf")
heatmap(1-cor.mat,hclustfun = hclust.ave,
main=expression("Heatmap of "~hat(Omega)~"("*d[ij]*"=1-"*rho[ij]*")"))
dev.off()

#Estimates of the XtX differentiation measures (using the calibrated XtXst estimator)
core.snp.res=read.table("core-model_summary_pi_xtx.out", h=T)
#check behavior of the p-values associated to the XtXst estimator
hist(10**(-1*core.snp.res$log10.1.pval.),freq=F,breaks=50)
abline(h=1)


core.snp.res.thresh <- core.snp.res$M_XtX
thresh=quantile(core.snp.res.thresh,probs=0.99)

layout(matrix(1:2,2,1))
pdf("core-model-pvalues.pdf")
plot(core.snp.res$XtXst)
plot(core.snp.res$log10.1.pval.,ylab="XtX P-value (-log10 scale)")
abline(h=thresh,lty=2)
dev.off()
#abline(h=3,lty=2) #0.001 p--value theshold
thresh
     99%
8.323625

write.table(core.snp.res.thresh, file="core-model.snp.scores.txt", sep="\t", quote=FALSE)
```

# 10 January 24
# re-ran all models cause I had d0ij as 10 not 20.

# c2 stats model
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


# ok, at the end I ran five models for each (one with no seed option and one with a seed option), so I can compare the output.
# not sure how to do that yet, but I will first print the qvalues for both C2 models, 
# I also gave up using the st cov model as this is not the best fit for our dataset. So it will be between the auxilary and the C2 models.
# okay there was something wrong with the file path (not sure why), so need to rerun all models. Will start the analysis tomorrow. 
# 15.01.2024

# visualising outputs
# 19.01.2024

 module load gcc12-env
 module load R/4.3.1

```R
#st model

covaux.snp.res=read.table("std-model_summary_betai.out",h=T)
covaux.snp.xtx=read.table("std-model_summary_pi_xtx.out",h=T)$M_XtX
graphics.off()
layout(matrix(1:3,3,1))
plot(covaux.snp.res$BF.dB.,xlab="SNP",ylab="BFmc (in dB)")
plot(covaux.snp.res$M_Beta,xlab="SNP",ylab=expression(beta~"coefficient"))
plot(covaux.snp.xtx,xlab="SNP",ylab="XtX corrected for SMS")

#aux model

bf=read.table("aux-model_summary_betai.out",h=T)
hist(bf$BF.dB)
# Filter rows where the 'BF.dB' column is greater than a certain value (let's say 20)
filtered_data <- bf[bf$BF.dB > 20, ]

# save data to check if they are associated to END or START
write.table(filtered_data, file = "aux-bf.txt", sep = "\t", row.names = FALSE)

```

# Jan 12

# calculating qvalue for pvalues from the constrast model 
# done this in my pc

 module load R/4.3.1
 ```R
library(devtools)
setwd("~/Jenny-Eurytemora_affinis/baypass")
#install_github("jdstorey/qvalue")
browseVignettes(package = "qvalue")
library(qvalue)
library(dplyr)
library(purrr)


#set working directory
setwd("~/Jenny-Eurytemora_affinis/baypass/AUX-model")
#auxilary model convergence
df1 <- read.table("aux-model_summary_betai.out", header= TRUE)
df2 <- read.table("aux-model2_summary_betai.out", header= TRUE)
df3 <- read.table("aux-model-3_summary_betai.out", header= TRUE)
df4 <- read.table("aux-model-4_summary_betai.out", header= TRUE)
df5 <- read.table("aux-model-5_summary_betai.out", header= TRUE)

head(df1)
head(df2)

par(2,1)
plot(df1$BF.dB)
plot(df2$BF.dB)

#n of significant SNPs (BF > 20)

bf1 <- df1[df1$BF.dB >= 20, ]
nrow(bf1)
#3270
bf2 <- df2[df2$BF.dB >= 20, ]
nrow(bf2)
#3227
bf3 <- df3[df3$BF.dB >= 20, ]
nrow(bf3)
#3354

#list_of_dfs <- list(bf1, bf2, bf3)
#overlap <- reduce(list_of_dfs, inner_join, by = "MRK")
#333
#write.table(overlap, "overlap.txt", sep = "\t", quote = FALSE, row.names = FALSE)

df_list <- list(df1, df2, df3, df4, df5)

# Use lapply to extract the BF.dB column from each data frame
result_df <- data.frame(lapply(df_list, function(df) df$BF.dB))

# Rename the columns if needed
colnames(result_df) <- paste0("bf", 1:5)

bf_combined <- cbind(marker = 1:nrow(result_df), result_df)

#check for convergence
bf_columns <- grep("^bf", names(bf_combined), value = TRUE)

# Create an empty matrix to store correlation values
cor_matrix <- matrix(NA, nrow = length(bf_columns), ncol = length(bf_columns),
                     dimnames = list(bf_columns, bf_columns))

# Loop through the pairs of columns
for (i in seq_along(bf_columns)) {
  for (j in seq_along(bf_columns)) {
    cor_matrix[i, j] <- cor(x = bf_combined[[bf_columns[i]]], 
                            y = bf_combined[[bf_columns[j]]], 
                            method = "pearson", 
                            use = "complete.obs")
  }
}

# Print the correlation matrix
print(cor_matrix)

ggscatter(data = bf_combined, x = "bf1", y = "bf2", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "BF-model1", ylab = "BF-model2")

#calculate median for all of them
head(bf_combined)
bf_combined$median_bf <- apply(bf_combined[, c("bf1", "bf2", "bf3", "bf4", "bf5")], 1, median)

plot(bf_combined$median_bf)

#n of significant SNPs (BF > 20)

bf_sig <- bf_combined[bf_combined$median_bf >= 20, ]
nrow(bf_sig)
#951
bf_sig2 <- bf_combined[bf_combined$median_bf >= 15, ]
nrow(bf_sig2)
#2161




#######################################################
#now I'll do the same for c2

setwd("~/Jenny-Eurytemora_affinis/baypass/C2-model")

df1 <- read.table("c2-model_summary_contrast.out", header= TRUE)
df2 <- read.table("c2-model2_summary_contrast.out", header= TRUE)
df3 <- read.table("c2.3-model_summary_contrast.out", header= TRUE)
df4 <- read.table("c2.4-model_summary_contrast.out", header= TRUE)
df5 <- read.table("c2.5-model_summary_contrast.out", header= TRUE)
head(df1)

pval1 <- df1$log10.1.pval
range(df1$log10.1.pval)
plot(data$log10.1.pval)

pval2 <- df2$log10.1.pval
pval3 <- df3$log10.1.pval
pval4 <- df4$log10.1.pval
pval5 <- df5$log10.1.pval
head(pval2)
range(pval2)
pval1.2 <- 10^(-pval1)
pval2.2 <- 10^(-pval2)
pval3.2 <- 10^(-pval3)
pval4.2 <- 10^(-pval4)
pval5.2 <- 10^(-pval5)
range(pval2.2)
head(pval1.2)

pdf("Pvalues-5runs.pdf", width = 8, height = 6)
par(mfrow = c(3,2))
hist(pval1.2)
hist(pval2.2)
hist(pval3.2)
hist(pval4.2)
hist(pval5.2)
dev.off()

qobj <- qvalue(p = pval1.2)
qobj2 <- qvalue(p = pval2.2)
qobj3 <- qvalue(p = pval3.2)
qobj4 <- qvalue(p = pval4.2)
qobj5 <- qvalue(p = pval5.2)
#check if they're alright
plot(qobj2$qvalue)

pdf("q-values.pdf", width=8, height=6)
par(mfrow = c(3,2))
plot(qobj$qvalues)
plot(qobj2$qvalues)
plot(qobj3$qvalues)
plot(qobj4$qvalues)
plot(qobj5$qvalues)
dev.off()

#seem fine

#now calculate the median and then backtransform to -log10
qvalues <- qobj$qvalues
qvalues2 <- qobj2$qvalues
qvalues3 <- qobj3$qvalues
qvalues4 <- qobj4$qvalues
qvalues5 <- qobj5$qvalues

combinedq <- data.frame(q1=qvalues, q2=qvalues2, q3=qvalues3, q4=qvalues4, q5=qvalues5)

head(qvalues)
head(qvalues2)

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
#0.4 for all

#visualise correlation
#correlation test
install.packages("ggpubr")
library("ggpubr")
str(combined)
ggscatter(data = combinedq, x = "q1", y = "q2", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "c2-model1", ylab = "c2-model2")
ggscatter(data = combinedq, x = "q1", y = "q3", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "c2-model1", ylab = "c2-model3")



#what is going ooon? :D
#no convergence at all ! \o/
#so do a meadian of the values? ASK REID IF HE'S OKAY WITH THIS! 

#save output model runs

#The object can be summarized and visualized by:
# Open a text file for writing
sink("output_summary-modelruns.txt")

# Run the summary functions and print the results
sum1 <- summary(qobj)
print(sum1)

sum2 <- summary(qobj2)
print(sum2)

sum3 <- summary(qobj3)
print(sum3)

sum4 <- summary(qobj4)
print(sum4)

sum5 <- summary(qobj5)
print(sum5)

# Close the text file
sink()

pdf("sumplots-model-runs.pdf", width=8, height=6)
par(mfrow= c(3,2))
plot(qobj)
plot(qobj2)
plot(qobj3)
plot(qobj4)
plot(qobj5)
dev.off()

pdf("histo-model-runs.pdf", width=8, height=6)
par(mfrow= c(3,2))
hist(qobj)
hist(qobj2)
hist(qobj3)
hist(qobj4)
hist(qobj5)
dev.off()
##############

head(combinedq)
combinedq$median_q <- apply(combinedq[, c("q1", "q2", "q3", "q4", "q5")], 1, median)
range(combinedq$log_transformed)
combinedq$log_transformed <- -log10(combinedq$median_q)

pdf("qvalues-final.pdf", width=6, height=8)
par(mfrow=c(2,1))
plot(combinedq$median_q, xlab="marker", ylab="qvalue")
abline(h = 0.05, col = "red", lty = 2)
plot(combinedq$log_transformed, xlab="marker", ylab= "q(-log10 scale)")
abline(h = 1.30103, col = "red", lty = 2)
dev.off()

```

# recalculate MAF for the new vcf file

module load htslib/1.10.2  bzip2/1.0.8  

GRENEDALF=/gxfs_home/geomar/smomw573/software/grenedalf/bin/grenedalf
VCF=/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/variants/final_20230810.vcf.gz 
GENOME=/gxfs_work1/geomar/smomw573/seasonal_adaptation/genome/affinis.Atlantic.long_read_draft.Mar22.fasta
OUTPUT=/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/variants/

$GRENEDALF frequency --vcf-path $VCF --write-total-frequency --allow-file-overwriting --file-suffix maf-new-vcf-blabla --out-dir $OUTPUT
$GRENEDALF frequency --vcf-path $VCF --write-total-counts --file-suffix -total-counts-new-vcf --out-dir $OUTPUT
$GRENEDALF frequency --vcf-path $VCF --write-sample-alt-freq --file-suffix -alt-freq-new-vcf --out-dir $OUTPUT
$GRENEDALF frequency --vcf-path $VCF --write-sample-ref-freq --file-suffix -ref-freq-new-vcf --out-dir $OUTPUT

# reinstall grenedalf
```bash
module load gcc12-env/12.3.0
module load gcc/12.3.0
module load bzip2/1.0.8
module load xz/5.4.1
module load cmake/3.27.4

rm -rf libdeflate
git clone https://github.com/ebiggers/libdeflate.git
cd libdeflate
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=/gxfs_home/geomar/smomw573/software/test_grenedalf_final/usr ..
make
make install
cd ../..

rm -rf grenedalf
export CMAKE_PREFIX_PATH=/gxfs_home/geomar/smomw573/software/test_grenedalf_final/usr:$CMAKE_PREFIX_PATH
export INCLUDE=/gxfs_home/geomar/smomw573/software/test_grenedalf_final/usr/include:$INCLUDE
git clone --recursive https://github.com/lczech/grenedalf.git
cd grenedalf
make clean
make
```

# 23 January 2024

# okay now I finally got this table, so I will recalculate MAF

```


# open file in R
```r
freq <- read.csv("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/variants/frequencymaf-new-vcf.csv", header=T)

freq$MAF <- ifelse(freq$TOTAL.FREQ > 0.5, 1-freq$TOTAL.FREQ, freq$TOTAL.FREQ)
freq$marker <- cbind(marker = 1:nrow(freq), freq)
write.table(freq, file = "maf-new-vcf-corrected.txt", sep="\t", quote=FALSE, row.names=FALSE)
hist(freq$TOTAL.FREQ,breaks=50)
hist(freq$MAF,breaks=50)

sum(freq$TOTAL.FREQ < 0.01)
sum(freq$TOTAL.FREQ > 0.99)
sum(freq$MAF < 0.01)
sum(freq$MAF > 0.99)

min(freq$MAF)
max(freq$TOTAL.FREQ)
sum(freq$TOTAL.FREQ < 0.01)
#0
sum(freq$TOTAL.FREQ > 0.99)
#0
sum(freq$MAF < 0.01)
#0
sum(freq$MAF > 0.99)
#0
min(freq$MAF)
#0.01
max(freq$MAF)
#0.5
max(freq$TOTAL.FREQ)
#0.99
sum(freq$MAF < 0.5)
#1258257
sum(freq$MAF < 0.05)
#848949




#filter vcf
```bash
python ~/seasonal_adaptation/scripts/filter_sync_by_snplist.py -i /gxfs_work1/geomar/smomw504/seasonal_adaptation/analysis/variants/all.final.sync \
  -snps /gxfs_work1/geomar/smomw504/seasonal_adaptation/analysis/variants/snpdet \
  -o /gxfs_work1/geomar/smomw504/seasonal_adaptation/analysis/variants/filtered.sync

```


# 29 January
# installing poolseq

conda load r_env
R
```R

#install dependencies
install.packages("data.table")
install.packages("foreach")
install.packages("stringi")
install.packages("matrixStats")
install.packages("Rcpp")

library(data.table)
library(foreach)
library(stringi)
library(matrixStats)
library(Rcpp)

#install poolSeq
#wget https://github.com/ThomasTaus/poolSeq/archive/refs/tags/v0.3.5.tar.gz

install.packages("/gxfs_home/geomar/smomw573/software/v0.3.5.tar.gz", repos=NULL, type="source")

library(poolSeq)

#estimated number of generations in each year
#assuming generational turnover of 2-weeks
#2009: 2 generations
#2011: 1 generations
#2015: 4 generations
#2022: 4 generations


```
# pull counts REF ALT alleles from vcf
pre-sync-file.R

```R
library(stringr)
library(R.utils)

gunzip("final_20230810.vcf.gz", remove=FALSE)

setwd("~/Jenny-Eurytemora_affinis/AF-GLM")

#zipF<-file.choose("final_20230810.vcf.gz") # lets you choose a file and save its file path in R (at least for windows)
#outDir<-("~/Jenny-Eurytemora_affinis/AF-GLM") # Define the folder where the zip file should be unzipped to 
#unzip(zipF,exdir=outDir) 

dat <- read.table("final_20230810.vcf", header=FALSE)
head(dat)

#create new dataset with the info reid needs

colnames(dat) <- c("chrom","pos","pres","ref","alt","v6","v7","v8","v9",
                   "2009_t1","2009_t2", "2009_t3", "2009_t4",
                   "2011_t1","2011_t2",
                   "2015_t1", "2015_t2", "2015_t3", "2015_t4",
                   "2022_t1", "2022_t2", "2022_t3", "2022_t4")
filtdat <- dat[c("chrom","pos","ref","alt",
               "2009_t1","2009_t2", "2009_t3", "2009_t4",
               "2011_t1","2011_t2",
               "2015_t1", "2015_t2", "2015_t3", "2015_t4",
               "2022_t1", "2022_t2", "2022_t3", "2022_t4")]
head(filtdat)
columns_to_process <- 5:23

# Loop through each column and update the values
for (col in columns_to_process) {
  # Extract the numbers following the first ":"
  filtdat[, col] <- str_extract(filtdat[, col], ":(\\d+,\\d+):")
  filtdat[, col] <- gsub(":", "", filtdat[, col])
}

filtdat <- dat[c("chrom","pos","ref","alt",
                 "2009_t2",  "2009_t4",
                 "2011_t1","2011_t2",
                 "2015_t1", "2015_t4",
                 "2022_t1", "2022_t4")]
# Print the updated dataset
head(filtdat)
write.table(filtdat, "count_alleles_finalvcf.txt", sep="\t", quote=FALSE, row.names=FALSE)
install.packages("stringr")
install.packages("R.utils")
library(stringr)
library(R.utils)

gunzip("final_20230810.vcf.gz", remove=FALSE)

setwd("~/Jenny-Eurytemora_affinis/AF-GLM")

#zipF<-file.choose("final_20230810.vcf.gz") # lets you choose a file and save its file path in R (at least for windows)
#outDir<-("~/Jenny-Eurytemora_affinis/AF-GLM") # Define the folder where the zip file should be unzipped to 
#unzip(zipF,exdir=outDir) 

dat <- read.table("final_20230810.vcf", header=FALSE)
head(dat)

#create new dataset with the info reid needs

colnames(dat) <- c("chrom","pos","pres","ref","alt","v6","v7","v8","v9",
                   "2009_t1","2009_t2", "2009_t3", "2009_t4",
                   "2011_t1","2011_t2",
                   "2015_t1", "2015_t2", "2015_t3", "2015_t4",
                   "2022_t1", "2022_t2", "2022_t3", "2022_t4")
filtdat <- dat[c("chrom","pos","ref","alt",
               "2009_t1","2009_t2", "2009_t3", "2009_t4",
               "2011_t1","2011_t2",
               "2015_t1", "2015_t2", "2015_t3", "2015_t4",
               "2022_t1", "2022_t2", "2022_t3", "2022_t4")]

#now fix the comma
# Columns 5 onwards
columns_to_process <- 5:ncol(filtdat)
new_dataset <- filtdat[, 1:4]
# Loop through each column and update the values
for (col in columns_to_process) {
  # Convert the column to character
  filtdat[, col] <- as.character(filtdat[, col])
  
  # Split the column based on the comma
  split_values <- strsplit(filtdat[, col], ",")
  
  # Create new columns with the split values
  new_col_name_ref <- paste0(col, "_ref")
  new_col_name_alt <- paste0(col, "_alt")
  
  new_dataset[, new_col_name_ref] <- sapply(split_values, function(x) as.numeric(x[1]))
  new_dataset[, new_col_name_alt] <- sapply(split_values, function(x) as.numeric(x[2]))
}
# Print the updated dataset
head(filtdat)
# Loop through each column and update the values
for (col in columns_to_process) {
  # Split the column based on the comma
  split_values <- strsplit(as.character(filtdat[, col]), ",")
  
  # Create two new columns with the split values
  filtdat[, paste0(col, "REF")] <- sapply(split_values, function(x) x[1])
  filtdat[, paste0(col, "ALT")] <- sapply(split_values, function(x) x[2])
  
  # Remove the original column
  filtdat[, col] <- NULL
}

head(filtdat)
header <- names(new_dataset)
header

colnames(new_dataset) <- c("chrom",  "pos" ,   "ref" ,   "alt"  ,  
                           "2009_t1_ref" , "2009_t1_alt" , "2009_t2_ref" , "2009_t2_alt",  "2009_t3_ref", 
                           "2009_t3_alt",  "2009_t4_ref",  "2009_t4_alt", 
                           "2011_t1_ref" , "2011_t1_alt" , "2011_t2_ref" ,"2011_t2_alt" ,
                           "2015_t1_ref", "2015_t1_alt",
                           "2015_t2_ref", "2015_t2_alt" ,"2015_t3_ref" ,"2015_t3_alt", "2015_t4_ref" ,"2015_t4_alt", 
                           "2022_t1_ref", "2022_t1_alt" ,"2022_t2_ref",
                           "2022_t2_alt", "2022_t3_ref" ,"2022_t3_alt", "2022_t4_ref" ,"2022_t4_alt")
#data I want

final_df <-new_dataset[c("chrom",  "pos" ,   "ref" ,   "alt"  ,  
                         "2009_t2_ref" , "2009_t2_alt",  "2009_t4_ref",  "2009_t4_alt", 
                         "2011_t1_ref" , "2011_t1_alt" , "2011_t2_ref" ,"2011_t2_alt" ,
                         "2015_t1_ref", "2015_t1_alt", "2015_t4_ref" ,"2015_t4_alt", 
                         "2022_t1_ref", "2022_t1_alt" ,"2022_t4_ref", "2022_t4_alt")]
write.table(final_df, "count_alleles.txt", sep="\t", row.names=FALSE, quote=FALSE)

```
# python script to convert to sync file

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

python3 to.py --input count_alleles_finalvcf.txt --output output_sync.txt

# trial poolSeq
#subset file
awk '{print $1 "\t" $2 "\t" $3 "\t" $5 "\t" $7 "\t" $8 "\t" $9 "\t" $10 "\t" $13 "\t" $14 "\t" $17}' output_sync.txt > subset.sync

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




mySync <- read.sync(file="/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/variants/output_sync.txt", polarization="reference", 
  gen=c(1,2,3,4,1,2,1,2,3,4,1,2,3,4),repl=c(1,1,1,1,2,2,3,3,3,3,4,4,4,4))

subSync <- read.sync(file="/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/variants/subset.sync", polarization="reference", 
  gen=c(2,4,1,2,1,4,1,3),repl=c(1,1,2,2,3,3,4,4)) #import only the start/end points of interest

listSync <- list(mySync, subSync)
names(listSync) <- c("mySync", "subSync")
#   af(sync, chr, pos, repl, gen)
# af(mySync, "LS387016.1", 58, 1, 0)
# coverage(mySync, "LS387016.1", 58, 1, 0)
# af.traj(mySync, "LS387016.1", 58, 1)


## estimate effective population size:

af <- as.data.frame(subSync@alleles)
estNe_2009.1 <- estimateNe(p0=af[,7], pt=af[,9], cov0=af[,8], covt=af[,10], ploidy=2, t=2, poolSize=c(50, 50), method="JR.planII")
estNe_2009.2 <- estimateNe(p0=af[,7], pt=af[,9], cov0=af[,8], covt=af[,10], ploidy=2, t=2, poolSize=c(50, 50), method="W.planII")
estNe_2009.3 <- estimateNe(p0=af[,7], pt=af[,9], cov0=af[,8], covt=af[,10], ploidy=2, t=2, poolSize=c(50, 50), method="P.planII")
estNe_2009.4 <- estimateNe(p0=af[,7], pt=af[,9], cov0=af[,8], covt=af[,10], ploidy=2, t=2, poolSize=c(50, 50), method="P.alt.1step.planII")

estNe_2011.1 <- estimateNe(p0=af[,11], pt=af[,13], cov0=af[,12], covt=af[,14],ploidy=2, t=1, poolSize=c(50, 50), method="JR.planII")
estNe_2011.2 <- estimateNe(p0=af[,11], pt=af[,13], cov0=af[,12], covt=af[,14],ploidy=2, t=1, poolSize=c(50, 50), method="P.planII")
estNe_2011.3 <- estimateNe(p0=af[,11], pt=af[,13], cov0=af[,12], covt=af[,14],ploidy=2, t=1, poolSize=c(50, 50), method="W.planII")
estNe_2011.4 <- estimateNe(p0=af[,11], pt=af[,13], cov0=af[,12], covt=af[,14],ploidy=2, t=1, poolSize=c(50, 50), method="P.alt.1step.planII")

estNe_2015.1 <- estimateNe(p0=af[,15], pt=af[,17], cov0=af[,16], covt=af[,18],ploidy=2, t=4, poolSize=c(50, 50), method="JR.planII")
estNe_2015.2 <- estimateNe(p0=af[,15], pt=af[,17], cov0=af[,16], covt=af[,18],ploidy=2, t=4, poolSize=c(50, 50), method="W.planII")
estNe_2015.3 <- estimateNe(p0=af[,15], pt=af[,17], cov0=af[,16], covt=af[,18],ploidy=2, t=4, poolSize=c(50, 50), method="P.planII")
estNe_2015.4 <- estimateNe(p0=af[,15], pt=af[,17], cov0=af[,16], covt=af[,18],ploidy=2, t=4, poolSize=c(50, 50), method="P.alt.1step.planII")

estNe_2022.1 <- estimateNe(p0=af[,19], pt=af[,21], cov0=af[,20], covt=af[,22], ploidy=2,t=3, poolSize=c(50, 50),method="JR.planII")
estNe_2022.2 <- estimateNe(p0=af[,19], pt=af[,21], cov0=af[,20], covt=af[,22], ploidy=2,t=3, poolSize=c(50, 50),method="W.planII")
estNe_2022.3 <- estimateNe(p0=af[,19], pt=af[,21], cov0=af[,20], covt=af[,22], ploidy=2,t=3, poolSize=c(50, 50),method="P.planII")
estNe_2022.4 <- estimateNe(p0=af[,19], pt=af[,21], cov0=af[,20], covt=af[,22], ploidy=2,t=3, poolSize=c(50, 50),method="P.alt.1step.planII")


#subset
af.sub <- af[c(1:30),]
af <- af.sub
estNe_2009.5 <- estimateNe(p0=af[,7], pt=af[,9], cov0=af[,8], covt=af[,10], ploidy=2, t=2, poolSize=c(50, 50), method="JR.planII")
estNe_2009.6 <- estimateNe(p0=af[,7], pt=af[,9], cov0=af[,8], covt=af[,10], ploidy=2, t=2, poolSize=c(50, 50), method="W.planII")
estNe_2009.7 <- estimateNe(p0=af[,7], pt=af[,9], cov0=af[,8], covt=af[,10], ploidy=2, t=2, poolSize=c(50, 50), method="P.planII")
estNe_2009.8 <- estimateNe(p0=af[,7], pt=af[,9], cov0=af[,8], covt=af[,10], ploidy=2, t=2, poolSize=c(50, 50), method="P.alt.1step.planII")

estNe_2011.5 <- estimateNe(p0=af[,11], pt=af[,13], cov0=af[,12], covt=af[,14],ploidy=2, t=1, poolSize=c(50, 50), method="JR.planII")
estNe_2011.6 <- estimateNe(p0=af[,11], pt=af[,13], cov0=af[,12], covt=af[,14],ploidy=2, t=1, poolSize=c(50, 50), method="P.planII")
estNe_2011.7 <- estimateNe(p0=af[,11], pt=af[,13], cov0=af[,12], covt=af[,14],ploidy=2, t=1, poolSize=c(50, 50), method="W.planII")
estNe_2011.8 <- estimateNe(p0=af[,11], pt=af[,13], cov0=af[,12], covt=af[,14],ploidy=2, t=1, poolSize=c(50, 50), method="P.alt.1step.planII")

estNe_2015.5 <- estimateNe(p0=af[,15], pt=af[,17], cov0=af[,16], covt=af[,18],ploidy=2, t=4, poolSize=c(50, 50), method="JR.planII")
estNe_2015.6 <- estimateNe(p0=af[,15], pt=af[,17], cov0=af[,16], covt=af[,18],ploidy=2, t=4, poolSize=c(50, 50), method="W.planII")
estNe_2015.7 <- estimateNe(p0=af[,15], pt=af[,17], cov0=af[,16], covt=af[,18],ploidy=2, t=4, poolSize=c(50, 50), method="P.planII")
estNe_2015.8 <- estimateNe(p0=af[,15], pt=af[,17], cov0=af[,16], covt=af[,18],ploidy=2, t=4, poolSize=c(50, 50), method="P.alt.1step.planII")

estNe_2022.5 <- estimateNe(p0=af[,19], pt=af[,21], cov0=af[,20], covt=af[,22], ploidy=2,t=3, poolSize=c(50, 50),method="JR.planII")
estNe_2022.6 <- estimateNe(p0=af[,19], pt=af[,21], cov0=af[,20], covt=af[,22], ploidy=2,t=3, poolSize=c(50, 50),method="W.planII")
estNe_2022.7 <- estimateNe(p0=af[,19], pt=af[,21], cov0=af[,20], covt=af[,22], ploidy=2,t=3, poolSize=c(50, 50),method="P.planII")
estNe_2022.8 <- estimateNe(p0=af[,19], pt=af[,21], cov0=af[,20], covt=af[,22], ploidy=2,t=3, poolSize=c(50, 50),method="P.alt.1step.planII")



#hand trail to try understand what is going on with the different Ne estimators.
#Not really used for the analysis, but to decided which one to proceed with

# remove SNPs for which Ne cannot/should not be estimated

truncAF=NA
ploidy=2

#terms
p0=c(af$F1.R1.freq)
pt=c(af$F3.R1.freq)
cov0 <- af[,8]
covt <- af[,10]

estimateNe <- function(p0, pt, cov0, covt, t, ploidy=2, truncAF=NA, method="P.planI", Ncensus=NA, poolSize=rep(Ncensus, times=2), asList=FALSE))
  xi <- p0[keep]
  yi <- pt[keep]
  zi = (xi + yi)/2
  n = length(xi)
  # coverage is divided by 2, because later 2*S0, 2*St correction term will be used
  s0 <- cov0[keep]
  st <- covt[keep]

#Waples method - didnt work

Fc <- ((xi-yi)^2)/(zi-xi*yi)
t=2
count = Fc - ((1/cov1) + (1/cov2))
total = sum(count)
keep <- checkSNP(p0, pt, cov0, covt, truncAF=truncAF)

sum <- Fc - (1/s0 + 1/st)
Fc_planII <- (1/n)*sum( Fc - (1/s0 + 1/st))
final <- -t/(ploidy*log(1-Fc_planII))


#Agnes Jonas method (also didnt work)
C0i <- 1/s0 + 1/(2*50) - 1/(s0*2*50)
Cti <- 1/st + 1/(2*50) - 1/(st*2*50)

Ft <- (xi - yi)^2 - (zi-xi*yi)*( C0i + Cti ) / (zi-xi*yi) * (1 - Cti)
      res <- c(res, Np.planII=-t/(ploidy*log(1-Ft)))

#what are the allele freq of positive Fts
sum(Ft > 0)
pospos <- which(Ft > 0)
values_xi <- xi[pospos]
values_yi <- yi[pospos]
#they all vary around 20%
#now I got to figure out why


##table with output
#methods P.alt.1step.planII, W, JR and P
# Initialize an empty data frame
myTable <- data.frame(Variable_Name = character(), Year = character(), Ne_2009 = numeric(), Ne_2011 = numeric(), Ne_2015 = numeric(), Ne_2022 = numeric())

# Loop through years and variables
for (year in c("2009", "2011","2015", "2022")) {
  for (i in 1:4) {
    var_name <- paste0(year, ".", i)
    value <- get(paste0("estNe_", var_name))
    
    # Add to the data frame
    myTable <- rbind(myTable, data.frame(Variable_Name = var_name, Year = year, Value = value))
  }
}
table2 <- myTable[c("Year", "Value")]
write.table(file = "/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/variants/table_output_ne_pooseq.txt", table2, sep = "\t", quote =FALSE)

#estimate variance of p0 and pt for the different years
> var(af$F1.R1.freq)
[1] 0.02710909
> var(af$F3.R1.freq)
[1] 0.02681922
> var(af$F1.R2.freq)
[1] 0.02685312
> var(af$F2.R2.freq)
[1] 0.02686902
> var(af$F1.R3.freq)
[1] 0.02683984
> var(af$F4.R3.freq)
Error in var(af$F4.R3.freq) : 'x' is NULL
> var(af$F3.R3.freq)
[1] 0.0268215
> var(af$F1.R4.freq)
[1] 0.02678547
> var(af$F4.R4.freq)
[1] 0.02678056

> mean(af$F1.R1.freq)
[1] 0.9116324
>
> mean(af$F3.R1.freq)
[1] 0.9115647
> mean(af$F1.R2.freq)
[1] 0.9120435
> mean(af$F2.R2.freq)
[1] 0.911593
> mean(af$F1.R3.freq)
[1] 0.9115756
> mean(af$F4.R3.freq)
[1] NA
Warning message:
In mean.default(af$F4.R3.freq) :
  argument is not numeric or logical: returning NA
> mean(af$F3.R3.freq)
[1] 0.9117126
> mean(af$F1.R4.freq)
[1] 0.9119215
> mean(af$F4.R4.freq)
[1] 0.9113805

 mean(af$F1.R1.cov)
[1] 203.402
> mean(af$F3.R1.cov)
[1] 341.0051
> mean(af$F1.R2.cov)
[1] 265.4084
> mean(af$F2.R2.cov)
[1] 340.0921
> mean(af$F1.R3.cov)
[1] 381.3563
> mean(af$F3.R3.cov)
[1] 397.2056
> mean(af$F1.R4.cov)
[1] 351.3534
> mean(af$F4.R4.cov)
[1] 383.9032


> var(af$F1.R1.cov)
[1] 19737.84
> var(af$F3.R1.cov)
[1] 55106.82
> var(af$F1.R2.cov)
[1] 34870.97
> var(af$F2.R2.cov)
[1] 54461.01
> var(af$F1.R3.cov)
[1] 69002.83
> var(af$F3.R3.cov)
[1] 73365.29
> var(af$F1.R4.cov)
[1] 59558.41
> var(af$F4.R4.cov)
[1] 69708.34


```
# 08.02.2024
## okay, so based on the results with poolSeq, trials and literature research, we decided to follow the analysis using the Jorde 2007 Ne estimator.
## these values will be used as an input for estimating drift, and consequently the significance values for our analysis. 

#estNe_2009
#Njr.planII
#2350.291
#estNe_2011
#Njr.planII 
#384.5897
#estNe_2015
#estNe_2015
#Njr.planII
#13258.17
#estNe_2022
#Njr.planII
#2055.779


Now I have to run the simulations using the Ne values

```R
library(poolSeq)
#make one file for each pop

subSync <- read.sync(file="/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/drift_simulation/subset.sync", polarization="reference", 
  gen=c(2,4,1,2,1,4,1,3),repl=c(1,1,2,2,3,3,4,4)) #import only the start/end points of interest

t1Sync <- read.sync(file="/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/drift_simulation/2009.sync", polarization="reference", 
  gen=c(2,3),repl=c(1,2)) #import only the start/end points of interest

t2Sync <- read.sync(file="/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/drift_simulation/2011.sync", polarization="reference", 
  gen=c(1,2),repl=c(1,2)) #import only the start/end points of interest

t3Sync <- read.sync(file="/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/drift_simulation/2015.sync", polarization="reference", 
  gen=c(1,4),repl=c(1,2)) #import only the start/end points of interest

t4Sync <- read.sync(file="/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/drift_simulation/2022.sync", polarization="reference", 
  gen=c(1,3),repl=c(1,2)) #import only the start/end points of interest


#checking why number of alleles is different across the datasets
nrow(df1)
nrow(df2)
df1 <- subSync@alleles
df2 <- t1Sync@alleles

names1 = set(df1['posID'])
names2 = set(df2['posID'])
head(df1$posID)
names1 <- df1$posID
names2 <- df2$posID
common_names <- intersect(names1, names2)
unique_namesdf1 <- setdiff(names1, names2)
unique_namesdf2 <- setdiff(names2, names1)

unique_namesdf1
str(unique_namesdf1)
name_to_search <- "Scz6wRH_1034;HRSCAF=1162.14997"

# Check if the name is present in df1 and get the corresponding rows
matching_rows <- df1[df1$posID == name_to_search, ]
matching_rows
name_to_search <- "Scz6wRH_9;HRSCAF=61.48321"
name_to_search <- "Scz6wRH_1034;HRSCAF=1162.14997"

# Check if the name is present in df1 and get the corresponding rows
matching_rows <- df1[df1$posID == name_to_search, ]
matching_rows

#fixed alleles

listSync <- list(subSync, t1Sync, t2Sync, t3Sync, t4Sync)
names(listSync) <- c("subSync", "t1Sync", "t2Sync", "t3Sync", "t4Sync")
#   af(sync, chr, pos, repl, gen)
# af(mySync, "LS387016.1", 58, 1, 0)
# coverage(mySync, "LS387016.1", 58, 1, 0)
# af.traj(mySync, "LS387016.1", 58, 1)


# run simulation based on the starting allele frequency of each year (check this correlation)
## Ne is based on above
# mean of AA af

af_2009 <- af(subSync,,, 1, c(2,4))
af_2011 <- af(subSync,,, 2, c(1,2))
af_2015 <- af(subSync,,, 3, c(1,4))
af_2022 <- af(subSync,,, 4, c(1,3))

##



###############################
###############################

# simulation for 2009
#first set the sim.R script in which I simulate the AF

#


```
# join single year simulations into one file
## 15-02-2024
```bash
for i in {1..500}

do {
    cut -f 4 simulated.2009.${i}.sync | cut -f 1,3 -d ":" | sed 's/:/\t/'> ./temp/tmp.2009.1
    cut -f 5 simulated.2009.${i}.sync | cut -f 1,3 -d ":" | sed 's/:/\t/'> ./temp/tmp.2009.2
    cut -f 4 simulated.2011.${i}.sync | cut -f 1,3 -d ":" | sed 's/:/\t/'> ./temp/tmp.2011.1
    cut -f 5 simulated.2011.${i}.sync | cut -f 1,3 -d ":" | sed 's/:/\t/'> ./temp/tmp.2011.2
    cut -f 4 simulated.2015.${i}.sync | cut -f 1,3 -d ":" | sed 's/:/\t/'> ./temp/tmp.2015.1
    cut -f 5 simulated.2015.${i}.sync | cut -f 1,3 -d ":" | sed 's/:/\t/'> ./temp/tmp.2015.2
    cut -f 4 simulated.2022.${i}.sync | cut -f 1,3 -d ":" | sed 's/:/\t/'> ./temp/tmp.2022.1
    cut -f 5 simulated.2022.${i}.sync | cut -f 1,3 -d ":" | sed 's/:/\t/'> ./temp/tmp.2022.2

paste ./temp/tmp.2009.1 ./temp/tmp.2009.2 ./temp/tmp.2011.1 ./temp/tmp.2011.2 ./temp/tmp.2015.1 ./temp/tmp.2015.2 ./temp/tmp.2022.1 ./temp/tmp.2022.2 > ./final_files/rep${i}.counts

}

done
```
# 16.02.2024
## Now I will run the baypass c2 model for all 500 simulations to get a p value
## this will be done with an array script

# baypass-simulations.sh

```bash 
#!/bin/bash
#SBATCH -D /gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/baypass-simulations
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


inputdir=/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/variants
files=/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/drift_simulation/final_files
omega=/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/baypass
outputdir=/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/baypass-simulations
baypassdir=/gxfs_home/geomar/smomw573/software/baypass_public-master/sources

$baypassdir/g_baypass \
-gfile $files/rep${SLURM_ARRAY_TASK_ID}.counts \
-efile $inputdir/cov-baypass.txt \
-omegafile $omega/core-model_mat_omega.out \
-poolsizefile $inputdir/haploid-size-baypass-txt \
-d0yij 20 \
-outprefix $outputdir/c2-model.${SLURM_ARRAY_TASK_ID} -nthreads 10 -contrastfile $inputdir/cov-baypass.txt

```

# 21.02.2024
# rerun simulations cause 10 snps were missing

## now all have all the simulations, so I need to create a new file which has all the p values for all snps in all 500 runs 
file ending with summary_pi_xtx.out are the ones we need
```bash 
#dir baypass-simulations
for i in {1..500}; do
    awk '{print $6}' c2-model.${i}_summary_contrast.out > ./temp/column_${i}.txt
done

#create combined file
cd temp
paste *.txt > trial.txt
#check if number of lines are 500
head -1 trial.txt | tr '\t' '\n' | wc -l
#if I need to know the order of the files
head *txt | grep "^==" > order_pvalues_combinedfile
#500 all good! 

#calculate percentils in R
srun --pty --x11 --nodes=1 --cpus-per-task=4 --mem=32000 --time=06:00:00 /bin/bash
conda activate r_env
R
```

```R
df <- read.table("trial.txt", header = TRUE)
#calculate percentil 2.5 and 97.5 for each row (dont really need this)
# Calculate the percentile 2.5 for each row across the 500 columns
# Calculate the percentiles for each row across the 500 columns
#lower_percentile <- apply(df, 1, function(row) quantile(row, 0.025))
#upper_percentile <- apply(df, 1, function(row) quantile(row, 0.975))

# Add the calculated percentiles as new columns to your data frame
#data$lower_2.5 <- lower_percentile
#data$upper_97.5 <- upper_percentile

# Write the updated data frame back to a file (if needed)
write.table(data, "u.txt", sep="\t", row.names=FALSE)

#now calculate the median for all original 5 runs of the c2 model

setwd("~/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/baypass/baypass/C2-model")

df1 <- read.table("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/baypass/c2-model_summary_contrast.out", header= TRUE)
df2 <- read.table("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/baypass/c2-model2_summary_contrast.out", header= TRUE)
df3 <- read.table("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/baypass/c2.3-model_summary_contrast.out", header= TRUE)
df4 <- read.table("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/baypass/c2.4-model_summary_contrast.out", header= TRUE)
df5 <- read.table("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/baypass/c2.5-model_summary_contrast.out", header= TRUE)
head(df1)

#fetch p values
pval1 <- df1$log10.1.pval
pval2 <- df2$log10.1.pval
pval3 <- df3$log10.1.pval
pval4 <- df4$log10.1.pval
pval5 <- df5$log10.1.pval


combinedp <- data.frame(p1=pval1, p2=pval2, p3=pval3, p4=pval4, p5=pval5)
combinedp$medianp <- apply(combinedp,1,median)
write.table(combinedp, "medianc2model-uncorrected-pvalues.txt", sep="\t", row.names=FALSE)
write.table(combinedp, "pvalues-with-simcounts.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
write.table(result, "final-list-sig-snps.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

#import data with 500 simluations per SNP
# Read the data from the file
dfsim <- read.table("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/baypass-simulations/temp/trial.txt", header = TRUE)
combinedp <- read.table("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/baypass-simulations/temp/medianc2model-uncorrected-pvalues.txt", sep="\t", header=TRUE)

#okay no just quickly try empVal function
library(devtools)
library(qvalue)

#create my list
# Create a list with dfsim as stat0 and combinedp$medianp as stat

# Use my_data in the empPvals function

dfsim <- as.matrix(dfsim)

pvalues <- empPvals(stat = combinedp$medianp, stat0 = dfsim, pool=TRUE)
qobj <- qvalue(p = pvalues)
qvalues <- qobj$qvalues
#1678
pi0 <- qobj$pi0
lfdr <- qobj$lfdr
summary(qobj)
#48705
hist(qobj)
plot(qobj)
write.table(qobj$pvalues,"sig-snps-empval.txt", row.names=FALSE, quote=FALSE)

pvalues2 <- empPvals(stat = combinedp$medianp, stat0 = dfsim, pool=FALSE)
qobj2 <- qvalue(p = pvalues2)
qvalues2 <- qobj2$qvalues2
#1678
pi02 <- qobj2$pi02
lfdr2 <- qobj2$lfdr2
summary(qobj2)
#50582
write.table(qobj2$pvalues,"sig-snps-empval-falsepool.txt", row.names=FALSE, quote=FALSE)

finalpvals <- data.frame(combinedp$medianp, qobj2$pvalues, qobj$pvalues)

write.table(finalpvals$pool, "/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/baypass-simulations/temp/pval_pool.txt", row.names=FALSE, quote=FALSE)
write.table(finalpvals$falsepool, "/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/baypass-simulations/temp/pval_falsepool.txt", row.names=FALSE, quote=FALSE)
write.table(finalpvals$qvaluesfalse, "/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/baypass-simulations/temp/qvalue_false.txt", row.names=FALSE, quote=FALSE)


#okay now I need to pull out the snps - restart
#we want table false pool
chrom <- read.table(file="/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/GLM/chroms.txt", header= T)
pvaluesfalsepool <- read.table("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/baypass-simulations/temp/pval_falsepool.txt", header=TRUE)
qvaluesfalsepool <- read.table("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/baypass-simulations/temp/qvalue_false.txt", header=TRUE)
pvaluespool <- read.table("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/baypass-simulations/temp/pval_pool.txt", header=TRUE)
freq <- read.table(file="\\\\helmholtz.geomar.de/Users$/jnascimento/Daten/Jenny-Eurytemora_affinis/AF-GLM/freq_new_vcf.txt", header= T)

combined <- cbind(chrom, pvaluesfalsepool, qvaluesfalsepool, pvaluespool)
names(combined) <- c("chrom", "pvalues_false", "qvalues_false", "pvalues_pool")
write.table(combined, "/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/baypass-simulations/temp/final-significance-simulations.txt", row.names=FALSE, quote=FALSE)
#significant snps
#sigp false
sigpfalse <- combined[combined$pvalues_false < 0.05, ]

#sigq false
sigqfalse <- combined[combined$qvalues_false < 0.05, ]

```
# for the section above I used figures-snps-baypass-simulations.R in my personal computer
# 06.03.2024
# compare p values from 5 models to the one generated by simulations
# personal computer again
# also ran script ploting-SNP-density.R


# cheking out what AF of snps that are significant in one set (simulations) and not in baypass runs, to try to figure out why

srun --pty --x11 --nodes=1 --cpus-per-task=4 --mem=32000 --time=06:00:00 /bin/bash
conda activate r_env

```R
library(dplyr)
library(ggplot2)
library(tidyr)

df <- read.table("blabla.txt", header=TRUE)
dfsig <- df %>% filter(df$pvalues_false < 0.05)
nrow(dfsig)
dfsig$logsimp <- -log10(dfsig$pvalues_false)
plot(dfsig$logsimp, dfsig$combinedp)
colors <- ifelse(dfsig$combinedp > 1.30103, "red", "black")
sum(dfsig$combinedp > 1.30103)
plot(dfsig$combinedp, dfsig$logsimp, col=colors)

colors <- ifelse(dfsig$combinedp > 1.30103, "red", "black")
sum(dfsig$combinedp > 1.30103)
plot(dfsig$logsimp, dfsig$combinedp, col=colors)
dfsig$pvalbaypass <- 10^(-dfsig$combinedp)
sum(dfsig$pvalbaypass < 0.05)
wonkysnps <- dfsig %>% filter(pvalbaypass > 0.05)
nrow(wonkysnps)
str(wonkysnps)

write.table(wonkysnps, file="/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/variants/incoherent-snp-list.txt", col.names=TRUE, quote=FALSE)

#now check the AF and coverage of these snps

cov <- read.table("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/GLM/cov_new_vcf.txt", header = TRUE, sep = "\t")
freq <- read.table("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/GLM/freq_new_vcf.txt", header = TRUE, sep = "\t")
pos <- read.table("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/GLM/chrom-poly-position-seasonal-data.txt", header = FALSE, sep = "\t")

covfinal <- cov[, c("EA_2009_T2", "EA_2009_T4", "EA_2011_T1", "EA_2011_T2", "EA_2015_T1", "EA_2015_T4", "EA_2022_T1", "EA_2022_T4")]
freqfinal <- freq[, c("EA_2009_T2", "EA_2009_T4", "EA_2011_T1", "EA_2011_T2", "EA_2015_T1", "EA_2015_T4", "EA_2022_T1", "EA_2022_T4")]
str(freqfinal)
freqfinal$chrom <- pos$V1
freqfinal$pos <- pos$V2
covfinal$chrom <- pos$V1
covfinal$pos <- pos$V2
str(freqfinal)
str(dfsig)

covfinal$chrompos <- paste(covfinal$chrom, covfinal$pos, sep = "_")
freqfinal$chrompos <- paste(freqfinal$chrom, freqfinal$pos, sep = "_")

str(covfinal)

dfsig$chrompos <- dfsig$chrom
str(wonkysnps)
wonkysnps$chrompos <- wonkysnps$chrom
wonkyfreq <- merge(freqfinal, wonkysnps, by = "chrompos", all.x = FALSE)

common_cols <- intersect(colnames(wonkyfreq), colnames(freqfinal))

# Subset wankyfreq using common columns
wonkyfreq <- wonkyfreq[, common_cols]
str(wonkyfreq)
write.table(wonkyfreq, file="/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/variants/incoherent-snp-freq.txt", col.names=TRUE, quote=FALSE)

#plot what was going on with the snps
changes <- wonkyfreq %>%
  mutate(
    '2009.af' = `EA_2009_T4` - `EA_2009_T2`,
    '2011.af' = `EA_2011_T2` - `EA_2011_T1`,
    '2015.af' = `EA_2015_T4` - `EA_2015_T1`,
    '2022.af' = `EA_2022_T1` - `EA_2022_T1`
  )
changes <- select(changes, chrompos, '2009.af', '2011.af', '2015.af', '2022.af')
long_changes <- changes %>%
  pivot_longer(cols = -chrompos, 
               names_to = "Year", 
               values_to = "AF")
plot(long_changes$AF)
par(mfcol = c(1, 2))
hist(long_changes$AF)
mean(long_changes$AF)
#-0.0008932972
sd(long_changes$AF)
#0.0338069

#check initial frequency
start_values <- select(wonkyfreq, chrompos, 'EA_2009_T2', 'EA_2011_T1', 'EA_2015_T1', 'EA_2022_T1')
long_start_values <- start_values %>%
  pivot_longer(cols = -chrompos, 
               names_to = "Yearhist", 
               values_to = "AF")

history()
long_wonkyfreq <- wonkyfreq %>%
  pivot_longer(cols = -chrompos, 
               names_to = "Year", 
               values_to = "AF")

#now I will check the variation in AF of the significant ones
nonwonkysnps <- dfsig %>% filter(pvalbaypass < 0.05)
nonwonkysnps$chrompos <- nonwankysnps$chrom
nonwonkyfreq <- merge(freqfinal, nonwonkysnps, by = "chrompos", all.x = FALSE)

common_cols <- intersect(colnames(nonwonkyfreq), colnames(freqfinal))

# Subset wankyfreq using common columns
nonwonkyfreq <- nonwonkyfreq[, common_cols]
str(wonkyfreq)
write.table(wonkyfreq, file="/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/variants/incoherent-snp-freq.txt", col.names=TRUE, quote=FALSE)

#plot what was going on with the snps
changes2 <- nonwonkyfreq %>%
  mutate(
    '2009.af' = `EA_2009_T4` - `EA_2009_T2`,
    '2011.af' = `EA_2011_T2` - `EA_2011_T1`,
    '2015.af' = `EA_2015_T4` - `EA_2015_T1`,
    '2022.af' = `EA_2022_T1` - `EA_2022_T1`
  )
changes2 <- select(changes2, chrompos, '2009.af', '2011.af', '2015.af', '2022.af')
long_changes2 <- changes2 %>%
  pivot_longer(cols = -chrompos, 
               names_to = "Year", 
               values_to = "AF")
hist(long_changes2$AF)
plot(long_changes2$AF)
mean(long_changes2$AF)
#-0.001005056
sd(long_changes2$AF)
#0.03758289

#check initial frequency
start_values2 <- select(nonwonkyfreq, chrompos, 'EA_2009_T2', 'EA_2011_T1', 'EA_2015_T1', 'EA_2022_T1')
long_start_values2 <- start_values2 %>%
  pivot_longer(cols = -chrompos, 
               names_to = "Year", 
               values_to = "AF")

history()
long_wonkyfreq <- wonkyfreq %>%
  pivot_longer(cols = -chrompos, 
               names_to = "Year", 
               values_to = "AF")

 hist(long_start_values$AF)
 hist(long_start_values2$AF)


##############now same for coverage

wonkycov<- merge(covfinal, wonkysnps, by = "chrompos", all.x = FALSE)

common_cols <- intersect(colnames(wankycov), colnames(covfinal))

# Subset wankyfreq using common columns
wonkycov <- wonkycov[, common_cols]
str(wonkycov)
write.table(wonkycov, file="/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/variants/incoherent-snp-cov.txt", col.names=TRUE, quote=FALSE)

#plot what was going on with the snps

long_cov <- wonkycov %>%
  pivot_longer(cols = -chrompos, 
               names_to = "Year", 
               values_to = "COV")

hist(long_cov$COV)
  mean(long_cov$COV)
#254.3453
sd(long_cov$COV)
#163.7068

#now I will check the variation in cov of the significant ones
nonwonkysnps <- dfsig %>% filter(pvalbaypass < 0.05)
nonwonkysnps$chrompos <- nonwonkysnps$chrom
nonwonkycov <- merge(covfinal, nonwonkysnps, by = "chrompos", all.x = FALSE)

common_cols <- intersect(colnames(nonwonkyfreq), colnames(covfinal))

# Subset wankyfreq using common columns
nonwonkycov <- nonwonkycov[, common_cols]
long_cov2 <- nonwonkycov %>%
  pivot_longer(cols = -chrompos, 
               names_to = "Year", 
               values_to = "COV")

hist(long_cov2$COV)

#check initial cov
start_values <- select(wonkycov, chrompos, 'EA_2009_T2', 'EA_2011_T1', 'EA_2015_T1', 'EA_2022_T1')
long_start_values <- start_values %>%
  pivot_longer(cols = -chrompos, 
               names_to = "Year", 
               values_to = "COV")
hist(long_start_values$COV)


#plot what was going on with the snps
start_values2  <- select(nonwonkycov, chrompos, 'EA_2009_T2', 'EA_2011_T1', 'EA_2015_T1', 'EA_2022_T1')
long_start_values2 <- start_values2 %>%
  pivot_longer(cols = -chrompos, 
               names_to = "Year", 
               values_to = "COV")
hist(long_start_values2$COV)

#okay so now I will select the 1K fluctuating snps, and check which ones have high p values in the fluctuating plot


fluctsnps <- read.table("chrompos_fluctuatingsnps.txt", header=TRUE)
colnames(fluctsnps) <- c("chrompos")
#set of highly significant in baypass
dfsigsim <- df %>% filter(df$pvalues_false < 0.05)
df$logsimp <- -log10(df$pvalues_false)

plot(df$logsimp, df$combinedp)
plot(dfsig$logsimp, dfsig$combinedp)
colors <- ifelse(dfsig$combinedp > 2, "red", "black")
sum(dfsig$combinedp > 2.5)
plot(dfsig$combinedp, dfsig$logsimp, col=colors)
plot(df$logsimp, df$combinedp, col=colors)
abline(h = 2, col = "red", lty = 2)
abline(v = 1.30103, col = "red", lty = 2)
#now filter for the ones that are highly significant in one but not in other

filtered_df <- subset(df, combinedp > 2 & logsimp < 1.30103)
#lets see what these snps are doing!
#merge with AF data
merged <- merge(filtered_df, freqfinal, by.x = "chrom", by.y = "chrompos")
common_cols <- intersect(colnames(merged), colnames(freqfinal))

# Subset wankyfreq using common columns
merged <- merged[, common_cols]

colnames(merged) <- c("chrompos","2009.start", "2009.end",
                  "2011.start", "2011.end",
                  "2015.start", "2015.end",
                  "2022.start", "2022.end"
)

trial <- gather(merged, key = "year", value = "af", -chrompos)
levels <- c("2009.start", "2009.end", "2011.start", "2011.end", "2015.start", "2015.end", "2022.start", "2022.end")
trial$year <- factor(trial$year, levels = levels)
factorsnames <- c("2009.start", "2009.end",
                  "2011.start", "2011.end",
                  "2015.start", "2015.end",
                  "2022.start", "2022.end"
)
d <- ggplot(data = trial, aes(x = year, y = af, color = chrompos)) +
  geom_line(aes(group = chrompos), size = 1.5, alpha = 0.5) + # Increased line size
  labs(x = "Time", y = "AF") +  # Removed color label
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    plot.title = element_text(size = 16)
  ) +
  scale_x_discrete(labels = factorsnames) +
  guides(color = "none")  # Remove legend
d


start_columns <- merged[, c("chrompos", grep("\\.start", names(merged), value = TRUE))]
long_start_values3 <- start_columns %>%
  pivot_longer(cols = -chrompos, 
               names_to = "Year", 
               values_to = "AF")
hist(long_start_values3$AF)

#ok so now I will merge the highly significant ones in baypass 

colors <- ifelse(df$combinedp > 2 & df$logsimp < 1.30103, "red", "black")
plot(df$logsimp, df$combinedp, col=colors)

#those are the ones I'll be looking at
merged <- merge(filtered_df,freqfinal, by="chrompos")

#okay, now the ones which are highly significant in baypass and slightly significant in simulations

filtered_df2 <- subset(df, combinedp > 2 & logsimp > 1.30103)
colors <- ifelse(df$combinedp > 3 & df$logsimp > 1.30103 & df$logsimp < 2.5, "red", "black")
sum(dfsig$combinedp > 2.5)
plot(df$logsimp,df$combinedp,  col=colors)




#ok from those I want to color them in nroAF 1K plot! 
filtered_df2 <- subset(df, combinedp > 3 & logsimp > 1.30103 & logsimp < 2.5)
#lets see what these snps are doing!
#merge with AF data
merged2 <- merge(filtered_df2, freqfinal, by.x = "chrom", by.y = "chrompos")
common_cols <- intersect(colnames(merged2), colnames(freqfinal))


merged2 <- merged2[, common_cols]

colnames(merged2) <- c("chrompos","2009.start", "2009.end",
                  "2011.start", "2011.end",
                  "2015.start", "2015.end",
                  "2022.start", "2022.end"
)


#ok now merge with the 1K snps to see how many are significant 
merged3 <- merge(merged2, fluctsnps, by="chrompos")
#ok, only 793 of these are in the 1K dataset
#those are the ones which did not merge! 
non_merged_snps <- anti_join(fluctsnps, merged3, by = "chrompos")
sum(df, combinedp > 3 & logsimp > 1.30103)
dfsig$pvalbaypass <- 10^(-dfsig$combinedp)
sum(dfsig$pvalbaypass < 0.05)



#forget all what ive tried so far. lets start fresh



#ok, now I am just looking at the ones which are highly significant
filtered_df2 <- subset(df, combinedp > 3)
#lets see what these snps are doing!
#merge with AF data
merged2 <- merge(filtered_df2, freqfinal, by.x = "chrom", by.y = "chrompos")
common_cols <- intersect(colnames(merged2), colnames(freqfinal))


merged2 <- merged2[, common_cols]

colnames(merged2) <- c("chrompos","2009.start", "2009.end",
                  "2011.start", "2011.end",
                  "2015.start", "2015.end",
                  "2022.start", "2022.end"
)

#901 snps with pvalues <0.001


#ok now merge with the 1K snps to see how many are significant 
merged3 <- merge(merged2, fluctsnps, by="chrompos")
#ok, only 254 something of these are in the 1K dataset
#those are the ones which did not merge! 
#meaning they are highly significant in baypass, but are not part of the 1K fluctuating snps!) 
non_merged_snps <- anti_join(merged2, fluctsnps, by = "chrompos")
#600something snps

#now I can plot the 1K snps, and colour these funky ones! 
#first I will plot the ones that did not merge
trial2 <- gather(non_merged_snps, key = "year", value = "af", -chrompos)
levels <- c("2009.start", "2009.end", "2011.start", "2011.end", "2015.start", "2015.end", "2022.start", "2022.end")
trial2$year <- factor(trial2$year, levels = levels)
factorsnames <- c("2009.start", "2009.end",
                  "2011.start", "2011.end",
                  "2015.start", "2015.end",
                  "2022.start", "2022.end"
)
d <- ggplot(data = trial2, aes(x = year, y = af, color = chrompos)) +
  geom_line(aes(group = chrompos), size = 1.5, alpha = 0.5) + # Increased line size
  labs(x = "Time", y = "AF") +  # Removed color label
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    plot.title = element_text(size = 16)
  ) +
  scale_x_discrete(labels = factorsnames) +
  guides(color = "none")  # Remove legend
d

trial3 <- gather(merged3, key = "year", value = "af", -chrompos)
levels <- c("2009.start", "2009.end", "2011.start", "2011.end", "2015.start", "2015.end", "2022.start", "2022.end")
trial3$year <- factor(trial3$year, levels = levels)
d <- ggplot(data = trial3, aes(x = year, y = af, color = chrompos)) +
  geom_line(aes(group = chrompos), size = 1.5, alpha = 0.5) + # Increased line size
  labs(x = "Time", y = "AF") +  # Removed color label
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    plot.title = element_text(size = 16)
  ) +
  scale_x_discrete(labels = factorsnames) +
  guides(color = "none")  # Remove legend
d

#proper plot
merged3_new <- merged3

# Exclude the "chrompos" column and apply the transformation
merged3_new[, -1] <- 1 - merged3_new[, -1]
norm <- merged3 %>%
  mutate(
    '2009.end' = `2009.end` - `2009.start`,
    '2011.end' = `2011.end` - `2011.start`,
    '2015.end' = `2015.end` - `2015.start`,
    '2022.end' = `2022.end` - `2022.start`
  )

#set start of the year to zero
norm <- norm %>%
  mutate(
    '2009.start' = 0,
    '2011.start' = 0,
    '2015.start' = 0,
    '2022.start' = 0
  )

#for neg snps
normdec <- merged3_new %>%
  mutate(
    '2009.end' = `2009.end` - `2009.start`,
    '2011.end' = `2011.end` - `2011.start`,
    '2015.end' = `2015.end` - `2015.start`,
    '2022.end' = `2022.end` - `2022.start`
  )

#set start of the year to zero
normdec <- normdec %>%
  mutate(
    '2009.start' = 0,
    '2011.start' = 0,
    '2015.start' = 0,
    '2022.start' = 0
  )
  #positive
result <- norm %>%
  mutate(across(starts_with("2009.") | starts_with("2011.") | starts_with("2015.") | starts_with("2022."), ~ . > 0))
head(result)

result2 <- result %>%
  filter_at(vars(ends_with(".end")), all_vars(. == TRUE))
nrow(result2)
#40 SNPs
noclue <- merge(result2, merged3, by = "chrompos", suffixes = c(".result", ".af"))
nrow(merged)
head(merged)
positivesnps <- noclue %>%
  select(chrompos, ends_with(".af"))
head(positivesnps)

#now I will correct them all to 2009
list <- c("2009.start.af", "2009.end.af", "2011.start.af", "2011.end.af", "2015.start.af", "2015.end.af", "2022.start.af", "2022.end.af")

#positive <- positivesnps %>%
#  mutate(across(all_of(list), ~ . - `2009.start.af`))
#head(positive)
#head(final)
positive <- positivesnps
finalpost <- gather(positive, key = "year", value = "af", -chrompos)
levels <- c("2009.start.af", "2009.end.af", "2011.start.af", "2011.end.af", "2015.start.af", "2015.end.af", "2022.start.af", "2022.end.af")
finalpost$year <- factor(finalpost$year, levels = levels)

factorsnames <- c("2009.start", "2009.end",
                  "2011.start", "2011.end",
                  "2015.start", "2015.end",
                  "2022.start", "2022.end"
)

d <- ggplot(data = finalpost, aes(x = year, y = af)) +
  geom_line(aes(group = chrompos), color = "black", size = 0.8, alpha = 0.1) +
  labs(x = "Time", y = "AF", color = "Year") +
  theme_classic()+
  theme(
    axis.title.x = element_text(size = 14),   # Adjust x-axis label size
    axis.title.y = element_text(size = 14),   # Adjust y-axis label size
    plot.title = element_text(size = 16)     # Adjust plot title size
  ) +
  scale_x_discrete(labels = factorsnames)
#stat_summary(fun=median, color="dodgerblue3", geom = "line", group=1, size=1.5)
d



#now negative ones

#select positive snps
#inverse AF with af3
resultdec <- normdec %>%
  mutate(across(starts_with("2009.") | starts_with("2011.") | starts_with("2015.") | starts_with("2022."), ~ . > 0))
head(resultdec)
result2dec <- resultdec %>%
  filter_at(vars(ends_with(".end")), all_vars(. == TRUE))
nrow(result2dec)
#3207 SNPs
#flip to positive

mergeddec <- merge(result2dec, merged2_new, by = "chrompos", suffixes = c(".result", ".af"))
nrow(mergeddec)
head(mergeddec)
positivesnpsdec <- mergeddec %>%
  select(chrompos, ends_with(".af"))
head(positivesnpsdec)

#now I will correct them all to 2009
list <- c("2009.start.af", "2009.end.af", "2011.start.af", "2011.end.af", "2015.start.af", "2015.end.af", "2022.start.af", "2022.end.af")

positivedec <- positivesnpsdec %>%
  mutate(across(all_of(list), ~ . - `2009.start.af`))
head(positivedec)


finalpostdec <- gather(positivedec, key = "year", value = "af", -chrompos)
levels <- c("2009.start.af", "2009.end.af", "2011.start.af", "2011.end.af", "2015.start.af", "2015.end.af", "2022.start.af", "2022.end.af")
finalpostdec$year <- factor(finalpostdec$year, levels = levels)

factorsnames <- c("2009.start", "2009.end",
                  "2011.start", "2011.end",
                  "2015.start", "2015.end",
                  "2022.start", "2022.end"
)

df_combined <- rbind(finalpostdec, finalpost)

d <- ggplot(data = df_combined, aes(x = year, y = af)) +
  geom_line(aes(group = chrompos), color = "black", size = 0.8, alpha = 0.1) +
  labs(x = "Time", y = "AF", color = "Year") +
  theme_classic()+
  theme(
    axis.title.x = element_text(size = 14),   # Adjust x-axis label size
    axis.title.y = element_text(size = 14),   # Adjust y-axis label size
    plot.title = element_text(size = 16)     # Adjust plot title size
  ) +
  scale_x_discrete(labels = factorsnames) +
  stat_summary(fun="mean", color="dodgerblue3", geom = "line", group=1, size=1.5)
d
ggsave("af-normalizedperyear-all-pvalues.pdf",d, w=7, h=5)
sum(df, combinedp > 3 & logsimp > 1.30103)
dfsig$pvalbaypass <- 10^(-dfsig$combinedp)
sum(dfsig$pvalbaypass < 0.05)
#ok now look at the behaviour of these snps! 


#create COV and AF datasets for those snps
covfluct <- merge(fluctsnps, sigpvals, by = "chrompos", suffixes = c(".result", ".af"))
affluct








df<- read.table("final-significance-simulations.txt", header=TRUE)
freq <- read.table(file="\\\\helmholtz.geomar.de/Users$/jnascimento/Daten/Jenny-Eurytemora_affinis/AF-GLM/freq_new_vcf.txt", header= T)
str(df)
str(freq)
head(df)
sum(df$pvalues_false <= 0.004)
#
af <- freq[, c("EA_2009_T2","EA_2009_T4" ,"EA_2011_T1", "EA_2011_T2","EA_2015_T1", "EA_2015_T4", "EA_2022_T1", "EA_2022_T4")]
af2 <- cbind(df$chrom, af)
af3 <- 1-af
afdec <- cbind(df$chrom, af3)

names(af2)[1] <- c("chrom")
names(afdec)[1] <- c("chrom")
str(af2)
#sigp false
sigpfalse <- df[df$pvalues_false < 0.05, ]

#sigq false
sigqfalse <- df[df$qvalues_false < 0.05, ]

#AF of pvalue significant snps
str(af2)
sigpvals <- merge(sigpfalse, af2, by = "chrom")
sigqvals <- merge(sigqfalse, af2, by = "chrom")

#for decreasing snps
sigpvalsdec <- merge(sigpfalse, afdec, by = "chrom")
sigqvalsdec <- merge(sigqfalse, afdec, by = "chrom")
write.table(sigpvals,"significant_pvalues_simulations.txt", sep="\t", quote=FALSE)
write.table(sigqvals,"significant_qvalues_simulations.txt", sep="\t", quote=FALSE)

#start with pvalues
sigpvals <- sigpvals[, names(af2)]
sigqvals <- sigqvals[, names(af2)]
sigpvalsdec <- sigpvalsdec[, names(afdec)]
sigqvalsdec <- sigqvalsdec[, names(afdec)]

colnames(sigpvals) <- c("chrompos", "2009.start", "2011.start",
                        "2015.start", "2022.start",
                        "2009.end", "2011.end",
                        "2015.end", "2022.end")
colnames(sigqvals) <- c("chrompos", "2009.start", "2011.start",
                        "2015.start", "2022.start",
                        "2009.end", "2011.end",
                        "2015.end", "2022.end")
colnames(sigpvalsdec) <- c("chrompos", "2009.start", "2011.start",
                           "2015.start", "2022.start",
                           "2009.end", "2011.end",
                           "2015.end", "2022.end")
colnames(sigqvalsdec) <- c("chrompos", "2009.start", "2011.start",
                           "2015.start", "2022.start",
                           "2009.end", "2011.end",
                           "2015.end", "2022.end")

#okay, fisrt visualise what each snp is doing normally
#say 2009.start is zero

#add row with numbers so I can ID SNPs


#first check how much they changed
norm <- sigqvals %>%
  mutate(
    '2009.end' = `2009.end` - `2009.start`,
    '2011.end' = `2011.end` - `2011.start`,
    '2015.end' = `2015.end` - `2015.start`,
    '2022.end' = `2022.end` - `2022.start`
  )

#set start of the year to zero
norm <- norm %>%
  mutate(
    '2009.start' = 0,
    '2011.start' = 0,
    '2015.start' = 0,
    '2022.start' = 0
  )

#first check how much they changed
normdec <- sigqvalsdec %>%
  mutate(
    '2009.end' = `2009.end` - `2009.start`,
    '2011.end' = `2011.end` - `2011.start`,
    '2015.end' = `2015.end` - `2015.start`,
    '2022.end' = `2022.end` - `2022.start`
  )

#set start of the year to zero
normdec <- normdec %>%
  mutate(
    '2009.start' = 0,
    '2011.start' = 0,
    '2015.start' = 0,
    '2022.start' = 0
  )
#select positive snps

result <- norm %>%
  mutate(across(starts_with("2009.") | starts_with("2011.") | starts_with("2015.") | starts_with("2022."), ~ . > 0))
head(result)

result2 <- result %>%
  filter_at(vars(ends_with(".end")), all_vars(. == TRUE))
nrow(result2)
#3278 SNPs
merged <- merge(result2, sigpvals, by = "chrompos", suffixes = c(".result", ".af"))
nrow(merged)
head(merged)
positivesnps <- merged %>%
  select(chrompos, ends_with(".af"))
head(positivesnps)

#now I will correct them all to 2009
list <- c("2009.start.af", "2009.end.af", "2011.start.af", "2011.end.af", "2015.start.af", "2015.end.af", "2022.start.af", "2022.end.af")

#positive <- positivesnps %>%
#  mutate(across(all_of(list), ~ . - `2009.start.af`))
#head(positive)
#head(final)
positive <- positivesnps
finalpost <- gather(positive, key = "year", value = "af", -chrompos)
levels <- c("2009.start.af", "2009.end.af", "2011.start.af", "2011.end.af", "2015.start.af", "2015.end.af", "2022.start.af", "2022.end.af")
finalpost$year <- factor(finalpost$year, levels = levels)

factorsnames <- c("2009.start", "2009.end",
                  "2011.start", "2011.end",
                  "2015.start", "2015.end",
                  "2022.start", "2022.end"
)

d <- ggplot(data = finalpost, aes(x = year, y = af)) +
  geom_line(aes(group = chrompos), color = "black", size = 0.8, alpha = 0.1) +
  labs(x = "Time", y = "AF", color = "Year") +
  theme_classic()+
  theme(
    axis.title.x = element_text(size = 14),   # Adjust x-axis label size
    axis.title.y = element_text(size = 14),   # Adjust y-axis label size
    plot.title = element_text(size = 16)     # Adjust plot title size
  ) +
  scale_x_discrete(labels = factorsnames)
#stat_summary(fun=median, color="dodgerblue3", geom = "line", group=1, size=1.5)
d
ggsave("af-normalizedperyear-pos.pdf",d, w=7, h=5)


#now negative ones

#select positive snps
#inverse AF with af3
resultdec <- normdec %>%
  mutate(across(starts_with("2009.") | starts_with("2011.") | starts_with("2015.") | starts_with("2022."), ~ . > 0))
head(resultdec)
result2dec <- resultdec %>%
  filter_at(vars(ends_with(".end")), all_vars(. == TRUE))
nrow(result2dec)

```

# 23.03.2024
# filtering vcf for MAF 0.05 to see how the SNPs look like!

```bash
conda create -n vcf

vcftools --gzvcf /gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/variants/final_20230810.vcf.gz  --maf 0.05 --recode --out /gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/variants/final_20230810_maf05.vcf.gz
```






# estimate pi with grenedalf

# run job 02.02.2024

#!/bin/bash
#SBATCH -D /gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/
#SBATCH --mail-type=END
#SBATCH --mail-user=jnascimento@geomar.de
#SBATCH --partition=base
#SBATCH --nodes=1
#SBATCH --tasks-per-node=2
#SBATCH --cpus-per-task=1
#SBATCH --mem=50G
#SBATCH --time=48:00:00
#SBATCH --job-name=bp-c2
#SBATCH --output=genome-diversity.out
#SBATCH --error=genome-diversity.err

module load gcc12-env/12.3.0
module load gcc/12.3.0
module load bzip2/1.0.8
module load xz/5.4.1
module load cmake/3.27.4

GRENEDALF=/gxfs_home/geomar/smomw573/software/grenedalf/bin/grenedalf
MPILEUP=/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/variants/all.mpileup-excluding20072011
GENOME=/gxfs_work1/geomar/smomw573/seasonal_adaptation/genome/Eaffinis.Baltic.PseudoRef.Mar22.fasta.gz 
OUTPUT=/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/variants/


$GRENEDALF diversity --pileup-path $MPILEUP --filter-sample-min-count 1 --reference-genome-fasta-file $GENOME --allow-file-overwriting --window-type genome --pool-sizes 50 --measure theta-pi --file-prefix diversity- --out-dir $OUTPUT



# run job 05.04.2024
# calculate pi with windows (10k)

#!/bin/bash
#SBATCH -D /gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/
#SBATCH --mail-type=END
#SBATCH --mail-user=jnascimento@geomar.de
#SBATCH --partition=base
#SBATCH --nodes=1
#SBATCH --tasks-per-node=2
#SBATCH --cpus-per-task=1
#SBATCH --mem=50G
#SBATCH --time=48:00:00
#SBATCH --job-name=bp-c2
#SBATCH --output=genome-diversity.out
#SBATCH --error=genome-diversity.err

module load gcc12-env/12.3.0
module load gcc/12.3.0
module load bzip2/1.0.8
module load xz/5.4.1
module load cmake/3.27.4

GRENEDALF=/gxfs_home/geomar/smomw573/software/grenedalf/bin/grenedalf
MPILEUP=/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/variants/all.mpileup-excluding20072011
GENOME=/gxfs_home/geomar/smomw573/work/seasonal_adaptation/genome/Eaffinis.Baltic.PseudoRef.Mar22.fasta.gz 
OUTPUT=/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/variants/


$GRENEDALF diversity --pileup-path $MPILEUP --filter-sample-min-count 1 --reference-genome-fasta-file $GENOME --allow-file-overwriting --window-type sliding  --window-sliding-width 10000 --pool-sizes 50 --measure theta-pi --file-prefix diversity-10kwindow- --out-dir $OUTPUT





























































# poolfstat for new vcf

```R
/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/variants/final_20230810.vcf.gz
cd 

#poolfstat for vcf

library(poolfstat)

#import vcf file
ea.readcount30X <- vcf2pooldata(vcf.file="depth-corrected.recode.vcf.gz",poolsizes=rep(50,14))
# summary of the resulting pooldata object
ea.readcount30X

#ha?
selected.snps.idx <- as.numeric(sub("rs","",rownames(ea.readcount30X@snp.info)))
head(selected.snps.idx)

#estimate genome wide Fst across all the popuatlions
ea.readcount30X.fst<-computeFST(ea.readcount30X)
ea.readcount30X.fst$FST
#genome wide Fst is -0.01960076, so 0

# Block-Jackknife estimation of FST standard-error and confidence intervals:
ea.readcount30X.fst<-computeFST(ea.readcount30X,nsnp.per.bjack.block = 1000, verbose=FALSE)
ea.readcount30X.fst$FST
#same value  -0.01960076

ea.readcount30X.fst$mean.fst #block-jacknife estimate of s.e.
#-0.02009452
ea.readcount30X.fst$se.fst #s.e. of the genome-wide Fst estimate
#2.90547e-05
ea.readcount30X.fst$mean.fst+c(-1.96,1.96)*ea.readcount30X.fst$se.fst
# -0.02015146 -0.02003757


#Computing multi-locus FST to scan the genome over sliding-windows of SNPs

ea.readcount30X.fst<-computeFST(ea.readcount30X,sliding.window.size=50)
ea.readcount30X.fst<-computeFST(ea.readcount30X,sliding.window.size=100)
ea.readcount30X.fst<-computeFST(ea.readcount30X,sliding.window.size=10)
#we have 42K scaffolds so loosing a lot of data with sliding window size of 50-100 chrom?


#6722 chromosomes scanned (with more than 50 SNPs)
#Average (min-max) Window Sizes 0.3 ( 0.1 - 4.9 ) kb


#931 chromosomes scanned (with more than 100 SNPs)
#Average (min-max) Window Sizes 0.8 ( 0.2 - 7.4 ) kb


#31251 chromosomes scanned (with more than 10 SNPs)
#Average (min-max) Window Sizes 0 ( 0 - 6.9 ) kb

#I will just play around with 100 as I dont really understand
ea.readcount30X.fst<-computeFST(ea.readcount30X,sliding.window.size=100)

plot(ea.readcount30X.fst$sliding.windows.fst$CumulatedPosition/1e6,
     ea.readcount30X.fst$sliding.windows.fst$MultiLocusFst,
     xlab="Cumulated Position (in Mb)",ylab="Muli-locus Fst")
     #col=as.numeric(ea.readcount30X.fst$sliding.windows.fst$Chr),pch=16) Doesnt work as we dont have chromossome numbers
abline(h=ea.readcount30X.fst$FST,lty=2)

head(ea.readcount30X.fst$sliding.windows.fst$CumulatedPosition/1e6)
head(ea.readcount30X.fst$sliding.windows.fst$MultiLocusFst)
head(ea.readcount30X.fst$sliding.windows.fst$Chr)

#Manhattan plot of the multi-locus FST computed over sliding-windows of 50 SNPs on the PoolSeq example data. The dashed line indicates the estimated overall genome-wide FST . The 20 simulated
#chromosomes are represented by alternate colors


#pairwise FST

ea.pairwisefst<-compute.pairwiseFST(ea.readcount30X,verbose=FALSE)
#heatmap
#Heatmap representing the pairwise-population FST matrix of the 14 populations of the 30XPool-Seq example data set
heatmap(ea.pairwisefst)
#it moves pops which are more similar to each other 

#Block-Jackknife estimation of FST standard-error and visualisation of confidence intervals
ea.pairwisefst@PairwiseFSTmatrix
plot(ea.pairwisefst)


head(ea.pairwisefst@values)
library(ggplot2)

# add mean to ggplot2 boxplot
ggplot(ds, aes(x = label, y = temperature, fill = label)) +
  geom_boxplot() +
  stat_summary(fun = "mean", geom = "point", shape = 8,
               size = 2, color = "white")

# Transpose the data using the gather function
data_cov_transposed <- gather(data_COV, key = "ID", value = "weight", -CHROM, -POS, -REF, -ALT)


####count_data####

data2 <- data
#count data
data_CNT <- data2[, c("CHROM", "POS", "REF", "ALT")]
# Identify columns with names containing "_CNT"
count_columnsref <- grep(".REF_CNT$", names(data2))
count_columnsalt <- grep(".ALT_CNT$", names(data2))
# Rename the columns to remove "_CNT"
new_colnames_ref <- gsub(".REF_CNT", "", names(data2)[count_columnsref])
new_colnames_alt <- gsub(".ALT_CNT", "", names(data2)[count_columnsalt])
colnames(data)[count_columnsref] <- new_colnames_ref
colnames(data)[count_columnsalt] <- new_colnames_alt

# Add the selected columns to the new dataset
data_cnt_ref <- cbind(data_CNT, data[, count_columnsref, drop = FALSE])
data_cnt_alt <- cbind(data_CNT, data[, count_columnsalt, drop = FALSE])

# Transpose the data using the gather function
data_cnt_ref_transposed <- gather(data_cnt_ref, key = "ID", value = "REF", -CHROM, -POS, -REF, -ALT)
data_cnt_alt_transposed <- gather(data_cnt_alt, key = "ID", value = "ALT", -CHROM, -POS, -REF, -ALT)

#merge files

sample_info <- data.frame(
  CHROM = data_cnt_alt_transposed$CHROM,
  POS = data_cnt_alt_transposed$POS,
  ID = data_cnt_alt_transposed$ID,
  ALT = data_cnt_alt_transposed$ALT,  # Final column from data_cnt_alt_transposed
  REF = data_cnt_ref_transposed$REF,  # Final column from data_cnt_ref_transposed
  COV = data_cov_transposed$weight  # Final column from data_cov_transposed
)

#create time and year factors

# Extract 'year' from the 'ID' column
sample_info$year <- as.numeric(sub("^EA_(\\d+)_.*", "\\1", sample_info$ID))

# Extract 'time' from the 'ID' column
sample_info$time <- sub("^EA_\\d+_(.*)", "\\1", sample_info$ID)

# Print the head of the updated 'sample_info' DataFrame
head(sample_info)


#allele frequency matrix

new_data <- data.frame(
  ID = sample_info$ID,
  ALT_REF = paste(sample_info$ALT, sample_info$REF, sep = ",")
)
 colnames(new_data) <- c("ID","freq")

 #pivot data
#spread_data <- new_data %>%
#  spread(ID, freq) %>%
#  t()

#save(list = ls(all.names = TRUE), file = "files_for_glm.Rdata")

#this is taking too long so whilst script runs I'll try it another way as I also need to estimate number of chromosomes per pooled sample

t12009 <- new_data[new_data$ID == "EA_2009_T1", ]
t22009 <- new_data[new_data$ID == "EA_2009_T2", ]
t32009 <- new_data[new_data$ID == "EA_2009_T3", ]
t42009 <- new_data[new_data$ID == "EA_2009_T4", ]

t12011 <- new_data[new_data$ID == "EA_2011_T1", ]
t22011 <- new_data[new_data$ID == "EA_2011_T2", ]

t12015 <- new_data[new_data$ID == "EA_2015_T1", ]
t22015 <- new_data[new_data$ID == "EA_2015_T2", ]
t32015 <- new_data[new_data$ID == "EA_2015_T3", ]
t42015 <- new_data[new_data$ID == "EA_2015_T4", ]

t12022 <- new_data[new_data$ID == "EA_2022_T1", ]
t22022 <- new_data[new_data$ID == "EA_2022_T2", ]
t32022 <- new_data[new_data$ID == "EA_2022_T3", ]
t42022 <- new_data[new_data$ID == "EA_2022_T4", ]

merged_data <- cbind(t12009$freq, t22009$freq, t32009$freq, t42009$freq,
                    t12011$freq, t22011$freq,
                    t12015$freq, t22015$freq, t32015$freq, t42015$freq,
                    t12022$freq, t22022$freq, t32022$freq, t42022$freq)
merged_data2 <- cbind(data$CHROM, data$POS, t12009$freq, t22009$freq, t32009$freq, t42009$freq,
                    t12011$freq, t22011$freq,
                    t12015$freq, t22015$freq, t32015$freq, t42015$freq,
                    t12022$freq, t22022$freq, t32022$freq, t42022$freq)

# Remove quotes from the matrix
merged_data <- gsub("\"", "", merged_data)
# Set column names
colnames(merged_data) <- c("EA_2009_T1", "EA_2009_T2", "EA_2009_T3", "EA_2009_T4", 
  "EA_2011_T1", "EA_2011_T2", 
  "EA_2015_T1", "EA_2015_T2", "EA_2015_T3", "EA_2015_T4", 
  "EA_2022_T1", "EA_2022_T2", "EA_2022_T3", "EA_2022_T4")

colnames(merged_data2) <- c("chrom", "pos", "EA_2009_T1", "EA_2009_T2", "EA_2009_T3", "EA_2009_T4", 
  "EA_2011_T1", "EA_2011_T2", 
  "EA_2015_T1", "EA_2015_T2", "EA_2015_T3", "EA_2015_T4", 
  "EA_2022_T1", "EA_2022_T2", "EA_2022_T3", "EA_2022_T4")
# Remove the quotes
merged_data <- gsub("\"", "", merged_data)
merged_data <- as.matrix(merged_data)

merged_data2 <- gsub("\"", "", merged_data2)
merged_data2 <- as.matrix(merged_data2)

# Specify the file path where you want to save the matrix
file_path <- "AF_matrix.txt"  # Change this to your desired file path
file_path2 <- "site_info.txt"  # Change this to your desired file path
file_path3 <- "AF_matrix_with_chrom.txt"  # Change this to your desired file path
# Save the matrix to a text file
write.table(merged_data, file = file_path, quote = FALSE, sep = "\t")
write.table(new_data, file = file_path2, quote = FALSE, sep = "\t")
write.table(new_data, file = file_path3, quote = FALSE, sep = "\t")

#count number of chromossomes per sample


a <- sample_info[sample_info$ID == "EA_2009_T1", ]
library(dplyr)
n_distinct(a$CHROM)


#####
#chrom pos data

chrom_poly_position_data1 <- a[,c(1:2)]
file_path4 <- "/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/variants/chrom-poly-position-seasonal-data.txt"  # Change this to your desired file path
# Save the matrix to a text file
write.table(chrom_poly_position_data1, file = file_path4, quote = FALSE, sep = "\t", row.names= FALSE)


###############GLM SCRIPT#######################


# Script to run seasonal analysis on mel seasonal data
#use early late (2011 t2=t4)
#GLM loop#
# Open the file for writing (creates a new file if it doesn't exist)
#fileout <- file("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/GLM/seasonal-glm-output.txt", "w")

# Open the file for writing (creates a new file if it doesn't exist)
fileout <- file("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/GLM/seasonal-glm-output-time-uncorrected-coverage.txt", "w")
fileout2 <- file("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/GLM/seasonal-glm-output-year-uncorrected-coverage.txt", "w")
# Loop through each line in freq_matrix
for (i in 1:nrow(freq_matrix)) {
  # Fit the GLM model for the current line

  
  out <- summary(glm(freq_matrix[i, ] ~ popinfo_glm$time + popinfo_glm$Y, family = binomial, weights = cov_matrix[i, ]))

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


```
