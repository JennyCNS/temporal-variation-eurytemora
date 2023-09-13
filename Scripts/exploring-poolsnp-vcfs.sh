### investigating output of vcf file trial6

#what I want to do is check how each sample behaved in terms of datamissingness etc.
#trial 5 contains only polymorphic variants whilst trial 6 contains all types of variants - read more on what that means
#aim is to visualise how these samples differ and then plot some stats

cd /gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/PoolSNP/



#check number of CHROM in analyisis

zgrep -v "^##" ../PoolSNP/output_maf0.01mincov50trial5.vcf.gz | awk '{print $1}' | uniq | wc -l
#output
#1178 chrom positiond in 21377 variants called



module load bcftools/1.10.2
conda install -c bioconda tabix


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



#####STATS

#now we have tab index (tbi) for each of the files

gunzip output_maf0.01mincov50trial5.vcf.gz
bgzip output_maf0.01mincov50trial5.vcf.gz
tabix -p vcf output_maf0.01mincov50trial5.vcf.gz

#running vcf tools to check data missingness (or trying to transform data into plink format) I get that we have an issue with ploidy. How to fix this?

#fix polyploidy

bcftools +fixploidy output_maf0.01mincov50trial5.vcf.gz > final_vcfs/output_maf0.01mincov50trial5fixed.vcf.gz

#explore data statistics with plink

#bash script 
conda install -c bioconda plink

# Output directory where processed files will be stored
output=/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/PoolSNP/stats

# Loop through all files in the current directory
for f in *; do
  if [ -f "$f" ]; then
    echo "Processing file: $f"

    # Run the plink command on the current file
    plink --vcf "$f" --allow-extra-chr --double-id --recode --make-bed --out "$output/${f%.*}"
  fi
done

plink --file final-maf0.01-mincov50-mincount10-allsites0-allsamples_b.vcf --allow-extra-chr --missing --hardy --freq --out /output/final-allsamples
plink --file final-maf0.01-mincov50-mincount10-allsites0-excluding20072011.vcf --allow-extra-chr --missing --hardy --freq --out /output/10mincount_filtered
plink --file final-maf0.01-mincov50-mincount5-allsites0-excluding20072011.vcf --allow-extra-chr --missing --hardy --freq --out /output/5mincount_filtered


module load vcftools/0.1.14
vcftools --vcf final-maf0.01-mincov50-mincount10-allsites0-allsamples_b.vcf.gz --missing-indv
 
#plot outputs in R
#R scripts are saved in computer


#check snps in the output files

module load bcftools/1.10.2
bcftools view -v snps input.vcf.gz -o output_snps.vcf.gzls