#plink PCA

module load plink/1.07
input=/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/grenedalf-freq/finalfile.vcf.gz
output=/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/grenedalf-freq/plink-files/
prefix=finalsnps

#prepare data for PCA
plink --vcf $input --double-id --make-bed --allow-extra-chr --out $output/finalsnps
plink --bfile $prefix --recode --allow-extra-chr --out $prefix


plink --file $prefix --double-id --allow-extra-chr --set-missing-var-ids @:# --indep-pairwise 50 10 0.1 --out $prefix
plink --file $prefix --double-id --allow-extra-chr --set-missing-var-ids @:# --extract $prefix.prune.in --make-bed --pca --out $prefix

module load R/4.0.2  
