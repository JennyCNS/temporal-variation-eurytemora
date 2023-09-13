#!/bin/bash

module load bcftools/1.10.2

# Path to the pooled VCF file
pooled_vcf_file="/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/grenedalf-freq/finalfile.vcf.gz"
output_file="maf_output.txt"


# Use bcftools to preprocess the VCF file and calculate MAF
bcftools view "$pooled_vcf_file" | \
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%AD]\n' | \
awk -F'\t' -v OFS='\t' '
    function get_alt_count(ad) {
        split(ad, depths, ",");
        return depths[2];
    }
    function get_total_alleles(ad) {
        split(ad, depths, ",");
        return depths[1] + depths[2];
    }
    {
        alt_count = 0;
        total_alleles = 0;
        for (i = 6; i <= NF; i++) {
            alt_count += get_alt_count($i);
            total_alleles += get_total_alleles($i);
        }
        maf = (total_alleles > 0) ? alt_count / total_alleles : 0;
        print $1, $2, $3, $4, maf;
    }' \
> "$output_file"

echo "MAF values have been calculated and saved to $output_file"



##do it manually to check values


bcftools view "$pooled_vcf_file" | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%AD]\n' > maf_output.txt  \


#!/bin/bash

# Input VCF file
input_vcf="maf_output.txt "

# Output file
output_file="trial.txt"

# Use awk to modify the VCF file
awk -F'\t' -v OFS='\t' '
{
    # Initialize the total and sum variables
    total = 0;
    sum = 0;

    # Loop through columns 5 to 18
    for (i = 5; i <= 21; i++) {
        # Replace "." with 0 in each column
        gsub(/\./, "0", $i);

        # Split the column using comma as delimiter
        n = split($i, arr, ",");

        # Loop through the elements after the comma and sum them up
        for (j = 2; j <= n; j++) {
            sum += arr[j];
        }

        # Add the total of each column to the overall total
        total += $i;
    }

    # Append the total and sum values as new columns
    print $0, total, sum;
}
' "$input_vcf" > "$output_file"

echo "Modified VCF file saved to $output_file"


###############################

####calculate MAF by hand


bcftools view "$pooled_vcf_file" | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%AD]\n' > noncalculatedfreqs.txt

awk -v OFS='\t' '{print $1, $2, $3, $4, $19, $20}' noncalculatedfreqs.txt > maf.txt
awk -F'\t' -v OFS='\t' '{ sum = $NF + $(NF-1); print $0, sum }' maf.txt > sum_output.txt
awk -v OFS='\t' '{print $1, $2, $3, $4, $19, $20, $6 / $7}' sum_output.txt > final_maf.txt


### pdocuing allele frequency table for pca
awk 'BEGIN{OFS=FS="\t"} {
    printf "TRINITY_%s\t%s", $1, $2;
    for(i=6; i<=18; i++) {
        split($i, counts, ",");
        sum = counts[1] + counts[2];
        result = (sum == 0) ? 0 : counts[2] / sum;
        printf "\t%.8f", result;
    }
    printf "\n";
}' noncalculatedfreqs.txt > output_file.txt


zgrep -v "^#" finalfile.vcf.gz | head -n 1
zgrep -v "^##" finalfile.vcf.gz | head -n 1 > header.txt
cat header.txt output_file.txt > frequencies_maf-notrounded.txt


#MAF fits with grenedalf estimations. 