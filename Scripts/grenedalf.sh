#grenedalf frequeny estimation

module load htslib/1.10.2  cmake/3.18.4  bzip2/1.0.8   autoconf/2.69   automake/1.16.3

GRENEDALF=/gxfs_home/geomar/smomw573/software/grenedalf/bin/grenedalf
VCF=/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/PoolSNP/finalfile.vcf.gz
GENOME=/gxfs_work1/geomar/smomw573/seasonal_adaptation/genome/Eaffinis.Baltic.PseudoRef.Mar22.fasta \
OUTPUT=/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/grenedalf-freq


$GRENEDALF frequency --vcf-path $VCF --allow-file-overwriting --reference-genome-fasta-file $GENOME --write-total-frequency> $OUTPUT/grenedalftrial.log
























#solution
#use mpileup
#that was not the solution


#trial with mpileup
#redo the mpileup file keeping only the snps that were kept after running poolsnp with the specified parameters


#fix mpileup to only keep the variants I want
#this will be done with python 


def extract_variants_from_mpileup(mpileup_file, variants_list_file, output_file):
    variants_to_extract = set()

    # Read the variants list file and store the CHROM and POS positions in a set
    with open(variants_list_file, 'r') as variants_file:
        for line in variants_file:
            chrom, pos = line.strip().split()
            variants_to_extract.add((chrom, pos))

    # Open the output file to write the extracted variants
    with open(output_file, 'w') as output:
        # Read the mpileup file and extract variants based on positions
        with open(mpileup_file, 'r') as mpileup:
            for line in mpileup:
                fields = line.strip().split('\t')
                chrom, pos = fields[0], fields[1]

                if (chrom, pos) in variants_to_extract:
                    output.write(line)

# Example usage:
mpileup_file = 'all.mpileup-excluding2011T3'
variants_list_file = 'final.bed'
output_file = 'extracted_variants.mpileup'
extract_variants_from_mpileup(mpileup_file, variants_list_file, output_file)




