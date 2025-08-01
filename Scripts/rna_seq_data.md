#RNA seq data
#copepods
#seasonal adaptation

#stats raw and trimmed files with seqkit
```bash
seqkit stats *.gz
```
#aligning raw files

#create file with samples for quality check
``` bash
awk 'BEGIN {print "sample,fastq_1,,fastq_2,strandedness"} 
NR%2==1 {
    f1=$0;
    getline;
    f2=$0;
    match(f1, /^[^_]+_[^_]+_[^_]+/);
    sample=substr(f1, RSTART, RLENGTH);
    print sample "," f1 ",," f2 ",auto";
}' samplesheet.csv > new_samplesheet.csv
```

# doesnt work
gff transform into gtf
cp Purple_Eaffinis_maker_annotation.gtf Purple_Eaffinis_maker_annotation.gtf.bak

#somethings wrong with the annotation

sed -i -r 's/(transcript_id "([A-Za-z0-9\-]+)";)/\1 gene_id "\2";/g' Purple_Eaffinis_maker_annotation.gtf


run nextflow to check the quality of RNA-samples
```bash
nextflow run nf-core/rnaseq \
-c /gxfs_home/geomar/smomw573/work/seasonal_adaptation/transcriptome/variants/nf_core_NEC.config \
-profile singularity \
--gtf /gxfs_home/geomar/smomw573/work/seasonal_adaptation/new_genome/Purple_Eaffinis_maker_annotation.gtf \
--igenomes_ignore true \
--aligner star_salmon \
--input samplesheet.csv \
--fasta /gxfs_home/geomar/smomw573/work/seasonal_adaptation/new_genome/Purple_Eaffinis_genome_assembly.fa \
--outdir /gxfs_home/geomar/smomw573/work/seasonal_adaptation/transcriptome/analysis/nextflow-rna-output \
--trimmer fastp \
-resume

#run 2

nextflow run nf-core/rnaseq \
-c /gxfs_home/geomar/smomw573/work/seasonal_adaptation/transcriptome/variants/nf_core_NEC.config \
-profile singularity \
--gtf /gxfs_home/geomar/smomw573/work/seasonal_adaptation/new_genome/Purple_Eaffinis_maker_annotation.gtf \
--igenomes_ignore true \
--aligner star_salmon \
--input samplesheet.csv \
--fasta /gxfs_home/geomar/smomw573/work/seasonal_adaptation/new_genome/Purple_Eaffinis_genome_assembly.fa \
--outdir /gxfs_home/geomar/smomw573/work/seasonal_adaptation/transcriptome/analysis/nextflow-rna-output.2 \
--trimmer fastp \
--forward_strand unstranded \
--save_trimmed \
-resume

#run 3
#with the annotation file till gave me

nextflow run nf-core/rnaseq \
-c /gxfs_home/geomar/smomw573/work/seasonal_adaptation/transcriptome/variants/nf_core_NEC.config \
-profile singularity \
--gtf /gxfs_home/geomar/smomw573/work/seasonal_adaptation/new_genome/Purple_Eaffinis_annotation.gtf \
--igenomes_ignore true \
--aligner star_salmon \
--input samplesheet.csv \
--fasta /gxfs_home/geomar/smomw573/work/seasonal_adaptation/new_genome/Purple_Eaffinis_genome_assembly.fa \
--outdir /gxfs_home/geomar/smomw573/work/seasonal_adaptation/transcriptome/analysis/nextflow-rna-output.2 \
--trimmer fastp \
--forward_strand unstranded \
--save_trimmed \
-resume

```

something isnt working so I will try out manually
```bash
gffread /gxfs_home/geomar/smomw573/work/seasonal_adaptation/new_genome/Purple_Eaffinis_maker_annotation.gff -g /gxfs_home/geomar/smomw573/work/seasonal_adaptation/new_genome/Purple_Eaffinis_genome_assembly.fa -w transcripts.fa
```
didnt work
till converted with agat
agat_convert_sp_gff2gtf.pl -gff Purple_Eaffinis_maker_annotation.gff -o Purple_Eaffinis_maker_annotation.gtf

# ok the pipeline worked!

# now I will start with differential gene expression analysis

# I will use deseq2

# create conda env rna-seq
install subread
conda install bioconda::subread

run slurm job

```
featureCounts \
-t exon  \
-g gene_id  \
-a /gxfs_home/geomar/smomw573/work/seasonal_adaptation/new_genome/Purple_Eaffinis_maker_annotation.gtf \
-T 16 -p -B -C \
-o /gxfs_home/geomar/smomw573/work/seasonal_adaptation/transcriptome/variants/raw_gene_counts.txt \
/gxfs_home/geomar/smomw573/work/seasonal_adaptation/transcriptome/analysis/nextflow-rna-output.2/star_salmon/*.bam
```

#new trial with salmon
#first produce the transcript fasta file

agat_convert_sp_gff2gtf.pl -gff Purple_Eaffinis_maker_annotation.gff -o Purple_Eaffinis_maker_annotation.gtf

agat_sp_extract_sequences.pl \
  --gtf annotation.gtf \
  -g Purple_Eaffinis_genome_assembly.fa \
  -o Purple_Eaffinis_transcripts.fa



# create transcript.fa with ggfread

gffread -w transcripts.fa -g Purple_Eaffinis_genome_assembly.fa Purple_Eaffinis_maker_annotation.gff

# index with salmon

salmon index -t transcripts.fa -i salmon_index
Elapsed time: 0.442899s

[2025-06-23 11:06:08.106] [jointLog] [warning] Removed 131 transcripts that were sequence duplicates of indexed transcripts.
[2025-06-23 11:06:08.106] [jointLog] [warning] If you wish to retain duplicate transcripts, please use the `--keepDuplicates` flag
[2025-06-23 11:06:08.106] [jointLog] [info] Replaced 4 non-ATCG nucleotides
[2025-06-23 11:06:08.106] [jointLog] [info] Clipped poly-A tails from 0 transcripts
[2025-06-23 11:06:08.107] [jointLog] [info] Building rank-select dictionary and saving to disk
[2025-06-23 11:06:08.108] [jointLog] [info] done

```bash
#!/bin/bash

input_dir="/gxfs_home/geomar/smomw573/work/seasonal_adaptation/transcriptome/analysis/nextflow-rna-output.2/fastp"
output_dir="/gxfs_home/geomar/smomw573/work/seasonal_adaptation/transcriptome/analysis/salmon/quants"

mkdir -p "${output_dir}"

for fq1 in "${input_dir}"/*_1.fastp.fastq.gz; do
    # Extract just the filename from the full path
    filename=$(basename "${fq1}")
    # Remove _1.fastp.fastq.gz to get sample prefix
    samp=${filename%_1.fastp.fastq.gz}

    echo "Processing sample ${samp}"

    salmon quant -i /gxfs_home/geomar/smomw573/work/seasonal_adaptation/new_genome/salmon_index -l A \
        -1 "${input_dir}/${samp}_1.fastp.fastq.gz" \
        -2 "${input_dir}/${samp}_2.fastp.fastq.gz" \
        -p 8 --validateMappings --gcBias --numGibbsSamples 20 -o "${output_dir}/${samp}_quant"
done
```

#got this warning for all samples
[2025-06-23 11:54:43.874] [jointLog] [warning] NOTE: Read Lib [[ /gxfs_home/geomar/smomw573/work/seasonal_adaptation/transcriptome/analysis/nextflow-rna-output.2/fastp/T2F3_cold_5_1.fastp.fastq.gz, /gxfs_home/geomar/smomw573/work/seasonal_adaptation/transcriptome/analysis/nextflow-rna-output.2/fastp/T2F3_cold_5_2.fastp.fastq.gz]] :

Detected a *potential* strand bias > 1% in an unstranded protocol check the file: /gxfs_home/geomar/smomw573/work/seasonal_adaptation/transcriptome/analysis/salmon/quants/T2F3_cold_5_quant/lib_format_counts.json for details


# producing metadata file
```bash
#!/bin/bash

# Set paths
input_dir="/gxfs_home/geomar/smomw573/work/seasonal_adaptation/transcriptome/analysis/nextflow-rna-output.2/fastp"
salmon_dir="/gxfs_home/geomar/smomw573/work/seasonal_adaptation/transcriptome/analysis/salmon/quants"

# Output CSV
output_csv="metadata_final.csv"

# Write header
echo "Sample,Forward,Reverse,Condition,Quants" > "$output_csv"

# Loop through all _1.fastp.fastq.gz files
for fq1 in "$input_dir"/*_1.fastp.fastq.gz; do
    # Extract sample name
    filename=$(basename "$fq1")
    sample=${filename%_1.fastp.fastq.gz}

    # Build paths
    fq2="${input_dir}/${sample}_2.fastp.fastq.gz"
    forward="$fq1"
    reverse="$fq2"
    condition=$(echo "$sample" | grep -oE 'cold|warm')
    quant="${salmon_dir}/${sample}_quant"

    # Write line
    echo "${sample},${forward},${reverse},${condition},${quant}" >> "$output_csv"
done
```


#something is wrong and I will try to use the old transcriptome (from 2014)

# create transcript.fa with ggfread

salmon index -t GBGO01.1.fsa_nt -i eaff_transcriptome_index

# index with salmon

```bash
salmon index -t transcripts_old_genome.fa -i salmon_index

input_dir="/gxfs_home/geomar/smomw573/work/seasonal_adaptation/transcriptome/analysis/nextflow-rna-output.2/fastp"
output_dir="/gxfs_home/geomar/smomw573/work/seasonal_adaptation/transcriptome/analysis/salmon/quants_2"

mkdir -p "${output_dir}"

for fq1 in "${input_dir}"/*_1.fastp.fastq.gz; do
    # Extract just the filename from the full path
    filename=$(basename "${fq1}")
    # Remove _1.fastp.fastq.gz to get sample prefix
    samp=${filename%_1.fastp.fastq.gz}

    echo "Processing sample ${samp}"

    salmon quant -i /gxfs_home/geomar/smomw573/work/seasonal_adaptation/new_genome/eaff_transcriptome_index -l A \
        -1 "${input_dir}/${samp}_1.fastp.fastq.gz" \
        -2 "${input_dir}/${samp}_2.fastp.fastq.gz" \
        -p 8 --validateMappings --gcBias --numGibbsSamples 20 -o "${output_dir}/${samp}_quant"
done
```

# produce tx2gene.csv file
grep "^>" GBGO01.1.fsa_nt | cut -d " " -f1 | sed 's/>//' | awk '{print $1","$1}' > tx2gene.csv

