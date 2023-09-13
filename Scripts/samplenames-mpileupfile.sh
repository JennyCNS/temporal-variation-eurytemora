#get sample names in mpileup file

module load samtools/1.10
module load gatk

#sample name in BAM file
gatk GetSampleName -I EA_2007_T1.bam -O EA_2007_T1.txt
 
 #sample name RG bam file
samtools view EA_2007_T1.bam | head
samtools view -H EA_2007_T1.bam | head | grep '@HD'
samtools view -H EA_2007_T1.bam | grep '@RG'
