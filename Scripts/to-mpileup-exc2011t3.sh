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
        -b $WORK/seasonal_adaptation/analysis/variants/bamfiles-exc2011T3.txt  \
        -o all.mpileup-excluding2011T3
