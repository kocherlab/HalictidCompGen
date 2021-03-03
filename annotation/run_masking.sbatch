#!/bin/bash
#SBATCH -n 1
#SBATCH --cpus-per-task=10
#SBATCH --time=23:59:59 --qos=1day
#SBATCH --mem=20000
#SBATCH --mail-user=berubin@princeton.edu
#SBATCH --job-name=HLIG_masking
#SBATCH --output=/Genomics/kocherlab/berubin/annotation/hic/HLIG/repeats/masking-%j.out
#SBATCH --error=/Genomics/kocherlab/berubin/annotation/hic/HLIG/repeats/masking-%j.err

module add RepeatMasker
export PATH=$PATH:/Genomics/kocherlab/berubin/local/src/ncbi-blast-2.5.0+/bin/

cat /Genomics/kocherlab/berubin/annotation/hic/HLIG/repeats/rm_noredun_noprot.fa /Genomics/kocherlab/berubin/annotation/arthropoda_repeatmasker_lib.fa > HLIG_arthropoda_consensi.fa

mkdir HLIG_rm_masked

cp /Genomics/kocherlab/berubin/assembly/hic/HLIG/HLIG_genome_v2.0.fasta HLIG_rm_masked/HLIG_genome_v2.0.fasta

chmod u+w HLIG_rm_masked/HLIG_genome_v2.0.fasta

RepeatMasker -e ncbi -pa 10 -lib HLIG_arthropoda_consensi.fa -gff -xsmall -source -dir HLIG_rm_masked HLIG_rm_masked/HLIG_genome_v2.0.fasta

