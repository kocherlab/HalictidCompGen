#!/bin/bash
#SBATCH -n 1
#SBATCH --cpus-per-task=10
#SBATCH --time=6-23:00 --qos=1wk
#SBATCH --mem=20000
#SBATCH --mail-user=berubin@princeton.edu
#SBATCH --job-name=HLIG_repeats
#SBATCH --output=/Genomics/kocherlab/berubin/annotation/hic/HLIG/repeats/repeats-%j.out
#SBATCH --error=/Genomics/kocherlab/berubin/annotation/hic/HLIG/repeats/repeats-%j.err

module add RepeatModeler
export PATH=$PATH:/Genomics/kocherlab/berubin/local/src/ncbi-blast-2.5.0+/bin/

BuildDatabase -name HLIG_db -engine ncbi /Genomics/kocherlab/berubin/assembly/hic/HLIG/HLIG_genome_v2.0.fasta

RepeatModeler -database HLIG_db -pa 3 -engine ncbi >& HLIG.out

/Genomics/kocherlab/berubin/local/src/cdhit/cd-hit-est -i ./RM_*/consensi.fa -o rm_noredun.fa -c 0.8 -M 20000 -n 5 -aS 0.8 -r 1 -T 10

blastx -query rm_noredun.fa -db /Genomics/kocherlab/berubin/annotation/reference_genomes/uniprot/uniprot_sprot_blastdb -outfmt '6 qseqid sseqid pident length evalue bitscore qlen slen' -out rm_noredun_against_uniprot.txt -num_threads 10 

blastx -query rm_noredun.fa -db /Genomics/kocherlab/berubin/annotation/reference_genomes/drosophila/dmel_blastdb -outfmt '6 qseqid sseqid pident length evalue bitscore qlen slen' -out rm_noredun_against_dmel.txt -num_threads 10

python /Genomics/kocherlab/berubin/annotation/filter_repeats.py rm_noredun.fa rm_noredun_noprot.fa rm_noredun_against_uniprot.txt rm_noredun_against_dmel.txt
