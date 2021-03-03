#!/bin/bash
#SBATCH -n 1
#SBATCH --cpus-per-task=10
#SBATCH --time=6-23:00 --qos=1wk
#SBATCH --mem=150000
#SBATCH --mail-user=berubin@princeton.edu
#SBATCH --job-name=tri_LZEP
#SBATCH --output=/Genomics/kocherlab/berubin/transcriptomes/outputs/denovo_LZEP-%j.out
#SBATCH --error=/Genomics/kocherlab/berubin/transcriptomes/errors/denovo_LZEP-%j.err

/Genomics/kocherlab/berubin/local/src/Trinity/Trinity --seqType fq --left /Genomics/kocherlab/berubin/data/rnaseq/all_demult_rnaseq/S77A_R1.fastq.gz,/Genomics/kocherlab/berubin/data/rnaseq/all_demult_rnaseq/S85A_R1.fastq.gz --right /Genomics/kocherlab/berubin/data/rnaseq/all_demult_rnaseq/S77A_R2.fastq.gz,/Genomics/kocherlab/berubin/data/rnaseq/all_demult_rnaseq/S85A_R2.fastq.gz --max_memory 100G --CPU 10 --trimmomatic --normalize_reads --jaccard_clip --output /scratch/tmp/berubin/transcriptomes/trinity_denovo/LZEP/LZEP_trinity_denovo_callum_jaccard
rm -r /scratch/tmp/berubin/transcriptomes/trinity_denovo/LZEP/LZEP_trinity_denovo_callum_jaccard/chrysalis/
rm -r /scratch/tmp/berubin/transcriptomes/trinity_denovo/LZEP/LZEP_trinity_denovo_callum_jaccard/insilico_read_normalization/
rm -r /scratch/tmp/berubin/transcriptomes/trinity_denovo/LZEP/LZEP_trinity_denovo_callum_jaccard/read_partitions/
rm /scratch/tmp/berubin/transcriptomes/trinity_denovo/LZEP/LZEP_trinity_denovo_callum_jaccard/LZEP_R1.fastq.gz.PwU.qtrim.fq
rm /scratch/tmp/berubin/transcriptomes/trinity_denovo/LZEP/LZEP_trinity_denovo_callum_jaccard/LZEP_R2.fastq.gz.PwU.qtrim.fq
rm /scratch/tmp/berubin/transcriptomes/trinity_denovo/LZEP/LZEP_trinity_denovo_callum_jaccard/jellyfish.kmers.fa
rm /scratch/tmp/berubin/transcriptomes/trinity_denovo/LZEP/LZEP_trinity_denovo_callum_jaccard/both.fa
rm /scratch/tmp/berubin/transcriptomes/trinity_denovo/LZEP/LZEP_trinity_denovo_callum_jaccard/*.qtrim.gz
cp /scratch/tmp/berubin/transcriptomes/trinity_denovo/LZEP/LZEP_trinity_denovo_callum_jaccard/Trinity.fasta /Genomics/kocherlab/berubin/transcriptomes/trinity_denovo/LZEP/LZEP_Trinity_callum_jaccard.fasta
