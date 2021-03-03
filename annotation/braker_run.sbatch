#!/bin/bash
#SBATCH -n 1
#SBATCH --cpus-per-task=10
#SBATCH --time=6-23:00 --qos=1wk
#SBATCH --mem=100000
#SBATCH --mail-user=berubin@princeton.edu
#SBATCH --job-name=brakLZEP
#SBATCH --output=/Genomics/kocherlab/berubin/annotation/hic/LZEP/braker/outputs/brak_LZEP-%j.out
#SBATCH --error=/Genomics/kocherlab/berubin/annotation/hic/LZEP/braker/errors/brak_LZEP-%j.err

module add BRAKER1

module add exonerate

module add ncbi-blast

/Genomics/kocherlab/berubin/local/src/hisat2-2.0.5/hisat2-build /Genomics/kocherlab/berubin/annotation/hic/LZEP/repeats/LZEP_rm_masked/LZEP_genome_v2.0.fasta.masked /scratch/tmp/berubin/annotation/braker/LZEP/LZEP_build

cd /scratch/tmp/berubin/annotation/braker/LZEP

/Genomics/kocherlab/berubin/local/src/hisat2-2.0.5/hisat2 --max-intronlen 100000 -p 10 -x /scratch/tmp/berubin/annotation/braker/LZEP/LZEP_build --phred33 -1 /Genomics/kocherlab/berubin/data/rnaseq/all_demult_rnaseq/BER1-E2_head_R1.fastq.gz,/Genomics/kocherlab/berubin/data/rnaseq/all_demult_rnaseq/BER1-E2_abd_R1.fastq.gz,/Genomics/kocherlab/berubin/data/rnaseq/all_demult_rnaseq/BER2-B10_ant_R1.fastq.gz -2 /Genomics/kocherlab/berubin/data/rnaseq/all_demult_rnaseq/BER1-E2_head_R2.fastq.gz,/Genomics/kocherlab/berubin/data/rnaseq/all_demult_rnaseq/BER1-E2_abd_R2.fastq.gz,/Genomics/kocherlab/berubin/data/rnaseq/all_demult_rnaseq/BER2-B10_ant_R2.fastq.gz | samtools view -bS - | samtools sort -o /scratch/tmp/berubin/annotation/braker/LZEP/LZEP_sorted.bam

/Genomics/kocherlab/berubin/local/src/BRAKER_v2.1.0/braker.pl --genome=/Genomics/kocherlab/berubin/annotation/hic/LZEP/repeats/LZEP_rm_masked/LZEP_genome_v2.0.fasta.masked --bam=/scratch/tmp/berubin/annotation/braker/LZEP/LZEP_sorted.bam --species=LZEP_brakerold --gff3 --cores=10 --AUGUSTUS_CONFIG_PATH=/Genomics/kocherlab/berubin/local/src/augustus/config/ --GENEMARK_PATH=/Genomics/kocherlab/berubin/local/src/gm_et_linux_64/gmes_petap/ --softmasking --overwrite --AUGUSTUS_ab_initio

cp /scratch/tmp/berubin/annotation/braker/LZEP/braker/LZEP_braker/augustus.hints.gff3 /Genomics/kocherlab/berubin/annotation/hic/LZEP/braker/LZEP_braker_augustus.gff3

cp /scratch/tmp/berubin/annotation/braker/LZEP/braker/LZEP_braker/augustus.ab_initio.gff3 /Genomics/kocherlab/berubin/annotation/hic/LZEP/braker/LZEP_braker_augustus_abinitio.gff3

