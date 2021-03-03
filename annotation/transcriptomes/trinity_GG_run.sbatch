#!/bin/bash
#SBATCH -n 1
#SBATCH --cpus-per-task=10
#SBATCH --time=6-23:00 --qos=1wk
#SBATCH --mem=100000
#SBATCH --mail-user=berubin@princeton.edu
#SBATCH --job-name=triggLLEU
#SBATCH --output=/Genomics/kocherlab/berubin/transcriptomes/outputs/GG_LLEU-%j.out
#SBATCH --error=/Genomics/kocherlab/berubin/transcriptomes/errors/GG_LLEU-%j.err

/Genomics/kocherlab/berubin/local/src/hisat2-2.0.5/hisat2-build /Genomics/kocherlab/berubin/assembly/hic/LLEU/LLEU_genome_v2.0.fasta /scratch/tmp/berubin/transcriptomes/trinity_GG/LLEU/LLEU_build

cd /scratch/tmp/berubin/transcriptomes/trinity_GG/LLEU

/Genomics/kocherlab/berubin/local/src/hisat2-2.0.5/hisat2 --max-intronlen 100000 -p 10 -x /scratch/tmp/berubin/transcriptomes/trinity_GG/LLEU/LLEU_build --phred33 -1 /Genomics/kocherlab/berubin/data/rnaseq/all_demult_rnaseq/BER2-B8_head_R1.fastq.gz,/Genomics/kocherlab/berubin/data/rnaseq/all_demult_rnaseq/BER2-B8_abd_R1.fastq.gz,/Genomics/kocherlab/berubin/data/rnaseq/all_demult_rnaseq/BER2-B8_ant_R1.fastq.gz -2 /Genomics/kocherlab/berubin/data/rnaseq/all_demult_rnaseq/BER2-B8_head_R2.fastq.gz,/Genomics/kocherlab/berubin/data/rnaseq/all_demult_rnaseq/BER2-B8_abd_R2.fastq.gz,/Genomics/kocherlab/berubin/data/rnaseq/all_demult_rnaseq/BER2-B8_ant_R2.fastq.gz | samtools view -bS - | samtools sort -o /scratch/tmp/berubin/transcriptomes/trinity_GG/LLEU/LLEU_sorted.bam

/Genomics/kocherlab/berubin/local/src/Trinity/Trinity --genome_guided_bam /scratch/tmp/berubin/transcriptomes/trinity_GG/LLEU/LLEU_sorted.bam --max_memory 100G --CPU 10 --SS_lib_type RF --genome_guided_max_intron 100000 --output /scratch/tmp/berubin/transcriptomes/trinity_GG/LLEU/LLEU_trinity_GG

cp /scratch/tmp/berubin/transcriptomes/trinity_GG/LLEU/LLEU_trinity_GG_jaccard/Trinity-GG.fasta /Genomics/kocherlab/berubin/transcriptomes/trinity_GG/LLEU/LLEU_Trinity-GG.fasta
cp /scratch/tmp/berubin/transcriptomes/trinity_GG/LLEU/LLEU_sorted.bam /Genomics/kocherlab/berubin/transcriptomes/trinity_GG/LLEU/LLEU_sorted.bam
rm -r /scratch/tmp/berubin/transcriptomes/trinity_GG/LLEU/LLEU_trinity_GG/

/Genomics/kocherlab/berubin/local/src/hisat2-2.0.5/hisat2 --max-intronlen 100000 -p 10 -x /scratch/tmp/berubin/transcriptomes/trinity_GG/LLEU/LLEU_build --phred33 -1 /Genomics/kocherlab/berubin/data/rnaseq/all_demult_rnaseq/S81A_R1.fastq.gz -2 /Genomics/kocherlab/berubin/data/rnaseq/all_demult_rnaseq/S81A_R2.fastq.gz | samtools view -bS - | samtools sort -o /scratch/tmp/berubin/transcriptomes/trinity_GG/LLEU/LLEU_sorted_callum.bam

/Genomics/kocherlab/berubin/local/src/Trinity/Trinity --genome_guided_bam /scratch/tmp/berubin/transcriptomes/trinity_GG/LLEU/LLEU_sorted_callum.bam --max_memory 100G --CPU 10 --genome_guided_max_intron 100000 --output /scratch/tmp/berubin/transcriptomes/trinity_GG/LLEU/LLEU_trinity_GG_callum

cp /scratch/tmp/berubin/transcriptomes/trinity_GG/LLEU/LLEU_trinity_GG_callum/Trinity-GG.fasta /Genomics/kocherlab/berubin/transcriptomes/trinity_GG/LLEU/LLEU_Trinity-GG_callum.fasta
cp /scratch/tmp/berubin/transcriptomes/trinity_GG/LLEU/LLEU_sorted_callum.bam /Genomics/kocherlab/berubin/transcriptomes/trinity_GG/LLEU/LLEU_sorted_callum.bam
rm -r /scratch/tmp/berubin/transcriptomes/trinity_GG/LLEU/LLEU_trinity_GG_callum/

rm -r /scratch/tmp/berubin/transcriptomes/trinity_GG/LLEU/

