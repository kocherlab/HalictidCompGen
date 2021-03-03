#!/bin/bash
#SBATCH -n 1
#SBATCH --cpus-per-task=8
#SBATCH --time=23:00:00 --qos=1day
#SBATCH --mem=100000
#SBATCH --mail-user=berubin@princeton.edu
#SBATCH --job-name=AAUR_trinotate
#SBATCH --output=/Genomics/kocherlab/berubin/annotation/hic/trinotate/AAUR/AAUR_trinotate-%j.out
#SBATCH --error=/Genomics/kocherlab/berubin/annotation/hic/trinotate/AAUR/AAUR_trinotate-%j.err

module add trinotate

/usr/local/Trinotate/3.0.1/admin/Build_Trinotate_Boilerplate_SQLite_db.pl Trinotate

#makeblastdb -in uniprot_sprot.pep -dbtype prot
#gunzip Pfam-A.hmm.gz
#hmmpress Pfam-A.hmm

blastx -query AAUR.fna -db /Genomics/kocherlab/berubin/annotation/trinotate/uniprot_sprot.pep -num_threads 8 -max_target_seqs 1 -outfmt 6 > blastx.outfmt6

blastp -query AAUR.faa -db /Genomics/kocherlab/berubin/annotation/trinotate/uniprot_sprot.pep -num_threads 8 -max_target_seqs 1 -outfmt 6 > blastp.outfmt6

hmmscan --cpu 8 --domtblout TrinotatePFAM.out /Genomics/kocherlab/berubin/annotation/trinotate/Pfam-A.hmm AAUR.faa > pfam.log

signalp -f short -n signalp.out AAUR.faa

tmhmm --short < AAUR.faa > tmhmm.out

python /Genomics/kocherlab/berubin/annotation/trinotate/make_map.py AAUR.faa AAUR.faa_lens map.txt

Trinotate Trinotate.sqlite init --transdecoder_pep AAUR.faa_lens --transcript_fasta AAUR.fna --gene_trans_map map.txt

Trinotate Trinotate.sqlite LOAD_swissprot_blastp blastp.outfmt6

Trinotate Trinotate.sqlite LOAD_pfam TrinotatePFAM.out

Trinotate Trinotate.sqlite LOAD_tmhmm tmhmm.out

Trinotate Trinotate.sqlite LOAD_signalp signalp.out

Trinotate Trinotate.sqlite report > trinotate_annotation_report.xls
