#!/bin/bash
#SBATCH -n 1
#SBATCH --cpus-per-task=20
#SBATCH --time=6-00:00 --qos=1wk
#SBATCH --mem=100000
#SBATCH --mail-user=berubin@princeton.edu
#SBATCH --job-name=aur_longranger
#SBATCH --output=/Genomics/kocherlab/berubin/assembly/aurata/fragscaff/longranger-%j.out
#SBATCH --error=/Genomics/kocherlab/berubin/assembly/aurata/fragscaff/longranger-%j.err


python /Genomics/kocherlab/berubin/assembly/rename_scafs.py /Genomics/kocherlab/berubin/assembly/aurata/supernova_both_libs/aurata_both_libs_gc_filter_v0.fasta aurata_both_libs_gc_filter_v0_sorted.fasta AAUR

/Genomics/kocherlab/berubin/local/src/longranger-2.1.2/longranger mkref aurata_both_libs_gc_filter_v0_sorted.fasta

/Genomics/kocherlab/berubin/local/src/longranger-2.1.2/longranger align --reference=refdata-aurata_both_libs_gc_filter_v0_sorted --fastqs=/Genomics/kocherlab/berubin/data/10x/february/aurata --id=aurata_both_libs_gc_filter_v0_align --localcores=20

python /Genomics/kocherlab/berubin/assembly/add_read_group.py aurata_both_libs_gc_filter_v0_align/outs/possorted_bam.bam possorted_rg.bam aurata_both_libs_gc_filter_v0_sorted.fasta

perl /Genomics/kocherlab/berubin/local/src/fragscaff/fragScaff.pl -B possorted_rg.bam -b 1 -E 3000 -G R -N n_coord.bed

perl /Genomics/kocherlab/berubin/local/src/fragscaff/fragScaff.pl -B possorted_rg.bam.E3000.o10000.N.bamParse -A -O E3k_o10k_out -t 1

perl /Genomics/kocherlab/berubin/local/src/fragscaff/fragScaff.pl -B possorted_rg.bam.E3000.o10000.N.bamParse -K possorted_rg.bam.E3000.o10000.N.r1.links.txt -F aurata_both_libs_gc_filter_v0_sorted.fasta -O aurata_E3k_o10k


