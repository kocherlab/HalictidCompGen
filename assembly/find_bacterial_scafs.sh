#!/bin/bash
#SBATCH -n 1
#SBATCH --cpus-per-task=12
#SBATCH --time=6-00:00 --qos=1wk
#SBATCH --mail-user=berubin@princeton.edu
#SBATCH --job-name=intros
#SBATCH --output=/Genomics/kocherlab/berubin/sodalis/introgressions/outputs/intro_blast-%j.out
#SBATCH --err=/Genomics/kocherlab/berubin/sodalis/introgressions/errors/intro_blast-%j.err


export PATH=$PATH:/Genomics/kocherlab/berubin/local/src/ncbi-blast-2.5.0+/bin/

makeblastdb -in animal_genomes.fa -out animal_genomes_db -dbtype nucl

makeblastdb -in bacterial_genomes.fa -out bacterial_genomes_db -dbtype nucl

species="pura"
blastn -query /Genomics/kocherlab/berubin/assembly/pura/Augochlora_pura_v1.0.fa -db bacterial_genomes_db -outfmt 6 -num_threads 3 -out pura_against_bacteria.txt 

python parse_blast.py /Genomics/kocherlab/berubin/assembly/${species}/Augochlora_${species}_v1.0.fa ${species}_against_bacteria.txt ${species}_against_bacteria.fa 

blastn -query ${species}_against_bacteria.fa -db bacterial_genomes_db -outfmt 6 -num_threads 2 -out ${species}_hits_against_bacteria.txt

blastn -query ${species}_against_bacteria.fa -db animal_genomes_db -outfmt 6 -num_threads 2 -out ${species}_hits_against_animals.txt 

python parse_blast.py /Genomics/kocherlab/berubin/assembly/${species}/Augochlora_${species}_v1.0.fa ${species}_against_bacteria.txt ${species}_against_bacteria_noanimals.fa ${species}_hits_against_animals.txt 

python parse_blast_id.py ${species}_hits_against_bacteria.txt ${species}_taxonomy.txt



