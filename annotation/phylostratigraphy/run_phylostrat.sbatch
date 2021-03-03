#!/bin/bash
#SBATCH -n 1
#SBATCH --cpus-per-task=18
#SBATCH --time=6-23:59 --qos=1wk
#SBATCH --mem=100000
#SBATCH --mail-user=berubin@princeton.edu
#SBATCH --job-name=MGEN_animals
#SBATCH --output=/Genomics/kocherlab/berubin/stratigraphy/outputs/MGEN_animals-%j.out
#SBATCH --error=/Genomics/kocherlab/berubin/stratigraphy/errors/MGEN_animals-%j.err

module add ncbi-blast

makeblastdb -dbtype prot -in protein_db.faa

/Genomics/kocherlab/berubin/local/src/Phylostratigraphy/createPSmap.pl -o individual_prots/MGEN.faa --database protein_db.faa --prefix MGEN_animals --seqOffset 50 --evalue 1e-5 --threads 18 --blastPlus