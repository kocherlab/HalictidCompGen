#!/bin/bash
#SBATCH -n 1
#SBATCH --cpus-per-task=16
#SBATCH --time=6-23:00 --qos=1wk
#SBATCH --mem=50000
#SBATCH --mail-user=berubin@princeton.edu
#SBATCH --job-name=ortho
#SBATCH --output=/Genomics/kocherlab/berubin/comparative/halictids/orthology/orthofinder/outputs/orthofinder-%j.out
#SBATCH --error=/Genomics/kocherlab/berubin/comparative/halictids/orthology/orthofinder/errors/orthofinder-%j.err

export PATH=$PATH:/Genomics/kocherlab/berubin/local/src/OrthoFinder/depedencies/mmseqs2/bin/
export PATH=$PATH:/Genomics/kocherlab/berubin/local/src/OrthoFinder/dependencies/

/Genomics/kocherlab/berubin/local/src/OrthoFinder/orthofinder/orthofinder.py -t 16 -a 16 -s /Genomics/kocherlab/berubin/comparative/halictids/orthology/orthofinder/species.tree -f /Genomics/kocherlab/berubin/comparative/halictids/orthology/orthofinder/protein_seqs 