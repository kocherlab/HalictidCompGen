#!/bin/bash
#SBATCH -n 1
#SBATCH --cpus-per-task=6
#SBATCH --time=6-23:00 --qos=1wk
#SBATCH --mem=50000
#SBATCH --mail-user=berubin@princeton.edu
#SBATCH --job-name=orthomcl
#SBATCH --output=/Genomics/kocherlab/berubin/annotation/hic/orthomcl/mcl-%j.out
#SBATCH --error=/Genomics/kocherlab/berubin/annotation/hic/orthomcl/mcl-%j.err


#To log in to mysql use
#mysql -u orthomcl -h knees.Princeton.EDU -p orthomcl
#Then can "create database orthomcl;"
#and "drop database orthomcl;" when you're done to keep space use low.

/Genomics/kocherlab/berubin/local/src/orthomclSoftware-v2.0.9/bin/orthomclInstallSchema hali.config hali.sqllog

cd compliant_fasta
for species in AAUR APUR AVIR HLIG HQUA HRUB LALB LCAL LFIG LLEU LMAL LMAR LOEN LPAU LVIE LZEP NMEL MGEN DNOV
do
    if [ ${species} == "DNOV" ] 
    then
	/Genomics/kocherlab/berubin/local/src/orthomclSoftware-v2.0.9/bin/orthomclAdjustFasta ${species} /Genomics/kocherlab/berubin/annotation/reference_genomes/Dufourea_novaeangliae_v1.1.pep.fa 1
    elif [ ${species} == "MGEN" ]
    then
	/Genomics/kocherlab/berubin/local/src/orthomclSoftware-v2.0.9/bin/orthomclAdjustFasta ${species} /Genomics/kocherlab/berubin/data/mgen/Mgen_v1.0.pep.fa 1
    else
	/Genomics/kocherlab/berubin/local/src/orthomclSoftware-v2.0.9/bin/orthomclAdjustFasta ${species} /Genomics/kocherlab/berubin/official_release_v2.1/${species}/${species}_OGS_v2.1_longest_isoform.pep.fasta 1
    fi
done

cd ..

/Genomics/kocherlab/berubin/local/src/orthomclSoftware-v2.0.9/bin/orthomclFilterFasta compliant_fasta/ 10 20

export PATH=$PATH:/Genomics/kocherlab/berubin/local/src/ncbi-blast-2.5.0+/bin/

makeblastdb -in goodProteins.fasta -out goodProteins.db -dbtype prot

blastp -query goodProteins.fasta -db goodProteins.db -dbsize 233411 -evalue 1e-5 -outfmt 6 -num_threads 16 -seg yes -num_descriptions 100000 -num_alignments 100000 -out hali_allvall.txt

/Genomics/kocherlab/berubin/local/src/orthomclSoftware-v2.0.9/bin/orthomclBlastParser hali_allvall.txt compliant_fasta/ >> similarSequences.txt

/Genomics/kocherlab/berubin/local/src/orthomclSoftware-v2.0.9/bin/orthomclLoadBlast hali.config similarSequences.txt
#for this step, the similarSequencesTable parameter in the config file has to correspond with the name of the database given in step one ("HALI"). So in this case it has to be "SimilarSequencesHALI".

/Genomics/kocherlab/berubin/local/src/orthomclSoftware-v2.0.9/bin/orthomclPairs hali.config orthomclPairs.log cleanup=yes

/Genomics/kocherlab/berubin/local/src/orthomclSoftware-v2.0.9/bin/orthomclDumpPairsFiles hali.config

mcl mclInput --abc -I 1.5 -o mclOutput

/Genomics/kocherlab/berubin/local/src/orthomclSoftware-v2.0.9/bin/orthomclMclToGroups OG_ 10000 < mclOutput > groups.txt