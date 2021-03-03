#!/bin/bash
#SBATCH -n 1
#SBATCH --cpus-per-task=14
#SBATCH --time=6-23:59 --qos=1wk
#SBATCH --mem=200000
#SBATCH --mail-user=berubin@princeton.edu
#SBATCH --job-name=APUR_maker
#SBATCH --output=/Genomics/kocherlab/berubin/annotation/hic/APUR/maker/maker-%j.out
#SBATCH --error=/Genomics/kocherlab/berubin/annotation/hic/APUR/maker/maker-%j.err

module add mpich-x86_64
module add RepeatMasker
module add rmblast
module add exonerate
module add bioperl
module add boost
module add hmmer
export PATH=$PATH:/Genomics/kocherlab/berubin/local/src/ncbi-blast-2.5.0+/bin/
export AUGUSTUS_CONFIG_PATH=/Genomics/kocherlab/berubin/local/src/augustus/config/
export PATH=$PATH:/Genomics/kocherlab/berubin/local/src/augustus/scripts/
export PATH=/Genomics/kocherlab/berubin/local/src/augustus/bin/:$PATH
export PATH=$PATH:/Genomics/kocherlab/berubin/local/src/snap/
export PATH=$PATH:/Genomics/kocherlab/berubin/local/src/EVidenceModeler-1.1.1

rm -r /tmp/maker_*
rm -r /tmp/berubin_*

python /Genomics/kocherlab/berubin/annotation/reformat_gff.py /Genomics/kocherlab/berubin/annotation/hic/APUR/APUR_DB.pasa_assemblies.gff3 /Genomics/kocherlab/berubin/annotation/hic/APUR/APUR_DB.pasa_assemblies_form.gff3 APUR

if [ "x$SLURM_JOB_ID" != "x" ] && [ ! -e "/tmp/$SLURM_JOB_ID" ]; then
    mkdir "/tmp/berubin_$SLURM_JOB_ID"
fi

/Genomics/kocherlab/berubin/local/src/maker3/maker/bin/maker -fix_nucleotides -base APUR -TMP "/tmp/berubin_$SLURM_JOB_ID" &

sleep 20m

for thread in {1..13}
do
     /Genomics/kocherlab/berubin/local/src/maker3/maker/bin/maker -fix_nucleotides -base APUR -TMP "/tmp/berubin_$SLURM_JOB_ID" &
done
wait

rm -rf "/tmp/berubin_$SLURM_JOB_ID"

mkdir round_1_output_brakerpred
cd round_1_output_brakerpred

/Genomics/kocherlab/berubin/local/src/maker3/maker/bin/fasta_merge -d ../APUR.maker.output/APUR_master_datastore_index.log -o APUR

/Genomics/kocherlab/berubin/local/src/maker3/maker/bin/gff3_merge -d ../APUR.maker.output/APUR_master_datastore_index.log -o APUR_all.gff3

/Genomics/kocherlab/berubin/local/src/maker3/maker/bin/gff3_merge -d ../APUR.maker.output/APUR_master_datastore_index.log -o APUR_genes.gff3 -g

/Genomics/kocherlab/berubin/local/src/busco/BUSCO.py -i APUR.all.maker.proteins.fasta -o busco_prots -l /Genomics/kocherlab/berubin/local/src/busco/hymenoptera_odb9 -m proteins -c 14 -sp APUR_braker

/Genomics/kocherlab/berubin/local/src/interproscan-5.21-60.0/interproscan.sh -i APUR.all.maker.non_overlapping_ab_initio.proteins.fasta --seqtype p -goterms -iprlookup -b APUR.all.maker.non_overlapping_ab_initio.proteins.fasta.iprscan -pa

grep "IPR" APUR.all.maker.non_overlapping_ab_initio.proteins.fasta.iprscan.tsv | awk '{print $1}' | sort | uniq > ipr_hits.txt

/Genomics/kocherlab/berubin/local/src/maker3/maker/bin/gff3_select APUR_all.gff3 ipr_hits.txt > ipr_hits.gff3
cd ..

if [ "x$SLURM_JOB_ID" != "x" ] && [ ! -e "/tmp/$SLURM_JOB_ID" ]; then
   mkdir "/tmp/berubin_$SLURM_JOB_ID"
fi

/Genomics/kocherlab/berubin/local/src/maker3/maker/bin/maker -fix_nucleotides -base APUR_ipradd -TMP "/tmp/berubin_$SLURM_JOB_ID" maker_opts_ipradd.ctl

rm -r /tmp/berubin_$SLURM_JOB_ID

mkdir round_2_output_ipradd
cd round_2_output_ipradd

/Genomics/kocherlab/berubin/local/src/maker3/maker/bin/fasta_merge -d ../APUR_ipradd.maker.output/APUR_ipradd_master_datastore_index.log -o pura

/Genomics/kocherlab/berubin/local/src/maker3/maker/bin/gff3_merge -d ../APUR_ipradd.maker.output/APUR_ipradd_master_datastore_index.log -o pura_all.gff3

/Genomics/kocherlab/berubin/local/src/maker3/maker/bin/gff3_merge -d ../APUR_ipradd.maker.output/APUR_ipradd_master_datastore_index.log -o pura_genes.gff3 -g

/Genomics/kocherlab/berubin/local/src/maker3/maker/bin/maker_map_ids --prefix APUR_ --justify 5 pura_genes.gff3 > APUR_genes.ids

cp pura_genes.gff3 APUR_genes.gff3

/Genomics/kocherlab/berubin/local/src/maker3/maker/bin/map_gff_ids APUR_genes.ids APUR_genes.gff3

python /Genomics/kocherlab/berubin/annotation/make_name_set.py pura.all.maker.proteins.fasta pura_genes.gff3 APUR_genes.ids APUR_prots.ids

cp pura.all.maker.proteins.fasta APUR_proteins.fasta

/Genomics/kocherlab/berubin/local/src/maker3/maker/bin/map_fasta_ids APUR_prots.ids APUR_proteins.fasta

cp pura.all.maker.transcripts.fasta APUR_transcripts.fasta

/Genomics/kocherlab/berubin/local/src/maker3/maker/bin/map_fasta_ids APUR_prots.ids APUR_transcripts.fasta


python /Genomics/kocherlab/berubin/annotation/make_ogs_hic.py /Genomics/kocherlab/berubin/assembly/hic/APUR/APUR_genome_v2.0.fasta APUR_proteins.fasta APUR_transcripts.fasta APUR_genes.gff3 /Genomics/kocherlab/berubin/annotation/acacias/PCON/maker/nobac_scafs.txt /Genomics/kocherlab/berubin/official_release APUR

python /Genomics/kocherlab/berubin/annotation/longest_isoform_gff.py /Genomics/kocherlab/berubin/official_release/APUR/APUR_OGS_v2.1.gff3 /Genomics/kocherlab/berubin/official_release/APUR/APUR_OGS_v2.1_prots_longest_isoform.fasta /Genomics/kocherlab/berubin/official_release/APUR/APUR_OGS_v2.1_longest_isoform.gff3

dummy_line="dummy\tdummy\tmRNA\t0\t0\t.\t+\t.\tID=dummy;Parent=dummy;Name=dummy;"
cat /Genomics/kocherlab/berubin/official_release/APUR/APUR_OGS_v2.1_longest_isoform.gff3 > /Genomics/kocherlab/berubin/official_release/APUR/tmp_longest_isoform.gff3
echo -e ${dummy_line} >> /Genomics/kocherlab/berubin/official_release/APUR/tmp_longest_isoform.gff3

perl /Genomics/kocherlab/berubin/local/src/common_scripts/gff2fasta.pl /Genomics/kocherlab/berubin/official_release/APUR/APUR_genome_v2.1.fasta /Genomics/kocherlab/berubin/official_release/APUR/tmp_longest_isoform.gff3 /Genomics/kocherlab/berubin/official_release/APUR/APUR_OGS_v2.1_longest_isoform

rm /Genomics/kocherlab/berubin/official_release/APUR/tmp_longest_isoform.gff3
rm /Genomics/kocherlab/berubin/official_release/APUR/APUR_genome_v2.1.fasta.index
rm /Genomics/kocherlab/berubin/official_release/APUR/APUR_OGS_v2.1_longest_isoform.cdna.fasta
rm /Genomics/kocherlab/berubin/official_release/APUR/APUR_OGS_v2.1_trans_longest_isoform.fasta
mv /Genomics/kocherlab/berubin/official_release/APUR/APUR_OGS_v2.1_longest_isoform.exon.fasta /Genomics/kocherlab/berubin/official_release/APUR/APUR_OGS_v2.1_longest_isoform.trans.fasta
mv /Genomics/kocherlab/berubin/official_release/APUR/APUR_OGS_v2.1_prots.fasta /Genomics/kocherlab/berubin/official_release/APUR/APUR_OGS_v2.1_pep.fasta

python /Genomics/kocherlab/berubin/assembly/hic/make_release_genomes.py APUR

/Genomics/kocherlab/berubin/local/src/busco/BUSCO.py -i /Genomics/kocherlab/berubin/official_release_v2.1/APUR/APUR_OGS_v2.1_longest_isoform.pep.fasta -o busco_prots -l /Genomics/kocherlab/berubin/local/src/busco/hymenoptera_odb9 -m proteins -c 14 -sp APUR_braker -f
