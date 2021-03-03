species=$1

export PATH=$PATH:/home/ben/local/src/blatSrc/bin
export PATH=$PATH:/home/ben/local/src/fasta-36.3.8e/bin

mkdir ${species}

cp alignAssembly.config ${species}/alignAssembly.config
sed -i "s/SPECIESNAME/${species}/g" ${species}/alignAssembly.config

python format_transcriptomes.py ${species} #make sure that transcript names are unique

cd ${species}

cat ${species}_denovo.fasta ${species}_gg.fasta > ${species}_transcripts.fasta

/home/ben/local/src/PASApipeline-2.0.2/seqclean/seqclean/seqclean ${species}_transcripts.fasta

/home/ben/local/src/PASApipeline-2.0.2/misc_utilities/accession_extractor.pl < ${species}_denovo.fasta > tdn.accs
    
/home/ben/local/src/PASApipeline-2.0.2/scripts/Launch_PASA_pipeline.pl -c alignAssembly.config -C -r -R -g ${species}_genome_v2.1.fasta -t ${species}_transcripts.fasta.clean -T -u ${species}_transcripts.fasta --ALIGNERS gmap,blat --CPU 6 -I 100000 --TDN tdn.accs --TRANSDECODER &> launch_pasa_pipe.log

/home/ben/local/src/PASApipeline-2.0.2/scripts/pasa_asmbls_to_training_set.dbi --pasa_transcripts_fasta ${species}_DB.assemblies.fasta --pasa_transcripts_gff3 ${species}_DB.pasa_assemblies.gff3
