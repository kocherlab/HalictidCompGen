#!/bin/bash
#SBATCH -n 1
#SBATCH --cpus-per-task=20
#SBATCH --time=6-23:00 --qos=1wk
#SBATCH --mem=360000
#SBATCH --mail-user=berubin@princeton.edu
#SBATCH --job-name=HQUA_1.0
#SBATCH --output=/Genomics/kocherlab/berubin/assembly/HQUA/supernova1.0/supernova-%j.out 
#SBATCH --error=/Genomics/kocherlab/berubin/assembly/HQUA/supernova1.0/supernova-%j.err

sample_num="22"
s_num="1"
barcode="AATAATGG"
chunk="1"
lane="4"
cur_species="HQUA"
reads_per_file=210000000

mkdir "first_trim_logs"
mkdir /scratch/tmp/berubin/assembly/${cur_species}/
mkdir /scratch/tmp/berubin/assembly/${cur_species}/supernova1.0/
outdir="/scratch/tmp/berubin/assembly/${cur_species}/supernova1.0/first_trim_data"
mkdir ${outdir}

/Genomics/grid/users/berubin/.local/bin/cutadapt --interleaved --discard-trimmed -a file:/Genomics/kocherlab/berubin/assembly/illumina_adapters_fifthtrim.fa -A file:/Genomics/kocherlab/berubin/assembly/illumina_adapters_fifthtrim.fa -o ${outdir}/read-RA_si-${barcode}_lane-00${lane}-chunk-00${chunk}.fastq.gz /Genomics/kocherlab/berubin/data/10x/december/3992-BER-00${sample_num}_S${s_num}_L00${lane}_R1_001.fastq.gz /Genomics/kocherlab/berubin/data/10x/december/3992-BER-00${sample_num}_S${s_num}_L00${lane}_R2_001.fastq.gz > ./first_trim_logs/infotrim_walker_everywhere_${barcode}_${lane}_${chunk}.log 

mkdir "second_trim_logs"
mkdir "second_trim_data"

/Genomics/grid/users/berubin/.local/bin/cutadapt --interleaved --discard-trimmed --overlap 4 -b file:/Genomics/kocherlab/berubin/assembly/illumina_adapters_fifthtrim.fa -B file:/Genomics/kocherlab/berubin/assembly/illumina_adapters_fifthtrim.fa -o ./second_trim_data/read-RA_si-${barcode}_lane-00${lane}-chunk-00${chunk}.fastq.gz ${outdir}/read-RA_si-${barcode}_lane-00${lane}-chunk-00${chunk}.fastq.gz > ./second_trim_logs/infotrim_walker_everywhere_${barcode}_${lane}_${chunk}.log

rm -r ${outdir}

for random in 1 2 3
do
    outdir="/scratch/tmp/berubin/assembly/${cur_species}/supernova1.0/second_trim_rand_${random}"
    mkdir ${outdir}
    zcat second_trim_data/read-RA_si-${barcode}_lane-00${lane}-chunk-00${chunk}.fastq.gz | awk '{printf("%s%s",$0,(NR%8==0)?"\n":"\0")}' | shuf | head -n ${reads_per_file} | tr "\0" "\n" | gzip > ${outdir}/read-RA_si-${barcode}_lane-00${lane}-chunk-00${chunk}.fastq.gz 
    ln -s /Genomics/kocherlab/berubin/data/10x/december/3992-BER-00${sample_num}_S${s_num}_L00${lane}_I1_001.fastq.gz ${outdir}/read-I1_si-${barcode}_lane-00${lane}-chunk-00${chunk}.fastq.gz
    cd /scratch/tmp/berubin/assembly/${cur_species}/supernova1.0/

    /Genomics/kocherlab/berubin/local/src/supernova-1.0.0/supernova run --id=${cur_species}_rand_${random} --fastqs=${outdir} --maxreads=420000000 &> supernova_${cur_species}_rand_${random}.log

    /Genomics/kocherlab/berubin/local/src/supernova-1.0.0/supernova mkfasta --asmdir=${cur_species}_rand_${random}/outs/assembly --outprefix=/Genomics/kocherlab/berubin/assembly/${cur_species}/supernova1.0/${cur_species}_rand_${random} --style pseudohap --minsize=1000 --headers=short

done