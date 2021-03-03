#!/bin/bash
#SBATCH -n 1
#SBATCH --cpus-per-task=10
#SBATCH --time=6-23:00 --qos=1wk
#SBATCH --mem=200000
#SBATCH --mail-user=berubin@princeton.edu
#SBATCH --job-name=oen_seal
#SBATCH --output=/Genomics/kocherlab/berubin/assembly/oenotherae/sealer/sealer-%j.out
#SBATCH --error=/Genomics/kocherlab/berubin/assembly/oenotherae/sealer/sealer-%j.err

module add boost
cur_species="oenotherae"
merged_genome="oenotherae_1h_3h_2r"

barcode="ACTTCATA"
chunk="1"

out_data_dir="/Genomics/kocherlab/berubin/assembly/${cur_species}/sealer/no10x_data"
second_trim_dir="/Genomics/kocherlab/berubin/assembly/${cur_species}/supernova/second_trim_data"

mkdir ${out_data_dir}

/Genomics/grid/users/berubin/.local/bin/cutadapt --interleaved --cut 24 -q 20,20 --trim-n -m 115 -o ${out_data_dir}/read-RA_si-${barcode}_lane-00${lane}-chunk-00${chunk}_no10x_1.fastq.gz -p ${out_data_dir}/read-RA_si-${barcode}_lane-00${lane}-chunk-00${chunk}_no10x_2.fastq.gz ${second_trim_dir}/read-RA_si-${barcode}_lane-00${lane}-chunk-00${chunk}.fastq.gz &

zcat /Genomics/kocherlab/berubin/assembly/${cur_species}/sealer/no10x_data/read-RA_si-*_lane-00*-chunk-00*_no10x_2.fastq.gz | gzip >> read_no10x_2.fastq.gz &
zcat /Genomics/kocherlab/berubin/assembly/${cur_species}/sealer/no10x_data/read-RA_si-*_lane-00*-chunk-00*_no10x_1.fastq.gz | gzip >> read_no10x_1.fastq.gz &
wait

/Genomics/kocherlab/berubin/local/abyss2/bin/abyss-sealer -S /Genomics/kocherlab/berubin/assembly/${cur_species}/sealer/${merged_genome}.fa -B 3000 -G 1000 -P 10 -L 125 -j 10 -k125 -k120 -k115 -k110 -k105 -k100 -k95 -k90 -k85 -k80 -k75 -k70 -k65 -k60 -k55 -k50 -b 164G -o ${merged_genome}_sealer read_no10x_1.fastq.gz read_no10x_2.fastq.gz



