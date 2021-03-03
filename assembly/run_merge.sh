#!/bin/bash
#SBATCH -n 1
#SBATCH --cpus-per-task=20
#SBATCH --time=3-00:00 --qos=1wk
#SBATCH --mem=20000
#SBATCH --mail-user=berubin@princeton.edu
#SBATCH --job-name=fig_merge
#SBATCH --output=/Genomics/kocherlab/berubin/assembly/figueresi/merging/outputs/merging-%j.out
#SBATCH --error=/Genomics/kocherlab/berubin/assembly/figueresi/merging/errors/merging-%j.err

module add MUMmer
module add gcc
module add hmmer
module add funannotate
export PATH=$PATH:/Genomics/kocherlab/berubin/local/src/quickmerge/merger/
#export PATH=$PATH:/Genomics/kocherlab/berubin/local/src/quickmerge/MUMmer3.23/
export AUGUSTUS_CONFIG_PATH=/Genomics/kocherlab/berubin/local/src/augustus-3.2.2/config/
export PATH=$PATH:/Genomics/kocherlab/berubin/local/src/augustus-3.2.2/scripts/
export PATH=$PATH:/Genomics/kocherlab/berubin/local/src/augustus-3.2.2/bin/
export PATH=$PATH:/Genomics/kocherlab/berubin/local/src/ncbi-blast-2.5.0+/bin/
cur_species="figueresi"
scaf_pre="LFIG"
hco=7
quickc=2
minlen=10000
seedlen=2000000

python /Genomics/kocherlab/berubin/assembly/rename_scafs.py /Genomics/kocherlab/berubin/assembly/${cur_species}/supernova/${cur_species}_rand_1.fasta /Genomics/kocherlab/berubin/assembly/${cur_species}/merging/${cur_species}_rand_1.fa ${scaf_pre}_1 &
python /Genomics/kocherlab/berubin/assembly/rename_scafs.py /Genomics/kocherlab/berubin/assembly/${cur_species}/supernova/${cur_species}_rand_2.fasta /Genomics/kocherlab/berubin/assembly/${cur_species}/merging/${cur_species}_rand_2.fa ${scaf_pre}_2 &
python /Genomics/kocherlab/berubin/assembly/rename_scafs.py /Genomics/kocherlab/berubin/assembly/${cur_species}/supernova/${cur_species}_rand_3.fasta /Genomics/kocherlab/berubin/assembly/${cur_species}/merging/${cur_species}_rand_3.fa ${scaf_pre}_3 &
wait

export PERL5LIB=/Genomics/grid/users/berubin/perl5/Statistics-Descriptive-3.0612/lib:$PERL5LIB
perl /Genomics/kocherlab/berubin/local/src/N50stats.pl -in ${cur_species}_rand_1.fa -overwrite
perl /Genomics/kocherlab/berubin/local/src/N50stats.pl -in ${cur_species}_rand_2.fa -overwrite
perl /Genomics/kocherlab/berubin/local/src/N50stats.pl -in ${cur_species}_rand_3.fa -overwrite

mkdir 1h_2r
cd 1h_2r
/Genomics/kocherlab/berubin/local/src/quickmerge/merge_wrapper.py /Genomics/kocherlab/berubin/assembly/${cur_species}/merging/${cur_species}_rand_1.fa /Genomics/kocherlab/berubin/assembly/${cur_species}/merging/${cur_species}_rand_2.fa -pre ${cur_species}_1h_2r -hco ${hco} -c ${quickc} -lm ${minlen} -l ${seedlen} &
cd ..

mkdir 2h_1r
cd 2h_1r
/Genomics/kocherlab/berubin/local/src/quickmerge/merge_wrapper.py /Genomics/kocherlab/berubin/assembly/${cur_species}/merging/${cur_species}_rand_2.fa /Genomics/kocherlab/berubin/assembly/${cur_species}/merging/${cur_species}_rand_1.fa -pre ${cur_species}_2h_1r -hco ${hco} -c ${quickc} -lm ${minlen} -l ${seedlen} &
cd ..

mkdir 1h_3r
cd 1h_3r
/Genomics/kocherlab/berubin/local/src/quickmerge/merge_wrapper.py /Genomics/kocherlab/berubin/assembly/${cur_species}/merging/${cur_species}_rand_1.fa /Genomics/kocherlab/berubin/assembly/${cur_species}/merging/${cur_species}_rand_3.fa -pre ${cur_species}_1h_3r -hco ${hco} -c ${quickc} -lm ${minlen} -l ${seedlen} &
cd ..

mkdir 3h_1r
cd 3h_1r
/Genomics/kocherlab/berubin/local/src/quickmerge/merge_wrapper.py /Genomics/kocherlab/berubin/assembly/${cur_species}/merging/${cur_species}_rand_3.fa /Genomics/kocherlab/berubin/assembly/${cur_species}/merging/${cur_species}_rand_1.fa -pre ${cur_species}_3h_1r -hco ${hco} -c ${quickc} -lm ${minlen} -l ${seedlen} &
cd ..

mkdir 3h_2r
cd 3h_2r
/Genomics/kocherlab/berubin/local/src/quickmerge/merge_wrapper.py /Genomics/kocherlab/berubin/assembly/${cur_species}/merging/${cur_species}_rand_3.fa /Genomics/kocherlab/berubin/assembly/${cur_species}/merging/${cur_species}_rand_2.fa -pre ${cur_species}_3h_2r -hco ${hco} -c ${quickc} -lm ${minlen} -l ${seedlen} &
cd ..

mkdir 2h_3r
cd 2h_3r
/Genomics/kocherlab/berubin/local/src/quickmerge/merge_wrapper.py /Genomics/kocherlab/berubin/assembly/${cur_species}/merging/${cur_species}_rand_2.fa /Genomics/kocherlab/berubin/assembly/${cur_species}/merging/${cur_species}_rand_3.fa -pre ${cur_species}_2h_3r -hco ${hco} -c ${quickc} -lm ${minlen} -l ${seedlen} &
cd ..
wait

for reference in 1h_2r 2h_1r 1h_3r 3h_1r 3h_2r 2h_3r
do
    cd ${reference}
    /Genomics/kocherlab/berubin/local/src/N50stats.pl -in merged.fasta -overwrite
    cd ..
done

for hybrid in 1 2 3
do
    for reference in 1h_2r 2h_1r 1h_3r 3h_1r 3h_2r 2h_3r
    do
	mkdir ${hybrid}h_${reference}
	cd ${hybrid}h_${reference}
	/Genomics/kocherlab/berubin/local/src/quickmerge/merge_wrapper.py /Genomics/kocherlab/berubin/assembly/${cur_species}/merging/${cur_species}_rand_${hybrid}.fa /Genomics/kocherlab/berubin/assembly/${cur_species}/merging/${reference}/merged.fasta -pre ${cur_species}_${hybrid}h_${reference} -hco ${hco} -c ${quickc} -lm ${minlen} -l ${seedlen} &
	cd ..
    done
done
wait

for merging in 1h_1h_2r 1h_2h_1r 1h_1h_3r 1h_3h_1r 1h_3h_2r 1h_2h_3r 2h_1h_2r 2h_2h_1r 2h_1h_3r 2h_3h_1r 2h_3h_2r 2h_2h_3r 3h_1h_2r 3h_2h_1r 3h_1h_3r 3h_3h_1r 3h_3h_2r 3h_2h_3r 
do
    cd ${merging}
    /Genomics/kocherlab/berubin/local/src/N50stats.pl -in merged.fasta -overwrite
    cd ..
done

for index in 1 2 3
do
    /Genomics/kocherlab/berubin/local/src/busco/BUSCO.py -i ${cur_species}_rand_${index}.fa -o ${cur_species}_rand_${index}_busco -m geno -l /Genomics/kocherlab/berubin/local/src/busco/hymenoptera_odb9 -c 20 -sp honeybee1 -f
done
END

for merging in 1h_1h_2r 1h_2h_1r 1h_1h_3r 1h_3h_1r 1h_3h_2r 1h_2h_3r 2h_1h_2r 2h_2h_1r 2h_1h_3r 2h_3h_1r 2h_3h_2r 2h_2h_3r 3h_1h_2r 3h_2h_1r 3h_1h_3r 3h_3h_1r 3h_3h_2r 3h_2h_3r 1h_2r 2h_1r 1h_3r 3h_1r 3h_2r 2h_3r
do
    if [ ! -e /Genomics/kocherlab/berubin/assembly/${cur_species}/merging/${merging}/run_${merging}_busco/short_summary_${merging}_busco.txt ]
    then
	cd ${merging}
	/Genomics/kocherlab/berubin/local/src/busco/BUSCO.py -i merged.fasta -o ${merging}_busco -m geno -l /Genomics/kocherlab/berubin/local/src/busco/hymenoptera_odb9 -c 20 -sp honeybee1 -f
	cd ..
    fi
done
