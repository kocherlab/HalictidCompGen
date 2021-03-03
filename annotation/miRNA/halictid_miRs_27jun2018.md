# Halictid miR discovery and analysis
### Karen Kapheim
### Summer 2018

## Script and data location

*Working directory:* `/uufs/chpc.utah.edu/common/home/kapheim-group1/Kocher_miRNA`
*Script directory:* `/uufs/chpc.utah.edu/common/home/kapheim-group1/Kocher_miRNA/scripts`
*Job reports directory:* `/uufs/chpc.utah.edu/common/home/kapheim-group1/Kocher_miRNA/jobreports`
*Small RNA sequences:* `/uufs/chpc.utah.edu/common/home/kapheim-group1/Kocher_miRNA/miRNA`
*Genome assemblies:* `/uufs/chpc.utah.edu/common/home/kapheim-group1/Kocher_miRNA/genomes`

## Notes on nomenclature
Stem-loop in miRBase is the pre-mir (Griffiths-Jones et al., 2006).
The pre-mir is excised from a longer primary transcript (pri-mir).
The mature sequences are designated ‘miR’ and the precursor hairpins are labelled ‘mir’.
Suffixes indicate paralogs that differ at only one or two positions (Griffiths-Jones et al., 2006).

## Initial run

### Acquire known mature miRNA sequences
Steps 1-4 were initially completed as part of an earlier project aimed at identifying miRs in 6 other bee species.
For the current project, I added the final set of miRs from these 6 bees to the original known mirs file.
Saved the final combined file in `known_miRs` directory within working directory.

1. miRBase (v21)  (all mature miRNA seqs) (Griffiths-Jones et al., 2006, Griffiths-Jones et al., 2008, Kozomara & Griffiths-Jones, 2014).
	* Apis mellifera, 29 July 2016
	* Drosophila melanogaster, 02 Aug 2016
	* Nasonia vitripennis, 02 Aug 2016
	* Tribolium castenum, 02 Aug 2016
	* Bombyx mori, 02 Aug 2016
2. Ashby et al. 2016 Scientific Reports
	* Downloaded Table S1
	* Used sequence converter website (http://sequenceconversion.bugaco.com) to convert to fasta
	* Appended this to the mirs from miRBase in a file called `amel_mature_miRs_combined.fa`
3. Needed to modify the identifying line in the miRBase files to get rid of white spaces:

```
cat tcas_mature_miRs_miRBase21_02aug2016.fa > tcas_mature_miRs_miRBase21_02aug2016_cp.fa
cat tcas_mature_miRs_miRBase21_02aug2016_cp.fa | tr " " _ > tcas_mature_miRs_miRBase21_02aug2016_mod.fa
```
4. Combined all of the above into a single file called `insect_mature_miRs.fa`

```
cat amel_mature_miRs_combined.fa > insect_mature_miRs.fa
cat bmor_mature_miRs_miRBase21_02aug2016_mod.fa >> insect_mature_miRs.fa
cat dmel_mature_miRs_miRBase21_02aug2016_mod.fa >> insect_mature_miRs.fa
cat nvit_mature_miRs_miRBase21_02aug2016_mod.fa >> insect_mature_miRs.fa
cat tcas_mature_miRs_miRBase21_02aug2016_mod.fa >> insect_mature_miRs.fa
```

5. Add discovered bee miRs from earlier project (unpubl data)

First add species identifiers to the headers.

```
cd /uufs/chpc.utah.edu/common/home/kapheim-group1/mirdeep2_0_0_8/knownMIRs
sed 's/^>/>mrot-/g' mrot_novel_mature_homologs_final_nowhite_rna.fa > mrot_novel_mature_homologs_final_nowhite_rna_id.fa
sed 's/^>/>mgen-/g' mgen_novel_mature_homologs_final_nowhite_rna.fa > mgen_novel_mature_homologs_final_nowhite_rna_id.fa
sed 's/^>/>nmel-/g' nmel_novel_mature_homologs_final_nowhite_rna.fa > nmel_novel_mature_homologs_final_nowhite_rna_id.fa
sed 's/^>/>amel-/g' amel_novel_mature_homologs_final_nowhite_rna.fa > amel_novel_mature_homologs_final_nowhite_rna_id.fa
sed 's/^>/>bimp-/g' bimp_novel_mature_homologs_final_nowhite_rna.fa > bimp_novel_mature_homologs_final_nowhite_rna_id.fa
sed 's/^>/>bter-/g' bter_novel_mature_mirs_final_rna_nowhite.fa > bter_novel_mature_final_nowhite_rna_id.fa
```

Now combine files and copy to working directory.

```
cd /uufs/chpc.utah.edu/common/home/kapheim-group1/mirdeep2_0_0_8/knownMIRs
ml mirdeep2
remove_white_space_in_id.pl insect_mature_miRs.fa > insect_mature_miRs_nowhite.fa
cat insect_mature_miRs_nowhite.fa > insect_bees_mature_miRs_nowhite.fa
cat mrot_novel_mature_homologs_final_nowhite_rna_id.fa >> insect_bees_mature_miRs_nowhite.fa
cat mgen_novel_mature_homologs_final_nowhite_rna_id.fa >> insect_bees_mature_miRs_nowhite.fa
cat nmel_novel_mature_homologs_final_nowhite_rna_id.fa >> insect_bees_mature_miRs_nowhite.fa
cat amel_novel_mature_homologs_final_nowhite_rna_id.fa >> insect_bees_mature_miRs_nowhite.fa
cat bimp_novel_mature_homologs_final_nowhite_rna_id.fa >> insect_bees_mature_miRs_nowhite.fa
cat bter_novel_mature_final_nowhite_rna_id.fa >> insect_bees_mature_miRs_nowhite.fa
cp insect_bees_mature_miRs_nowhite.fa /uufs/chpc.utah.edu/common/home/kapheim-group1/Kocher_miRNA/known_miRs/
```

### Prepare genomes
Modify headers of the Mgen genome assembly to remove white spaces.
All other genomes seem to have clean headers.

```
cp MGEN_genome_v2.0.fasta cpMGEN_genome_v2.0.fasta
cat cpMGEN_genome_v2.0.fasta | tr " " _ > MGEN_genome_v2.0_modhead.fasta
rm cpMGEN_genome_v2.0.fasta
```

### Prepare sequences

1. Concatenate sequences
Each library was run split and run across two lanes.
After QC check (fastqc), decided to merge these prior to any other trimming or mapping.

```
gunzip *R1.fastq.gz
mkdir catseqs
cat AAUR_1_R1.fastq AAUR_2_R1.fastq > ./catseqs/AAUR_cat_R1.fastq
cat APUR_1_R1.fastq APUR_2_R1.fastq > ./catseqs/APUR_cat_R1.fastq
cat AVIR_1_R1.fastq AVIR_2_R1.fastq > ./catseqs/AVIR_cat_R1.fastq
cat HLIG_1_R1.fastq HLIG_2_R1.fastq > ./catseqs/HLIG_cat_R1.fastq
cat HRUB_1_R1.fastq HRUB_2_R1.fastq > ./catseqs/HRUB_cat_R1.fastq
cat LALB_SOC_1_R1.fastq LALB_SOC_2_R1.fastq > ./catseqs/LALB_SOC_cat_R1.fastq
cat LALB_SOL_1_R1.fastq LALB_SOL_2_R1.fastq > ./catseqs/LALB_SOL_cat_R1.fastq
cat LCAL_1_R1.fastq LCAL_2_R1.fastq > ./catseqs/LCAL_cat_R1.fastq
cat LFIG_1_R1.fastq LFIG_2_R1.fastq > ./catseqs/LFIG_cat_R1.fastq
cat LLEU_1_R1.fastq LLEU_2_R1.fastq > ./catseqs/LLEU_cat_R1.fastq
cat LMAL_1_R1.fastq LMAL_2_R1.fastq > ./catseqs/LMAL_cat_R1.fastq
cat LMAR_1_R1.fastq LMAR_2_R1.fastq > ./catseqs/LMAR_cat_R1.fastq
cat LOEN_1_R1.fastq LOEN_2_R1.fastq > ./catseqs/LOEN_cat_R1.fastq
cat LPAU_1_R1.fastq LPAU_2_R1.fastq > ./catseqs/LPAU_cat_R1.fastq
cat LVIE_1_R1.fastq LVIE_2_R1.fastq > ./catseqs/LVIE_cat_R1.fastq
cat LZEP_1_R1.fastq LZEP_2_R1.fastq > ./catseqs/LZEP_cat_R1.fastq

```



### Build the index files for each genome

example script: `bowtiebuild_AAUR.slurm`

```
#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH -c 4
#SBATCH --account=kapheim-kp
#SBATCH --partition=kapheim-kp
#SBATCH -o slurm-%j.out-%N
#SBATCH -e slurm-%j.err-%N
#
#SET VARS
GENOME=AAUR_genome_v2.0.fasta
INDEX=AAUR_v2.0
DIR=/uufs/chpc.utah.edu/common/home/kapheim-group1/Kocher_miRNA/genomes
#
#LOAD MODULES
ml mirdeep2
#
cd $DIR
#
echo "running MirDeep2 bowtie-build for $INDEX"
#
#run
bowtie-build $GENOME $INDEX
#
echo "complete"
```


| species | script | jobreport |
| --- | --- | --- |
| AAUR | bowtiebuild_AAUR.slurm | slurm-5302006.out-kp292 |
| APUR | bowtiebuild_APUR.slurm | slurm-5302025.out-kp031 |
| AVIR | bowtiebuild_AVIR.slurm | slurm-5302026.out-kp031 |
| HLIG | bowtiebuild_HLIG.slurm | slurm-5302027.out-kp166 |
| HRUB | bowtiebuild_HRUB.slurm | slurm-5302028.out-kp110 |
| LALB | bowtiebuild_LALB.slurm | slurm-5302029.out-kp166 |
| LCAL | bowtiebuild_LCAL.slurm | slurm-5302030.out-kp196 |
| LFIG | bowtiebuild_LFIG.slurm | slurm-5302031.out-kp110 |
| LLEU | bowtiebuild_LLEU.slurm | slurm-5302032.out-kp031 |
| LMAL | bowtiebuild_LMAL.slurm | slurm-5302033.out-kp196 |
| LMAR | bowtiebuild_LMAR.slurm | slurm-5302034.out-kp166 |
| LOEN | bowtiebuild_LOEN.slurm | slurm-5302035.out-kp022 |
| LPAU | bowtiebuild_LPAU.slurm | slurm-5302036.out-kp110 |
| LVIE | bowtiebuild_LVIE.slurm | slurm-5302037.out-kp031 |
| LZEP | bowtiebuild_LZEP.slurm | slurm-5302038.out-kp196 |

### Process and map microRNA reads to genome

Quality filtering and trimming can be done in one step.

example script: `mapper_AAUR.slurm`

```
#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH -c 6
#SBATCH --account=kapheim
#SBATCH --partition=kingspeak
#SBATCH -o slurm-%j.out-%N
#SBATCH -e slurm-%j.err-%N
#
#LOAD MODULES
ml mirdeep2
#
#SET VARS
SPECIES=AAUR
INDEX=AAUR_v2.0
SEQS=AAUR_cat_R1.fastq
GENDIR=/uufs/chpc.utah.edu/common/home/kapheim-group1/Kocher_miRNA/genomes
SEQDIR=/uufs/chpc.utah.edu/common/home/kapheim-group1/Kocher_miRNA/miRNA/catseqs
OUTDIR1=/uufs/chpc.utah.edu/common/home/kapheim-group1/Kocher_miRNA/miRNA/processed_reads
OUTDIR2=/uufs/chpc.utah.edu/common/home/kapheim-group1/Kocher_miRNA/mapped_reads
#
#
cd $GENDIR
#
echo "running MirDeep2 mapper.pl - trimming and removing adapters for ${SPECIES}"
#
###NOTES
# -e denotes the input files are .fastq
# -j removes reads with noncanonical letters
# -k removes adapter sequence (taken from https://www.neb.com/faqs/2017/07/17/how-should-my-nebnext-small-rna-library-be-trimmed)
# -l removes reads shorter than 18 nt
# -m collapse reads
#-h parse to fasta format
# -s output the processed reads
# -t output the read mappings
# -v output progress
###Run
mapper.pl ${SEQDIR}/$SEQS -e -j -k AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC  -l 18 -m -p $INDEX -h -s ${OUTDIR1}/${SPECIES}_reads_collapsed.fa -t ${OUTDIR2}/${SPECIES}_reads_collapsed_vs_genome.arf -v
#
echo "complete"
```

Mapping results:

| Species | Total reads | Mapped |Unmapped | % mapped | % unmapped | Logfile | Jobreport |
| --- | --- | --- | --- | --- | --- | --- | --- |
| AAUR | 7302517 | 5475057 | 1827460 | 0.750 | 0.250 | mapper.log_23315 | slurm-5303224.err-kp021 |
| APUR | 9564928 | 5351218 | 4213710 | 0.559 | 0.441 | mapper.log_29235 | slurm-5303251.err-kp006 |
| AVIR | 16365653 | 12560170  | 3805483 | 0.767 | 0.233 | mapper.log_7780 | slurm-5303165.err-kp199 |
| HLIG | 7014763 | 4942540 | 2072223 | 0.705 | 0.295 | mapper.log_8792 | slurm-5303166.err-kp002 |
| HRUB | 16916305 | 11221470 | 5694835 | 0.663 | 0.337 | mapper.log_8607 | slurm-5303167.err-kp199 |
| LALB_SOC | 12203318 | 9258676 | 2944642 | 0.759 | 0.241 | mapper.log_22210 | slurm-5303168.err-kp021 |
| LALB_SOL | 14926162 | 11179935 | 3746227 | 0.749 | 0.251 | mapper.log_28149 | slurm-5303169.err-kp006 |
| LCAL | 11358364 | 8452184 | 2906180 | 0.744 | 0.256 | mapper.log_9169 | slurm-5303170.err-kp002 |
| LFIG | 6516928 | 5079523 | 1437405 | 0.779 | 0.221 | mapper.log_108677 | slurm-5303171.err-kp159 |
| LLEU | 12232074 | 9488941 | 2743133 | 0.776 | 0.224 | mapper.log_22759 | slurm-5303172.err-kp021 |
| LMAL | 19297078 | 14125594 | 5171484 | 0.732 | 0.268 | mapper.log_9717 | slurm-5303173.err-kp002 |
| LMAR | 15044076 | 9080923 | 5963153 | 0.604 | 0.396 | mapper.log_109088 | slurm-5303174.err-kp159 |
| LOEN | 8100471 | 6595739 | 1504732 | 0.814 | 0.186 | mapper.log_28797 | slurm-5303175.err-kp006 |
| LPAU | 6662136 | 4289164 | 2372972 | 0.644 | 0.356 | mapper.log_9472 | slurm-5303176.err-kp199 |
| LVIE | 7368310 | 4967314 | 2400996 | 0.674 | 0.326 | mapper.log_28666 | slurm-5303177.err-kp001 |
| LZEP | 4952116 | 3557570 | 1394546 | 0.718 | 0.282 | mapper.log_7880 | slurm-5302997.err-kp292 |

### miRNA detection

example script: `md2_APUR.slurm`

```
#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH -c 6
#SBATCH --account=kapheim
#SBATCH --partition=kingspeak
#SBATCH -o slurm-%j.out-%N
#SBATCH -e slurm-%j.err-%N
#
#LOAD MODULES
ml mirdeep2
#
#SET VARS
SPECIES=APUR
WORKDIR=/uufs/chpc.utah.edu/common/home/kapheim-group1/Kocher_miRNA
#
#
cd $WORKDIR
#
echo "running MirDeep2 miRDeep2.pl to identify miRNAs for ${SPECIES}"
#
#run
miRDeep2.pl ./miRNA/processed_reads/${SPECIES}_reads_collapsed.fa ./genomes/${SPECIES}_genome_v2.0.fasta ./mapped_reads/${SPECIES}_reads_collapsed_vs_genome.arf none ./known_miRs/insect_bees_mature_miRs_nowhite.fa none -r _${SPECIES} -z _${SPECIES}  2>${SPECIES}_report.log
#
echo "complete"
```

| species | jobreport | error |
| --- | --- | --- |
| AAUR  | slurm-5305608.out-kp292 |
| APUR    | slurm-5305616.out-kp158 | genome fasta contains not allowed characters in sequences |
| AVIR    | slurm-5305617.out-kp158  | none |
| HLIG    |  slurm-5305621.out-kp160 | genome fasta contains not allowed characters in sequences |
| HRUB    | slurm-5305623.out-kp160 | none |
| LALB_SOC    | slurm-5305624.out-kp021 | first line of genome fasta does not start with '>identifier' |
| LALB_SOL    | slurm-5305625.out-kp021 | first line of genomes fasta does not start with '>identifier' |
| LCAL    | slurm-5305626.out-kp021 | genome fasta contains not allowed characters in sequences |
| LFIG    | slurm-5305627.out-kp021 | none |
| LLEU    | slurm-5305628.out-kp111 | genome fasta contains not allowed characters in sequences |
| LMAL    | slurm-5305629.out-kp111 | genome fasta contains not allowed characters in sequences |
| LMAR    | slurm-5305630.out-kp111 | genome fasta contains not allowed characters in sequences |
| LOEN    | slurm-5305631.out-kp111 | genome fasta contains not allowed characters in sequences |
| LPAU    | slurm-5305632.out-kp111 | genome fasta contains not allowed characters in sequences |
| LVIE    | slurm-5305633.out-kp111 | genome fasta contains not allowed characters in sequences |
| LZEP    | slurm-5305634.out-kp111 | genome fasta contains not allowed characters in sequences |

Fix the genome files and rerun.

For each species:
```
sed -e '/^[^>]/s/[^ACGTNacgtn]/N/g' APUR_genome_v2.0.fasta > APUR_genome_v2.0_clean.fasta
sed -e '/^[^>]/s/[^ACGTNacgtn]/N/g' HLIG_genome_v2.0.fasta > HLIG_genome_v2.0_clean.fasta
sed -e '/^[^>]/s/[^ACGTNacgtn]/N/g' LALB_genome_v2.0.fasta > LALB_genome_v2.0_clean.fasta
sed -e '/^[^>]/s/[^ACGTNacgtn]/N/g' LCAL_genome_v2.0.fasta > LCAL_genome_v2.0_clean.fasta
sed -e '/^[^>]/s/[^ACGTNacgtn]/N/g' LLEU_genome_v2.0.fasta > LLEU_genome_v2.0_clean.fasta
sed -e '/^[^>]/s/[^ACGTNacgtn]/N/g' LMAL_genome_v2.0.fasta > LMAL_genome_v2.0_clean.fasta
sed -e '/^[^>]/s/[^ACGTNacgtn]/N/g' LMAR_genome_v2.0.fasta > LMAR_genome_v2.0_clean.fasta
sed -e '/^[^>]/s/[^ACGTNacgtn]/N/g' LOEN_genome_v2.0.fasta > LOEN_genome_v2.0_clean.fasta
sed -e '/^[^>]/s/[^ACGTNacgtn]/N/g' LPAU_genome_v2.0.fasta > LPAU_genome_v2.0_clean.fasta
sed -e '/^[^>]/s/[^ACGTNacgtn]/N/g' LVIE_genome_v2.0.fasta > LVIE_genome_v2.0_clean.fasta
sed -e '/^[^>]/s/[^ACGTNacgtn]/N/g' LZEP_genome_v2.0.fasta > LZEP_genome_v2.0_clean.fasta
sanity_check_genome.pl LALB_genome_v2.0_clean.fasta
```
Reruns are in `md2_species_r2.slurm`

| species | jobreport | error |
| --- | --- | --- |
| APUR | slurm-5314321.out-kp292 |
| HLIG    | slurm-5314392.out-kp158 |
| LCAL    | slurm-5314393.out-kp013  |
| LLEU    | slurm-5314406.out-kp166 |
| LMAL    | slurm-5314407.out-kp012 |
| LMAR    | slurm-5314410.out-kp032  |
| LOEN    | slurm-5314411.out-kp005  |
| LPAU    | slurm-5314412.out-kp004 |
| LVIE    | slurm-5314413.out-kp167  |
| LZEP    | slurm-5314414.out-kp165  |
| LALB_SOC   | slurm-5314987.out-kp011  |
| LALB_SOL   | slurm-5314988.out-kp022 |

Mapping mature, star, and loop seqs against RFAM:

| Species | Reads processed | Reads with at least one reported alignment | % reads with at least one reported alignment | Reads that failed to align | % reads that failed to align | Reported alignment | subfolder |
| ---- | ---- | --- | --- | --- | --- | --- | --- |
| AAUR |  4961  |  21 | 0.42  | 4940  | 99.58  | 58  | mirna_results_27_06_2018_t_21_54_10_AAUR |
| AVIR   |  5721 | 25  | 0.44  | 5696  | 99.56  | 139  |  mirna_results_27_06_2018_t_22_03_58_AVIR |
| HRUB   |  6511 |  26 | 0.40  | 6485  | 99.60  | 91  | mirna_results_27_06_2018_t_22_06_15_HRUB  |
| LFIG   | 4012  | 15  | 0.37  | 3997  | 99.63  | 77  | mirna_results_27_06_2018_t_22_11_06_LFIG  |
| APUR   | 4129  | 17  | 0.41  | 4112  | 99.59  | 2660  | mirna_results_28_06_2018_t_15_00_19_APUR  |
| HLIG   | 5958  |  24 | 0.40  | 5934  | 99.60  | 107  |  mirna_results_28_06_2018_t_15_36_34_HLIG |
| LALB_SOC   | 6401  | 30  | 0.47  | 6371  | 99.53  | 347  | mirna_results_28_06_2018_t_17_50_18_LALB_SOC  |
| LALB_SOL    | 6634  | 31  | 0.47  | 6603  | 99.53  | 348  | mirna_results_28_06_2018_t_17_50_19_LALB_SOL  |
| LCAL    | 6515  | 19  | 0.29  | 6496  | 99.71  | 58  | mirna_results_28_06_2018_t_15_37_21_LCAL  |
| LFIG   | 4012  | 15  | 0.37  | 3997  | 99.63  | 77  | mirna_results_27_06_2018_t_22_11_06_LFIG  |
| LLEU   | 5039  | 28  | 0.56  | 5011  | 99.44  | 162  | mirna_results_28_06_2018_t_15_42_17_LLEU  |
| LMAL    | 9317  | 32  | 0.34  | 9285  | 99.66  | 136  | mirna_results_28_06_2018_t_15_48_39_LMAL  |
| LMAR    | 7176  | 29  | 0.40  | 7147  | 99.60  | 408  | mirna_results_28_06_2018_t_15_51_30_LMAR  |
| LOEN   | 3568  | 22  | 0.62  | 3546  | 99.38  | 58  | mirna_results_28_06_2018_t_15_54_52_LOEN  |
| LPAU   | 11606  | 68  | 0.59  | 11538  | 99.41  | 352  | mirna_results_28_06_2018_t_15_59_26_LPAU  |
| LVIE    | 13031  | 50  | 0.38  | 12981  | 99.62  | 1352  | mirna_results_28_06_2018_t_16_04_54_LVIE  |
| LZEP    | 10912  | 41  | 0.38  | 10871  | 99.62  | 366  | mirna_results_28_06_2018_t_16_10_17_LZEP  |

### Quantify expression

example script: `quant_LFIG.slurm`

```
#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH -c 6
#SBATCH --account=kapheim
#SBATCH --partition=kingspeak
#SBATCH -o slurm-%j.out-%N
#SBATCH -e slurm-%j.err-%N
#
#LOAD MODULES
ml mirdeep2
#
#SET VARS
SPECIES=LFIG
WORKDIR=/uufs/chpc.utah.edu/common/home/kapheim-group1/Kocher_miRNA
PRE=./mirna_results_27_06_2018_t_22_11_06_LFIG/novel_pres_27_06_2018_t_22_06_15_HRUB_score-50_to_na.fa
MAT=./mirna_results_27_06_2018_t_22_11_06_LFIG/novel_mature_27_06_2018_t_22_06_15_HRUB_score-50_to_na.fa
STAR=./mirna_results_27_06_2018_t_22_11_06_LFIG/novel_star_27_06_2018_t_22_06_15_HRUB_score-50_to_na.fa
#
#
cd $WORKDIR
#
echo "running MirDeep2 quantifier.pl to quantify expression of  miRNAs in ${SPECIES}"
#
#run
quantifier.pl -p $PRE -m $MAT -r ./miRNA/processed_reads/${SPECIES}_reads_collapsed.fa -s $STAR
#
echo "complete"
```

| species | reads | mapped | unmapped | % mapped | % unmapped | file | jobreport |
| --- | --- | --- | --- | --- | --- | --- | --- |
| AAUR | 7219044 | 4489156 | 2729888 | 0.622 | 0.378 | miRNAs_expressed_all_samples_1530201408.csv  | slurm-5308736.err-kp292 |
| AVIR | 16187328 | 6208909 | 9978419 | 0.384 | 0.616 | miRNAs_expressed_all_samples_1530223843.csv | slurm-5314532.err-kp027 |
| HRUB | 16727989 | 6785996 | 9941993 | 0.406 | 0.594 | miRNAs_expressed_all_samples_1530223864.csv | slurm-5314533.err-kp027 |
| APUR | 9454164 | 4373159 | 5081005 | 0.463 | 0.537 | miRNAs_expressed_all_samples_1530295135.csv | slurm-5318792.err-kp199 |
| HLIG   | 6930609  | 4019878  | 2910731  |  0.580 |  0.420 | miRNAs_expressed_all_samples_1530295153.csv  | slurm-5318793.err-kp199 |
| LALB_SOC  | 12064339  | 6918808  | 5145531  | 0.573  | 0.427  | miRNAs_expressed_all_samples_1530295170.csv  | slurm-5318794.err-kp199  |
| LALB_SOL   | 14759668  | 6531633  | 8228035  | 0.443  | 0.557  |  miRNAs_expressed_all_samples_1530295189.csv | slurm-5318795.err-kp199 |
| LCAL   | 11230133  | 6003544  | 5226589  |  0.535 | 0.465  | miRNAs_expressed_all_samples_1530295212.csv  | slurm-5318796.err-kp199 |
| LFIG   | 6442402  | 4027109  | 2415293  |  0.625 | 0.375  | miRNAs_expressed_all_samples_1530295234.csv  | slurm-5318797.err-kp199  |
| LLEU   | 12096673  | 7171286  | 4925387  | 0.593  | 0.407  | miRNAs_expressed_all_samples_1530295249.csv  | slurm-5318798.err-kp199  |
| LMAL   | 19075315  | 9561455  | 9513860  |  0.501 |  0.499 | miRNAs_expressed_all_samples_1530295277.csv  | slurm-5318799.err-kp199  |
| LMAR   | 14870831  |  6597102 | 8273729  | 0.444  | 0.556  | miRNAs_expressed_all_samples_1530295303.csv  | slurm-5318800.err-kp199 |
| LOEN   | 8014396  | 5746811  | 2267585  | 0.717  | 0.283  | miRNAs_expressed_all_samples_1530295315.csv  | slurm-5318801.err-kp025 |
| LPAU   | 6579103  | 2551590  | 4027513  |  0.388 | 0.612  | miRNAs_expressed_all_samples_1530295324.csv  | slurm-5318802.err-kp199 |
| LVIE   | 7280436  | 3736893  | 3543543  | 0.513  | 0.487  | miRNAs_expressed_all_samples_1530295333.csv  | slurm-5318803.err-kp025  |
| LZEP   | 4891192  | 2399312  | 2491880  |  0.491 |  0.509 | miRNAs_expressed_all_samples_1530295343.csv  | slurm-5318804.err-kp199  |

### Filtering novel miRs

I used the following criteria to filter novel miRs from results of md2 run (e.g. `mirna_results_27_06_2018_t_21_54_10_AAUR.html`)
1. No rRNA/tRNA suspects
2. Minimum of 5 reads each on the mature and star strands of the hairpin sequence
3. Has significant randfold p-value

Filtering results:

| Criterion | AAUR  | AVIR  |  HRUB |  LFIG | APUR  | HLIG  | LCAL   | LLEU  | LMAL  | LMAR  | LOEN  | LPAU  | LVIE  | LZEP  | LALB_SOC | LALB_SOL |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Predicted miRs  |  198  | 247   |   255 | 232   | 181   | 233   | 286    | 243   | 356   | 252   | 259   | 336   | 233   | 202   | 309      | 319      |
| flagged as rRNA or tRNA |  0    | 1     |     0 | 0     | 0     | 0     | 0      | 0     | 0     | 0     | 0     | 1     | 0     | 0     | 0   | 0        |
| < 5 reads on either mature or star |  127  | 106   |  163  |  120  | 108   | 131   | 206 | 151   252   | 177   | 156   | 248   | 159   | 150   | 212 | 207 |
| non-significant randfold p |  2  | 9     | 9     |  7    | 7     | 2     |  6     | 3     | 6     | 5     | 6     | 7     | 8     | 4     | 8  | 5   |
| Novel miRs retained after filtering  |   69  | 77    | 83    | 76    | 66  | 100   | 74     | 89    | 98    | 70  | 97  | 81 | 66    | 48 | 89 | 107      |
| Have identical seed match to known miR  | 56  | 53 |  61  |  50   | 51  | 57   |  52  |  53   | 65    | 51    | 55  | 43    | 50  | 37  | 53  | 60  |


## Identify homologs across species
Now run through the pipeline a second time, but with all of the halictid mirs included as known miRs. This allows identifiaction of homologs across species.

### Acquire known miRNAs
Collect the final set of miRs from each halictid species, and combine with the first set.

1. Pull out fasta for final set of miRs from eachs species.

```
mkdir initial_miRs
ll mirna_results_27_06_2018_t_21_54_10_AAUR
cp mirna_results_27_06_2018_t_21_54_10_AAUR/novel_mature_27_06_2018_t_21_54_10_AAUR_score-50_to_na.fa initial_miRs/
```

2. Use blast to pull out the seqs for the filtered set of miRs.

```
ml blast/2.3.0+
 cd ../kapheim-group1/Kocher_miRNA/initial_miRs/
 for i in *.fa; do makeblastdb -in $i -dbtype nucl -parse_seqids ; done
 blastdbcmd -db novel_mature_27_06_2018_t_21_54_10_AAUR_score-50_to_na.fa -entry_batch AAUR_initial_mirs.txt -out AAUR_initial_mirs.fa
grep -c '>' AAUR_initial_mirs.fa
 blastdbcmd -db novel_mature_28_06_2018_t_15_00_19_APUR_score-50_to_na.fa -entry_batch APUR_initial_mirs.txt -out APUR_initial_mirs.fa
grep -c '>' APUR_initial_mirs.fa
blastdbcmd -db novel_mature_27_06_2018_t_22_03_58_AVIR_score-50_to_na.fa -entry_batch AVIR_initial_mirs.txt -out AVIR_initial_mirs.fa
grep -c '>' AVIR_initial_mirs.fa
blastdbcmd -db novel_mature_28_06_2018_t_15_36_34_HLIG_score-50_to_na.fa -entry_batch HLIG_initial_mirs.txt -out HLIG_initial_mirs.fa
grep -c '>' HLIG_initial_mirs.fa
blastdbcmd -db novel_mature_27_06_2018_t_22_06_15_HRUB_score-50_to_na.fa -entry_batch HRUB_initial_mirs.txt -out HRUB_initial_mirs.fa
grep -c '>' HRUB_initial_mirs.fa
blastdbcmd -db novel_mature_28_06_2018_t_17_50_18_LALB_SOC_score-50_to_na.fa -entry_batch LALB_SOC_initial_mirs.txt -out LALB_SOC_initial_mirs.fa
grep -c '>' LALB_SOC_initial_mirs.fa
blastdbcmd -db novel_mature_28_06_2018_t_17_50_19_LALB_SOL_score-50_to_na.fa -entry_batch LALB_SOL_initial_mirs.txt -out LALB_SOL_initial_mirs.fa
grep -c '>' LALB_SOL_initial_mirs.fa
blastdbcmd -db novel_mature_28_06_2018_t_15_37_21_LCAL_score-50_to_na.fa -entry_batch LCAL_initial_mirs.txt -out LCAL_initial_mirs.fa
grep -c '>' LCAL_initial_mirs.fa
blastdbcmd -db novel_mature_27_06_2018_t_22_11_06_LFIG_score-50_to_na.fa -entry_batch LFIG_initial_mirs.txt -out LFIG_initial_mirs.fa
grep -c '>' LFIG_initial_mirs.fa
blastdbcmd -db novel_mature_28_06_2018_t_15_42_17_LLEU_score-50_to_na.fa -entry_batch LLEU_initial_mirs.txt -out LLEU_initial_mirs.fa
grep -c '>' LLEU_initial_mirs.fa
blastdbcmd -db novel_mature_28_06_2018_t_15_48_39_LMAL_score-50_to_na.fa -entry_batch LMAL_initial_mirs.txt -out LMAL_initial_mirs.fa
grep -c '>' LMAL_initial_mirs.fa
blastdbcmd -db novel_mature_28_06_2018_t_15_51_30_LMAR_score-50_to_na.fa -entry_batch LMAR_initial_mirs.txt -out LMAR_initial_mirs.fa
grep -c '>' LMAR_initial_mirs.fa
blastdbcmd -db novel_mature_28_06_2018_t_15_54_52_LOEN_score-50_to_na.fa -entry_batch LOEN_initial_mirs.txt -out LOEN_initial_mirs.fa
grep -c '>' LOEN_initial_mirs.fa
blastdbcmd -db novel_mature_28_06_2018_t_15_59_26_LPAU_score-50_to_na.fa -entry_batch LPAU_initial_mirs.txt -out LPAU_initial_mirs.fa
grep -c '>' LPAU_initial_mirs.fa
blastdbcmd -db novel_mature_28_06_2018_t_16_04_54_LVIE_score-50_to_na.fa -entry_batch LVIE_initial_mirs.txt -out LVIE_initial_mirs.fa
grep -c '>' LVIE_initial_mirs.fa
blastdbcmd -db novel_mature_28_06_2018_t_16_10_17_LZEP_score-50_to_na.fa -entry_batch LZEP_initial_mirs.txt -out LZEP_initial_mirs.fa
grep -c '>' LZEP_initial_mirs.fa
```

3. Now add to known miRNAs

```
ml mirdeep2
for i in *initial_mirs.fa; do remove_white_space_in_id.pl $i > no_white_${i}; done
for i in no_white*; do sed '/^[^>]/ y/tT/uU/' $i > rna_${i}; done
cd ../known_miRs/
cat insect_bees_mature_miRs_nowhite.fa > insect_bees_initialHalictids_mature_miRs.fa
for species in AAUR APUR AVIR HLIG HRUB LALB_SOC LALB_SOL LCAL LFIG LLEU LMAL LMAR LOEN LPAU LVIE LZEP; do cat ../initial_miRs/rna_no_white_${species}_initial_mirs.fa >> insect_bees_initialHalictids_mature_miRs.fa; done
grep -c '>' insect_bees_initialHalictids_mature_miRs.fa
```

### Rerun miR detection

Replace known miRs with the revised set of known miRNAs

Example script: `md2_homologs_AAUR.slurm`

```
#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH -c 6
#SBATCH --account=kapheim-kp
#SBATCH --partition=kapheim-kp
#SBATCH -o slurm-%j.out-%N
#SBATCH -e slurm-%j.err-%N
#
#LOAD MODULES
ml mirdeep2
#
#SET VARS
SPECIES=AAUR
WORKDIR=/uufs/chpc.utah.edu/common/home/kapheim-group1/Kocher_miRNA
#
#
cd $WORKDIR
#
echo "running MirDeep2 miRDeep2.pl a second time with updated known miRs to identify homologs in ${SPECIES}"
#
#run
miRDeep2.pl ./miRNA/processed_reads/${SPECIES}_reads_collapsed.fa ./genomes/${SPECIES}_genome_v2.0.fasta ./mapped_reads/${SPECIES}_reads_collapsed_vs_genome.arf none ./known_miRs/insect_bees_initialHalictids_mature_miRs.fa none -r _${SPECIES}_homologs -z _${SPECIES}_homologs  2>${SPECIES}_homologs_report.log
#
echo "complete"
```


| species | jobreport | error |
| --- | --- | --- |
| AAUR  | slurm-5354491.out-kp292 |
| APUR    | slurm-5354551.out-kp029 |
| AVIR    | slurm-5354572.out-kp163 |
| HLIG    | slurm-5354573.out-kp025 |
| HRUB    | slurm-5354574.out-kp020 |
| LALB_SOC    | slurm-5354576.out-kp018 |
| LALB_SOL    | slurm-5354577.out-kp021 |
| LCAL    | slurm-5354578.out-kp023 |
| LFIG    | slurm-5354579.out-kp029 |
| LLEU    | slurm-5354580.out-kp025 |
| LMAL    | slurm-5354581.out-kp029 |
| LMAR    | slurm-5354582.out-kp163 |
| LOEN    | slurm-5354583.out-kp018 |
| LPAU    | slurm-5354584.out-kp030 |
| LVIE    | slurm-5354585.out-kp001 |
| LZEP    | slurm-5354586.out-kp032 |

Glean mapping stats from *_homologs_report.log
Mapping mature, star, and loop seqs against RFAM:

| Species | Reads processed | Reads with at least one reported alignment | % reads with at least one reported alignment | Reads that failed to align | % reads that failed to align | Reported alignment | subfolder |
| ---- | ---- | --- | :---: | --- | :---: | --- | ---|
| AAUR |  4961 | 21  | 0.42  | 4940  | 99.58  | 58  | mirna_results_03_07_2018_t_20_03_54_AAUR_homologs  |
|  APUR  | 4129  | 17  | 0.41  | 4112  | 99.59  | 2660  |  mirna_results_03_07_2018_t_20_21_48_APUR_homologs |
| AVIR   |  5721 | 25  | 0.44  | 5696  | 99.56  | 139  | mirna_results_03_07_2018_t_20_43_22_AVIR_homologs  |
| HLIG   | 5958  | 24  | 0.40  | 5934  | 99.60  | 107  | mirna_results_03_07_2018_t_20_45_37_HLIG_homologs  |
| HRUB   | 6511  |  26 | 0.40  | 6485  | 99.60  | 91  | mirna_results_03_07_2018_t_20_49_40_HRUB_homologs  |
| LALB_SOC   | 6401  | 30  | 0.47  | 6371  | 99.53  | 347  |  mirna_results_03_07_2018_t_21_10_39_LALB_SOC_homologs |
| LALB_SOL   | 6634  | 31  | 0.47  | 6603  | 99.53  | 348  | mirna_results_03_07_2018_t_21_17_05_LALB_SOL_homologs  |
| LCAL   | 6515  | 19  | 0.29  | 6496  | 99.71  | 58  | mirna_results_03_07_2018_t_21_31_14_LCAL_homologs  |
| LFIG   | 4012  | 15  | 0.37  | 3997  | 99.63  | 77  |  mirna_results_03_07_2018_t_21_37_16_LFIG_homologs |
| LLEU   | 5039  | 28  | 0.56  | 5011  | 99.44  | 162  | mirna_results_03_07_2018_t_22_07_11_LLEU_homologs  |
| LMAL  |  9317 | 32  | 0.34  | 9285  | 99.66  | 136  | mirna_results_03_07_2018_t_22_32_45_LMAL_homologs |
| LMAR   | 7176  | 29  | 0.40  | 7147  | 99.60  | 408  | mirna_results_03_07_2018_t_22_39_12_LMAR_homologs  |
|  LOEN | 3568  | 22  | 0.62  | 3546  | 99.38  | 58  |  mirna_results_03_07_2018_t_22_54_00_LOEN_homologs |
| LPAU   | 11606  | 68  | 0.59  | 11538  |99.41   | 352  | mirna_results_03_07_2018_t_22_59_12_LPAU_homologs  |
| LVIE    | 13037  | 50  | 0.38  | 12987  | 99.62  | 1352  | mirna_results_03_07_2018_t_23_06_34_LVIE_homologs  |
| LZEP   | 10916  | 41  | 0.38  | 10875  | 99.62  | 366  | mirna_results_03_07_2018_t_23_07_05_LZEP_homologs  |



### Filtering novel miRs

I used the following criteria to filter novel miRs from results of md2 run (e.g. `result_03_07_2018_t_20_03_54_AAUR_homologs.html`)
1. No rRNA/tRNA suspects
2. Minimum of 5 reads each on the mature and star strands of the hairpin sequence
3. Has significant randfold p-value

Filtering results:
| Criterion                                   | AAUR  | AVIR  |  HRUB |  LFIG | APUR  | HLIG  | LCAL   | LLEU  | LMAL  | LMAR  | LOEN  | LPAU  | LVIE  | LZEP  | LALB_SOC | LALB_SOL |
| ---------------------------------------     | :---: | :---: | :---: | :---: | :---: | :---: | :--- : | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :------: | :------: |
| Predicted miRs                              |  198  | 252   |   269 | 234   | 188   | 237   | 294    | 250   | 367   | 266   | 269   | 359   | 235   | 218   | 322      | 330      |
| flagged as rRNA or tRNA                     |  0    | 1     |     0 | 0     | 0     | 0     | 0      | 0     | 0     | 0     | 0     | 1     | 0     | 0     | 0        | 0        |
| < 5 reads on either mature or star          |  127  | 166   |  177  |  151  | 114   | 134   | 213    | 158   | 262   | 189   | 166   | 270   | 161   | 165   | 225      | 216      |
| non-significant randfold p                  |    2  | 9     | 9     |  9    | 5     | 1     |  5     | 4     | 7     | 8     | 7     | 8     | 9     | 5     | 8        | 7        |
| Novel miRs retained after filtering         |   69  | 76    | 83    | 74    | 69    | 102   | 76     | 88    | 98    | 69    | 96    | 80    | 65    | 48    | 89       | 107      |
| Previously detected miRs from same species  |  13   |  18   |   14  |   17  |   10  |  40   |  11    |  22   | 17    |  12   |  13   | 12    | 6     |  4    | 28       | 40       |
| Have identical seed match to known miR      |   56  | 58    |  69   |  57   | 58    | 62    |  65    |  66   | 81    | 57    | 83    | 68    | 59    | 44    | 61       | 67       |
| No match to a known miR                     |   0   |  0    |  0    |    0  |   1   |  0    |   0    |    0  |  0    |   0   |  0    |   0   |   0   |   0   | 0        | 0        |

The final miRs are saved as `miRs_AAUR.txt`.


### Acquire sequences of final miRNA sets
Collect the final set of miRs from each halictid species, and combine with the first set.

1. Pull out fasta for final set of miRs from eachs species.

```
mkdir final_miRs_fasta
ll mirna_results_03_07_2018_t_20_03_54_AAUR_homologs
cp mirna_results_03_07_2018_t_20_03_54_AAUR_homologs/novel_mature_03_07_2018_t_20_03_54_AAUR_homologs_score-50_to_na.fa final_miRs_fasta/
```

2. Use blast to pull out the seqs for the filtered set of miRs.

```
ml blast/2.3.0+
cd final_miRs_fasta/
for i in *.fa; do makeblastdb -in $i -dbtype nucl -parse_seqids ; done
blastdbcmd -db novel_mature_03_07_2018_t_20_03_54_AAUR_homologs_score-50_to_na.fa -entry_batch final_miRs_AAUR.txt -out AAUR_final_mirs.txt
grep -c '>' AAUR_final_mirs.txt
```

## Redo homolog search
B. Rubin identified a set of scaffolds in each species that are duplicates. He created a list of these for filtering.
Need to redo the filtering from the initial run to remove any miRs on these scaffolds.   Then redo the homolog identification run.

### Additional filtering

The list of problem scaffolds are in `https://www.dropbox.com/sh/3p2ux491ee4lwvq/AADrS1DuqtF7GZ9wi7571K4Na?dl=0`.

Copied each to `/uufs/chpc.utah.edu/common/home/kapheim-group1/Kocher_miRNA/problem_scaffolds`

Count number of problem scaffolds

```
for i in *.txt; do wc -l $i ; done
```

| species | # problem scaffolds |
| --- | --- |
| AAUR | 4 |
| APUR | 0 |
| AVIR | 3 |
| HLIG | 4 |
| HRUB | 1 |
| LALB | 176 |
| LCAL | 3 |
| LFIG | 3 |
| LLEU | 7 |
| LMAL | 44 |
| LMAR | 6  |
| LOEN | 11 |
| LPAU | 58 |
| LVIE | 15 |
| LZEP | 7  |

Search against the initial mir lists

```
for i in *.txt; do sed 's/[[:blank:]]*$//' $i > no_white_${i} ; done
cd ../initial_miRs/
fgrep -f ../problem_scaffolds/no_white_AAUR_filter_these.txt ./no_white_AAUR_initial_mirs.fa
fgrep -f ../problem_scaffolds/no_white_APUR_filter_these.txt ./no_white_APUR_initial_mirs.fa
fgrep -f ../problem_scaffolds/no_white_AVIR_filter_these.txt ./no_white_AVIR_initial_mirs.fa
fgrep -f ../problem_scaffolds/no_white_HLIG_filter_these.txt ./no_white_HLIG_initial_mirs.fa
fgrep -f ../problem_scaffolds/no_white_HRUB_filter_these.txt ./no_white_HRUB_initial_mirs.fa
fgrep -f ../problem_scaffolds/no_white_LALB_filter_these.txt ./no_white_LALB_SOL_initial_mirs.fa
fgrep -f ../problem_scaffolds/no_white_LALB_filter_these.txt ./no_white_LALB_SOC_initial_mirs.fa
fgrep -f ../problem_scaffolds/no_white_LCAL_filter_these.txt ./no_white_LCAL_initial_mirs.fa
fgrep -f ../problem_scaffolds/no_white_LFIG_filter_these.txt ./no_white_LFIG_initial_mirs.fa
fgrep -f ../problem_scaffolds/no_white_LLEU_filter_these.txt ./no_white_LLEU_initial_mirs.fa
fgrep -f ../problem_scaffolds/no_white_LMAL_filter_these.txt ./no_white_LMAL_initial_mirs.fa
fgrep -f ../problem_scaffolds/no_white_LMAR_filter_these.txt ./no_white_LMAR_initial_mirs.fa
fgrep -f ../problem_scaffolds/no_white_LOEN_filter_these.txt ./no_white_LOEN_initial_mirs.fa
fgrep -f ../problem_scaffolds/no_white_LPAU_filter_these.txt ./no_white_LPAU_initial_mirs.fa
fgrep -f ../problem_scaffolds/no_white_LVIE_filter_these.txt ./no_white_LVIE_initial_mirs.fa
fgrep -f ../problem_scaffolds/no_white_LZEP_filter_these.txt ./no_white_LZEP_initial_mirs.fa
```

Manually inspected list of hits, because matched scaffolds with additional numerals (e.g. LPAU_scaf_31 hits LPAU_scaf_3122).


Filtering results:
| Criterion                                   | AAUR  | AVIR  |  HRUB |  LFIG | APUR  | HLIG  | LCAL   | LLEU  | LMAL  | LMAR  | LOEN  | LPAU  | LVIE  | LZEP  | LALB_SOC | LALB_SOL |
| ------------------------------------------- | :---: | :---: | :---: | :---: | :---: | :---: | :--- : | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :------: | :------: |
| Predicted miRs                              |  198  | 247   |   255 | 232   | 181   | 233   | 286    | 243   | 356   | 252   | 259   | 336   | 233   | 202   | 309      | 319      |
| flagged as rRNA or tRNA                     |  0    | 1     |     0 | 0     | 0     | 0     | 0      | 0     | 0     | 0     | 0     | 1     | 0     | 0     | 0        | 0        |
| < 5 reads on either mature or star          |  127  | 106   |  163  |  120  | 108   | 131   | 206    | 151   | 252   | 177   | 156   | 248   | 159   | 150   | 212      | 207      |
| non-significant randfold p                  |    2  | 9     | 9     |  7    | 7     | 2     |  6     | 3     | 6     | 5     | 6     | 7     | 8     | 4     | 8        | 5        |
| Novel miRs retained after initial filtering | 69    | 77    | 83    | 76    | 66    | 100   | 74     | 89    | 98    | 70    | 97    | 81    | 66    | 48    | 89       | 107      |
| Located on a problem scaffold               |    0  |  0    |   0   |   0   | 0     |   0   |  0     | 0     |  11   |   0   |  5    |  1    |   0   |  0    |  0       | 0        |
| miRs retain after remove problem scaffolds  |   69  |  77   | 83    | 76    | 66    | 100   | 74     |  89   |   87  | 70    |   92  |   80  |   66  |  48   |  89      |  107     |
| Have identical seed match to known miR      |   56  | 53    |  61   |  50   | 51    | 57    |  52    |  53   | 55    | 51    | 51    | 42    | 50    | 37    | 53       | 60       |

miRs that need to be filtered out:
LMAL
>LMAL_chr_57_34103
>LMAL_chr_57_33552
>LMAL_scaf_36_33133
>LMAL_scaf_36_33127
>LMAL_scaf_4229_39554
>LMAL_chr_57_34102
>LMAL_chr_57_33553
>LMAL_chr_57_34113
>LMAL_chr_57_33542
>LMAL_scaf_4233_39610
>LMAL_scaf_4233_39633
LOEN
>LOEN_chr_4661_14324
>LOEN_chr_4661_14489
>LOEN_chr_4661_14488
>LOEN_chr_4661_14505
>LOEN_chr_4659_14296
LPAU
>LPAU_chr_133_30413

Saved new list of miRs as 'LMAL_initial_mirs_noprobs.txt'

### Redo homolog identification for all species with updated miR list

1. Edit known miRs file
Do this manually with text editor, since there are so few.

```
cd ../known_miRs/
grep -c '>' insect_bees_initialHalictids_mature_miRs.fa
emacs insect_bees_initialHalictids_mature_miRs.fa
grep -c '>' insect_bees_initialHalictids_mature_miRs_noprobs.fa
```
saved as `/uufs/chpc.utah.edu/common/home/kapheim-group1/Kocher_miRNA/known_miRs/insect_bees_initialHalictids_mature_miRs_noprobs.fa`

### Rerun miR detection

Replace known miRs with the revised set of known miRNAs

Example script: `md2_homologs_r2_LFIG.slurm`

```
#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH -c 6
#SBATCH --account=kapheim
#SBATCH --partition=kingspeak
#SBATCH -o slurm-%j.out-%N
#SBATCH -e slurm-%j.err-%N
#
#LOAD MODULES
ml mirdeep2
#
#SET VARS
SPECIES=LFIG
WORKDIR=/uufs/chpc.utah.edu/common/home/kapheim-group1/Kocher_miRNA
#
#
cd $WORKDIR
#
echo "running MirDeep2 miRDeep2.pl a third time with updated known miRs to identify homologs in ${SPECIES}"
#
#run
miRDeep2.pl ./miRNA/processed_reads/${SPECIES}_reads_collapsed.fa ./genomes/${SPECIES}_genome_v2.0.fasta ./mapped_reads/${SPECIES}_reads_collapsed_vs_genome.arf none ./known_miRs/insect_bees_initialHalictids_mature_miRs_noprobs.fa none -r _${SPECIES}_homologs_r2 -z _${SPECIES}_homologs_r2  2>${SPECIES}_homologs_r2_report.log
#
echo "complete"
```


| species | jobreport |
| AAUR | slurm-5508392.out-kp292 |
| APUR  | slurm-5508486.out-kp016 |
| AVIR   | slurm-5508487.out-kp016  |
| HLIG   | slurm-5508498.out-kp026 |
| HRUB   | slurm-5508499.out-kp016 |
| LALB_SOC   | slurm-5508500.out-kp026  |
| LALB_SOL   | slurm-5508501.out-kp111  |
| LCAL   |  slurm-5508502.out-kp009 |
| LFIG   |  slurm-5508503.out-kp016 |
| LLEU   | slurm-5508504.out-kp026  |
| LMAL   |  slurm-5508505.out-kp016 |
| LMAR   |  slurm-5508506.out-kp111 | error - don't understand the problem, reran - slurm-5515740.out-kp002
| LOEN   | slurm-5508507.out-kp009 |
|  LPAU  | slurm-5508508.out-kp111 |
| LVIE   | slurm-5508510.out-kp011  |
| LZEP   | slurm-5508511.out-kp031 |
| MGEN   | slurm-5625897.out-kp292 |
| NMEL   | slurm-5626015.out-kp028 |


Filtering results:
| Criterion                                   | AAUR  | AVIR  |  HRUB |  LFIG | APUR  | HLIG  | LCAL   | LLEU  | LMAL  | LMAR  | LOEN  | LPAU  | LVIE  | LZEP  | LALB_SOC | LALB_SOL |  MGEN |  NMEL |
| ---------------------------------------     | :---: | :---: | :---: | :---: | :---: | :---: | :--- : | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :------: | :------: | :---: | :---: |
| Predicted miRs                              |  203  | 253   |   268 | 236   | 186   | 238   | 289    | 248   | 366   | 260   | 266   | 362   | 238   | 220   | 320      | 330      | 335   | 291 |
| flagged as rRNA or tRNA                     |  0    | 2     |     0 | 0     | 0     | 0     | 0      | 0     | 0     | 0     | 0     | 1     | 0     | 0     | 0        | 0        | 1     | 1 |
| < 5 reads on either mature or star          |  132  | 166   |  176  |  153  | 111   | 135   | 208    | 156   | 265   | 183   | 163   | 273   | 164   | 167   | 223      | 216      | 220   | 175 |
| non-significant randfold p                  |    2  | 10    | 8     |  9    | 6     | 1     |  4     | 5     | 7     | 8     | 4     | 9     | 8     | 5     | 10       | 7        | 8     | 17 |
| Novel miRs retained after initial filtering |   69  | 75    | 84    | 74    | 69    | 102   | 77     | 87    | 94    | 69    | 99    | 79    | 66    | 48    | 87       | 107      | 105   | 97 |
| Located on a problem scaffold               | 0     |   0   |   0   |  0    |  0    |   0   |   0    |  0    |  11   | 0     | 5     |  1    |  0    | 0     |  0       |   0      | na    | na |
| miRs retain after remove problem scaffolds  | 69    |    75 |    84 |   74  |  69   |  102  | 77     |   87  |  83   |  69   | 94    |  78   | 66    | 48    |  87      |    107   | 105   | 97 |
| Previously detected miRs from same species  |  13   |  18   |   15  |   17  |   10  |  40   |  11    |  22   | 17    |  12   |  12   | 12    | 6     |  4    | 29       | 40       | 37    | 23 |
| Have identical seed match to known miR      |   56  | 57    |  69   |  57   | 57    | 62    |  66    |  65   | 66    | 57    | 82    | 66    | 60    | 44    | 58       | 67       | 66    | 73 |
| No match to a known miR                     |   0   |  0    |  0    |    0  |   2   |  0    |   0    |    0  |  0    |   0   |  0    |   0   |   0   |   0   | 0        | 0        | 2     | 1 |


## Characterize genome location of miRs

1. Pull out bed files for pre-miRs

```
mkdir final_miRs_bed
ll mirna_results_23_07_2018_t_17_04_23_AAUR_homologs_r2
cp mirna_results_23_07_2018_t_17_04_23_AAUR_homologs_r2/novel_pres_23_07_2018_t_17_04_23_AAUR_homologs_r2_score-50_to_na.bed final_miRs_bed/AAUR_pre.bed
```

2. Prepare bed files
The bed files have headers. Need to remove them.

```
for i in *; do head -n 5 $i; done
for i in *; do sed '1,4d' $i > clean_${i}; done
for i in *; do dos2unix $i; done
for i in $(ls clean*); do cat $i | awk '{if($3>$2) print $0}' > fix_${i}; done
```

3. Get GFF files

Downloaded gff files from https://www.dropbox.com/sh/3p2ux491ee4lwvq/AADrS1DuqtF7GZ9wi7571K4Na?dl=0.

stored in `/uufs/chpc.utah.edu/common/home/kapheim-group1/Kocher_miRNA/gffs`


4. Use bedtools to find intersectin with gene models

```
mkdir intersect
ml bedtools
bedtools intersect -wao -loj -a final_miRs_bed/fix_clean_AAUR_pre.bed -b gffs/AAUR_OGS_v2.1_longest_isoform.gff3.gz > intersect/AAUR_genes_premiRs.txt
bedtools intersect -wao -loj -a final_miRs_bed/fix_clean_APUR_pre.bed -b gffs/APUR_OGS_v2.1_longest_isoform.gff3.gz > intersect/APUR_genes_premiRs.txt
bedtools intersect -wao -loj -a final_miRs_bed/fix_clean_AVIR_pre.bed -b gffs/AVIR_OGS_v2.1_longest_isoform.gff3.gz > intersect/AVIR_genes_premiRs.txt
bedtools intersect -wao -loj -a final_miRs_bed/fix_clean_HLIG_pre.bed -b gffs/HLIG_OGS_v2.1_longest_isoform.gff3.gz > intersect/HLIG_genes_premiRs.txt
bedtools intersect -wao -loj -a final_miRs_bed/fix_clean_HRUB_pre.bed -b gffs/HRUB_OGS_v2.1_longest_isoform.gff3.gz > intersect/HRUB_genes_premiRs.txt
bedtools intersect -wao -loj -a final_miRs_bed/fix_clean_LALB_SOC_pre.bed -b gffs/LALB_OGS_v2.1_longest_isoform.gff3.gz > intersect/LALB_SOC_genes_premiRs.txt
bedtools intersect -wao -loj -a final_miRs_bed/fix_clean_LALB_SOL_pre.bed -b gffs/LALB_OGS_v2.1_longest_isoform.gff3.gz > intersect/LALB_SOL_genes_premiRs.txt
bedtools intersect -wao -loj -a final_miRs_bed/fix_clean_LCAL_pre.bed -b gffs/LCAL_OGS_v2.1_longest_isoform.gff3.gz > intersect/LCAL_genes_premiRs.txt
bedtools intersect -wao -loj -a final_miRs_bed/fix_clean_LFIG_pre.bed -b gffs/LFIG_OGS_v2.1_longest_isoform.gff3.gz > intersect/LFIG_genes_premiRs.txt
bedtools intersect -wao -loj -a final_miRs_bed/fix_clean_LLEU_pre.bed -b gffs/LLEU_OGS_v2.1_longest_isoform.gff3.gz > intersect/LLEU_genes_premiRs.txt
bedtools intersect -wao -loj -a final_miRs_bed/fix_clean_LMAL_pre.bed -b gffs/LMAL_OGS_v2.1_longest_isoform.gff3.gz > intersect/LMAL_genes_premiRs.txt
bedtools intersect -wao -loj -a final_miRs_bed/fix_clean_LMAR_pre.bed -b gffs/LMAR_OGS_v2.1_longest_isoform.gff3.gz > intersect/LMAR_genes_premiRs.txt
bedtools intersect -wao -loj -a final_miRs_bed/fix_clean_LOEN_pre.bed -b gffs/LOEN_OGS_v2.1_longest_isoform.gff3.gz > intersect/LOEN_genes_premiRs.txt
bedtools intersect -wao -loj -a final_miRs_bed/fix_clean_LPAU_pre.bed -b gffs/LPAU_OGS_v2.1_longest_isoform.gff3.gz > intersect/LPAU_genes_premiRs.txt
bedtools intersect -wao -loj -a final_miRs_bed/fix_clean_LVIE_pre.bed -b gffs/LVIE_OGS_v2.1_longest_isoform.gff3.gz > intersect/LVIE_genes_premiRs.txt
bedtools intersect -wao -loj -a final_miRs_bed/fix_clean_LZEP_pre.bed -b gffs/LZEP_OGS_v2.1_longest_isoform.gff3.gz > intersect/LZEP_genes_premiRs.txt
```

5. Filter out excluded miRs

Filtered based on the final set of miRs (species_homologs_r2.txt).

```
awk 'BEGIN{i=0}
FNR==NR {a[i++]=$1; next}
{for(j=0;j<i;j++)
if(index($0,a[j]))
{print $0; break}
}'  ../final_miRs_textfiles/AAUR_homologs_r2.txt AAUR_genes_premiRs.txt  > AAUR_genes_premiRs_filtered.txt
for i in *_filtered.txt; do wc -l $i; done
```

Filtered output is located in files name 'species_genes_premiRs_filtered.txt'.


The table below reports the number of miRs that fall within each type of genome feature.
Numbers not in parentheses represent features on the same strand as the pre-miR.
Numbers in parentheses indicate strand mismatch.

| species |total miRs | intergenic | exon | intron | 5'UTR | 3'UTR |
| ------- | :-------: | :--------: | :---: | :---: | :---: | :---: |
| AAUR    |   69      |   41       |  6       |   18 (4)  | 0 | 0 |
| APUR    |   69      |   41       |  1 (1)   |    21 (5) | 0 | 0 |
| AVIR    |  75       |   54       |  1       |  16 (4)   | 0 | 0 |
| HLIG    | 102       |   77       |  0       |  19 (6)   | 0 | 0 |
| HRUB    |  84       |   55       |  2       |   23 (4)   | 0 | 0 |
| LALB_SOC |  87      |   58       |  1 (1)   |  24 (3)    | 0 | 0 |
| LALB_SOL |  107     |   68       |  2 (2)   |  30 (5)    | 0 | 0 |
| LCAL    |   77      |   49       |  2 (3)   |  22 (1)    | 0 | 0 |
| LFIG    |   74      |   53       |  0 (1)   |  18 (2)    | 0 | 0 |
| LLEU    |   87      |   57       |  1       |   25 (4)   | 0 | 0 |
| LMAL    |  83       |   53       |  1 (2)   |   21 (6)   | 0 | 0 |
| LMAR    |   69      |   39       |  2 (1)   |   24 (3)   | 0 | 0 |
| LOEN    |   94      |   57       |  2 (1)   |  28 (6)    | 0 | 0 |
| LPAU    |  78       |   55       |  0 (2)   |  17 (4)    | 0 | 0 |
| LVIE    |  66       |   46       |  1       |  19        | 0 | 0 |
| LZEP    |  48       |   32       |  2       |  12 (1)    | 1 | (1) |

* For LZEP, one premiR overlapped with a genes on both the same and opposite strands, and is thus counted twice.

## Computational target prediction

1. Gather the final miR sequences for each species

Extract final set of miR IDs from each species.

```
awk  'NR > 1 {print $1}' AAUR_homologs_r2.txt > AAUR_homologs_r2_IDs.txt
awk  'NR > 1 {print $1}' APUR_homologs_r2.txt > APUR_homologs_r2_IDs.txt
awk  'NR > 1 {print $1}' HLIG_homologs_r2.txt > HLIG_homologs_r2_IDs.txt
awk  'NR > 1 {print $1}' HRUB_homologs_r2.txt > HRUB_homologs_r2_IDs.txt
awk  'NR > 1 {print $1}' LALB_SOC_homologs_r2.txt > LALB_SOC_homologs_r2_IDs.txt
awk  'NR > 1 {print $1}' LALB_SOL_homologs_r2.txt > LALB_SOL_homologs_r2_IDs.txt
awk  'NR > 1 {print $1}' LCAL_homologs_r2.txt > LCAL_homologs_r2_IDs.txt
awk  'NR > 1 {print $1}' LFIG_homologs_r2.txt > LFIG_homologs_r2_IDs.txt
awk  'NR > 1 {print $1}' LLEU_homologs_r2.txt > LLEU_homologs_r2_IDs.txt
awk  'NR > 1 {print $1}' LMAL_homologs_r2.txt > LMAL_homologs_r2_IDs.txt
awk  'NR > 1 {print $1}' LMAR_homologs_r2.txt > LMAR_homologs_r2_IDs.txt
awk  'NR > 1 {print $1}' LOEN_homologs_r2.txt > LOEN_homologs_r2_IDs.txt
awk  'NR > 1 {print $1}' LPAU_homologs_r2.txt > LPAU_homologs_r2_IDs.txt
awk  'NR > 1 {print $1}' LVIE_homologs_r2.txt > LVIE_homologs_r2_IDs.txt
awk  'NR > 1 {print $1}' LZEP_homologs_r2.txt > LZEP_homologs_r2_IDs.txt
for i in *IDs.txt;do wc -l $i;done
```

Make blast database from fasta output and filter based on these IDs.

```
ml blast/2.3.0+
cd mirna_results_23_07_2018_t_17_04_23_AAUR_homologs_r2/
makeblastdb -in novel_mature_23_07_2018_t_17_04_23_AAUR_homologs_r2_score-50_to_na.fa -dbtype nucl -parse_seqids
blastdbcmd -db novel_mature_23_07_2018_t_17_04_23_AAUR_homologs_r2_score-50_to_na.fa -entry_batch ../final_miRs_textfiles/AAUR_homologs_r2_IDs.txt -out ../final_miRs_fasta/AAUR_homologs_r2_mature.fa
cd mirna_results_23_07_2018_t_20_23_07_APUR_homologs_r2/
makeblastdb -in novel_mature_23_07_2018_t_20_23_07_APUR_homologs_r2_score-50_to_na.fa -dbtype nucl -parse_seqids
 blastdbcmd -db novel_mature_23_07_2018_t_20_23_07_APUR_homologs_r2_score-50_to_na.fa -entry_batch ../final_miRs_textfiles/APUR_homologs_r2_IDs.txt -out ../final_miRs_fasta/APUR_homologs_r2_mature.fa
 cd ../mirna_results_23_07_2018_t_21_41_12_AVIR_homologs_r2/
 makeblastdb -in novel_mature_23_07_2018_t_21_41_12_AVIR_homologs_r2_score-50_to_na.fa -dbtype nucl -parse_seqids
 blastdbcmd -db novel_mature_23_07_2018_t_21_41_12_AVIR_homologs_r2_score-50_to_na.fa -entry_batch ../final_miRs_textfiles/AVIR_homologs_r2_IDs.txt -out ../final_miRs_fasta/AVIR_homologs_r2_mature.fa
 cd ../mirna_results_23_07_2018_t_22_44_31_HLIG_homologs_r2/
 makeblastdb -in novel_mature_23_07_2018_t_22_44_31_HLIG_homologs_r2_score-50_to_na.fa -dbtype nucl -parse_seqids
blastdbcmd -db novel_mature_23_07_2018_t_22_44_31_HLIG_homologs_r2_score-50_to_na.fa -entry_batch ../final_miRs_textfiles/HLIG_homologs_r2_IDs.txt -out ../final_miRs_fasta/HLIG_homologs_r2_mature.fa
cd ../mirna_results_23_07_2018_t_23_16_43_HRUB_homologs_r2/
makeblastdb -in novel_mature_23_07_2018_t_23_16_43_HRUB_homologs_r2_score-50_to_na.fa -dbtype nucl -parse_seqids
blastdbcmd -db novel_mature_23_07_2018_t_23_16_43_HRUB_homologs_r2_score-50_to_na.fa -entry_batch ../final_miRs_textfiles/HRUB_homologs_r2_IDs.txt -out ../final_miRs_fasta/HRUB_homologs_r2_mature.fa
cd ../mirna_results_24_07_2018_t_00_07_06_LALB_SOC_homologs_r2/
makeblastdb -in novel_mature_24_07_2018_t_00_07_06_LALB_SOC_homologs_r2_score-50_to_na.fa -dbtype nucl -parse_seqids
blastdbcmd -db novel_mature_24_07_2018_t_00_07_06_LALB_SOC_homologs_r2_score-50_to_na.fa -entry_batch ../final_miRs_textfiles/LALB_SOC_homologs_r2_IDs.txt -out ../final_miRs_fasta/LALB_SOC_homologs_r2_mature.fa
cd ../mirna_results_24_07_2018_t_00_39_28_LALB_SOL_homologs_r2/
makeblastdb -in novel_mature_24_07_2018_t_00_39_28_LALB_SOL_homologs_r2_score-50_to_na.fa -dbtype nucl -parse_seqids
blastdbcmd -db novel_mature_24_07_2018_t_00_39_28_LALB_SOL_homologs_r2_score-50_to_na.fa -entry_batch ../final_miRs_textfiles/LALB_SOL_homologs_r2_IDs.txt -out ../final_miRs_fasta/LALB_SOL_homologs_r2_mature.fa
cd ../mirna_results_24_07_2018_t_00_46_28_LCAL_homologs_r2/
makeblastdb -in novel_mature_24_07_2018_t_00_46_28_LCAL_homologs_r2_score-50_to_na.fa -dbtype nucl -parse_seqids
blastdbcmd -db novel_mature_24_07_2018_t_00_46_28_LCAL_homologs_r2_score-50_to_na.fa -entry_batch ../final_miRs_textfiles/LCAL_homologs_r2_IDs.txt -out ../final_miRs_fasta/LCAL_homologs_r2_mature.fa
cd ../mirna_results_24_07_2018_t_01_30_32_LFIG_homologs_r2/
makeblastdb -in novel_mature_24_07_2018_t_01_30_32_LFIG_homologs_r2_score-50_to_na.fa -dbtype nucl -parse_seqids
 blastdbcmd -db novel_mature_24_07_2018_t_01_30_32_LFIG_homologs_r2_score-50_to_na.fa -entry_batch ../final_miRs_textfiles/LFIG_homologs_r2_IDs.txt -out ../final_miRs_fasta/LFIG_homologs_r2_mature.fa
 cd ../mirna_results_24_07_2018_t_01_54_56_LLEU_homologs_r2/
 makeblastdb -in novel_mature_24_07_2018_t_01_54_56_LLEU_homologs_r2_score-50_to_na.fa -dbtype nucl -parse_seqids
 blastdbcmd -db novel_mature_24_07_2018_t_01_54_56_LLEU_homologs_r2_score-50_to_na.fa -entry_batch ../final_miRs_textfiles/LLEU_homologs_r2_IDs.txt -out ../final_miRs_fasta/LLEU_homologs_r2_mature.fa
 cd ../mirna_results_24_07_2018_t_02_27_34_LMAL_homologs_r2/
makeblastdb -in novel_mature_24_07_2018_t_02_27_34_LMAL_homologs_r2_score-50_to_na.fa -dbtype nucl -parse_seqids
blastdbcmd -db novel_mature_24_07_2018_t_02_27_34_LMAL_homologs_r2_score-50_to_na.fa -entry_batch ../final_miRs_textfiles/LMAL_homologs_r2_IDs.txt -out ../final_miRs_fasta/LMAL_homologs_r2_mature.fa
cd ../mirna_results_24_07_2018_t_21_21_35_LMAR_homologs_r2/
makeblastdb -in novel_mature_24_07_2018_t_21_21_35_LMAR_homologs_r2_score-50_to_na.fa -dbtype nucl -parse_seqids
blastdbcmd -db novel_mature_24_07_2018_t_21_21_35_LMAR_homologs_r2_score-50_to_na.fa -entry_batch ../final_miRs_textfiles/LMAR_homologs_r2_IDs.txt -out ../final_miRs_fasta/LMAR_homologs_r2_mature.fa
cd ../mirna_results_24_07_2018_t_02_28_27_LOEN_homologs_r2/
makeblastdb -in novel_mature_24_07_2018_t_02_28_27_LOEN_homologs_r2_score-50_to_na.fa -dbtype nucl -parse_seqids
blastdbcmd -db novel_mature_24_07_2018_t_02_28_27_LOEN_homologs_r2_score-50_to_na.fa -entry_batch ../final_miRs_textfiles/LOEN_homologs_r2_IDs.txt -out ../final_miRs_fasta/LOEN_homologs_r2_mature.fa
cd ../mirna_results_24_07_2018_t_02_32_42_LPAU_homologs_r2/
makeblastdb -in novel_mature_24_07_2018_t_02_32_42_LPAU_homologs_r2_score-50_to_na.fa -dbtype nucl -parse_seqids
blastdbcmd -db novel_mature_24_07_2018_t_02_32_42_LPAU_homologs_r2_score-50_to_na.fa -entry_batch ../final_miRs_textfiles/LPAU_homologs_r2_IDs.txt -out ../final_miRs_fasta/LPAU_homologs_r2_mature.fa
cd ../mirna_results_24_07_2018_t_02_46_21_LVIE_homologs_r2/
makeblastdb -in novel_mature_24_07_2018_t_02_46_21_LVIE_homologs_r2_score-50_to_na.fa -dbtype nucl -parse_seqids
blastdbcmd -db novel_mature_24_07_2018_t_02_46_21_LVIE_homologs_r2_score-50_to_na.fa -entry_batch ../final_miRs_textfiles/LVIE_homologs_r2_IDs.txt -out ../final_miRs_fasta/LVIE_homologs_r2_mature.fa
cd ../mirna_results_24_07_2018_t_03_10_13_LZEP_homologs_r2/
makeblastdb -in novel_mature_24_07_2018_t_03_10_13_LZEP_homologs_r2_score-50_to_na.fa -dbtype nucl -parse_seqids
blastdbcmd -db novel_mature_24_07_2018_t_03_10_13_LZEP_homologs_r2_score-50_to_na.fa -entry_batch ../final_miRs_textfiles/LZEP_homologs_r2_IDs.txt -out ../final_miRs_fasta/LZEP_homologs_r2_mature.fa
cd ../final_miRs_fasta/
for i in *.fa; do grep -c '>' $i ; done
```

Clean and convert to RNA

```
ml perl
for i in AAUR APUR AVIR HLIG HRUB LALB_SOC LALB_SOL LCAL LFIG LLEU LMAL LMAR LOEN LPAU LVIE LZEP ; do remove_white_space_in_id.pl ${i}_homologs_r2_mature.fa  > ${i}_homologs_r2_mature_clean.fa; done
for i in AAUR APUR AVIR HLIG HRUB LALB_SOC LALB_SOL LCAL LFIG LLEU LMAL LMAR LOEN LPAU LVIE LZEP ; do sed '/^[^>]/ y/tT/uU/'   ${i}_homologs_r2_mature_clean.fa  > ${i}_homologs_r2_mature_clean_rna.fa; done
```

2. Extract 500 bp downstream of each genes

     * Chose 500 bp following Ashby et al. 2016 and the average length of 3' UTRs in Dmel is 442 nt

For each species, make a genome file of the assembly for bedtools

```
ml samtools
cd ../genomes/
for i in *_v2.0.fasta; do samtools faidx $i; done
for i in AAUR APUR AVIR HLIG HRUB LALB LCAL LFIG LLEU LMAL LMAR LOEN LPAU LVIE LZEP ; do cut -f 1-2 ${i}_genome_v2.0.fasta.fai  > ${i}_genome_v2.0.bed; done
```

Use bedtools flank to extract the sequences

```
ml bedtools
# check if should pull out gene or mRNA
for i in *.gff3; do grep -c 'mRNA' $i ; grep -c 'gene' $i ; done
# no difference
for i in *.gff3; do grep 'gene' $i > geneOnly_${i}; done
for i in geneOnly*; do wc -l $i; done
#
for i in AAUR APUR AVIR HLIG HRUB LALB LCAL LFIG LLEU LMAL LMAR LOEN LPAU LVIE LZEP ; do bedtools flank -i geneOnly_${i}_OGS_v2.1_longest_isoform.gff3 -g ../genomes/${i}_genome_v2.0.bed -l 0 -r 500 -s > ${i}_3p500.gff; done
```

Use bedtools to get fasta for the flanking downstream region


_Modify gff files_

```
cd gffs/
for i in AAUR APUR AVIR HLIG HRUB LALB LCAL LFIG LLEU LMAL LMAR LOEN LPAU LVIE LZEP;
do
awk 'BEGIN{FS="['\t'=;]"; OFS="\t"} {print $1,$2,$10,$4,$5,$6,$7,$8}' ${i}_3p500.gff > ${i}_3p500_format_1.gff;
awk 'BEGIN{OFS="\t"} {print $1,$2,$10,$4,$5,$6,$7,$8}' ${i}_3p500_format_1.gff > ${i}_3p500_format_2.gff;
done
for i in AAUR APUR AVIR HLIG HRUB LALB LCAL LFIG LLEU LMAL LMAR LOEN LPAU LVIE LZEP;
do
awk '{if (($5 - $4) > 1) print $0}' ${i}_3p500_format_2.gff > ${i}_3p500_format_3.gff;
wc -l ${i}_3p500_format_2.gff;
wc -l  ${i}_3p500_format_3.gff ;
done
```

13856 AAUR_3p500_format_2.gff
13796 AAUR_3p500_format_3.gff
11660 APUR_3p500_format_2.gff
11630 APUR_3p500_format_3.gff
13670 AVIR_3p500_format_2.gff
13589 AVIR_3p500_format_3.gff
11697 HLIG_3p500_format_2.gff
11658 HLIG_3p500_format_3.gff
11994 HRUB_3p500_format_2.gff
11938 HRUB_3p500_format_3.gff
12455 LALB_3p500_format_2.gff
12394 LALB_3p500_format_3.gff
11982 LCAL_3p500_format_2.gff
11937 LCAL_3p500_format_3.gff
12537 LFIG_3p500_format_2.gff
12497 LFIG_3p500_format_3.gff
11402 LLEU_3p500_format_2.gff
11370 LLEU_3p500_format_3.gff
12491 LMAL_3p500_format_2.gff
12446 LMAL_3p500_format_3.gff
12451 LMAR_3p500_format_2.gff
12405 LMAR_3p500_format_3.gff
12172 LOEN_3p500_format_2.gff
12150 LOEN_3p500_format_3.gff
19868 LPAU_3p500_format_2.gff
19665 LPAU_3p500_format_3.gff
12549 LVIE_3p500_format_2.gff
12513 LVIE_3p500_format_3.gff
12610 LZEP_3p500_format_2.gff
12525 LZEP_3p500_format_3.gff

_Pull out flanking sequence_

```
mkdir targets
cd targets/
ml bedtools
for i in AAUR APUR AVIR HLIG HRUB LALB LCAL LFIG LLEU LMAL LMAR LOEN LPAU LVIE LZEP; do bedtools getfasta -fi ../genomes/${i}_genome_v2.0.fasta -bed ../gffs/${i}_3p500_format_3.gff -fo ${i}_3p500.fa -s -name; grep -c '>' ${i}_3p500.fa; done
```

13796
11630
13589
11658
11938
12394
11937
12497
11370
12446
12405
12150
19665
12513
12525

3. Run miranda

```
miranda ../final_miRs_fasta/AAUR_homologs_r2_mature_clean_rna.fa AAUR_3p500.fa -en -20 -sc 140 -strict -out AAUR_miranda_out
grep '>>' AAUR_miranda_out > AAUR_miranda_hits
wc -l AAUR_miranda_hits
cut -f 1 AAUR_miranda_hits | sort | uniq | wc -l
miranda ../final_miRs_fasta/APUR_homologs_r2_mature_clean_rna.fa APUR_3p500.fa -en -20 -sc 140 -strict -out APUR_miranda_out
grep '>>' APUR_miranda_out > APUR_miranda_hits
wc -l APUR_miranda_hits
cut -f 1 APUR_miranda_hits | sort | uniq | wc -l
miranda ../final_miRs_fasta/AVIR_homologs_r2_mature_clean_rna.fa AVIR_3p500.fa -en -20 -sc 140 -strict -out AVIR_miranda_out
grep '>>' AVIR_miranda_out > AVIR_miranda_hits
wc -l AVIR_miranda_hits
cut -f 1 AVIR_miranda_hits | sort | uniq | wc -l
miranda ../final_miRs_fasta/HLIG_homologs_r2_mature_clean_rna.fa HLIG_3p500.fa -en -20 -sc 140 -strict -out HLIG_miranda_out
grep '>>' HLIG_miranda_out > HLIG_miranda_hits
wc -l HLIG_miranda_hits
cut -f 1 HLIG_miranda_hits | sort | uniq | wc -l
miranda ../final_miRs_fasta/HRUB_homologs_r2_mature_clean_rna.fa HRUB_3p500.fa -en -20 -sc 140 -strict -out HRUB_miranda_out
grep '>>' HRUB_miranda_out > HRUB_miranda_hits
wc -l HRUB_miranda_hits
cut -f 1 HRUB_miranda_hits | sort | uniq | wc -l
miranda ../final_miRs_fasta/LALB_SOC_homologs_r2_mature_clean_rna.fa LALB_3p500.fa -en -20 -sc 140 -strict -out LALB_SOC_miranda_out
grep '>>' LALB_SOC_miranda_out > LALB_SOC_miranda_hits
wc -l LALB_SOC_miranda_hits
cut -f 1 LALB_SOC_miranda_hits | sort | uniq | wc -l
miranda ../final_miRs_fasta/LALB_SOL_homologs_r2_mature_clean_rna.fa LALB_3p500.fa -en -20 -sc 140 -strict -out LALB_SOL_miranda_out
grep '>>' LALB_SOL_miranda_out > LALB_SOL_miranda_hits
wc -l LALB_SOL_miranda_hits
cut -f 1 LALB_SOL_miranda_hits | sort | uniq | wc -l
miranda ../final_miRs_fasta/LCAL_homologs_r2_mature_clean_rna.fa LCAL_3p500.fa -en -20 -sc 140 -strict -out LCAL_miranda_out
grep '>>' LCAL_miranda_out > LCAL_miranda_hits
wc -l LCAL_miranda_hits
cut -f 1 LCAL_miranda_hits | sort | uniq | wc -l
miranda ../final_miRs_fasta/LFIG_homologs_r2_mature_clean_rna.fa LFIG_3p500.fa -en -20 -sc 140 -strict -out LFIG_miranda_out
grep '>>' LFIG_miranda_out > LFIG_miranda_hits
wc -l LFIG_miranda_hits
cut -f 1 LFIG_miranda_hits | sort | uniq | wc -l
miranda ../final_miRs_fasta/LLEU_homologs_r2_mature_clean_rna.fa LLEU_3p500.fa -en -20 -sc 140 -strict -out LLEU_miranda_out
grep '>>' LLEU_miranda_out > LLEU_miranda_hits
wc -l LLEU_miranda_hits
cut -f 1 LLEU_miranda_hits | sort | uniq | wc -l
miranda ../final_miRs_fasta/LMAL_homologs_r2_mature_clean_rna.fa LMAL_3p500.fa -en -20 -sc 140 -strict -out LMAL_miranda_out
grep '>>' LMAL_miranda_out > LMAL_miranda_hits
wc -l LMAL_miranda_hits
cut -f 1 LMAL_miranda_hits | sort | uniq | wc -l
grep '>>' LMAR_miranda_out > LMAR_miranda_hits
wc -l LMAR_miranda_hits
cut -f 1 LMAR_miranda_hits | sort | uniq | wc -l
miranda ../final_miRs_fasta/LOEN_homologs_r2_mature_clean_rna.fa LOEN_3p500.fa -en -20 -sc 140 -strict -out LOEN_miranda_out
grep '>>' LOEN_miranda_out > LOEN_miranda_hits
wc -l LOEN_miranda_hits
cut -f 1 LOEN_miranda_hits | sort | uniq | wc -l
miranda ../final_miRs_fasta/LPAU_homologs_r2_mature_clean_rna.fa LPAU_3p500.fa -en -20 -sc 140 -strict -out LPAU_miranda_out
grep '>>' LPAU_miranda_out > LPAU_miranda_hits
wc -l LPAU_miranda_hits
cut -f 1 LPAU_miranda_hits | sort | uniq | wc -l
miranda ../final_miRs_fasta/LVIE_homologs_r2_mature_clean_rna.fa LVIE_3p500.fa -en -20 -sc 140 -strict -out LVIE_miranda_out
grep '>>' LVIE_miranda_out > LVIE_miranda_hits
wc -l LVIE_miranda_hits
cut -f 1 LVIE_miranda_hits | sort | uniq | wc -l
miranda ../final_miRs_fasta/LZEP_homologs_r2_mature_clean_rna.fa LZEP_3p500.fa -en -20 -sc 140 -strict -out LZEP_miranda_out
grep '>>' LZEP_miranda_out > LZEP_miranda_hits
wc -l LZEP_miranda_hits
cut -f 1 LZEP_miranda_hits | sort | uniq | wc -l
```

| species | # hits | # uniq miRs with hits |
| --- | --- | --- |
| AAUR   | 2730  |  69  |
| APUR  | 2398  |  67 |
| AVIR | 3418  |  75 |
| HLIG   |  4466 | 102 |
| HRUB  | 3027  | 84  |
| LALB_SOC   | 3493  |  86 |
| LALB_SOL   |  4157 | 106  |
| LCAL   | 2813  |  77 |
| LFIG   | 2552  |  73 |
| LLEU   | 3208  | 87  |
| LMAL   | 3338  | 82  |
| LMAR   | 2901  | 69  |
| LOEN   |  3206 | 93  |
| LPAU   | 5155  |  78 |
| LVIE   | 2536  |  66 |
| LZEP   | 1947  | 48  |


4. Run RNAHybrid

Setting minimum free energy to -20, using the fly set (3utr_fly) to set xi and theta values - the position and shape parameters of the extreme value distribution from which p-values are calculated.

example script: `rnahybrid_AAUR.slurm`

```
#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH -c 12
#SBATCH --account=kapheim-kp
#SBATCH --partition=kapheim-kp
#SBATCH -o slurm-%j.out-%N
#SBATCH -e slurm-%j.err-%N
#
#LOAD MODULES
ml rnahybrid
#
#SET VARS
SPECIES=AAUR
WORKDIR=/uufs/chpc.utah.edu/common/home/kapheim-group1/Kocher_miRNA/targets
#
#
cd $WORKDIR
#
echo "running rnahybrid for target prediction in ${SPECIES}"
#
#run
RNAhybrid -e -20 -c -s 3utr_fly -t ${SPECIES}_3p500.fa -q ../final_miRs_fasta/${SPECIES}_homologs_r2_mature_clean_rna.fa  > ${SPECIES}_rnahybrid_out
#
echo "complete"
```

Job reports:

| Species | jobreport |
| --- | --- |
| AAUR  | slurm-5617556.out-kp292  |
| APUR   | slurm-5617557.out-kp292  |
| AVIR   |  slurm-5617558.out-kp292 |
| HLIG  | slurm-5617559.out-kp292  |
| HRUB   |  slurm-5617560.out-kp292 |
| LALB_SOC   | slurm-5625690.out-kp010  |
| LALB_SOL   |  slurm-5625691.out-kp017 |
| LCAL   |  slurm-5625692.out-kp028 |
| LFIG   | slurm-5625693.out-kp161  |
| LLEU   | slurm-5625694.out-kp002  |
| LMAL   |  slurm-5625695.out-kp006 |
| LMAR   | slurm-5625696.out-kp013  |
| LOEN   |  slurm-5625697.out-kp014 |
| LPAU   | slurm-5625698.out-kp160  |
| LVIE   | slurm-5625700.out-kp164  |
| LZEP   | slurm-5625701.out-kp158  |


Notes on output from http://manpages.ubuntu.com/manpages/bionic/man1/RNAhybrid.1.html

> For each target/query pair one line of output is generated. Each line is a colon (:) separated list of the following fields: target name, query name, minimum free energy, position in target, alignment line 1, line  2,  line  3, line  4.  If  a  target  or  a  query  is  given  on  command line (ie. no -t or -q respectively), its name in the output will be "command line".


But this does not seem to be the case. The minimum free energy is in column 8 and the p-value is in column 9.
The query miR is in column 6.

Filter based on p-value < 0.001

```
for i in AAUR APUR AVIR HLIG HRUB LALB_SOC LALB_SOL LCAL LFIG LLEU LMAL LMAR LOEN LPAU LVIE LZEP; do awk -F ":" '$9 < 0.001 {print $0}' ${i}_rnahybrid_out > ${i}_rnahybrid_filtered001; done
```

Summarize output

```
for i in *rnahybrid_filtered001; do wc -l $i; cut -f 6 -d : $i | sort | uniq | wc -l; done
```

| species | # hits | # uniq miRs with hits |
| --- | :---: | :---: |
| AAUR   | 2195  | 69   |
| APUR   | 1684  | 68  |
| AVIR   | 2749  | 73  |
| HLIG   | 3001  | 102   |
| HRUB   | 2762  | 84  |
| LALB_SOC   | 2643  | 86  |
| LALB_SOL   | 2549  | 106  |
| LCAL   | 1875  | 77  |
| LFIG   | 2064  | 74  |
| LLEU   | 1570  | 87  |
| LMAL   | 2001  | 83  |
| LMAR   | 1792  | 69  |
| LOEN   | 1919  | 94  |
| LPAU   | 3148  | 78  |
| LVIE   | 1724  | 66  |
| LZEP   | 1218  | 48  |
