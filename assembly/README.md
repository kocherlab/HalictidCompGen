# Assembly

### 10x Genomics Supernova assembly
For each species, we constructed a single 10x Genomics library from a single individual and sequenced the resulting library on half a lane of HiSeqX10. We noticed substantial adapter contamination within the resulting data and aggressively trimmed these adapter sequences using `cutadapt` and the `illumina_adapters.fa` file provided here. `cutadapt` was run twice to thoroughly remove all adapters. 

We obtained more sequencing coverage than absolutely necessary and found that including 210 million reads in the Supernova assembly pipeline yielded the most contiguous assemblies. For each genome, we assembled three random subsets of 210 million reads separately. The code for the adapter trimming, tandom selection of reads, Supernova assembly, and production of haploid assembly sequences is given in `run_supernova.sbatch`.

The reads from *Augochlorella aurata* (AAUR) showed overrepresentation of particular sequences that prevented the successful assembly of this genome. We removed 99% of the reads that included these sequences in order to obtain assemblies for this taxon. `filter_AAUR_reads.py` was used to perform the filtering. We performed subsequent adapter trimming as we did for all other species but only assembled a single version of this genome and did not attempt to merge multiple versions. However, we were able to substantially improve the contiguity of the assembly by running fragscaff. The code for this is given in `run_AAUR_fragscaff.sbatch`.

We found a similar issue with the *Halictus rubicundus* (HRUB). This time, the overrepresented sequences used for filtering were "TATCCTATCCTATCC" and its reverse complement, "GGATAGGATAGGATA", otherwise using the same approach as for AAUR. Again, no merging was performed and `fragscaff` was run in the same way for HRUB to yield the full assembly.

For *Halictus quadricinctus*, we found two slightly different sequences overrepresented and filtered both these and the reverse complements: "TATCCTATCCTATCC", "GGATAGGATAGGATA", "TATCCTAACCTATCC", "GGATAGGTTAGGATA". No merging or other post-processing was performed for this genome.

For *Agapostemon virescens* (AVIR), the sequences "GGTCTGGGTTAGGTC" and its reverse complement "TAACCCAGACCTAAC" were filtered. We also were only able to obtain a usable assembly by using Supernova v1.1.5 on this genome. For all other genomes we used Supernova v1.0.0. The commands used for assembly were the same as for the other taxa. 

### Merging and gap-filling
We then merged the resulting assemblies using `quickmerge` in order to obtain the most contiguous and complete assemblies possible. `quickmerge` only works in a pairwise fashion and the order in which genomes are merged impacts the final result. Thus, we merged the three genome assemblies in every possible order and examined the results for the best assembly as in `run_merge.sh`. This was assessed by genome contiguity and BUSCO content, working to optimize contiguity and the percent of complete BUSCOs while not increasing the numbers of BUSCOs present in multiple copies. 

Finally, we used the `sealer` module that is a part of the `ABySS` assembly package to fille gaps. Note that we first trimmed off the part of the first read in each read pair that contains the 10x Genomics barcode. The commands for this are included in `run_sealer.sbatch`.

### Contaminant removal

As expected, some of the assembled scaffolds in the resulting assemblies were derived from bacteria or other microorganisms. In order to filter as many of these out as possible, we used `BLASTN` to compare the similarity of each scaffold sequence to NCBI's databases of animal genomes and bacterial genomes. These two datasets (all animal genome sequences and all bacterial genome sequences) were downloaded on Feb. 1, 2017. The procedure for running this filtering is given in `find_bacterial_scafs.sh` for APUR but the same procedure was used for all taxa. The basic procedure is first to use `BLASTN` to compare the genome of interest to the database of bacterial genomes. If more than 25% of the length of a given scaffold has a significant hit (e-value < 10) to any bacterial genome or if more than 10% of the length of a given scaffold has a significant hit to any bacterial genome *and* there are more than 10 HSPs on that scaffold, then these scaffolds are taken as putative bacterial sequences. These scaffolds are then blasted against both the bacterial genomes database as well as the animal genomes database. If scaffolds have a closer hit to the animal genome database than to the bacterial genome database, they are not considered to be bacterial. Otherwise, they are marked as likely bacterial in origin. `find_bacterial_scafs.sh` makes use of the scripts `parse_blast.py` and `parse_blast_id.py`. The output of the pipeline is a file named `[species]_taxonomy.txt` which lists the sequence regions of likely bacterial origin in the genome of interest and the bacterial taxon with the closest match to that region. Although only scaffolds that passed the previously described requirements were removed from the assembly, this file also lists regions of scaffolds thare similar to bacterial genomes and the most similar bacterial taxon. These may represent introgressions into the host genomes.

### Initial draft genome quality

The resulting draft genomes with scaffolds of likely bacterial origin removed are relatively contiguous. These results are also present in n50_stats.txt.

Species	| Tot. length | Largest scaffold | N10 length | N10 number | N50 length | N50 number | N90 length | N90 number |
--------|-------------|------------------|------------|------------|------------|------------|------------|------------|
AAUR|	352,077,334|	12,073,071|	8,959,723|	4|	1,247,309|	54|	5,743|	4,504|
APUR | 304,283,628 | 15,436,641 | 10,672,764 | 3 | 3,127,818 | 24 | 146,456 | 134 |
AVIR|	469,330,841|	2,137,112|	901,313|	39|	127,334|	670|	2,342|	25,197|
HLIG|	312,298,520|	17,419,600|	15,635,795|	2|	4,779,453|	15|	10,436|	581|
HRUB|	327,565,881|	4,597,427|	2,491,684|	10|	581,289|	125|	16,000|	1,640|
HQUA|	289,685,969|	7,185,806|	4,240,294|	6|	564,288|	89|	4,021|	6,788|	
LLEU|	283,650,161|	27,371,161|	20,289,172|	2|	12,125,545|	9|	355,488|	36|
LMAR|	306,433,948|	11,039,667|	10,188,095|	3|	6,467,210|	19|	21,236|	215|
LFIG|	318,857,596|	15,247,764|	8,422,116|	3|	2,462,014|	34|	5,292|	1,700|
LZEP|	276,047,489|	23,298,098|	19,308,951|	2|	4,523,099|	13|	5,387|	1,972|
LVIE|	287,167,640|	23,343,689|	16,942,357|	2|	5,119,798|	14|	5,381|	1,200|
LPAU|	378,389,323|	22,645,087|	12,528,215|	3|	1,548,541|	41|	3,907|	5,575|
LOEN|	290,802,534|	23,012,633|	16,235,044|	2|	11,127,264|	10|	967,202|	32|
LMAL|	276,356,024|	24,117,772|	17,745,610|	2|	13,453,691|	9|	975,305|	23|
LCAL|	283,842,111|	15,852,017|	15,566,962|	2|	11,106,183|	11|	7,641|	672|
LALB|	337,527,399|	18,476,510|	15,962,212|	2|	4,813,235|	15|	141,786|	204|

And are also relatively complete as determined by BUSCO comparing to the set of 4,415 highly conserved single-copy genes expected to be present in all Hymenoptera species based on OrthoDB v9. These data are also present in `busco_summary.txt`.

Species|	Complete|	Complete single|	Complete duplicated|	Fragmented|	Missing|
-------|----------------|----------------------|---------------------------|--------------|------------|
AAUR|	95.9|	95.2|	0.7|	2|	2.1|
APUR|	97.7|	97.3|	0.4|	1.1|	1.2|
AVIR|	90.0|	89.6|	0.4|	4.4|	5.6|
HLIG|	96.1|	95.8|	0.3|	1.6|	2.3|
HRUB|	95.7|	95.2|	0.5|	2.2|	2.1|
HQUA|	95|	94.7|	0.3|	2.7|	2.3|
LLEU|	97.5|	97.2|	0.3|	1.1|	1.4|
LMAR|	96.6|	96|	0.6|	1.4|	2|
LFIG|	96.4|	95.8|	0.6|	1.5|	2.1|
LZEP|	96.3|	95.8|	0.5|	1.7|	2|
LVIE|	96.3|	95.4|	0.9|	1.7|	2|
LPAU|	94.9|	93.5|	1.4|	2|	3.1|
LOEN|	96.8|	96.3|	0.5|	1.4|	1.8|
LMAL|	95.9|	95.5|	0.4|	1.6|	2.5|
LCAL|	96.4|	96.1|	0.3|	1.7|	1.9|
LALB|	97.5|	96.5|	1.0|	1.2|	1.3|


### Hi-C scaffolding

These draft genomes were then scaffolded using Hi-C sequencing, improving contiguity as shown below. This information is also in `hic_n50_stats.txt`. These assemblies are official release v2.1.

Species|Tot. length|Largest scaffold|N10 length|N10 number|N50 length|N50 number|N90 length|N90 number|
-------|-----------|----------------|----------|----------|----------|----------|----------|----------|
AAUR|	350,330,594|	22,690,499|	21,065,757|	2|	13,236,083|	11|	5,499|	2,902|
APUR|	304,429,207|	41,380,077|	41,380,077|	1|	26,696,440|	5|	50,819|	14|
AVIR|	470,107,380|	25,019,036|	24,601,988|	2|	17,608,045|	11|	2,344|	22,446|
HLIG|	311,043,946|	26,879,809|	21,141,814|	2|	16,590,429|	8|	9,930|	214|
HRUB|	327,929,407|	29,488,006|	25,357,459|	2|	17,162,222|	8|	17,280|	228|
HQUA|	290,059,369|	21,219,412|	15,646,358|	2|	13,249,058|	10|	4,022|	6,054|
LLEU|	283,111,487|	27,342,619|	21,848,561|	2|	13,877,945|	8|	41,975|	19|
LMAR|	302,926,204|	20,544,363|	15,377,833|	2|	7,692,625|	15|	31,725|	51|
LFIG|	318,834,509|	34,263,369|	34,263,369|	1|	17,020,604|	8|	5,280|	1,284|
LVIE|	286,575,882|	35,722,830|	35,722,830|	1|	24,830,099|	5|	5,333|	984|
LZEP|	275,725,122|	25,729,413|	23,461,232|	2|	14,651,701|	7|	5,338|	1,795|
LPAU|	365,716,076|	41,377,013|	41,377,013|	1|	31,106,527|	5|	3,791|	5,277|
LOEN|	281,109,366|	29,315,781|	29,315,781|	1|	13,884,473|	9|	4,230,112|	22|
LMAL|	267,046,309|	34,372,178|	34,372,178|	1|	28,945,805|	5|	57,000|	13|
LCAL|	283,837,535|	21,743,198|	20,511,341|	2|	14,689,314|	9|	7,580|	531|
LALB|	289,555,301|	27,694,846|	26,518,496|	2|	14,869,788|	7|	4,694|	1,749|
NMEL|	298,732,867|	35,637,317|	35,637,317|	1|	19,732,370|	6|	59,726|	16|

The Hi-C sequencing revealed the presence of large duplicated parts of a number of genomes. These sequences remain in the v2.1 assemblies but the names of scaffolds that appear to be extraneous have been marked with "_altto_" followed by the index of the scaffold that they mirror. These scaffolds are not included in any further analyses.