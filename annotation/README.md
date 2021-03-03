# Annotation

## Repeat masking

In order to obtain high quality genome annotations, we first characterized the repetitive elements present in each individual genome using RepeatModeler. We removed redundancy from the resulting set of repetitive elements by clustering those with 80% or greater similarity using `CD-HIT`. We compared the resulting repeat elements to all uniprot proteins and *Drosophila melanogaster* proteins using `BLASTX`. Putative repetitive elements with bitscores of at least 100 or more than 50% similar over 50% of the length of either the protein or repeat sequence were removed from the set of repetitive elements so as to avoid masking protein sequences in the genome. Only repetitive elements longer than 80 bases were retained. This filtering was conducted with `filter_repeats.py`.

The resulting set of repeat elements was combined with the repetitive elements for Arthropods available from RepBase (downloaded from https://www.girinst.org/ on March 8, 2017) and then used to mask each genome as in `run_masking.sbatch`.

## BRAKER gene prediction

We generated gene predictions for each genome using BRAKER v2.1.0. First, all RNAseq reads for each genome were mapped to the repeat masked genomes using hisat2 version 2.0.5 with a maximum intron length of 100000. BRAKER was run on these mapped reads and the repeat masked genome, indicating that the genome was softmasked. The commands used are given in `braker_run.sbatch`. For *Halictus quadricinctus* (HQUA), for which we were unable to obtain RNAseq data, we used the RNAseq data from *Halictus rubicundus* (HRUB) for BRAKER prediction.

## MAKER annotation

MAKER v3.00.0 was run on the repeat-masked genomes. The GFF files of aligned ESTs from PASA were used as EST evidence. All high quality protein predictions from Transdecoder from all species were combined and used as protein evidence for each genome. In addition, OGSv3.2 from *Apis mellifera*, OGSv5.42 from *Lasioglossum albipes*, and all uniprot proteins downloaded on December 2, 2016. 

All gene predictions from BRAKER were fed to MAKER as `pred_gff` and we also used the MAKER implementation of EvidenceModeler. We used both the `always_complete` and `correct_est_fusion` options in MAKER.

Included in the output from MAKER is a set of gene predictions for which protein and transcript evidence did not overlap (`*.non_overlapping_ab_initio.*`). We used InterProScan v5.21-60.0 to identify those sequences with conserved protein domains. We considered anything that InterProScan assigned an InterPro family to as a true real gene and reran MAKER, incorporating these predictions into the annotation. The resulting annotations were used as the basis for the OGS. The only subsequent processing that occurred was the removal of any annotations on scaffolds deemed to be spurious dubplications in the assembly step. The full set of commands used to run MAKER is in `run_maker.sbatch` and all accessory scripts are included here.

Details of the final Official Gene Set v2.1 are given below. These gene sets are relatively complete as measured by BUSCO when comparing with the 4,415 genes expected to be present in all Hymenoptera species based on OrthoDB v9. The average percent of complete BUSCOs present in these Official Gene Sets is 93.6%.

Species|Genes|Isoforms|AA in longest<br>isoforms|Complete<br>BUSCOs|Complete single copy BUSCOs|Complete duplicate BUSCOs|Fragmented BUSCOs|Missing BUSCOs
-------|-----|--------|----------------------|--------------|---------------|------------------|----------|-------|
NMEL|11060|12075|6149908|95.9|95.8|0.1|2.4|1.7
AAUR|12511|13562|6384037|93.7|92.8|0.9|3.6|2.7
APUR|11660|12955|6416770|96.2|95.7|0.5|2.1|1.7
AVIR|13487|14526|6219784|87.5|87|0.5|5.9|6.6
HLIG|11669|12705|6213958|93.8|93.5|0.3|3.1|3.1
HRUB|11987|12802|6353270|94.2|93.9|0.3|3|2.8
HQUA|11969|12822|6163676|92.8|92.4|0.4|4.2|3
LLEU|11386|12471|6377109|94.9|94.5|0.4|3.4|1.7
LMAR|12321|13508|6458729|94.5|94|0.5|3.2|2.3
LFIG|12526|13605|6418041|94.5|93.8|0.7|3.2|2.3
LZEP|12589|14033|6445512|92.6|92|0.6|4|3.4
LVIE|12492|13711|6426298|94.2|93.7|0.5|3.3|2.5
LPAU|14982|16065|7051539|90.9|90.2|0.7|5.9|3.2
LOEN|11827|12793|6364080|94.1|93.5|0.6|3.4|2.5
LMAL|11802|12893|6385877|93.8|93.2|0.6|3.1|3.1
LCAL|11970|13075|6404538|94.2|93.7|0.5|3.2|2.6
LALB|11939|12900|6434887|93.6|92.7|0.9|3.4|3
