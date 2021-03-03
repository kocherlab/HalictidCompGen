# Transcriptome assembly

For almost every species included in this study, we obtained RNAseq data in order to assist in genome annotation. In general, we constructed and sequenced three RNAseq libraries from each species. These were built from RNA extracted from abdomens, heads (without antennae), and antennae. In addition, for several species we were able to sequence developmental timecourses included several larval and pupal stages. Sizes given for larval samples are subjective assessments of the age of the larva used. Smaller sizes are younger.

Sample|Barcode|Species|Tissue|Library type|#Reads
------|-------|-------|------|------------|------
BER1-E1_abd|GGTAGC|AAUR|abdomen|directional|49075436
BER3-C3_ant|TAGCTT|AAUR|antennae|directional|31204484
S79A|NA|AAUR|dufour’s|non-directional|82632360
BER1-E1_head|GTAGAG|AAUR|head|directional|55706252
BER1-A1_abd|ACTTGA|APUR|abdomen|directional|101643040
BER1-A1_ant|GATCAG|APUR|antennae|directional|84016268
S75|NA|APUR|antennae|non-directional|103273852
S80A|NA|APUR|dufour’s|non-directional|76298412
BER1-A1_head|CAGATC|APUR|head|directional|81086244
BER7-B11_abd|CGATGT|AVIR|abdomen|directional|84010808
BER7-B11_ant|TTAGGC|AVIR|antennae|directional|89824556
BER7-B11_head|ATCACG|AVIR|head|directional|66555324
BER1-E2_abd|CGATGT|HLIG|abdomen|directional|43224976
BER2-B10_ant|GGCTAC|HLIG|antennae|directional|19360380
BER1-E2_head|ATCACG|HLIG|head|directional|34597000
S87A|NA|LALB|dufour’s|non-directional|69420592
PD_nest12_abd|CAACTA|LCAL|abdomen|directional|51879264
PD_nest12_ant|CCAACA|LCAL|antennae|directional|38182012
S86A|NA|LCAL|dufour’s|non-directional|102780968
PD_nest12_head|ACTGAT|LCAL|head|directional|52930776
PD_nest12_L3|GTCCGC|LCAL|larva (size 3)|directional|39953972
PD_nest12_P|GTGGCC|LCAL|pupa|directional|26275152
Lfig_06_abd|ACAGTG|LFIG|abdomen|directional|84015776
Lfig_06_ant|GCCAAT|LFIG|antennae|directional|101029664
S98TA|NA|LFIG|dufour’s|non-directional|88371240
Lfig_06_head|TGACCA|LFIG|head|directional|84726804
BER2-B8_abd|CAAAAG|LLEU|abdomen|directional|72081244
BER2-B8_ant|CTTGTA|LLEU|antennae|directional|32867680
S81A|NA|LLEU|dufour’s|non-directional|75428312
BER2-B8_head|ATGAGC|LLEU|head|directional|35558624
ESW2-H4_abd|CATTTT|LMAL|abdomen|directional|42081884
S88A|NA|LMAL|dufour’s|non-directional|64685508
ESW2-H4_head|GTGAAA|LMAL|head|directional|41307496
ESW7-H3_L1.1|GTTTCG|LMAL|larva (size 1)|directional|21486688
ESW7-H3_L2|CAGATC|LMAL|larva (size 2)|directional|39169472
ESW7-H3_L3|ACTTGA|LMAL|larva (size 3)|directional|43660248
ESW7-H3_L4|GATCAG|LMAL|larva (size 4)|directional|38205888
ESW7-H3_L5|AGTCAA|LMAL|larva (size 5)|directional|25942192
KMT5L-H10_abd|GGCTAC|LMAR|abdomen|directional|87225056
KMT5L-H10_ant|CTTGTA|LMAR|antennae|directional|83685396
S76A|NA|LMAR|antennae|non-directional|90321872
S83A|NA|LMAR|dufour’s|non-directional|87850772
KMT5L-H10_head|TAGCTT|LMAR|head|directional|96708224
S78A|NA|LOEN|antennae|non-directional|145711020
S89|NA|LOEN|dufour’s|non-directional|107192404
ESW2-B6_abd|TGACCA|LPAU|abdomen|directional|46447104
ESW8-E10_ant|GAGTGG|LPAU|antennae|directional|43223608
S90|NA|LPAU|dufour’s|non-directional|119620408
ESW2-B6_head|TTAGGC|LPAU|head|directional|47437612
ESW8-E10_L1.1|CGTACG|LPAU|larva (size 1)|directional|39484528
ESW8-E10_L3|CCGTCC|LPAU|larva (size 3)|directional|23844444
ESW8-E10_PO|ATGTCA|LPAU|pupa (old)|directional|36280116
ESW8-E10_PY|AGTTCC|LPAU|pupa (young)|directional|32532892
BER5-F6_abd|GCCAAT|LVIE|abdomen|directional|31418112
BER3-E2_ant|ATTCCT|LVIE|antennae|directional|124113272
S84A|NA|LVIE|dufour’s|non-directional|87556076
BER5-F6_head|ACAGTG|LVIE|head|directional|32082344
S77A|NA|LZEP|antennae|non-directional|102098204
S85A|NA|LZEP|dufour’s|non-directional|104209136


These RNAseq libraries were assembled into transcriptomes. All data was combined from individual species and we conducted separate *de novo* and genome-guided assemblies using `Trinity`. Example commands used are in `trinity_GG_run.sbatch` and `trinity_denovo_run.sbatch`. Similar commands were used for all species. Note that we conducted separate assemblies for directional and non-directional libraries. We used previously published sequencing libraries for LALB, NMEL, and HRUB to assemble transcriptomes for these species.

The resulting assemblies (directional *de novo*, non-directional *de novo*, directional genome-guided, non-directional genome-guided) were combined into a single set of high-confidence transcripts using PASA and Transdecoder. The commands used for this are included in `run_pasa.sh`.