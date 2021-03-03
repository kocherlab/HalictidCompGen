# Evolutionary gains of CNEEs

It is difficult to choose the best way to count numbers of CNEEs in different clades. When extracting the alignment from cactus for input into phastcons, NMEL was used as a reference. This means that only sequences that included NMEL were included in the phastcons analysis so this can't really tell us about where CNEEs originated. Therefore, I reran the phastcons pipeline using individual clades. Note that I didn't explore parameters at all for each clade, I just used the same paramenters as I used for the full CNEE set from all data.

The total bases found for each clade are in "CNEE bases". "novel bases" are the bases in CNEEs that do not overlap with CNEEs from the full dataset including all taxa. And "novel bases (NMEL space)" are these bases transferred to NMEL coordinates. ("ingroup_lasioglossum" is all Lasioglossum excluding LLEU.)

clade|CNEE bases|novel bases|novel bases (NMEL space)
-----|----------|-----------|------------------------
augochlorini|49909016|13940196|9554497
halictini|34927318|9065552|6300234
halictus|46407897|10522381|6938386
halilasi|51471998|14599852|10186520
ingroup_lasioglossum|61567919|17294404|11813141
lasioglossum|59601588|16960592|11628794
lmarfigzepvie|60874561|14901127|10114429
lmaltoloen|59342183|12637796|8507156
lmalcalalb|13631120|2100339|1195925
nonmel|36358628|8186466|5366277

I don't see any particularly clear patterns here. Augochlorini and the Halictus/Lasioglossum ("halilasi") node are at the upper end as we would expect if new CNEEs are needed for social evolution. However, both Lasioglossum and ingroup_lasioglossum are higher so that doesn't seem particularly meaningful.

I used NMEL coordinates so that the CNEEs could be compared across clades. I especially wanted to see if the same CNEEs developed independently in different clades. The way that I did this was to first remove CNEEs that overlapped with CNEEs identified in the parent node. This is so that we know we are just looking at those CNEEs novel to the current node, rather than those that developed earlier. Here are the numbers of bases in CNEEs in each clade with the parent CNEEs removed.

clade|no parents|num genes|all CNEEs genes
-----|----------|---------|---------------
augochlorini|5015023|4909|6080
halictini|3563214|3958|5327
halictus|1970892|4004|5738
halilasi|5190759|4859|6456
lasioglossum|3262785|4270|6709

Then I looked for overlap between five pairs of clades; augochlorini with the other four and lasioglossum with halictus.

overlapping clades|num bases|num genes
------------------|---------|---------
augochlorini_noparents_halictini_noparents|413239|1199
augochlorini_noparents_halilasi_noparents|529451|1359
augochlorini_noparents_halictus_noparents|243007|865
augochlorini_noparents_lasioglossum_noparents|412375|1211
lasioglossum_noparents_halictus_noparents|110510|567

And tested these for GO enrichment. For the background set, I used the genes associated with all of the CNEEs identified in one of the taxon sets (Augochlorini for all tests with overlap in Augochlorini and Lasioglossum for overlaps with Halictus). I did permutations to make sure that these enrichments were significant.

Mushroom body development is enriched in the CNEEs that overlap the Augochlorini branch and the branch leading to Halictus/Lasioglossum, i.e., the two origin branches. There are 21/1357 genes associated with this term in the overlap as opposed to 43/4860 in the overall set. The FDR value is 0.06, the fraction of permutations with an FDR < 0.1 is 0.032 and the fraction of permutations with FDR < 0.06 is 0.016, so pretty strong signal. "axon extension", "dendrite self-avoidance", and "neuron projection extension" are also enriched but with weaker signals. 

The overlap between Lasioglossum and Halictus is reasonably large and has some interesting neural stuff. There is a huge amoung of functional overlap between Lasioglossum and Augochlorini including lots of synapse function.

A summary of all go results is in `enrichment_summary.txt`. To be included in this file the FDR value on the original test had to be <= 0.1 and the fraction of permutations with FDR < the test FDR had to be < 0.05. The second column is the FDR value from the GO enrichment test, the third is the fraction of permutations with an FDR < 0.1 and the fourth is the fraction of permutations with FDR < the FDR value from the GO enrichment test on real data.

Original GO enrichment files which list the associated genes are in `go_enrichment/` and counts of significant permutations are in `permutation_counts/`.