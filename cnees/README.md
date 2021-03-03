### Basic methods for CNEEs

The number of bases aligned by cactus for which there are at least 10 taxa represented in the alignment and with at least half the length of the alignment not as gaps is:  
NMEL_chr_10.final.maf	7044150  
NMEL_chr_11.final.maf	6725681  
NMEL_chr_12.final.maf	8552731  
NMEL_chr_13.final.maf	7824954  
NMEL_chr_14.final.maf	11828257  
NMEL_chr_1.final.maf	6388928  
NMEL_chr_2.final.maf	8565199  
NMEL_chr_3.final.maf	9593073  
NMEL_chr_4.final.maf	8543905  
NMEL_chr_5.final.maf	6639217  
NMEL_chr_6.final.maf	13436594  
NMEL_chr_7.final.maf	7716552  
NMEL_chr_8.final.maf	5166438  
NMEL_chr_9.final.maf	8736928  

That's approximately 117 Mb aligned.

Iterating parameters for phastCons:  45 0.3  

CNEEs  
66834781  
CNEEs min100  
45325096  
CNEEs - OGS  
53811684  
CNEEs - OGS min 100  
34861287  
CNEEs min100 intersect NCARs  
11033872  


### Methodological notes

I removed coding sequences from the CNEEs of all species individually. This, occasionally, split a CNEE in multiple parts. The last segment in that CNEE for that species was kept and the other segments discarded. This did not affect many CNEEs. For example, for HLIG, there are only 806 split CNEEs out of 167258 total CNEEs.

## Social/solitary exclusive CNEEs

There are a few possible hypotheses. First, social species might require changes in sequences that are already regulatory, allowing for "fine tuning" of gene expression. Alternatively, regulatory elements required for social behavior might no longer be needed in solitary losses of social behavior. These would lead to opposite predictions. For social "fine tuning", we expect there to be more changes in social branches whereas for relaxed selection in solitary we expect more changes in these branches. These hypotheses are both pretty well addressed by inferring the presence of CNEEs across all 19 species and then looking for differences in the 5 pairs of social and solitary species.

However, it is possible that the CNEEs hypothetically lost in solitary species would inhibit the detection of CNEEs when including all species. Therefore, I also inferred CNEEs from just the social species and just the solitary loss species separately. Here, we expect that there should be more total CNEEs in social species than solitary species. In particular we might see the presence of more CNEEs around genes that require more regulation in social species. This is perhaps less likely than changing already existing regulatory regions, since the complete loss of a highly conserved region would require very powerful selection.

CNEEs inferred from the 6 social species:  
171,666 total elements  
38,249,131 bp  
37,681,599 bp in elements > 100bp  
8,188,093 bp uniqe to social  
4,817,269 bp in elements that do not overlap solitary elements by more than 25%  

CNEEs inferred from the 6 solitary loss species:  
187,344 total elements  
39,263,837 bp  
38,430,745 bp in elements > 100bp  
9,203,152 bp unique to solitary  
7,088,866 bp in elements that do not overlap social elements by more than 25%  

So it actually appears that there are a few more CNEEs in solitary species than social species. On average, solitary elements are a bit shorter (209bp on average) than social elements (221bp) so they less often overlap. It could also be possible that multiple social elements overlap the same solitary element sometimes so they both get removed when only one solitary element gets removed. Finally, it's possible it may be easier to identify CNEEs in the solitary genomes because the overall rate of substitution in those species is higher, which means conserved elements will appear more conserved relative to background.

## TF enrichment in CNEEs

Overall, CNEEs are enriched for proximity to TFs relative to all genes:

There are 295 OGs with TF orthologs and 8291 OGs with HAv3.1 orthologs of any sort.

`odds_rat, pval = fisher_exact([[579, 7125], [232, 8059]])`


p: 1.8e-49, odds ratio: 2.287

## Files
CNEE_liftover_beds/ includes bed files for each halictid species with the CNEE regions lifted over from the NMEL coordinates used as a reference when identifying CNEEs. 

CNEE_gffs/ is the same information in gff format (which is what my analysis pipeline uses. 

cne_proximity.txt indicates which gene, if any, the CNEE has been associated with. Note that this gene association is in NMEL space. 
