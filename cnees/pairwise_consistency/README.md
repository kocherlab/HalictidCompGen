Phylogenies with branch lengths inferred for each locus are in `cnees_baseml.txt`.

## Pairwise branch lengths

There are 5 pairs of closely related social/solitary species. This means the odds of the social member of all pairs having a longer tip branch than the solitary member is (1/2)^5 or 0.03125. And the same is true of the solitary tips being longer. There are 49,223 CNEEs for which data is available for tip lengths for all 5 pairs of species and there are no ties in branch lengths for any of the pairs. So we expect there to be 1,538.2 loci with consistently greater branch lengths in social species and another 1,538.2 loci with consistently greater branch lengths in solitary species. What we actually observe is 606 CNEEs with greater branch lengths in all social taxa and 2,173 CNEEs with greater branch lengths in all solitary taxa. 

However, we don't want to work directly on those loci counts. The Solitary loss species (almost) always have longer branches than the most closely related social species using either protein sequences or 4-fold degenerate sites. The one exception to this is in *Halictus*. So if there is a faster rate of substitution in solitary losses overall then that will certainly bias our results. Therefore, I repeated the analysis except this time I divided the terminal CNEE branches by the 4-fold degenerate branch lengths (in "RAxML_bestTree.halictid_fourfold.tree").

After standardizing, there are 1,255 CNEEs faster-evolving in all social lineages and 1,876 CNEEs faster-evolving in all solitary lineages. 

binom.test(1255, 52421, 0.03125) p < 2.2e-16
binom.test(1876, 52421, 0.03125) p = 5.27e-09

Note that the total loci (52,421) includes loci with ties -- this makes little/no difference.

## Pairwise comparisons

`standardized_pairs_counts_min0.txt` has all of the data for these analyses. gives the differences (social - solitary) for each of the five pairs of species. "SocFast" is the number of pairs (out of 5) for which the social species has a longer branch. "SolFast" is the opposite. "Equal" means that the two branches are the same length (typically 0). The median of the standardized branch length differences is given in "MedSoc-Sol". The differences between pair of sequences (social & solitary) are given in the "AAUR_APUR", "HLIG_HQUA", "LMAR_LFIG", "LZEP_LVIE", and "LPAU_LOEN" columns. These come in tuples with 3 values separated by semicolons. The first value is the number of non-gap bases that are identical, the second is the number of non-gap bases that are different, and the third is the number of indels. "ProxOG" is the proximal orthogroup. "Fly" is the corresponding Drosophila gene. "HAv3.1" is the corresponding honey bee gene. "Age" is the phylostratigraphy age.

`soc_faster_cnees_intrprom.txt` are all of the CNEEs faster-evolving in all social lineages located in introns or promoters. `sol_faster_cnees_intrprom.txt` is the same for faster solitary lineages. These are sorted by median difference which may be a good starting point for choosing loci.

The "min0" in the file name indicates that any difference in branch lengths was counted. I did try using a minimum branch length difference of 0.006 as well and the number of loci was pretty minimally changed. Using this, there were 1,232 loci faster-evolving in social lineages and 1,781 in solitary. I chose 0.006 because if I combined all of the differences from all pairs into a single dataset, the quartiles were -0.0072 and 0.00495. Note that, of the 592,526 values, 107,592 were zero.

Four PGLS tests for differences in expression between social and solitary bees were conducted on heads, abdomens, antennae, and dufour's glands. The raw p-values for these are now in the `standardized_pairs_counts_min0.txt` files under "HedP", "AbdP", "AntP", and "DufP", respectively. If any of these p-values are less than 0.05, and the locus is consistently faster in social or solitary taxa and in an intron or promoter, then it is listed in `so[c,l]_faster_cnees_intrprom_sigexpress.txt`. There are 58 loci associated with 46 different genes that are faster in solitary and 43 loci associated with 34 different genes faster in social. I added a column in these files, "S2Expression", that lists the expression level of the gene in S2 cells from FlyBase. There are 30 loci with expression in loci evolving faster in solitary.
