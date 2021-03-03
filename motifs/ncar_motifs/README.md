# Motif (stubb scores) in NCARs

I recommend skipping right to the "Genome-wide normalization" section below. The results should be more dependable.

## Not normalized

I looked for correlations between TFBS binding strength and social behavior in the 616 NCARs with est_fast_pos signatures of selection. So, for each of the 223 motifs, I tested for a correlation with social behavior in every NCAR. The results of those tests are in `motif_pgls_nopoly_anc/`. Note that stubb failed to score all of the loci for all motifs so the numbers of NCARs tested are a bit less than 616. Also note that these correlations are based on raw stubb scores. This is different than the other analysis which ranks the stubb scores within each genome. 

Summaries of the numbers of positive correlations are in `sig_ncar_motif_p0.0[1,5]_summary.txt`. `soc_sig` is the number of NCARs with positive correlations with social behavior and `sol_sig` is the number of NCARs with negative correlations with social behavior. `genome-wide` indicates whether each motif was among those 94 found to be associated with social behavior (`soc`) or the 29 associated with solitary behavior (`sol`) in the other, genome-wide analysis. 

It looks like, regardless of whether the motifs were associated with social behavior in the other analysis, there tend to be more positive correlations with solitary species in the 616 loci. This actually makes a lot of sense from the perspective that the signal in the 616 loci is one of more change in social species and the motifs are based on Drosophila so more change should, in general, reduce binding strength.

However, within those constraints, there seem to be differences between the different sets of motifs. Among those motifs that are not associated with a behavior in the previous analysis, the proportion of significant associations (p<0.01) that is positively correlated with social behavior is 195/(195+235)=0.45. Among those motifs associated with social behavior in the previous analysis, it is similar at 177/(207+177)=0.46. Among the loci associated with solitary behavior in the previous analysis, it is 52/(52+78)=0.4. That does seem to be going in the expected direction but there are a lot of different types of selection being considered here so I'm having some trouble understanding this clearly. In those motifs more important for solitary behavior, the est_fast_pos changes in social bees could weaken binding sites? 

## Normalized

If we quantile normalize the stubb scores, however, then we do really start to see something that makes sense. In `sig_ncar_motif_p0.01_normal_summary.txt`, is the summary of results based on these data. Among the motifs previously associated with solitary behavior, about the same number have positive correlations with social (81) as solitary (79) behavior. The same is true for the motifs previously associated with either behavior (259 for social and 245 for solitary). However, among those 94 previously associated with social behavior, there are 274 positive correlations with social behavior and 206 positive correlations with solitary behavior. I am hopeful that this much more straightforward calculation is similar to what we would get with a rank-based test.

## Binding strength

If we make no requirement for a minimum binding strength, there are 130,603 total tests across NCARs and motifs. However, if we require that at least one species have a minimum binding strength of 2.5 (similar to the requirement for the genome-wide test), then a total of only 3,046 tests are done. This suggests that at least some of the correlations may not indicate true binding sites. So results are pretty different if we look at just those loci with a minimum binding strength. These results are in `sig_ncar_motif_p0.01_min2.5_summary.txt` and `sig_ncar_motif_p0.01_normal_min2.5_summary.txt`.

# Genome-wide normalization

The above analyses all just examined the 616 NCARs identified by the MK tests. So the normalization there is probably not ideal. Therefore, I obtained stubb scores for all NCARs across all species. This way we can use the whole genome to standardize scores which is likely better.

The way that the other motif test worked was that loci were ranked in individual genomes and the rank-scores were used for PGLSs. This is to avoid biases associated with individual genomes such as GC-content or mutation rate that might change the motifs at different rates in different species.

Data is all in the google drive at `/KocherLab/ComparativeGenomics/motifs_stubb/NCAR_stubb/`. Raw stubb scores are in `NCAR_compiled_stubb_scores/`. Scores quantile normalized are in `NCAR_compiled_stubb_normal/`. And scores rank-normalized are in `NCAR_compiled_stubb_ranks/`. PGLS results with no minimum value requirement are in `NCAR_stubb_[ranks,normal]_pgls_min0/` and  with minimum value requirements are in `NCAR_stubb_[ranks,normal]_pgls_min[0.95,1.5]/`.

## Genome-wide rank-normalization

### All loci

Here, we test all NCAR loci with at least 10 species with stubb scores. Among the 616 NCARs significantly associated with social behavior (PGLS p < 0.01) and across all 223 motifs, there are 564 positively correlated with social behavior and 669 correlated with solitary behavior. There are 100 motifs that were not associated with either social or solitary behavior previously, and if we just look at those motifs, there are 259 with positive correlations with social behavior and 305 with positive correlations with solitary behavior. Just looking at the 94 motifs previously associated with social behavior, there are 238 positive correlations with sociality and 259 positively associated with solitary behavior. And just looking at the 29 motifs previously associated with solitary behavior, there are 72 positively correlated with social behavior and 109 positively correlated with solitary behavior. 

To determine whether the 616 loci are different than the rest of the genome, I looked at 1,000 random sets of 616 loci. Just looking at the 94 motifs previously associated with social behavior, there are 970 permutations with 238 or more positive correlations with social behavior (p=0.97). Similarly, among the motifs associated with solitary behavior, there are 528 permutations with 109 or more positive correlations with solitary behavior. So these 616 loci may be different, but not in the way "expected" when compared with the genome as a whole. Amongst the loci previously associated with social behavior, there are not more positive correlations with social behavior than what is present in the genome as a whole. We should keep in mind that these loci are detected as faster-evolving in social species so it would be rather remarkable if we did find this pattern.

However, there are some loci that stand out. I also did 1,000-iteration permutation tests for each individual motif. In the 616 loci, there are 7 motifs that have significantly more positive correlations with social behavior than in random samplings (p < 0.05). These are foxo_SOLEXA, retn_SOLEXA, lola_PD_SOLEXA_5, Hmx_SOLEXA, CG7368_SOLEXA_2.5, grh_FlyReg, Hsf_compiled. lola is identified in Karen's 10 genomes paper, though it is a different motif than found here (lola_PQ_SOLEXA). If we increase the p-value cutoff for the permutation test to 0.1, then 9 motifs are significant, including ttk_PA_SOLEXA_5, which was identied in Karen's 10 genomes paper.

There are 4 motifs with more correlations with solitary behavior than amongst random permutations (p < 0.05). These are Ind_Cell, Hr83_SOLEXA, dl_NBT, and gt_SOLEXA. If we increase the p-value cutoff for the permutation test to 0.1, then 10 motifs are significant, including gsb_n_SOLEXA, which was identied in Karen's 10 genomes paper.

Data for this is in `sig_ncar_motif_p0.01_min0_genomewide_rank_summary_sorted.txt`, where:  
`soc_sig`: number of positive correlations with social behavior (p<0.01)  
`sol_sig`: number of positive correlations with solitary behavior (p<0.01)  
`num_tests`: number of tests done on that motif  
`previous`: whether the motif was previously associated with social behavior, solitary behavior, or neither  
`soc_perms`: the number of random sets of NCARs with more than or equal to `soc_sig`.  
`sol_perms`: the number of random sets of NCARs with more than or equal to `sol_sig`.  
`soc-sol`: the number of positive correlations with social behavior minus the number of positive correlations with solitary behavior.  

There are also some interesting things if we just look at the differences between social and solitary correlations for each motif within these loci. Hr39, which I know is of interest from some other studies, *is* significant at a p cutoff of 0.1 and has 6 positive correlations with solitary behavior and 0 correlations with social behavior.

### Only loci with max rank-normalized stubb score > 0.95

In the previous motif analysis, we only considered loci where at least one sequence had a rank stubb score of at least 0.95. This gives us more confidence that some of the motifs are actually present. If it was a correlation amongs a bunch of very low stubb scores, that seems a bit spurious. On the other hand, the previous analysis considered the maximum score across the k5b upstream and 2kb downstream of the TSS whereas this is only in a single 500 base NCAR. So this requirement is far more strict here.

Limiting the loci tested to just those for which there is at least one species with a ranked stubb score of at least 0.95, there is a single motif that has significantly more correlations with social species: vfl_SOLEXA_5. However, br_PL_SOLEXA_5 has a p-value of 0.099. This is the br that is identified in Karen's paper. 

There is also a single motif that has significantly (p < 0.05) more correlations with solitary species: jim_F1_9_SOLEXA_2.5. 

Data for these more strictly filtered results in `sig_ncar_motif_p0.01_min0.95_genomewide_rank_summary_sorted.txt`.

## Genome-wide quantile normalization

We can also just try a quantile normalization and see how different that is. Overall patterns are similar, though motifs of interest are pretty different.

### All loci

There are 492 significant (p<0.01) positive correlations with social behavior across all loci and motifs. There are 511 significant positive correlations with solitary behavior. If we just look at the 94 motifs previously associated with social behavior, there are a total of 208 positive correlations with social behavior and 201 positive correlations with solitary behavior. Looking at just the 29 motifs associated with solitary behavior, there are 68 positive correlations with social behavior and 78 positive correlations with solitary behavior. Among the 100 motifs not previously associated with either social or solitary behavior, there are 221 positive correlations with social behavior and 236 positive correlations with solitary behavior. 

To determine whether the 616 loci are different than the rest of the genome, I looked at 1,000 random sets of 616 loci. Just looking at the 94 motifs previously associated with social behavior, there are 960 permutations with 208 or more positive correlations with social behavior (p=0.96). Similarly, among the motifs associated with solitary behavior, there are 495 permutations with 78 or more positive correlations with solitary behavior. So these 616 loci may be different, but not in the way "expected" when compared with the genome as a whole. Amongst the loci previously associated with social behavior, there are not more positive correlations with social behavior than what is present in the genome as a whole. We should keep in mind that these loci are detected as faster-evolving in social species so it would be rather remarkable if we did find this pattern.

However, there are some loci that stand out. I also did 1,000-iteration permutation tests for each individual motif. In the 616 loci, there are seven motifs that have significantly more positive correlations with social behavior than in random samplings (p < 0.1). These are sr_SOLEXA_5, lola_PD_SOLEXA_5, CrebA_SOLEXA, Ara_SOLEXA, apt_SELEX, Hmx_SOLEXA, and Eip74EF_SOLEXA. They all also have at least 3 more positive correlations with sociality than with solitary behavior. Four of these were previously associated with social behavior and one was previously associated with solitary behavior. At a p cutoff of 0.05, there are only three: sr_SOLEXA_5, **lola_PD_SOLEXA_5**, and apt_SELEX. 

There are 11 loci with significantly more positive correlations with solitary behavior than in random samplings (p < 0.1). Hey_SOLEXA_5, Trl_FlyReg, esg_F3_5_SOLEXA, dl_NBT, Ind_Cell, CG8765_SOLEXA, CG12029_SOLEXA_5, CG11723_SOLEXA, Hr39_SOLEXA, CG4404_SOLEXA, Awh_SOLEXA. Dropping the p cutoff to 0.05 reduces this to just two (CG11723_SOLEXA and CG4404_SOLEXA).

So, overall, I think that we are seeing potentially interesting changes in a few motifs of interest. Perhaps lola is the most exciting one. Data for this is in `sig_ncar_motif_p0.01_min0_genomewide_normal_summary_sorted.txt`.

### Only loci with max normalized stubb score > 1.5

Limiting the loci tested to just those for which there is at least one species with a normalized stubb score of at least 1.5, there 2 motifs that have significantly more correlations with social species: Rel_SANGER and CG12236_PA_SOLEXA_5. There are 4 that have more correlations with solitary species: HLH4C_SOLEXA_5, Lbe_SOLEXA, Awh_SOLEXA, and Adf1_SOLEXA. Note that this restriction drastically limits the number of loci tested and the number of significant loci. So some of these "significantly more correlations" actually just have a single correlation. So this may not be the way to go in this case. Data are in `sig_ncar_motif_p0.01_min1.5_genomewide_normal_summary_sorted.txt` 