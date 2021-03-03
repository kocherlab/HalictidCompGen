# RELAX

HyPhy RELAX was run using two commands:

`python /Genomics/kocherlab/berubin/local/developing/selection_pipeline/selection_pipeline.py -b /Genomics/kocherlab/berubin/comparative/halictids/halictid_selection -o halictid -r /Genomics/kocherlab/berubin/comparative/halictids/orthology/orthofinder/protein_seqs/OrthoFinder/Results_Jan29_4/Orthogroups/Orthogroups.txt -p 16 -t 12 -a hyphy_relax -d halictids.params --taxa_inclusion /Genomics/kocherlab/berubin/comparative/halictids/halictid_selection/solloss_poly_anc.txt -e /Genomics/kocherlab/berubin/comparative/halictids/halictid_selection/RAxML_bestTree.halictid_fourfold.tree -c HQUA,APUR,LFIG,LVIE,LOEN,LLEU`

This takes all loci that fit the requirements and, using the Jarvis-filtered alignments, runs 

`HYPHYMP CPU=1 /usr/local/hyphy/2.3.11/lib/hyphy/TemplateBatchFiles/SelectionAnalyses/RELAX.bf Universal [working_dir]/og_cds_[OG].afa [working_dir]/og_[OG].tree foreground All > [working_dir]/og_[OG]_relax_unlabeledback.txt` 

and then compiles the results from all loci, performs FDR-correction, and outputs a single file. All files for each locus are also maintained. The only somewhat annoying part of this is that the tips of the tree and the sequence names have to exactly match, so trees have to be pruned correctly. In addition, the "foreground" in that command indicates that the taxa labeled with "{foreground}" in the tree file are those to be tested for relaxation. These are the ones listed with the `-c` parameter in the `selection_pipeline.py` command above. So, for example, the tree file might have this as its content:

"(Dnov:0.051781,(((APUR:0.0310698,AAUR{foreground}:0.0178756):0.00944955,Mgen:0.0255212):0.0432092,(AVIR:0.0640128,(((((LMAL{foreground}:0.0085019,(LALB:0.00213167,LCAL:0.00224266):0.00524523):0.00400913,(LOEN:0.0104539,LPAU{foreground}:0.0076911):0.00383719):0.00461634,((LZEP{foreground}:0.00552056,LVIE:0.00755234):0.021513,(LMAR{foreground}:0.0125433,LFIG:0.0182806):0.00133433):0.00106827):0.0134497,LLEU:0.0269634):0.0072458,(HLIG{foreground}:0.0166524,(HRUB:0.00682309,HQUA:0.00626899):0.00217839):0.0182128):0.0116671):0.0235487):0.0376657,NMEL:0.0727997);"

The equivalent command for extant social taxa was also run:

`python /Genomics/kocherlab/berubin/local/developing/selection_pipeline/selection_pipeline.py -b /Genomics/kocherlab/berubin/comparative/halictids/halictid_selection -o halictid -r /Genomics/kocherlab/berubin/comparative/halictids/orthology/orthofinder/protein_seqs/OrthoFinder/Results_Jan29_4/Orthogroups/Orthogroups.txt -p 16 -t 12 -a hyphy_relax -d halictids.params --taxa_inclusion /Genomics/kocherlab/berubin/comparative/halictids/halictid_selection/solloss_poly_anc.txt -e /Genomics/kocherlab/berubin/comparative/halictids/halictid_selection/RAxML_bestTree.halictid_fourfold.tree -c HLIG,AAUR,LMAR,LZEP,LPAU,LMAL`

The "solloss_poly_anc" taxonomy requires both HLIG and HQUA, both AAUR and APUR, and at least one pair of LMAR and LFIG, LZEP and LVIE, or LPAU and LOEN be present. Otherwise, the only requirement is the 12-taxon minimum.

HyPhy RELAX identifies relaxed selection by comparing dN/dS ratios between the background phylogeny and the lineages of interest. Simplistically, if dN/dS is closer to 1 on the focal lineages, that may indicate relaxed selection. Alternatively, if dN/dS is less than 1 in the background and even closer to 0 on the foreground, that is indicative of intensified selection. Similarly, if dN/dS is greater than 1 on the background and even larger on the foreground, that is also indicative of intensified selection.

We hypothesized that extant solitary taxa derived from lineages that were previously social would show convergent signals of relaxed selection on genes essential for sociality. Therefore, we ran HyPhy RELAX on the 6,904 loci that included at least 12 taxa, including HLIG, HQUA, AAUR, APUR, and at least 1 pair of LMAR & LFIG, LZEP & LVIE, and LPAU & LOEN. All polymorphic and ancestrally solitary taxa were left in the analysis (equivalent to "solloss_poly_anc"). We ran two tests on these loci. First, we examined the signal for relaxed or intensified selection on the ancestrally social tax that subsequently lost sociality (HQUA, APUR, LLEU, LFIG, LVIE, LOEN). In order to establish a baseline for comparison, we compared these results to tests for relaxed or intensified selection in the social taxa (HLIG, AAUR, LMAR, LZEP, LPAU, LMAL).

For tests of extant social taxa (HLIG, AAUR, LMAR, LZEP, LPAU, and LMAL), 305 loci showed signatures of relaxed selection (FDR p < 0.1), 393 showed signatures of intensification of selection. 6,887 were successfully tested.

For tests of extant solitary taxa that were ancestrally social (HQUA, APUR, LLEU, LFIG, LVIE, LOEN), 443 loci showed signatures of relaxed selection (FDR p < 0.1), 330 showed signatures of intensified selection and 6,890 were successfully tested. 

## Frequency of relaxed selection

We used Fisher exact tests to determine whether there is a greater frequency of relaxed selection in solitary taxa. This was implemented using `fisher_exact` in the `scipy.stats` package, e.g.,:

`import scipy.stats as stats  
oddsratio, pvalue = stats.fisher_exact([[443, 305], [6447, 6582]])`

Comparing the numbers of loci with signatures of relaxation based on FDR corrected p-values, there is a significantly greater frequency in solitary taxa (p = 2.42x10<sup>-7</sup>, odds-ratio = 1.48). 

Therefore, significantly more genes are experiencing relaxed selection in solitary losses of social behavior than social taxa.

### GO enrichment

GO enrichments for the genes experiencing relaxed selection are given in Table S4.

### Dates

There are 9 taxonomic levels used in the phylostratigraphy analysis that are then used for looking for correlations between gene ages and numbers of genes experiencing relaxed selection. I also attached ages to each taxonomic category. The papers where I found these dates and the dates found are below. These are crown ages.

Taxon|Prop. RELAX|Mya
-----|-----------|---
Bilateria|0.0633|684
Protostomia|0.03788|632
Neoptera|0.07159|373
Holometabola|0.06154|345
Hymenoptera|0.0679|240
Apocrita|0.11905|192
Aculeata|0.09091|162
Apoidea|0.16667|134
Halictidae|0.2|71

Pearson's correlation between these ages and the proportion of genes experiencing relaxed selection are significant (p = 0.018, r = 0.76).

Bilateria 684 Ma (https://www.nature.com/articles/s41598-017-03791-w)
or ~555 Ma (https://dev.biologists.org/content/129/13/3021)

Protostomia 632 (https://www.nature.com/articles/s41598-017-03791-w)

Neoptera (winged insects excluding Odonata) 373 (https://science.sciencemag.org/content/346/6210/763)

Holometabola 345 (https://science.sciencemag.org/content/346/6210/763)

Hymenoptera 240 (https://science.sciencemag.org/content/346/6210/763)
or maybe 256 from (https://www.sciencedirect.com/science/article/pii/S0960982217303251)

Apocrita 192 (https://www.sciencedirect.com/science/article/pii/S0960982217303251)

Aculeata 162 (https://www.sciencedirect.com/science/article/pii/S0960982217303251)

Apoidea 134 (https://www.sciencedirect.com/science/article/pii/S0960982217303251)

Halictidae 71 (https://www.sciencedirect.com/science/article/pii/S0960982217303251)

