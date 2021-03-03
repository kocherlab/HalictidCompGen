# GO term count evolution

The numbers of genes with different functions may differ in social and solitary species. We, therefore, examined the correlation between counts of GO terms and sociality across species. Counts of GO terms were assigned to species as is explained in the `annotation/` section, including Trinotate and OrthoDB information as well as incorporating information from orthologous groups of genes. These GO counts were then correlated with sociality using PGLS tests. 

We used the phylogeny and branch lengths generated from the matrix of all concatenated protein sequences for these tests. GO terms with fewer than 100 total genes annotated to them were not considered for those tests of all taxa and with fewer than 50 genes for those tests of taxa subsets.

We performed four sets of tests (in `go_counts_pgls.r`):

1. Using four categories of behavior (social (soc), polymorphic (poly), ancestrally solitary (solanc), solitary losses of social behavior (solloss)). These results are in `alltaxa_pgls_info.txt`.

2. Including all taxa but classifying all species as either social (soc) or solitary (sol). Polymorphic species are classified as social. Results are in `alltaxa_binary_pgls_info.txt`.

3. Excluding polymorphic species. Results are in `nopoly_anc_pgls_info.txt`.

4. Excluding polymorphic and ancestrally solitary species. Results are in `nopoly_noanc_pgls_info.txt`.

Very few terms pass multiple test correction across any of the tests. A few highlights:

> * In alltaxa_binary, "region of cytosol" is the only term to pass correction. It is more abundance in solitary taxa. "dopamine receptor signaling pathway" is more abundance in solitary taxa (FDR-p = 0.28 but uncorrected p = 0.001). "postsynaptic cytosol" (FDR-p = 0.12) is more abundant in solitary species.

> * In nopoly_anc, three terms pass correction, all of which are more abundant in social taxa: "nuclear localization sequence binding", "deoxyribose phosphate metabolic process", "2'-deoxyribonucleotide metabolic process". "positive regulation of synaptic plasticity" is more abundant in solitary taxa (FDR-p = 0.15), as is "postsynaptic cytosol" (FDR-p = 0.14) and "presynaptic cytosol" (FDR-p = 0.15). "GABA receptor binding" is more abundant in social (FDR-p = 0.16). 

> * Nothing passes correction or is additionally insightful in alltaxa or nopoly_noanc.