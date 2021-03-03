# Gene tree discordance

The species tree for the taxa included in our study is well-established and supported. However, some discordance between this species tree and gene trees is expected. Given that several of our analyses of protein-coding sequences measure evolution on the species tree given, this discordance could cause inaccurate estimates of evolutionary rates. Therefore, we assessed the likelihood of gene-tree species-tree discordance for every locus. This was done using `FastTree` and `consel` following the protocol here: http://www.microbesonline.org/fasttree/. This is also implemented in `selection_pipeline.py`. First gene trees are built from the filtered nucleotide alignments using RAxML. Note that paralogous sequences are first removed because they can't be easily compared to the species tree topology.

`python /Genomics/kocherlab/berubin/local/developing/selection_pipeline/selection_pipeline.py -b /Genomics/kocherlab/berubin/comparative/halictids/halictid_selection -o halictid -r /Genomics/kocherlab/berubin/comparative/halictids/orthology/orthofinder/protein_seqs/OrthoFinder/Results_Jan29_4/Orthogroups/Orthogroups.txt -p 18 -t 4 -a nopara_gene_trees -d halictids.params`

Then discordance with the species tree is measured using:

`python /Genomics/kocherlab/berubin/local/developing/selection_pipeline/selection_pipeline.py -b /Genomics/kocherlab/berubin/comparative/halictids/halictid_selection -o halictid -r /Genomics/kocherlab/berubin/comparative/halictids/orthology/orthofinder/protein_seqs/OrthoFinder/Results_Jan29_4/Orthogroups/Orthogroups.txt -p 18 -t 10 -a check_discordance -d halictids.params -e /Genomics/kocherlab/berubin/comparative/halictids/halictid_selection/RAxML_bestTree.halictid.tree`

`consel_consistency.txt` includes the raw p-values and FDR corrected p-values for each of the 8,990 loci for which gene tree discordance was checked. 545 loci (6.1%) have raw p-values < 0.01 and 277 (3.1%) have corrected p-values < 0.05. 