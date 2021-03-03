# Expression atlas

We created a rough reference of where genes are expressed across species using the RNAseq data collected for gene annotation (see https://github.com/berrubin/HalictidGenomics/tree/master/annotation/transcriptomes). For most species we collected RNAseq data from a single individual for whole abdomens, heads, antennae, and Dufour's glands. We mapped these data to the longest isoforms from the transcriptomes of each species using Salmon v0.9.1. The commands used were:

`salmon index -t [species]_cds.fna -i [species]_trans_index --type quasi -k 31`

`salmon quant -i [species]_trans_index -l [ISR or IU] -p 6 --gcBias -1 [sample]_R1.fastq.gz -2 [sample]_R2.fastq.gz -o [sample]_salmon`

ISR or IU was used depending on whether the library was directional or non-directional, as indicated on the transcriptome page. These commands were run and the results compiled using a python script (`atlas.py`). All compiled data is in `expression_atlas.txt`. Both TPM and the number of reads found for each gene from each tissue are given as well as the orthogologous group to which that gene has been assigned. `Sociality` indicates the behavior of the species.

Normalized TPMs for each tissue are in files ending with `_expression_table_normal.txt`.

### PGLS

I also ran PGLS correlating social behavior with TPM for each tissue. Results for these tests are in the files ending with `_express_min8_pgls_normal.txt`.

### Specificity

We also calculated tissue specificity for each orthogroup using the approach of Yanai et al 2005 (https://doi.org/10.1093/bioinformatics/bti042). These values were calculated in `ExpSpecHalictid.R` (written by Brendan Hunt) from normalized expression values and the results are in `normal_tauTissueSpec_randoUniqueOG_8sp4tissues_4soc_4sol.txt`. The results include specificity values based on data from all species (`tau.8sp4tissues`) as well as from just the social species (`tau.4soc4tissues`) and just the solitary species (`tau.4soc4tissues`). All analyses were based on the metric calculated from all species.