# Phylostratigraphy

In order to determine the evolutionary history of the genes in all halictid genomes annotated here, we used the Phylostratigraphy pipeline (https://github.com/AlexGa/Phylostratigraphy). We included the proteomes from a total of 36 species, including the 17 sepcies of halictids annotated in this study. The others span a diversity of animals. All taxa included and their taxonomy are given in `taxonomy.txt`. 

The results of the phylostratigraphic analysis for each species are provided in the files `[species]_animals_final_ps_map.csv.gz`. These files contain the taxonomic level at which each gene originated represented by a single number. These numbers correspond to:

For Halictidae:  
3 is Bilateria  
4 is Protostomia  
5 is Ecdysozoa  
6 is Panarthropoda  
7 is Arthropoda  
8 is Mandibulata  
9 is Pancrustacea  
10 is Hexapoda  
11 is Insecta  
12 is Dicondylia  
13 is Pterygota  
14 is Neoptera (winged insects)  
15 is Holometabola  
16 is Hymenoptera  
17 is Apocrita  
18 is Aculeata  
19 is Apoidea  
21 is Halictidae  
22 is Halictini/Augochlorini/Nomiinae/Rophitinae  
23 is genus-specific  
24 is species-specific

`OG_ages.txt` has ages by orthogroup. Any orthogroup with single-copy genes in at least five species with a majority (>50%) of those genes assigned to the same taxonomic level of origin was assigned that level (`PS` column). The `Highest_expression` column gives the highest RNAseq read count obtained from any tissue from any species in the orthogroup. This is just meant to be a way of validating whether a particular gene is real or not. The `all_PSs` and `genes` columns list the levels and genes inferred for each single-copy gene from that OG.

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

