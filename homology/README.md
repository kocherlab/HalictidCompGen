Raw orthogroups from OrthoFinder that are used in all analyses are in `Orthogroups_plus10k.txt`. `halictid_filtered_orthogroups.index` are the orthogroups used in analyses after sequence alignment and filtering. This file does still include some aligned paralogs that were removed for most analyses.

## OrthoFinder choice

OMA, OrthoMCL, and OrthoFinder were all run on the 19 genomes included in the study. Commands to run these and the resulting orthogroups are in the respective directories. These results were compared using the ```assess_orthos.py``` script. 

```exact_orthogroup_match_counts.txt``` lists the number of orthogroups represented by at least 10 species that were identical in composition between the three approaches. In addition, we mapped the genes from each genome to OrthoDBv10 and extracted the Apoidea orthogroup to which every gene mapped. These mappings were done using the web interface and are included in the Trinotate output files. These mappings were used to make a high confidence set of orthogroups for comparison with the other methods. The three methods were also compared to these OrthoDB (```odb```) orthogroups. The large number of orthogroups that are exactly the same between these methods is encouraging. The comparison with OrthoDB suggests that OrthoFinder may be the best approach for these data.

Identical orthogroups between OrthoMCL and OMA: 4628
Identical orthogroups between OrthoMCL and OrthoFinder: 6248
Identical orthogroups between OMA and OrthoFinder: 5032
Identical orthogroups between OrthoMCL and OrthoDBv10: 5090
Identical orthogroups between OMA and OrthoDBv10: 5242
Identical orthogroups between OrthoFinder and OrthoDB: 5641

The table below and ```methods_stats.txt``` is an assessment of the methods as a whole, examining the composition of orthogroups.

Method|	minspec|medsiz|	meansiz|totgen|	perc10|	perc20|	perc30|	perc40|	perc50|	numogs|	ogsw/para|
------|--------|------|--------|------|-------|-------|-------|-------|-------|-------|----------|
MCL|	4|	19.0|	15.725|	203491|	6.0|	11.0|	16.0|	18.0|	19.0|	11258|	3018|
OrthoF|	4|	19.0|	16.035|	213587|	6.0|	13.0|	17.0|	18.0|	19.0|	11796|	3324|
OMA|	4|	18.0|	15.104|	190266|	6.0|	10.0|	15.0|	17.0|	18.0|	12125|	1377|
MCL|	10|	19.0|	17.914|	190821|	15.0|	17.0|	18.0|	19.0|	19.0|	9231|	2678|
OrthoF|	10|	19.0|	17.925|	200997|	15.0|	17.0|	18.0|	19.0|	19.0|	9957|	2973|
OMA|	10|	18.0|	17.253|	174820|	14.0|	16.0|	17.0|	18.0|	18.0|	9846|	1069|
MCL|	15|	19.0|	18.453|	180259|	17.0|	18.0|	18.0|	19.0|	19.0|	8452|	2480|
OrthoF|	15|	19.0|	18.405|	188484|	17.0|	18.0|	18.0|	19.0|	19.0|	9190|	2731|
OMA|	15|	18.0|	17.982|	158510|	16.0|	17.0|	18.0|	18.0|	18.0|	8594|	874|


```Method```: the software used to calculate orthogroups. ```MCL``` is ```orthoMCL``` and ```OrthoF``` is ```OrthoFinder```.  
```minspec```: the information in the rest of the row is based on orthogroups with representatives from at least this number of species.  
```medsiz```: the median number of species in orthogroups.  
```meansiz```: the median number of species in orthogroups.  
```totgen```: the total number of genes represented by orthogroups.  
```perc```: these columns are percentiles of the number of species in orthogroups.  
```numogs```: the number of orthogroups with at least this number of species.  
```ogsw/para```: the number of orthogroups where at least on species is represented by multiple sequences.  

These results show that ```OrthoFinder``` generally produces orthogroups with at least slightly greater numbers of species represented and encompases a larger total number of genes. Although this appears to come at a cost of including additional paralogous sequences in individual orthogroups, there is not a drastically larger number of paralogs in these orthogroups so this should be an acceptable cost.

```assess_methods_counts.txt``` includes details of orthogroups from the three methods as they represent individual species.  

```#genes```: the number of genes in the genome of the current species.  
```#odb```: the number of genes that mapped to the Apoidea node in OrthoDBv10.  
```#mcls```: the number of OrthoMCL orthogroups that include at least one gene from this species. Note that the smallest size of OrthoMCL orthogroups is two.  
```#mcl10```: the number of OrthoMCL orthogroups with at least 10 sequences that include at least one gene from this species.   
```mclsp10```: the number of OrthoMCL orthogroups with at least 10 species that include at least one gene from this species.  
```mclsp10thispara```: the number of OrthoMCL orthogroups with at least 10 species that include at least two genes from this species.  
```mclsp10anypara```: the number of OrthoMCL orthogroups with at least 10 species that include at least two genes from any species.  
```mcl4odb```: the number of OrthoMCL orthogroups with at least 4 species that include at least one gene from this species.  
```mcl10odb```: the number of OrthoMCL orthogroups with at least 10 species that include at least one gene from this species.  

There are equivalent sets of columns for ```OrthoFinder``` results (```orf```) and ```OMA``` results (```oma```). The smallest group size for ```OrthoFinder``` is 1, hence why all genes present in the genome are assigned to orthogroups by this software. ```OMA``` uses a smallest group size of 2.

These results show that ```OrthoFinder``` includes the largest number of OrthoDB mappable genes in large orthogroups and seems less apt to overcluster than ```OrthoMCL``` based on the number of paralogs present. ```OrthoFinder``` also appears to include the largest number of genes present in large (>10 species) orthogroups.

Overall, results across methods are not drastically different. However, the several small improvements in orthogroup quality that appear to be present in ```OrthoFinder``` results led us to use these results moving forward.