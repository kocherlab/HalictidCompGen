## Regional divergence

I inferred branch lengths from 100 AA windows. Then, I looked at the total distance between each tip and DNOV. These numbers are in `region_branch_lengths/`. I also did the same thing on the sites predicted to be in pockets by `fpocket2` and those predicted to be involved in protein-protein interactions using `protindb`. I then looked to see if there were any regions where social species consistently either had more or less change relative to DNOV than ancestrally solitary species (NMEL and AVIR). I assume that these would represent a change in rate associated with the evolution of sociality. Note that this is completely agnostic to solloss and polymorphic lineages. I indicate if social species are faster or slower than ancestrally solitary species in those tables (`shift` column). (I just labeled shifts as "NA" if NMEL or AVIR was missing data.)

Some of the most interesting findings are in apolpp (OG_13369). Most importantly, the AAs in the predicted pocket have faster rates of change. There are also several windows in the beginning of the protein that show this pattern. The pocket sites are all between sites 193 and 932 so this is concordant. The later few thousand AAs of the protein do not show this pattern (though AVIR is missing quite a bit) and the only pattern shown is one window with the opposite trend -- slower evolution in social species. So this could be an interesting finding. I believe, from my notes below, that the beginning part of the protein is more directly associated with function and is what binds to the cell receptor.

There are also a few regions in Hex110 (OG_14622) that are faster in social species, though not in the predicted pocket or interaction sites. They include regions of all three Hemocyanin domains (N,M,C), though they are limited to the first half of the protein which I think is more conserved -- before the big insertion that splits Hemocyanin_C into two regions.

No regions are different in OG_14553 (Pat1).

## MEME analysis

I ran MEME on each of the four genes that have signals of selection on both of the social origin branches (using aBSREL) and have relaxed selection on solitary losses of social behavior (using RELAX). I ran the test three ways. First, including only the 2 origin branches as test lineages. Then also on the solitary loss branches. And finally on the extant social species. Results files are all in `meme_results/`.

### OG_14622 (Hex110)

The Hemocyanin N domain is residues 34-148, Hemocyanin M is 157-447, and Hemocyanin C is split at 456-656 with a big insertion followed by more of the domain at 963-1000. The signal peptide is about the first 22 residues. 23-31 are Danty's conserved hexamerin motif. There are also two glycosylation sites found in the honey bee sequence. The first is AAs 239-242 and is highly conserved. AVIR is the only one with any variation (NHS versus NYS). The second one is AAs 1039-1042 and is highly variable. The sequences for all of these regions are in `hex110_regions.txt`. 

Note that those coordinates are not comparable to those returned by MEME because the alignment was filtered before running MEME. The first glycosylation site is 196-198. The second site is not completely present because it was filtered. The signal peptid is AAs 1-22 and Danty's region is 23-32. None of these overlap with sites identified by MEME.

In the filtered alignment, Hemocyanin_N is 34-148, Hemocyanin_M is 157-444, and Hemocyanin_C is 453-651 and 954-991.

Information about hexamerin110 came from here: https://bmcmolbiol.biomedcentral.com/articles/10.1186/1471-2199-11-23

I matched these up with the 2origins MEME results. For Hex110, I get the following at the sites predicted by MEME. The last column is the AA present at that site in each species in the order: DNOV, NMEL, AAUR, APUR, MGEN, AVIR, HLIG, HQUA, HRUB, LLEU, LMAR, LFIG, LZEP, LVIE, LPAU, LOEN, LMAL, LCAL, LALB.

site|#brs|SuSPect|fpocket|protindb|AAs
----|----|-------|-------|--------|---|
80|2|1|False|False|KQQQQNDDDNDDDDDDDDD
94|1|4|False|False|SSSSSSSSSSSSSSSSSSS
294|1|6|False|False|NNSSCNNNNNNNNSNNNNN
509|0|1|False|False|IVQQQQAAAAAAAAAAATA
532|1|0|False|False|HHIIIYYHYYHHHHHHHHH
562|0|1|False|False|NSNNNKSGGATAVVAAAAA
563|1|4|False|False|KKQHHRRRRRRRRRRKRRR
577|1|1|False|False|QLAAAQEQQEEEEEEEEEE
595|0|1|True|False|FHHHHFFFFFFFFFFFFFF
621|2|2|False|False|EEKKKEKKKKKKKKKKKKK
728|1|1|False|False|--KKKKGGGGYGGGGGSSS
798|1|1|False|False|-----GVVVVVVVVVVVVV
898|0|0|False|False|-----KSSSSSSSTSSGSS
927|1|0|False|False|S-SNNSGGGGGGGGGGGGG
964|1|2|False|False|SSQQQSASSSSSSSSSSSS

Site 621 looks particularly compelling because of the AA pattern. However, the mutation from E to K at this site in particular has the lowest predicted functional consequence. On the other hand, it is immediately adjacent to a predicted pocket AA (site 620). So it could potentially be meaningful.

Site 294 is much more likely to have functional consequences.

The AA pattern at 563 is also pretty interesting and it is four AAs away from being in a pocket.

### OG_14553 (Pat1)

No ProteinDB predictions made so conservation values are below instead.

site|#brs|SuSPect|fpocket|conserve|AAs
----|----|-------|-------|--------|---|
50|0|3|False|0.70|XQAAAKEEEKKKKKKKKKK
60|1|4|False|0.70|XFGEECCCCCCCCCCCCCC
78|1|2|False|0.70|VVIIIVVVVVVVVVVVVVV
82|1|3|False|0.70|RREEERRRRRRRQQRQRRR
95|0|3|False|0.70|GSLXLENNNDSNNNNNNNN
99|1|3|False|0.66|IFLLLIIIIIIIIIIIIII
126|1|4|False|0.70|MNKKKMDDDNNNNNNNNNN
192|1|6|False|0.63|QQFYLYYYYYYYYYYYYYY
211|2|3|False|0.70|YYFFFYFFFFFFFFFFFFF
267|2|3|False|0.30|KKEEEKEEEEEEEEEEEEE
285|0|1|False|0.70|VEKKKVKKKKKKKKKKKKK
410|1|2|False|0.70|KALLMQCCCCCCCCXCCCC
416|1|1|False|0.29|LRLLLLNNNNNNNTXNNNN
510|1|2|False|0.70|LLTTTLLLSLLLLLLLLLL
514|1|1|True|0.69|SSSSSSHHYHHHHHHHHHH
565|1|2|False|0.49|DDDDDELLLLLPLLLLLLL

Site 50 looks interesting given the pattern of amino acid evolution. Mutations from the E of HLIG to K, have a consequence of 2 on the mutation scale. So not a lot, but some.

Site 126 also has a very interesting pattern. Mutations from D to N have consequences of 1 but D to M is 4 and D to K is 5. So this could be a meaningful site.

Site 211 looks promising but mutations from F to Y are of 0 consequence.

Site 267 mutations from E to K are only of consequence 1

Site 285 K to V are of consequence 0 and K to E are of consequence 1.

Site 514 H to S has a consequence of 0.

The predicted pocket is composed of many sites between 432 and 521.

### OG_13369

From what I can tell, the Vitellogenin-N Domain (residues 41-573) is likely the most functionally important part of OG_13369 (apolpp) (https://onlinelibrary.wiley.com/doi/full/10.1002/cbic.201300152). This part of the protein is what binds to the receptor on the cell. However, residues 605-918 (maybe? not clear) and 939-1033 seem to be where lipids get bound and carried. 

There are 52 sites identified on the two origin branches, including a number within the receptor-binding domain. There is no overlap in this domain with the other two versions of the test. However, there is a single site that overlaps in both the origins and solloss tests: 2554. Note that this coordinate is from the filtered alignment. Looking at the alignment at this position, all species have the same amino acid excep DNOV, NMEL, LLEU, and LFIG (AVIR has a gap). LLEU and LFIG have the same amino acid but DNOV and NMEL both have different ones. So that's not super compelling. Nor is the fact that the number of branches inferred to have selection by MEME at that site is zero in the origins test so may or may not be a real signal. (https://github.com/veg/hyphy/issues/666)

The crystal structure of vitellogenin, which is very related, can be seen here: http://www.rcsb.org/structure/1LSH and is from this paper: https://www.sciencedirect.com/science/article/pii/S0969212698000914. There are some specific active sites known in bacteria apolipoproteins (https://jb.asm.org/content/189/12/4456). However, the best blast hit to bacteria is only 32% similar over 15% of the sequence, so we probably can't infer much from this.

See also this summary: https://prosite.expasy.org/PDOC51211

site|#brs|SuSPect|fpocket|protindb|AAs
----|----|-------|-------|--------|---|
45|1|2|False|False|HHKXKHHXHHHHHHHHHHH
99|0|2|False|False|SSNXNQEXEEEEEEEEEEE
281|2|2|False|False|EEEXDEDDDDDDDDDDDDD
330|1|0|False|False|NDAXPKNNNNNNNNNNNNN
668|1|1|False|False|RRRRRRKKKKKKRRKRKKK
703|1|1|False|False|GGGGXGAVVAAAAAAAAAA
934|1|3|False|False|KKKKKNFFFHFFFFFFFFF
935|1|2|False|False|FFFFFFTTTATTTTTTTTT
1030|1|1|False|False|TNTTTTKKKKRRKKRRRKK
1128|0|2|False|False|SSAAAAIIITIIIIIIIII
1158|1|1|False|False|IAKKKGAAAAAAAAAAAAA
1167|1|1|False|False|DDDDDDAAAAAAAAAAAAA
1176|0|2|False|False|SSHHHQLLLFFFFFFFFFF
1297|0|1|False|False|GGAAXGTTTSTTTTTTTTT
1341|1|2|False|False|VETTTVTTTTTTTTTTTTT
1502|1|1|False|False|KKSSSKKKKKKKKKKKKKK
1505|0|2|False|False|LVIVIVLLLLLLLLLLLLL
1551|0|1|False|False|YYGWTGSSSSSSTTSSSSS
1564|0|1|False|False|FCLITKYYYYYYYYYYYYY
1598|1|1|False|False|AXQQXEHHHHNNNNNNNKK
1614|1|3|False|False|GGNNNGGGGGGGGGGGGGG
1679|1|1|False|False|XXXXXRKKKKKKKKKKKKK
1704|0|2|False|False|SSXXXXFFFFFCLLFFFFF
1708|0|2|False|False|SSXXXXKKKKKKKKKKKKK
1792|1|2|False|False|KKSASXKKKKKKKKKKTKK
1816|1|1|False|False|LLEEKXKKKKKKEKKKKKK
1884|1|1|False|False|VXTTTXVVVVVVLLVVVVV
1920|1|1|False|False|VVTTTXVVVVVVVVIVVVV
2125|1|1|False|False|SSSSSSTTTTTTTTTTTTT
2167|1|1|False|False|HYGAAQQQQQQQQQQQQQQ
2267|1|2|False|False|NNQQQYYYYYYYYYYYYYY
2365|1|1|False|False|NEIIIEEEEEEEDDEEEEE
2436|1|1|False|False|SSSSSSTTTSSSSSSSSSS
2479|0|3|False|False|LTTTTFWWWWWWSSWWWWW
2551|1|2|False|False|QDSSTXEEEEEEVVEEEEE
2554|0|0|False|False|TNDDDXSSSASASSSSSSS
2595|1|2|False|False|XAAAAXLLLLLLVVLLLLL
2660|1|2|False|False|YYYYYXLLLLLLLLLLLLL
2746|0|1|False|False|XXVAVXTTTTTTXXTTTTT
2756|1|2|False|False|ILFFFXIIIIIIIIIIIII
2763|0|1|False|False|KKEEEXTTTTTTXXTNTTT
2787|0|2|False|False|LLMTMXSSSSSSSSSSSSS
3024|0|1|False|False|EGYYYXNNNNNNNNNNNNN
3043|0|1|False|False|MLKKKXSSLLLLLLLLLLL
3050|1|1|False|False|SAPPPXSSSAAVAAAAAAA
3156|0|1|False|False|PTQKKXTTTAAAPPAAAPP
3163|0|1|False|False|KEVVVXTTTTTTTTTTTTT
3171|0|1|False|False|MINNNXPPPQPLLLLLLLL
3305|1|0|False|False|TTKKKXTTTTTTTTTTAPP
3352|0|1|False|False|LLTTTXLLLLLLLLLLLLL
3381|1|1|False|False|EEKKXXSSSNSSNNSSSSS

Site 99 is the only one with an interesting pattern that is clearly in a functional domain. Mutational consequences of E to Q/N/S are all 2.

The pattern at site 1176 is interesting. L to Q/H/S is 3, F is 1

Site 1341 T to V/E is 1.

Site 1505 L to V/I is 1.

Site 2267 Y to Q or N has a consequence of 3.

Site 2479 W to T is 3, S is 4, F is 2.

2551 E to S, t, and V is 2, and to D and Q is 1.

2595 L to V is 1 nd A is 2.

2746 T to A/V is 1.

2756 I to F/L is 1.

2763 T to K/E/N is 1.

2787 S to M/T/L is 2.

### OG_11519

MEME found some sites in this gene too. Since it's new to Hymenoptera Phyre couldn't really predict it so there isn't much other info about it.

## "Conservation" in social taxa in regions

Below is just comparing the five social/solitary pairs. So "higher rates of evolution in solitary" means that in each of the five pairs, the solitary member has more change than the social member.

### Apolpp

Regions 1100-1200, 1150-1250, 1700-1800, 1750-1850, and 2950-3050 all have higher rates of evolution in solitary species. None of these regions appear to be particularly important structurally as far as I recall. But they do exist. I think that what we would "want" to see is that the same regions that evolve faster in social species relative to AVIR and NMEL would also evolve slower relative to secondarily solitary species. We don't see that pattern anywhere.

No regions show the opposite pattern (i.e., lower rates in solitary species).

### Hex110

There is also one region in Hex11 that has higher rates (100-200). This does not overlap with regions of interest identified by comparison to AVIR and NMEL. 

Again, no regions show the opposite pattern (i.e., lower rates in solitary species).