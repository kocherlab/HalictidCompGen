#Purpose: calculate gene expression specificity indices for halictid tpm gene expression measures,
#within and across species.
#Ref for tissue specificity eqn: Yanai et al 2005 https://doi.org/10.1093/bioinformatics/bti042
#using some tidyverse stuff, particularly 'purrr' to apply functions

library(tidyverse)

#import halictid "expression atlas" (June 2019) data
atlas <- read.table(file="normal_atlas.txt", sep="\t", header=T)

#TPM tissue headers for analysis
TPMheaders <- c("antTPM", "headTPM", "dufTPM", "abdTPM")

#loop through rows and count number TPM measures with values (no longer needed)
#atlas$TPMcount <- apply(atlas[TPMheaders], 1, function(x) sum(!is.na(x)))

#loop through rows and record max TPM value
atlas$TPMmax <- apply(atlas[TPMheaders], 1, function(x) max(x, na.rm = TRUE))
atlas$TPMmax <- round(atlas$TPMmax, digits=4)

# #Initial test calculation of tau tissue specificity on small subset w/ 3 tissues
# #toy data set
# atlasHead <- head(atlas,15)
# 
# #turn columns into lists to apply functions to each gene for each tissue 
#toy data only has data for 3 tissues
# a <- as.list(atlasHead$antTPM)
# b <- as.list(atlasHead$headTPM)
# c <- as.list(atlasHead$abdTPM)
# d <- as.list(atlasHead$TPMmax)
# 
# #this calculates tau!
# #handles zero value TPMs, returning NaN, 
# #but likely need to filter by total TPM for minimal expression to avoid spurious values
# tau <- pmap(list(a,b,c,d),
#   function(a,b,c,d) (1-(a/d))/2 + (1-(b/d))/2 + (1-(c/d))/2)
# 
# #add tau calculations to data frame
# atlasHead$tau <- unlist(tau)
# atlasHead$tau <- round(atlasHead$tau, digits=4)

# #now that test has worked, need to split out contexts and run on each separately
# #looking to see which species have TPM data for which tissues
# atlas %>%
#   group_by(Species) %>%
#   summarise(antennae = any(!is.na(antTPM)), head = any(!is.na(headTPM)), dufour = any(!is.na(dufTPM)), abdomen = any(!is.na(abdTPM)))
# #will this work to further group by social designation?
# atlas %>%
#   group_by(Species, Sociality) %>%
#   summarise(antennae = any(!is.na(antTPM)), head = any(!is.na(headTPM)), dufour = any(!is.na(dufTPM)), abdomen = any(!is.na(abdTPM)))
# atlas %>%
#   group_by(Species) %>%
#   summarise(n_distinct(Sociality))
# #yes, but there is only one social designation per species, which simplifies analysis
# #can do inverse to see how many species per social designation
# atlas %>%
#   group_by(Sociality) %>%
#   summarise(n_distinct(Species))

#4 social species with all 4 tissues: AAUR, LCAL, LMAR, LPAU
#4 solitary species with all 4 tissues: APUR, LFIG, LLEU, LVIE
#makes for balanced comparison

#now to calculate tau tissue specificity index for all species with data for all 4 tissues
#filtering out rows with NAs for any tissue 
#and filtering for max TPM value > 1 to ensure bare minimum expression in one tissue
allTissuesAtlas  <- filter(atlas, !is.na(antTPM) & !is.na(headTPM) & !is.na(dufTPM) & !is.na(abdTPM) & TPMmax > 1)
summary(allTissuesAtlas)

#turn columns into lists to then apply tissues specificity function to each gene for each tissue
a <- as.list(allTissuesAtlas$antTPM)
b <- as.list(allTissuesAtlas$headTPM)
c <- as.list(allTissuesAtlas$dufTPM)
d <- as.list(allTissuesAtlas$abdTPM)
e <- as.list(allTissuesAtlas$TPMmax)

#this calculates tau tissue specificity with Yanai et al 2005 eqn for our 4 tissues
#have already filtered maximum TPM to be greater than 1
tau4tissues <- pmap(list(a,b,c,d,e),
            function(a,b,c,d,e) (1-(a/e))/3 + (1-(b/e))/3 + (1-(c/e))/3 + (1-(d/e))/3)

#add tau calculations to data frame
allTissuesAtlas$tau4tissues <- unlist(tau4tissues)
allTissuesAtlas$tau4tissues <- round(allTissuesAtlas$tau4tissues, digits=4)

#yay, it worked! write to file
write_tsv(allTissuesAtlas,"tauSpecificity_eachSpecies_4tissues_minMaxTPMof1.txt")

#later - can calculate for species with 3 tissues and append those to this table

#not sure why, but splitting included all rows that were previously filtered, so reimporting to exclude filtered data
allTissuesAtlasFilt <- read_tsv("tauSpecificity_eachSpecies_4tissues_minMaxTPMof1.txt")

#select only relevant TPM info prior to joining to reduce size
allTissuesAtlasFilt <- select(allTissuesAtlasFilt,Name,Sociality,Species,OG,antTPM,headTPM,dufTPM,abdTPM)
#drop rows where "OG" is NA
allTissuesAtlasFilt  <- filter(allTissuesAtlasFilt, !is.na(OG))
#round all numbers to 4 digits
allTissuesAtlasFilt <- allTissuesAtlasFilt %>% mutate_if(is.numeric, round, digits=4)

#split data by species
splitDF <- split(allTissuesAtlasFilt,allTissuesAtlasFilt$Species)
aaur <- splitDF$AAUR
apur <- splitDF$APUR
lcal <- splitDF$LCAL
lfig <- splitDF$LFIG
lleu <- splitDF$LLEU
lmar <- splitDF$LMAR
lpau <- splitDF$LPAU
lvie <- splitDF$LVIE

#having trouble with joining multicopy orthologs or completely eliminating IDs with paralogs
#will keep first entry of each ortholog by species
#IS THIS PROBLEMATIC? RUN BY SARAH AND BEN
aaur_uniqueOG <- distinct(aaur, OG, .keep_all = TRUE)
apur_uniqueOG <- distinct(apur, OG, .keep_all = TRUE)
lcal_uniqueOG <- distinct(lcal, OG, .keep_all = TRUE)
lfig_uniqueOG <- distinct(lfig, OG, .keep_all = TRUE)
lleu_uniqueOG <- distinct(lleu, OG, .keep_all = TRUE)
lmar_uniqueOG <- distinct(lmar, OG, .keep_all = TRUE)
lpau_uniqueOG <- distinct(lpau, OG, .keep_all = TRUE)
lvie_uniqueOG <- distinct(lvie, OG, .keep_all = TRUE)

#join AAUR and APUR by OG; "inner" only returns rows w/ matching values
#join1 <- inner_join(aaur_uniqueOG, apur_uniqueOG, by= "OG", suffix=c(".AAUR",".APUR"))

library(plyr); 
library(dplyr)

#join all data frames for one ortholog copy per species, species that have 4 tissue TPM only (8 species)
#"inner" only returns rows w/ matching values in all
join8s4t <- join_all(list(aaur_uniqueOG,apur_uniqueOG,lcal_uniqueOG,lfig_uniqueOG,
                          lleu_uniqueOG,lmar_uniqueOG,lpau_uniqueOG,lvie_uniqueOG),
                          by = 'OG', type = "inner")                     

write_tsv(join8s4t,"TPM_randoUniqueOG_8sp4tissues.txt")

#frustratingly, this does not give suffixes to column headers as an option, 
#will quickly rename in EXCEL! :) 
join8s4t<-read_tsv("TPM_randoUniqueOG_8sp4tissues.txt")


#now to calculate tau tissue specificity across species
#4 social species with all 4 tissues: AAUR, LCAL, LMAR, LPAU
#4 solitary species with all 4 tissues: APUR, LFIG, LLEU, LVIE

#turn columns into lists to then apply tissues specificity function to each gene for each tissue
#first all social species
a <- as.list(join8s4t$antTPM.AAUR)
b <- as.list(join8s4t$headTPM.AAUR)
c <- as.list(join8s4t$dufTPM.AAUR)
d <- as.list(join8s4t$abdTPM.AAUR)

e <- as.list(join8s4t$antTPM.LCAL)
f <- as.list(join8s4t$headTPM.LCAL)
g <- as.list(join8s4t$dufTPM.LCAL)
h <- as.list(join8s4t$abdTPM.LCAL)

i <- as.list(join8s4t$antTPM.LMAR)
j <- as.list(join8s4t$headTPM.LMAR)
k <- as.list(join8s4t$dufTPM.LMAR)
l <- as.list(join8s4t$abdTPM.LMAR)

m <- as.list(join8s4t$antTPM.LPAU)
n <- as.list(join8s4t$headTPM.LPAU)
o <- as.list(join8s4t$dufTPM.LPAU)
p <- as.list(join8s4t$abdTPM.LPAU)

#then all solitary species (whew, I am not good at automating)
q <- as.list(join8s4t$antTPM.APUR)
r <- as.list(join8s4t$headTPM.APUR)
s <- as.list(join8s4t$dufTPM.APUR)
t <- as.list(join8s4t$abdTPM.APUR)

u <- as.list(join8s4t$antTPM.LFIG)
v <- as.list(join8s4t$headTPM.LFIG)
w <- as.list(join8s4t$dufTPM.LFIG)
x <- as.list(join8s4t$abdTPM.LFIG)

y <- as.list(join8s4t$antTPM.LLEU)
z <- as.list(join8s4t$headTPM.LLEU)
aa <- as.list(join8s4t$dufTPM.LLEU)
ab <- as.list(join8s4t$abdTPM.LLEU)

ac <- as.list(join8s4t$antTPM.LVIE)
ad <- as.list(join8s4t$headTPM.LVIE)
ae <- as.list(join8s4t$dufTPM.LVIE)
af <- as.list(join8s4t$abdTPM.LVIE)

#view(colnames(join8s4t))
#4 social species with all 4 tissues: AAUR, LCAL, LMAR, LPAU
#4 solitary species with all 4 tissues: APUR, LFIG, LLEU, LVIE
#summary(join8s4t[, c(4:7, 10:13, 16:19, 22:25, 28:31, 34:37, 40:43, 46:49)])
#social <- join8s4t[, c(4:7, 16:19, 34:37, 40:43)] 
#solitary <- join8s4t[, c(10:13, 22:25, 28:31, 46:49)] 

#loop through rows and record max TPM value
join8s4t$TPMmax.8sp <- apply(join8s4t[, c(4:7, 10:13, 16:19, 22:25, 28:31, 34:37, 40:43, 46:49)], 1, function(x) max(x, na.rm = TRUE))
join8s4t$TPMmax.4social <- apply(join8s4t[, c(4:7, 16:19, 34:37, 40:43)], 1, function(x) max(x, na.rm = TRUE))
join8s4t$TPMmax.4solitary <- apply(join8s4t[, c(10:13, 22:25, 28:31, 46:49)], 1, function(x) max(x, na.rm = TRUE))

TPMmax8sp <- as.list(join8s4t$TPMmax.8sp)
TPMmax4soc <- as.list(join8s4t$TPMmax.4social)
TPMmax4sol <- as.list(join8s4t$TPMmax.4solitary)


####NOW TO CALCULATE TAU FOR 1) ALL SPECIES, 2) SOCIAL, 3) SOLITARY
#this calculates tau tissue specificity with Yanai et al 2005 eqn for our 4 tissues
#have already filtered maximum TPM to be greater than 1 AMONG 4 TPM VALUES FOR ONE SINGLE SPECIES


#8 species with 4 TPM tissues
tau8sp4tissues <- pmap(list(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,aa,ab,ac,ad,ae,af,TPMmax8sp),
                    function(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,aa,ab,ac,ad,ae,af,TPMmax8sp) 
                      (1-(a/TPMmax8sp))/31 + (1-(b/TPMmax8sp))/31 + (1-(c/TPMmax8sp))/31 + (1-(d/TPMmax8sp))/31
                    + (1-(e/TPMmax8sp))/31 + (1-(f/TPMmax8sp))/31 + (1-(g/TPMmax8sp))/31 + (1-(h/TPMmax8sp))/31
                    + (1-(i/TPMmax8sp))/31 + (1-(j/TPMmax8sp))/31 + (1-(k/TPMmax8sp))/31 + (1-(l/TPMmax8sp))/31
                    + (1-(m/TPMmax8sp))/31 + (1-(n/TPMmax8sp))/31 + (1-(o/TPMmax8sp))/31 + (1-(p/TPMmax8sp))/31
                    + (1-(q/TPMmax8sp))/31 + (1-(r/TPMmax8sp))/31 + (1-(s/TPMmax8sp))/31 + (1-(t/TPMmax8sp))/31
                    + (1-(u/TPMmax8sp))/31 + (1-(v/TPMmax8sp))/31 + (1-(w/TPMmax8sp))/31 + (1-(x/TPMmax8sp))/31
                    + (1-(y/TPMmax8sp))/31 + (1-(z/TPMmax8sp))/31 + (1-(aa/TPMmax8sp))/31 + (1-(ab/TPMmax8sp))/31
                    + (1-(ac/TPMmax8sp))/31 + (1-(ad/TPMmax8sp))/31 + (1-(ae/TPMmax8sp))/31 + (1-(af/TPMmax8sp))/31)

#add tau calculations to data frame
join8s4t$tau.8sp4tissues <- unlist(tau8sp4tissues)
join8s4t$tau.8sp4tissues <- round(join8s4t$tau.8sp4tissues, digits=4)


#4 social species with 4 TPM tissues
tau4soc4tissues <- pmap(list(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,TPMmax4soc),
                       function(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,TPMmax4soc) 
                         (1-(a/TPMmax4soc))/15 + (1-(b/TPMmax4soc))/15 + (1-(c/TPMmax4soc))/15 + (1-(d/TPMmax4soc))/15
                       + (1-(e/TPMmax4soc))/15 + (1-(f/TPMmax4soc))/15 + (1-(g/TPMmax4soc))/15 + (1-(h/TPMmax4soc))/15
                       + (1-(i/TPMmax4soc))/15 + (1-(j/TPMmax4soc))/15 + (1-(k/TPMmax4soc))/15 + (1-(l/TPMmax4soc))/15
                       + (1-(m/TPMmax4soc))/15 + (1-(n/TPMmax4soc))/15 + (1-(o/TPMmax4soc))/15 + (1-(p/TPMmax4soc))/15
                       )

#add tau calculations to data frame
join8s4t$tau.4soc4tissues <- unlist(tau4soc4tissues)
join8s4t$tau.4soc4tissues <- round(join8s4t$tau.4soc4tissues, digits=4)


#4 solitary species with 4 TPM tissues
tau4sol4tissues <- pmap(list(q,r,s,t,u,v,w,x,y,z,aa,ab,ac,ad,ae,af,TPMmax4sol),
                       function(q,r,s,t,u,v,w,x,y,z,aa,ab,ac,ad,ae,af,TPMmax4sol) 
                         (1-(q/TPMmax4sol))/15 + (1-(r/TPMmax4sol))/15 + (1-(s/TPMmax4sol))/15 + (1-(t/TPMmax4sol))/15
                       + (1-(u/TPMmax4sol))/15 + (1-(v/TPMmax4sol))/15 + (1-(w/TPMmax4sol))/15 + (1-(x/TPMmax4sol))/15
                       + (1-(y/TPMmax4sol))/15 + (1-(z/TPMmax4sol))/15 + (1-(aa/TPMmax4sol))/15 + (1-(ab/TPMmax4sol))/15
                       + (1-(ac/TPMmax4sol))/15 + (1-(ad/TPMmax4sol))/15 + (1-(ae/TPMmax4sol))/15 + (1-(af/TPMmax4sol))/15)

#add tau calculations to data frame
join8s4t$tau.4sol4tissues <- unlist(tau4sol4tissues)
join8s4t$tau.4sol4tissues <- round(join8s4t$tau.4sol4tissues, digits=4)


#write cross-species ortholog-based tissue specificity metrics
#Question: has gene expression normalization been conducted in a way that facilitates direct comparison of TPM across species?
write_tsv(join8s4t,"tauTissueSpec_randoUniqueOG_8sp4tissues_4soc_4sol.txt")
