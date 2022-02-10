## Trying combat-seq on Nest rather than Nest as block in DESeq2
#Load packages----
library(ggplot2); theme_set(theme_bw())
library(edgeR)
library(ggrepel)
library(ggfortify)
library(sva)
library(gridExtra)
library(DESeq2)

#Extract metadata file----
samples<-read.csv("LZEP_metadata.csv")
counts<-read.table("LZEP_featureCounts.tsv",header=TRUE)

#check samples match
samples$Sample == colnames(counts)

# separate tissues:
brain_samples<-samples[which(samples$Tissue=="Brain"),]
brain_counts<-counts[,brain_samples$Sample]
dim(brain_samples)
dim(brain_counts)
# 16 samples, 12589 genes

#Need to remove samples that don't have a nestmate: "Zeph_37_brain", "Zeph_91_brain"
brain_samples_filt<-brain_samples[which(!(brain_samples$Sample  %in% c("Zeph_37_brain","Zeph_91_brain"))),]
brain_counts_filt<-brain_counts[,brain_samples_filt$Sample]
dim(brain_samples_filt)
dim(brain_counts_filt)

# filter out genes with low coverage
dge.0 <- DGEList(brain_counts_filt)
cutoff <- 1
keep <- rowSums(cpm(dge.0)>=cutoff)>=4
dge.cpmfilt <- dge.0[keep,] 
dim(dge.cpmfilt) 
# 10126 genes pass low-coverage filtering

counts <- dge.cpmfilt$counts
coldata <- brain_samples_filt

#check samples match
colnames(counts) == coldata$Sample

# remove Nest effects with combat-seq
Nest=brain_samples_filt$Nest
adjusted <- ComBat_seq(as.matrix(counts), batch=Nest, group=NULL)
dge <- DGEList(counts=adjusted)
d<-calcNormFactors(dge, method="TMM")
logCPM <- cpm(d, log = T, prior.count = 3)
counts.transposed <- t(logCPM)
p <- prcomp(counts.transposed)
autoplot(p, data=brain_samples_filt, colour = "Group", shape = "Group", size=3, label=FALSE, main="Brain") + 
  geom_text_repel(label=brain_samples_filt$Group) +theme_bw()+theme_classic()
autoplot(p, data=brain_samples_filt, colour = "Group", shape = "Group", size=3, label=FALSE, main="Brain, Post-combat") + 
  geom_text_repel(label=brain_samples_filt$Nest) +theme_bw()+theme_classic()


dds <- DESeqDataSetFromMatrix(adjusted, coldata, design= ~Group)
dds <- estimateSizeFactors(dds)

### Transform counts for data visualization
rld <- rlog(dds, blind=TRUE)

### Plot PCA 
plotPCA(rld, intgroup="Group")

# Filter low count genes
dds <- dds[ rowSums(counts(dds)) >= 1, ]

## Check sample coverage
barplot(colSums(counts(dds))*1e-6, names=colnames(dds), ylab="Library size (millions)")

## Make DESeq table
dds.results <- DESeq(dds)
head(results(dds.results, tidy=TRUE))
resultsNames(dds.results)

## Plot Dispersion Estimates
plotDispEsts(dds.results)

## Get DE results
WvQResults <- results(dds.results, name="Group_W_vs_Q", alpha = 0.99, pAdjustMethod="BH")
WvQResults <- WvQResults[order(WvQResults$padj),]
plotMA(WvQResults)
summary(WvQResults, na.rm=TRUE)

# Output normalized counts
norm.counts<-counts(dds, normalized=T)
#write.csv(norm.counts,"LZEP_CombatNest_DESeq2_norm_counts_BRAIN.csv")

#output results
write.csv(WvQResults,"LZEP_CombatNest_brain_DESeq2.csv")

##########################################################
# fat body:
samples<-read.csv("LZEP_metadata_filtered.csv")
counts<-read.table("LZEP_featureCounts.tsv",header=TRUE)

fb_samples<-samples[which(samples$Tissue=="Fatbody"),]
fb_counts<-counts[,fb_samples$Sample]
dim(fb_samples)
dim(fb_counts)
# 14 samples, 12589 genes

#Need to remove samples that don't have a nestmate: "Zeph_37_FB", ""Zeph_36_FB"
fb_samples_filt<-fb_samples[which(!(fb_samples$Sample  %in% c("Zeph_37_FB","Zeph_36_FB"))),]
fb_counts_filt<-fb_counts[,fb_samples_filt$Sample]
dim(fb_samples_filt)
dim(fb_counts_filt)

# filter out genes with low coverage
dge.0 <- DGEList(fb_counts_filt)
cutoff <- 1
keep <- rowSums(cpm(dge.0)>=cutoff)>=4
dge.cpmfilt <- dge.0[keep,] 
dim(dge.cpmfilt) 
# 9383 genes pass low-coverage filtering

counts <- dge.cpmfilt$counts
coldata <- fb_samples_filt

#check samples match
colnames(counts) == coldata$Sample

# remove Nest effects with combat-seq
Nest=fb_samples_filt$Nest
adjusted <- ComBat_seq(as.matrix(counts), batch=Nest, group=NULL)
dge <- DGEList(counts=adjusted)
d<-calcNormFactors(dge, method="TMM")
logCPM <- cpm(d, log = T, prior.count = 3)
counts.transposed <- t(logCPM)
p <- prcomp(counts.transposed)
autoplot(p, data=fb_samples_filt, colour = "Group", shape = "Group", size=3, label=FALSE, main="Fatbody") + 
  geom_text_repel(label=fb_samples_filt$Group) +theme_bw()+theme_classic()
autoplot(p, data=fb_samples_filt, colour = "Group", shape = "Group", size=3, label=FALSE, main="Fatbody, Post-combat") + 
  geom_text_repel(label=fb_samples_filt$Nest) +theme_bw()+theme_classic()

dds <- DESeqDataSetFromMatrix(adjusted, coldata, design= ~Group)
dds <- estimateSizeFactors(dds)

### Transform counts for data visualization
rld <- rlog(dds, blind=TRUE)

### Plot PCA 
plotPCA(rld, intgroup="Group")

# Filter low count genes
dds <- dds[ rowSums(counts(dds)) >= 1, ]

## Check sample coverage
barplot(colSums(counts(dds))*1e-6, names=colnames(dds), ylab="Library size (millions)")

## Make DESeq table
dds.results <- DESeq(dds)
head(results(dds.results, tidy=TRUE))
resultsNames(dds.results)

## Plot Dispersion Estimates
plotDispEsts(dds.results)

## Get DE results
WvQResults <- results(dds.results, name="Group_W_vs_Q", alpha = 0.99, pAdjustMethod="BH",independentFiltering = FALSE)
WvQResults <- WvQResults[order(WvQResults$padj),]
plotMA(WvQResults)
summary(WvQResults, na.rm=TRUE)

# Output normalized counts
norm.counts<-counts(dds, normalized=T)
#write.csv(norm.counts,"LZEP_CombatNest_DESeq2_norm_counts_FATBODY.csv")

#output results
write.csv(WvQResults,"LZEP_CombatNest_fatbody_DESeq2.csv")