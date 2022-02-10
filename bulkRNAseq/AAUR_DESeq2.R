#Load packages----
library(ggplot2); theme_set(theme_bw())
library(edgeR)
library(ggrepel)
library(ggfortify)
library(sva)
library(gridExtra)
library(DESeq2)

#Extract metadata file----
samples<-read.csv("AAUR_metadata.csv")
counts<-read.table("AAUR_featureCounts.tsv",header=TRUE)

#check samples match
samples$Sample == colnames(counts)

# separate tissues:
brain_samples<-samples[which(samples$Tissue=="Brain"),]
brain_counts<-counts[,brain_samples$Sample]
dim(brain_samples)
dim(brain_counts)
# 18 samples, 12511 genes

# filter out genes with low coverage
dge.0 <- DGEList(brain_counts)
cutoff <- 1
keep <- rowSums(cpm(dge.0)>=cutoff)>=4
dge.cpmfilt <- dge.0[keep,] 
dim(dge.cpmfilt)
# 10137 genes pass low-coverage filtering

counts <- dge.cpmfilt$counts
coldata <- brain_samples

#check samples match
colnames(counts) == coldata$Sample

dds <- DESeqDataSetFromMatrix(counts, coldata, design= ~Group)
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

# Output normalized counts
norm.counts<-counts(dds, normalized=T)
# write.csv(norm.counts,"AAUR_DESeq2_norm_counts_BRAIN.csv")

## Get DE results
WvQResults <- results(dds.results, name="Group_W_vs_Q", alpha = 0.99, pAdjustMethod="BH")
WvQResults <- WvQResults[order(WvQResults$padj),]
plotMA(WvQResults)
summary(WvQResults, na.rm=TRUE)

#output results
write.csv(WvQResults,"AAUR_brain_DESeq2.csv")

# fat body:
samples<-read.csv("AAUR_metadata.csv")
counts<-read.table("AAUR_featureCounts.tsv",header=TRUE)

fb_samples<-samples[which(samples$Tissue=="Fatbody"),]
fb_counts<-counts[,fb_samples$Sample]
dim(fb_samples)
dim(fb_counts)
# 17 samples, 12511 genes

# filter out genes with low coverage
dge.0 <- DGEList(fb_counts)
cutoff <- 1
keep <- rowSums(cpm(dge.0)>=cutoff)>=4
dge.cpmfilt <- dge.0[keep,] 
dim(dge.cpmfilt) 
# 9311 genes pass low-coverage filtering

counts <- dge.cpmfilt$counts
coldata <- fb_samples

#check samples match
colnames(counts) == coldata$Sample

dds <- DESeqDataSetFromMatrix(counts, coldata, design= ~Group)
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

# Output normalized counts
norm.counts<-counts(dds, normalized=T)
#write.csv(norm.counts,"AAUR_filtered_DESeq2_norm_counts_FATBODY.csv")

## Get DE results
WvQResults <- results(dds.results, name="Group_W_vs_Q", alpha = 0.99, pAdjustMethod="BH")
WvQResults <- WvQResults[order(WvQResults$padj),]
plotMA(WvQResults)
summary(WvQResults, na.rm=TRUE)

#output results
write.csv(WvQResults,"AAUR_fatbody_DESeq2.csv")