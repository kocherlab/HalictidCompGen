#------------------------------------------------------------------------------- 
# Load analysis options and libraries
#------------------------------------------------------------------------------- 

#libraries
pkgs <- c("Seurat", "dplyr", "MAST", "RColorBrewer", "data.table", "corrplot",
          "ggplot2", "rstudioapi", "SingleCellExperiment", "cowplot", "MASS", "GeneOverlap",
          "AnnotationDbi", "plyr", "dittoSeq", "multtest")
invisible(lapply(pkgs, library, character.only = T))
rm(pkgs); gc()

library("org.Dm.eg.db") # We will be using AnnotationDbi to extract Flybase gene names and symbols

#options
options(stringsAsFactors = F, scipen = 9999)
set.seed(12345)

setwd(dirname(getActiveDocumentContext()$path)) # Set working directory to wherever you have this script.
# You will need to load in the contents of "filtered_feature_bc_matrix," which is outputted by Cell Ranger V6.
# The contents of this directory include three files: barcodes.tsv.gz, features.tsv.gz, matrix.mtx.gz. 
# Raw and processed data are accessible via the SRA link associated with this manuscript.
# Remember: LALB IDs must be converted to LZEP IDs before this analysis is run.

### Quick-start, all the relevant objects created in this analysis. ###

hal.combined.sct  <- readRDS(file = "../../out/hal.combined.sct_10PCs_DIM10_res0.6.rds") # Seurat object including all cells from LZEP and LALB bees.
cluster_markers_hal_integrated_Dims10_res0.6_FDR0.05 <- read.csv("../../out/cluster_markers_hal_integrated_Dims10_res0.6_FDR0.05.csv") #Cluster-defining genes; MAST algorithm with default settings.
hal_glia <- readRDS(file = "../../out/hal_glia_6PCs_res0.3.rds") # Seurat object consisting of glial cell subset (Cluster 7 from total cell analysis0. )
hal_glia_markers_res0.3_FDR0.05 <- read.csv("../../out/hal_glia_markers_res_FDR0.05.csv") #Cluster-defining genes from glia subsetting; MAST algorithm with default settings.

### End quick-start ###

lookup <- read.csv("../../src/lzep_hb_dmel.csv", fileEncoding = 'UTF-8-BOM')
lookup$gene <- gsub("_", "-", lookup$gene) # For converting LZEP genes / OGs to honey bee or fly orthologs. 

# Read in Cell Ranger V6 output for two LZEP and two LALB samples.

lzep_queen.data <- Read10X(data.dir = "../../data/LZEP_QW/lzep_queen_trimmed/filtered_feature_bc_matrix/")
lzep_queen <- CreateSeuratObject(counts = lzep_queen.data, project = "twoway_halictid_merge", min.cells = 3, min.features = 200) 
lzep_queen@meta.data$species <- "lzep"
lzep_queen@meta.data$sample <- "queen"

lzep_worker.data <- Read10X(data.dir = "../../data/LZEP_QW/lzep_worker_trimmed/filtered_feature_bc_matrix/")
lzep_worker <- CreateSeuratObject(counts = lzep_worker.data, project = "twoway_halictid_merge", min.cells = 3, min.features = 200) 
lzep_worker@meta.data$species <- "lzep"
lzep_worker@meta.data$sample <- "worker"

lzep_merge <- merge(lzep_queen, lzep_worker, project = "lzep")
lzep_merge
# An object of class Seurat 
# 9354 features across 1389 samples within 1 assay 
# Active assay: RNA (9354 features, 0 variable features)

lalb_A.data <- Read10X(data.dir = "../../data/LALB_converted_to_zep/A_trimmed/filtered_feature_bc_matrix/")
lalb_A <- CreateSeuratObject(counts = lalb_A.data, project = "twoway_halictid_merge", min.cells = 3, min.features = 200) 
lalb_A@meta.data$species <- "lalb"
lalb_A@meta.data$sample <- "A"

lalb_B.data <- Read10X(data.dir = "../../data/LALB_converted_to_zep/B_trimmed/filtered_feature_bc_matrix/")
lalb_B <- CreateSeuratObject(counts = lalb_B.data, project = "twoway_halictid_merge", min.cells = 3, min.features = 200) 
lalb_B@meta.data$species <- "lalb"
lalb_B@meta.data$sample <- "B"

lalb_merge <- merge(lalb_A, lalb_B, project = "lalb")
lzep_merge
# An object of class Seurat 
# 9354 features across 1389 samples within 1 assay 
# Active assay: RNA (9354 features, 0 variable features)

lzep_lalb_merge <- merge(lalb_merge, lzep_merge, add.cell.ids = c("lalb","lalb"), project = "twoway_halictid_merge")  
lzep_lalb_merge
# An object of class Seurat 
# 12249 features across 6716 samples within 1 assay 
# Active assay: RNA (12249 features, 0 variable features)

table(lzep_lalb_merge@meta.data$sample)
# A      B  queen worker 
# 2320   3007    692    697 

split_hal <- SplitObject(lzep_lalb_merge, split.by = "species")  # "hal" = "halictid," dear reader.

split_hal
# $lalb
# An object of class Seurat 
# 12249 features across 5327 samples within 1 assay 
# Active assay: RNA (12249 features, 0 variable features)
# 
# $lzep
# An object of class Seurat 
# 12249 features across 1389 samples within 1 assay 
# Active assay: RNA (12249 features, 0 variable features)

split_hal <- lapply(X = split_hal, FUN = SCTransform) # We will integrate by species, not by sample. Either is acceptable but I fear that the latter will overfit, and the former is more in-line with the biology. 
# We are using the SCTransform pipeline for integration and normalization.

features <- SelectIntegrationFeatures(object.list = split_hal, nfeatures = 3000)
sct.list <- PrepSCTIntegration(object.list = split_hal, anchor.features = features)

hal.anchors <- FindIntegrationAnchors(object.list = sct.list, normalization.method = "SCT",
                                      anchor.features = features)

hal.combined.sct <- IntegrateData(anchorset = hal.anchors, normalization.method = "SCT")

hal.combined.sct <- RunPCA(hal.combined.sct, verbose = FALSE)

x11()
ElbowPlot(hal.combined.sct)
# Go with 10 PCs. Up to 15 could probably be justified but I don't want to overcluster with so few cells. 

hal.combined.sct <- RunUMAP(hal.combined.sct, reduction = "pca", dims = 1:10)
hal.combined.sct <- FindNeighbors(hal.combined.sct, reduction = "pca", dims = 1:10)
hal.combined.sct <- FindClusters(hal.combined.sct, resolution = 0.6)
hal.combined.sct # Results are going to slightly change with different seeds, version of R/Seurat, etc. But here is a good place to double check your results.
# An object of class Seurat 
# 26623 features across 6716 samples within 3 assays 
# Active assay: integrated (2733 features, 2733 variable features)
# 2 other assays present: RNA, SCT
# 2 dimensional reductions calculated: pca, umap

# saveRDS(hal.combined.sct, file = "../../out/hal.combined.sct_10PCs_DIM10_res0.6.rds") # Creates the first object in the quick-start.

# pdf("../../out/hal_UMAP.pdf")
DimPlot(hal.combined.sct, reduction = "umap", label = FALSE) & NoAxes()
# dev.off()
# 
# These are some plots that help negotiate the data. Use if you like!
# p1 <- DimPlot(hal.combined.sct, reduction = "umap", label = TRUE) & NoAxes()
# p2 <- DimPlot(hal.combined.sct, reduction = "umap", group.by = "species")& NoAxes()
# p3 <- DimPlot(hal.combined.sct, reduction = "umap", group.by = "sample")& NoAxes()
# p4 <- dittoBarPlot(hal.combined.sct, "species", group.by = "integrated_snn_res.0.6")
# p5 <- dittoBarPlot(hal.combined.sct, "sample", group.by = "integrated_snn_res.0.6")


# Visualize some cell-type-specific markers.
DefaultAssay(hal.combined.sct) <- "SCT" #Do NOT use integrated object here, and especially do not use it for calling DEGs. See Seurat Github -- this has been a hot topic of discussion.

#Identify neuronal and glial markers
# pdf("../../out/neuron_glia_markers_UMAP.pdf")
FeaturePlot(hal.combined.sct, features = c("LZEP-02492", "LZEP-09329", "LZEP-06452", "LZEP-00406"), pt.size = 1, label = F) & NoLegend() & NoAxes()# fne and Syt1, neuron markers; repo and borderless, glial markers
# dev.off()

# Identify the localization of apolpp in the sc data.
# pdf("../../out/all_cells_dotplot.pdf")
DotPlot(hal.combined.sct, 
        features = c( 'LZEP-00406',#borderless
        'LZEP-08903', #apoltp
        'LZEP-01422', #apolpp
        'LZEP-07620', #lpr1
        'LZEP-05759', #lpr2
        'LZEP-01788'), #Jhe2
        scale = T,
        scale.by = "radius",
        col.min = 0) + coord_flip()
# dev.off()

cluster_markers_hal_integrated_Dims10_res0.6 <- FindAllMarkers(hal.combined.sct, only.pos = T, test.use = "MAST", assay = "SCT", slot = "data")
cluster_markers_hal_integrated_Dims10_res0.6_FDR0.05 <- subset(cluster_markers_hal_integrated_Dims10_res0.6, p_val_adj < 0.05)
cluster_markers_hal_integrated_Dims10_res0.6_FDR0.05 <- join(cluster_markers_hal_integrated_Dims10_res0.6_FDR0.05, lookup, by = "gene", type = "left")
cluster_markers_hal_integrated_Dims10_res0.6_FDR0.05 <- cluster_markers_hal_integrated_Dims10_res0.6_FDR0.05[order(cluster_markers_hal_integrated_Dims10_res0.6_FDR0.05$cluster, 
                                                                                                                   -cluster_markers_hal_integrated_Dims10_res0.6_FDR0.05$avg_log2FC),]
# Add Dmel information.
geneNames <- mapIds(org.Dm.eg.db,
  keys = cluster_markers_hal_integrated_Dims10_res0.6_FDR0.05$FB,
  column="GENENAME",
  keytype="ENSEMBL",
  multiVals="first")

geneNames[sapply(geneNames, is.null)] <- NA
geneNames <- as.data.frame(unlist(geneNames))
cluster_markers_hal_integrated_Dims10_res0.6_FDR0.05$name <- geneNames$`unlist(geneNames)`

geneSym <- mapIds(org.Dm.eg.db,
                  keys = cluster_markers_hal_integrated_Dims10_res0.6_FDR0.05$FB,
                  column="SYMBOL",
                  keytype="ENSEMBL",
                  multiVals="first")
geneSym[sapply(geneSym, is.null)] <- NA
geneSym <- as.data.frame(unlist(geneSym))
cluster_markers_hal_integrated_Dims10_res0.6_FDR0.05$sym <- geneSym$`unlist(geneSym)`

# write.csv(cluster_markers_hal_integrated_Dims10_res0.6_FDR0.05, file = "../../out/cluster_markers_hal_integrated_Dims10_res0.6_FDR0.05.csv") # This creates the second quick-start object.

# Subset glia
# This will follow a similar pipeline to the total cell analysis, now using fewer PCs (6) and a smaller resolution (0.3) to account for the fewer cells.

hal_glia <- subset(hal.combined.sct, idents = 7)
hal_glia
# An object of class Seurat 
# 26623 features across 523 samples within 3 assays 
# Active assay: integrated (2733 features, 2733 variable features)
# 2 other assays present: RNA, SCT
# 2 dimensional reductions calculated: pca, umap

DefaultAssay(hal_glia) <- "integrated"

hal_glia <- RunPCA(hal_glia, verbose = FALSE)
hal_glia <- RunUMAP(hal_glia, dims = 1:6)
hal_glia <- FindNeighbors(hal_glia, dims = 1:6)
hal_glia <- FindClusters(hal_glia, resolution = 0.3)
hal_glia
# An object of class Seurat 
# 26623 features across 523 samples within 3 assays 
# Active assay: integrated (2733 features, 2733 variable features)
# 2 other assays present: RNA, SCT
# 2 dimensional reductions calculated: pca, umap

# saveRDS(hal_glia, file = "../../out/hal_glia_6PCs_res0.3.rds") # This creates the third quick-start object.
# 

# pdf("../../out/glia_recluster_UMAP.pdf")
DimPlot(hal_glia, reduction = "umap", label = FALSE) & NoAxes()
# dev.off()

DefaultAssay(hal_glia) <- "SCT"

hal_glia_markers <- FindAllMarkers(hal_glia, only.pos = T, test.use = "MAST", assay = "SCT", slot = "data") 
hal_glia_markers_res0.3_FDR0.05 <- hal_glia_markers[hal_glia_markers$p_val_adj < 0.05,]
names(hal_glia_markers_res0.3_FDR0.05)[names(hal_glia_markers_res0.3_FDR0.05) == "X"] <- "gene"
hal_glia_markers_res0.3_FDR0.05 <- join(hal_glia_markers_res0.3_FDR0.05, lookup, by = "gene", type = "left")
geneName_glia <- mapIds(org.Dm.eg.db,
                           keys = hal_glia_markers_res0.3_FDR0.05$FB,
                           column="GENENAME",
                           keytype="ENSEMBL",
                           multiVals="first")

geneName_glia[sapply(geneName_glia, is.null)] <- NA
geneName_glia <- as.data.frame(unlist(geneName_glia))
hal_glia_markers_res0.3_FDR0.05$name <- geneName_glia$`unlist(geneName_glia)`

geneSym_glia <- mapIds(org.Dm.eg.db,
                          keys = hal_glia_markers_res0.3_FDR0.05$FB,
                          column="SYMBOL",
                          keytype="ENSEMBL",
                          multiVals="first")

geneSym_glia[sapply(geneSym_glia, is.null)] <- NA
geneSym_glia <- as.data.frame(unlist(geneSym_glia))
hal_glia_markers_res0.3_FDR0.05$sym <- geneSym_glia$`unlist(geneSym_glia)`

# 
# write.csv(hal_glia_markers_res0.3_FDR0.05, file = "../../out/hal_glia_markers_res_FDR0.05.csv") # This creates the fourth quick-start object.

# pdf("../../out/glial_cells_dotplot.pdf")
DotPlot(hal_glia, features = c('LZEP-08903', #apoltp
                               'LZEP-01422', #apolpp
                               'LZEP-07647', #Moody
                               "LZEP-02440", #zyd
                               'LZEP-01788'),#jheh2
        scale = T,
        scale.by = "radius",
        col.min = 0) + coord_flip()
# dev.off()


### Overlap analysis ###

# We will use the LZEP and AAUR bulk brain RNA-Sequencing data to test for overlap of OGs in the bulk and sc DEGs. The former are queen/worker-defining; the latter are cell-type-defining.
# This gives us an idea of where the queen and/or worker-based genes may be localized in the brain.

# First, we determine the OG universe. Only remove NAs for this calculation. Leave them in for everything else: they must count toward the list sizes for the hypergeometric test or the test will be incorrectly performed.
# Read in conversion table that contains halictid genes and OGs for all species analysed in this study.
coding_info <- read.table("../../data/coding_allinfo_synthesized.txt", header = T)

DefaultAssay(hal.combined.sct) <- "RNA" # Make sure you do this! Or else you will miss the majority of tested genes in the sc universe.
sc_gene_universe <- as.data.frame(rownames(hal.combined.sct))
colnames(sc_gene_universe)[1] <- "LZEP_gene" # Because we converted LALB to LZEP, if the conversion failed and the gene is still a LALB gene, then it cannot have an OG.
sc_gene_universe$LZEP_gene <- gsub("-", "_", sc_gene_universe$LZEP_gene)
sc_ogs <- na.omit(join(sc_gene_universe, coding_info, by = "LZEP_gene", type = "left")$OG)
length(sc_ogs)
# [1] 7540
# In the sc data, 7540 genes have OGs.

conversion <- data.frame(
  OG = coding_info$OG,
  LZEP_gene = coding_info$LZEP_gene,
  AAUR_gene = coding_info$AAUR_gene
)

# Read in bulk brain DEGs in the queen vs. worker pairwise comparison. 

lzep_WvQ <- read.csv("../../../../Shared drives/KocherLab/2021_CompGen_Bulk_RNAseq/DE_analysis/LZEP_filtered_CombatNest_brain_WvQ.csv")
colnames(lzep_WvQ)[1] <- "LZEP_gene"
lzep_WvQ_ogs <- join(lzep_WvQ, conversion, by = "LZEP_gene", type = "left")

aaur_WvQ <- read.csv("../../../../Shared drives/KocherLab/2021_CompGen_Bulk_RNAseq/DE_analysis/AAUR_filtered_brain_WvQ.csv")
colnames(aaur_WvQ)[1] <- "AAUR_gene"
aaur_WvQ_ogs <- join(aaur_WvQ, conversion, by = "AAUR_gene", type = "left")

table(duplicated(intersect(sc_ogs, intersect(lzep_WvQ_ogs$OG, aaur_WvQ_ogs$OG))))
# FALSE 
# 6668
# So, 6668 is our universe size, shared between the three lists, after NAs are removed.

# We should *not* remove NAs for this next step.

# OGs that define each of the 11 sc clusters
cluster_markers_hal_integrated_Dims10_res0.6_FDR0.05 <- read.csv("../../out/cluster_markers_hal_integrated_Dims10_res0.6_FDR0.05.csv")
cluster_markers_hal_integrated_Dims10_res0.6_FDR0.05[is.na(cluster_markers_hal_integrated_Dims10_res0.6_FDR0.05)]="NA_sc" # NA entries must be counted in list size, but must be renamed in to avoid inflating overlaps by 1 (as there are NAs in other lists, too).
cluster_ogs <- split(cluster_markers_hal_integrated_Dims10_res0.6_FDR0.05, f = cluster_markers_hal_integrated_Dims10_res0.6_FDR0.05$cluster)
cluster_ogs <- lapply(cluster_ogs, "[", , "OG")
# We now have a list of 11 items, each containing the cluster-specific OGs. 

# Establish differential expression cutoffs and separate queen- and worker-biased genes (from the queen vs. worker pairwise comparison). Use these data to make lists.
lzep_WvQ_FDR05 <- subset(lzep_WvQ_ogs, padj < 0.05)
lzep_deg_list <- list(subset(lzep_WvQ_FDR05, log2FoldChange > 0)$OG, subset(lzep_WvQ_FDR05, log2FoldChange < 0)$OG)
names(lzep_deg_list) <- c("worker", "queen")

aaur_WvQ_ogs_FDR05 <- subset(aaur_WvQ_ogs, padj < 0.05)
aaur_deg_list <- list(subset(aaur_WvQ_ogs_FDR05, log2FoldChange > 0)$OG, subset(aaur_WvQ_ogs_FDR05, log2FoldChange < 0)$OG)
names(aaur_deg_list) <- c("worker", "queen")  

worker_degs <- intersect(lzep_deg_list$worker, aaur_deg_list$worker)
queen_degs <- intersect(lzep_deg_list$queen, aaur_deg_list$queen)

deg_intersect <- list(worker_degs, queen_degs)
names(deg_intersect) <- c("worker", "queen")

scClusters_bothSpecies_vs_bulkDEGs_gom_obj <- newGOM(cluster_ogs, deg_intersect, genome.size = 6668)

print(scClusters_bothSpecies_vs_bulkDEGs_gom_obj)

or <- as.matrix(getMatrix(scClusters_bothSpecies_vs_bulkDEGs_gom_obj, "odds.ratio"))
p <- as.matrix(getMatrix(scClusters_bothSpecies_vs_bulkDEGs_gom_obj, "pval"))
p.adjust(p, method = "BH")
# overlap between worker DEGs and cluter 7 FDR = 0.001519925 

# Make overlap figure.
# pdf("../../out/odds_ratio.pdf", width = 5, height = 5)
par(mfrow = c(2,1))
m <- or[,1:2]
colfunc <- colorRampPalette(c("grey", "blue"))
colfunc(14)

pal <- plot(rep(0,20),col=colfunc(20),pch=15,cex=3)

corrplot(m, is.corr = FALSE, 
         method = c("color"), 
         cl.lim = c(0,20),
         col=colfunc(20), 
         p.mat = as.matrix(p),
         insig = "label_sig")
# dev.off()

### End of analysis ###

sessionInfo()
# R version 4.1.2 (2021-11-01)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 10 x64 (build 19044)
# 
# Matrix products: default
# 
# locale:
#   [1] LC_COLLATE=English_United States.1252  LC_CTYPE=English_United States.1252    LC_MONETARY=English_United States.1252 LC_NUMERIC=C                           LC_TIME=English_United States.1252    
# 
# attached base packages:
#   [1] stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] multtest_2.50.0             dittoSeq_1.6.0              plyr_1.8.6                  AnnotationDbi_1.56.2        GeneOverlap_1.30.0          MASS_7.3-54                 cowplot_1.1.1              
# [8] rstudioapi_0.13             ggplot2_3.3.5               corrplot_0.92               data.table_1.14.2           RColorBrewer_1.1-2          MAST_1.20.0                 SingleCellExperiment_1.16.0
# [15] SummarizedExperiment_1.24.0 Biobase_2.54.0              GenomicRanges_1.46.1        GenomeInfoDb_1.30.1         IRanges_2.28.0              S4Vectors_0.32.3            BiocGenerics_0.40.0        
# [22] MatrixGenerics_1.6.0        matrixStats_0.61.0          dplyr_1.0.7                 SeuratObject_4.0.4          Seurat_4.1.0               
# 
# loaded via a namespace (and not attached):
#   [1] Rtsne_0.15             colorspace_2.0-2       deldir_1.0-6           ellipsis_0.3.2         ggridges_0.5.3         XVector_0.34.0         spatstat.data_2.1-2    leiden_0.3.9          
# [9] listenv_0.8.0          bit64_4.0.5            ggrepel_0.9.1          fansi_1.0.2            codetools_0.2-18       splines_4.1.2          cachem_1.0.6           polyclip_1.10-0       
# [17] jsonlite_1.7.3         ica_1.0-2              cluster_2.1.2          png_0.1-7              pheatmap_1.0.12        uwot_0.1.11            shiny_1.7.1            sctransform_0.3.3     
# [25] spatstat.sparse_2.1-0  compiler_4.1.2         httr_1.4.2             Matrix_1.3-4           fastmap_1.1.0          lazyeval_0.2.2         cli_3.1.1              later_1.3.0           
# [33] htmltools_0.5.2        tools_4.1.2            igraph_1.2.11          gtable_0.3.0           glue_1.6.1             GenomeInfoDbData_1.2.7 RANN_2.6.1             reshape2_1.4.4        
# [41] Rcpp_1.0.8             scattermore_0.7        Biostrings_2.62.0      vctrs_0.3.8            nlme_3.1-153           lmtest_0.9-39          stringr_1.4.0          globals_0.14.0        
# [49] mime_0.12              miniUI_0.1.1.1         lifecycle_1.0.1        irlba_2.3.5            gtools_3.9.2           goftest_1.2-3          future_1.23.0          zlibbioc_1.40.0       
# [57] zoo_1.8-9              scales_1.1.1           spatstat.core_2.3-2    promises_1.2.0.1       spatstat.utils_2.3-0   parallel_4.1.2         memoise_2.0.1          reticulate_1.24       
# [65] pbapply_1.5-0          gridExtra_2.3          rpart_4.1-15           RSQLite_2.2.9          stringi_1.7.6          caTools_1.18.2         rlang_1.0.0            pkgconfig_2.0.3       
# [73] bitops_1.0-7           lattice_0.20-45        ROCR_1.0-11            purrr_0.3.4            tensor_1.5             patchwork_1.1.1        htmlwidgets_1.5.4      bit_4.0.4             
# [81] tidyselect_1.1.1       parallelly_1.30.0      RcppAnnoy_0.0.19       magrittr_2.0.2         R6_2.5.1               gplots_3.1.1           generics_0.1.2         DBI_1.1.2             
# [89] DelayedArray_0.20.0    withr_2.4.3            mgcv_1.8-38            pillar_1.7.0           fitdistrplus_1.1-6     KEGGREST_1.34.0        survival_3.2-13        abind_1.4-5           
# [97] RCurl_1.98-1.5         tibble_3.1.6           future.apply_1.8.1     crayon_1.4.2           KernSmooth_2.23-20     utf8_1.2.2             spatstat.geom_2.3-1    plotly_4.10.0         
# [105] grid_4.1.2             blob_1.2.2             digest_0.6.29          xtable_1.8-4           tidyr_1.2.0            httpuv_1.6.5           munsell_0.5.0          viridisLite_0.4.0   
