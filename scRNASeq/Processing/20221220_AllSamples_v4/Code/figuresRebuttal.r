## AJ - 03062024:
# See mail Yosta Vegting 24052024

#### Set default for a good start ####
source("initAnalysis.r")
init()

# Create a directory to hold results and plots
outDir <- "finalFigures"
dir.create(outDir)

dir.create(paste0(outDir,"/20240522_Rebuttal"))
figDir <- paste0(outDir,"/20240522_Rebuttal")

library(Seurat)
library(dplyr)
library(tidyr)
library(cowplot)
library(ggplot2)
library(alluvial)
library(clustree)
library(patchwork)

# Color function
# Also get the correct colors
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

#### Read the data ####
# Before reclustering
seurat <- readRDS("9_kidney_reference_mapping/Mito10.integrated.clustered.annot.kidney_immune_reference.rds")
DefaultAssay(seurat) <- "RNA"

# The monocytes and macrophages are in clusters 5, 7, 9 and 11 at resolution 0.6
# These clusters are assigned after integration, make sure you have selected the right assay (see line 31)
selClusters.seurat <- subset(x = seurat, subset = integrated_snn_res.0.6 %in% c("5","7","9","11"))
dim(selClusters.seurat)
# [1] 22769  4194

# How are they spread over the samples?
table(selClusters.seurat$orig.ident)
# MOMA17 MOMA302  MOMA52  MOMA57  MOMA67  MOMA68  MOMA72 
#    107     201     276    1517     292    1538     263

# Drop all superfluous factor levels
selClusters.seurat@meta.data <- droplevels(selClusters.seurat@meta.data)

# Plot
# https://stackoverflow.com/questions/60016390/set-axes-limits-in-patchwork-when-combining-ggplot2-objects

pdf(paste(figDir, "dimPlot_selectedMyeloidClusters_beforeAndAfter.pdf", sep="/"), width = 14, height = 7)
  Idents(seurat) <- "integrated_snn_res.0.6"
  p1 <- DimPlot(seurat, label=TRUE)
  p1 <- p1 + ggtitle("Before selection")

  Idents(selClusters.seurat) <- "integrated_snn_res.0.6"
  p2 <- DimPlot(selClusters.seurat, label=TRUE)
  p2 <- p2 + ggtitle("After selection") + scale_x_continuous(limits = c(-10, 10))
  
  p <- p1 | p2
  
  p_ranges_x <- c(ggplot_build(p[[1]])$layout$panel_scales_x[[1]]$range$range,
                  ggplot_build(p[[2]])$layout$panel_scales_x[[1]]$range$range)
  
  p_ranges_y <- c(ggplot_build(p[[1]])$layout$panel_scales_y[[1]]$range$range,
                  ggplot_build(p[[2]])$layout$panel_scales_y[[1]]$range$range)
  
  p & 
    xlim(min(p_ranges_x), max(p_ranges_x)) & 
    ylim(min(p_ranges_y), max(p_ranges_y))
  
dev.off()

# Clean the object
seurat.beforeReclustering <- DietSeurat(selClusters.seurat)

seurat.beforeReclustering
# An object of class Seurat 
# 46974 features across 4194 samples within 5 assays 
# Active assay: RNA (22769 features, 0 variable features)
# 2 layers present: counts, data
# 4 other assays present: ADT, SCT, integrated, prediction.score.celltype

rm(seurat, selClusters.seurat)

## After reclustering
dataDir <- "10_Reclustering_MP_and_MC"
seurat.afterReclustering <- readRDS(paste(dataDir, "selectedClusters.subjectIntegrated.reclustered.rds", sep="/"))

seurat.afterReclustering
# An object of class Seurat 
# 41653 features across 4194 samples within 5 assays 
# Active assay: integrated (3000 features, 3000 variable features)
# 2 layers present: data, scale.data
# 4 other assays present: RNA, ADT, SCT, prediction.score.celltype
# 2 dimensional reductions calculated: pca, umap

## OBSERVATION
# - lost around 5300 features ...
# - sample numbers match

# Can I now merge all the clustering results?
# - do the cellnames match?
all(rownames(seurat.beforeReclustering@meta.data) == rownames(seurat.afterReclustering@meta.data))
# [1] TRUE   - YES, they match!

# Get the cluster results BEFORE reclustering
clustResults.beforeReclustering <- seurat.beforeReclustering@meta.data[ , grep("integrated_snn",colnames(seurat.beforeReclustering@meta.data))]

# We used the clusters @0.6 to recluster
# - remove cluster @1
# - rename columns
clustResults.beforeReclustering$integrated_snn_res.1 <- NULL
colnames(clustResults.beforeReclustering) <- gsub("integrated_snn","before", colnames(clustResults.beforeReclustering))

# Add the cluster classification AFTER reclustering
clustResults.afterReclustering <- seurat.afterReclustering@meta.data[ , grep("integrated_snn",colnames(seurat.afterReclustering@meta.data))]

# Rename
colnames(clustResults.afterReclustering) <- gsub("integrated_snn","after", colnames(clustResults.afterReclustering))

# We only used the cluster @res.1 after reclustering
clustResults <- clustResults.beforeReclustering
clustResults$after_res.1 <- clustResults.afterReclustering$after_res.1

#### Flow diagram using alluvial ####
# Column 'Freq' to store the total count of all cluster combinations
clustResults$Freq = 1
clustResults2D = aggregate(Freq ~ before_res.0.3 + before_res.0.6 + after_res.1, data = clustResults, sum)

cols = rev(rainbow(nrow(clustResults2D), start = 0.1, end = 0.9))
pdf(paste0(figDir, "/Seurat_reclustering_flowdiagram.pdf"), height=14, width=28)
alluvial(
  select(clustResults2D, before_res.0.3, before_res.0.6, after_res.1),
  freq = clustResults2D$Freq,
  col = cols,
  alpha = 0.8,
  gap.width = 0.3,
  cw = 0.05,
  blocks = FALSE,
  axis_labels = colnames(clustResults2D)[1:3],
  cex = 3,
  cex.axis = 3)  
dev.off()


#### Clustering tree ####
# First put the cluster information into the seurat object
colnames(clustResults) <- paste0("cluster_", gsub("before|after", "", colnames(clustResults)))
clustResults$cluster_Freq <- NULL

pdf(paste0(figDir, "/Seurat_reclustering_tree.pdf"), height=14, width=14)
  p <- clustree(clustResults, prefix = "cluster_res.", node_text_size = 6, return = "plot",
                edge_width = 2.5, node_size_range = c(10,24))
  p <- p + theme(legend.text = element_text(size=16), legend.key.size = unit(0.6, "cm"),
                 legend.title = element_text(size=20)) +
  guides(color = guide_legend(override.aes = list(size=5)))
  print(p)
dev.off()

#### Table before vs after reclustering ####
table(clustResults$cluster_res.0.6, clustResults$cluster_res.1)

#      0   1   2   3   4   5   6   7   8   9  10  11  12
# 5   40 695 207   4  22  23   6 164  54  49   4  12   1
# 7    0   2  11 366 291 157   1  70 130  35   2   0  56
# 9  713  75  56   2  49   5   0   3   3  31   0  46   1
# 11  46  11 220   9   2 119 266   9   0  47  76   3   0

write.table(table(clustResults$cluster_res.0.6, clustResults$cluster_res.1), 
            file = paste0(figDir, "/Clusters_before2after.tab"), sep="\t", col.names = TRUE, row.names = TRUE, quote=FALSE)


#### CytoSig - Cytokine signatures ####
# "In Figure 2d, the authors present a dot plot showing the top 15 chemokine markers for each macrophage cluster. 
# The results show that the classical MDMs express higher levels of CCR1 and CCR2, as well as mediators for neutrophil recruitment. 
# To further strengthen this argument, I would also suggest running CytoSig (Jiang et al. 2021), which is a package that outputs 
# predicted cytokine signaling cascades that are active in cell states. It contains a curated set of pathway information and 
# uses both the mRNA levels of cytokines and receptors, as well as the downstream genes and other associated signaling partners. 
# This would make the analysis more robust and provide another level of support for this hypothesis."

library(Seurat)
library(scaper)
library(pheatmap)

dataDir <- "10_Reclustering_MP_and_MC"
seurat <- readRDS(paste(dataDir, "selectedClusters.subjectIntegrated.reclustered.rds", sep="/"))
DefaultAssay(seurat) <-"SCT"

seurat
# An object of class Seurat 
# 41653 features across 4194 samples within 5 assays 
# Active assay: SCT (15735 features, 0 variable features)
# 2 layers present: counts, data
# 4 other assays present: ADT, SCT, integrated, prediction.score.celltype
# 2 dimensional reductions calculated: pca, umap

# Why 41653 features:
#   RNA has 22769 features
#   SCT     15735                -> 38504 features
#   ADT       137                -> 37641 features
#   integrated 3000              -> 40641 features
#   prediction.score.celltype 12 -> 41653 features

# Make a new Ident
seurat@meta.data$CellTypes_res1 <- seurat@meta.data$integrated_snn_res.1
Idents(seurat) <- "CellTypes_res1"
levels(seurat)

# And rename
seurat <- RenameIdents(seurat, "0"="0: Non-classical MDMs", "1"="1: Classical MDMs", "2"="2: T cells",
                       "3"="3: cDCs","4"="4: Res-like C1q Mac","5"="5: T cells","6"="6: T & NK cells", "7"="7: cDCs & T cells",
                       "8"="8: SPP1 LAMs","9"="9: Overlap Mac", "10"="10: T cells","11"="11: Stress response non-classical MDMs",
                       "12"="12: cDCs")

supportedCytokines(database = "cytosig")
# [1] "ActivinA" "BDNF"     "BMP2"     "BMP4"     "BMP6"     "CD40LG"   "CXCL12"   "EGF"      "FGF2"     "GCSF"    
# [11] "GDF11"    "GMCSF"    "HGF"      "IFNG"     "IFNL"     "IL10"     "IL12"     "IL13"     "IL15"     "IL17A"   
# [21] "IL1A"     "IL1B"     "IL2"      "IL21"     "IL22"     "IL27"     "IL3"      "IL4"      "IL6"      "LIF"     
# [31] "LTA"      "MCSF"     "NO"       "OSM"      "TGFB1"    "TGFB3"    "TNFA"     "TNFSF12"  "TRAIL"    "VEGFA"   
# [41] "WNT3A" 

supportedCytokines(database = "reactome")
# [1] "BDNF"    "BMP2"    "BMP4"    "BMP6"    "CD40LG"  "CXCL12"  "EGF"     "FGF2"   
# [9] "HGF"     "IFNG"    "IL10"    "IL13"    "IL15"    "IL17A"   "IL1A"    "IL1B"   
# [17] "IL2"     "IL21"    "IL22"    "IL27"    "IL3"     "IL4"     "IL6"     "LIF"    
# [25] "OSM"     "TGFB1"   "TGFB3"   "TNFSF12" "VEGFA"   "WNT3A"  

# NOTE:
# No CCR1 or CCR2 in the databases??

CytoSig.score.output <- scapeForSeurat(seurat.object = seurat,
                                       database = "cytosig", cytokine = "all", normalize = TRUE)
# Error in vamForCollection(gene.expr = Matrix::t(normalized.counts), gene.set.collection = gene.set.collection,  : 
#   Length of tech.var.prop 10898 does not match the number of genes in the expression matrix 15735
# In addition: Warning message:
# In getTechVarPropForSCT(seurat.data) :
#   Multiple SCTransform models, the first model will be used.

# AJ - 26082024: As it uses the first SCT model it only gets the features from the first MOMA sample in the object:
# i.e. 10898 features from 276 cells...

DefaultAssay(seurat) <- "RNA"
CytoSig.score.output <- scapeForSeurat(seurat.object = seurat,
                                       database = "cytosig", cytokine = "all", normalize = TRUE)
# Did not find vst variance decomposition, setting technical variance proportion to 1
# Computing VAM distances for 41 gene sets, 4194 cells and 22769 genes.
# Min set size: 283, median size: 880
# Did not find vst variance decomposition, setting technical variance proportion to 1
# Computing VAM distances for 41 gene sets, 4194 cells and 22769 genes.
# Min set size: 80, median size: 613
# There were 50 or more warnings (use warnings() to see the first 50)

# This returns a assay with the cell-level cytokine activity scores ...

class(CytoSig.score.output)
# [1] "Seurat"
# attr(,"package")
# [1] "SeuratObject"

GetAssay(object = CytoSig.score.output, assay = "scape")
# Assay data with 41 features for 4194 cells
# First 10 features:
#  ActivinA, BDNF, BMP2, BMP4, BMP6, CD40LG, CXCL12, EGF, FGF2, GCSF

cytosig_mat <- as.data.frame(t(as.matrix(CytoSig.score.output@assays$scape@data)))
dim(cytosig_mat)
#[1] 4194   41

pdf(paste(figDir,"CytoSig_initialAttempt.pdf", sep="/"))
  pheatmap(cytosig_mat, fontsize_row = 4, fontsize_col = 7,
           cluster_rows = FALSE, cluster_cols = FALSE)

  pheatmap(cytosig_mat, fontsize_row = 4, fontsize_col = 7)
dev.off()

## OBSERVATION
# - Hmmmm, pretty big, should I not summarize first per cluster?

# Here, I calculate the median of the cell-level cytokine activity scores, but this might 
# not be the best approach.
# I could fx. also calculate the median score of the genes per cluster, first.
# And then score using CytoSig?

for (cluster in levels(CytoSig.score.output@meta.data$CellTypes_res1 )){
  subSet <- CytoSig.score.output@assays$scape@data[, which(CytoSig.score.output@meta.data$CellTypes_res1  == cluster)]
  df <- as.data.frame(apply(subSet, 1, function(x) median(x)))
  colnames(df) <- cluster
  if (!exists("df_all")){
    df_all <- df
  } else {
    df_all <- cbind(df_all, df)
  }
}

pdf(paste(figDir,"CytoSig_medianClusterActivity.pdf", sep="/"))
  pheatmap(df_all, fontsize_row = 4, fontsize_col = 7)
dev.off()

## Run CytoSig only for clusters of interest:
# Only the main clusters: 0,1,4,8
seurat.subset <- subset(x = seurat, subset = CellTypes_res1 %in% c("0","1","4","8"))
seurat.subset@meta.data$CellTypes_res1 <- factor(seurat.subset@meta.data$CellTypes_res1, levels = unique(seurat.subset@meta.data$CellTypes_res1))

CytoSig.score.output <- scapeForSeurat(seurat.object = seurat.subset,
                                       database = "cytosig", cytokine = "all", normalize = TRUE)

cytosig_mat <- as.data.frame(t(as.matrix(CytoSig.score.output@assays$scape@data)))
dim(cytosig_mat)
#[1] 2133   41

pdf(paste(figDir,"CytoSig_initialAttempt_selClusters.pdf", sep="/"))
  pheatmap(cytosig_mat, fontsize_row = 4, fontsize_col = 7,
           cluster_rows = FALSE, cluster_cols = FALSE)

  pheatmap(cytosig_mat, fontsize_row = 4, fontsize_col = 7)
dev.off()

for (cluster in levels(CytoSig.score.output@meta.data$CellTypes_res1 )){
  subSet <- CytoSig.score.output@assays$scape@data[, which(CytoSig.score.output@meta.data$CellTypes_res1  == cluster)]
  df <- as.data.frame(apply(subSet, 1, function(x) median(x)))
  colnames(df) <- cluster
  if (!exists("df_all")){
    df_all <- df
  } else {
    df_all <- cbind(df_all, df)
  }
}

pdf(paste(figDir,"CytoSig_medianClusterActivity_selClusters.pdf", sep="/"))
  pheatmap(df_all, fontsize_row = 4, fontsize_col = 7)
dev.off()

#### Doublets in cluster 9? ####
# -	Nalopen cluster 9 doublet removal, plaatje/tabel/cijfers waarmee duidelijk wordt dat daar geen doubletten inzitten

## Read data after reclustering
dataDir <- "10_Reclustering_MP_and_MC"
seurat.afterReclustering <- readRDS(paste(dataDir, "selectedClusters.subjectIntegrated.reclustered.rds", sep="/"))

# There are a couple of columns in the meta data that hold the results of the doublet detection, After integration they have not been merged into one column?
rm(df_d, df_s)  # Make sure to remove previous created objects!

for (i in grep("DF.classification", colnames(seurat.afterReclustering@meta.data))){
  print(colnames(seurat.afterReclustering@meta.data)[i])
  print(table(seurat.afterReclustering@meta.data$orig.ident, seurat.afterReclustering@meta.data[,i]))
  print(table(seurat.afterReclustering@meta.data$integrated_snn_res.1, seurat.afterReclustering@meta.data[,i]))
  
  a <- as.data.frame(table(seurat.afterReclustering@meta.data[,i], seurat.afterReclustering@meta.data$integrated_snn_res.1))
  # First column, 'Var1', contains state (Doublet/Singlet), second column, 'Var2', the cluster number. The third column contains the frequency
  if (!(exists("df_d")) ){ 
    if ( "Doublet" %in% a$Var1){
      df_d <- a
    } 
  }
  if (!(exists("df_s")) ){
    if ( !("Doublet" %in% a$Var1)){
      df_s <- a
    } 
  }
  if (exists("df_d") ){ 
    if ( "Doublet" %in% a$Var1){
      df_d <- merge(df_d, a, by.x=c('Var1', 'Var2'), by.y=c('Var1', 'Var2'), all = TRUE)
      # Replace all NA's by zero
      df_d[is.na(df_d)] <- 0
      # And sum all 'Freq' columns
      df_d$Total <- rowSums(df_d[,3:ncol(df_d)])
      df_d$Freq.x <- NULL
      df_d$Freq.y <- NULL
    }
  if (exists("df_s") ){
    if ( !("Doublet" %in% a$Var1)){
      df_s <- merge(df_s, a, by.x=c('Var1', 'Var2'), by.y=c('Var1', 'Var2'), all = TRUE)
      # Replace all NA's by zero
      df_s[is.na(df_s)] <- 0
      # And sum all 'Freq' columns
      df_s$Total <- rowSums(df_s[,3:ncol(df_s)])
      df_s$Freq.x <- NULL
      df_s$Freq.y <- NULL
    }
    }
  }
}

colnames(df_s) <- c("State", "Cluster", "Total")
colnames(df_d) <- c("State", "Cluster", "Total")
# Hmm, numbers do not match with the original object, probably I should only take certain annotations into account??
sum(df_s[,3])
# [1] 7726

# I could also do it the other way around ...
# i.e. look which cells belong to cluster 9 and then trace back in the different individual MOMA samples how they were assigned.
Idents(seurat.afterReclustering) <- "seurat_clusters"
cellCluster9 <- unlist(CellsByIdentities(seurat.afterReclustering, idents = "9", cells = NULL, return.null = FALSE))

# Sample list
sampleList <- c("MOMA52","MOMA57","MOMA17" ,"MOMA302","MOMA67" ,"MOMA68" ,"MOMA72")
doubletDir <- "3_doublet_detection"

for (i in 1:length(sampleList)){
  sc <- sampleList[i]
  print(paste0(i, " - ", sc))
  seurat <- readRDS(paste0(doubletDir,"/", sc, ".Mito10.doubletsRemoved.rds"))
  cellsToUse <- gsub(paste0("_",i), "", cellCluster9[which(grepl(paste0("_",i), cellCluster9))])
  # find columns holding doublet/singlet information
  idx <- grep("DF.classification", colnames(seurat@meta.data))
  print(seurat@meta.data[cellsToUse, c("orig.ident", colnames(seurat@meta.data)[idx])])
}

# OBSERVATIONS #
# - Strange, that some cells are not found back in the original dataset... makes you wonder...
# - Nevertheless, it shows that cluster 9 does not contain only doublets!!

# [1] "1 - MOMA52"
# orig.ident DF.classifications_0.25_0.09_458 DF.classifications_0.25_0.09_401
# NA                       <NA>                             <NA>                             <NA>
# NA.1                     <NA>                             <NA>                             <NA>
# ACTGAACGTTGGAGGT-1     MOMA52                          Singlet                          Singlet
# ACTGTCCCACAACGTT-1     MOMA52                          Singlet                          Singlet
# [1] "2 - MOMA57"
# orig.ident DF.classifications_0.25_0.09_688 DF.classifications_0.25_0.09_609
# NA                       <NA>                             <NA>                             <NA>
# AAGTCTGGTTGCGCAC-1     MOMA57                          Singlet                          Singlet
# ACCTTTATCGCGTTTC-1     MOMA57                          Singlet                          Singlet
# CATGGCGAGAGTAATC-1     MOMA57                          Singlet                          Singlet
# CGATGTATCAGATAAG-1     MOMA57                          Singlet                          Singlet
# CTCGTCAAGCTTTGGT-1     MOMA57                          Singlet                          Singlet
# GGCGTGTTCACCAGGC-1     MOMA57                          Singlet                          Singlet
# GGGAGATCAGTAAGAT-1     MOMA57                          Doublet                          Singlet
# NA.1                     <NA>                             <NA>                             <NA>
# TCGCGTTGTTAAGGGC-1     MOMA57                          Singlet                          Singlet
# TTTGTCAGTGACTACT-1     MOMA57                          Singlet                          Singlet
# [1] "3 - MOMA17"
# orig.ident DF.classifications_0.25_0.09_98 DF.classifications_0.25_0.09_81
# ACCATTTTCGCGCTGA-1     MOMA17                         Doublet                         Singlet
# CTCCATGGTCAGGTGA-1     MOMA17                         Singlet                         Singlet
# GATGATCCAAGACTGG-1     MOMA17                         Singlet                         Singlet
# TCGGGACAGGTAGTCG-1     MOMA17                         Singlet                         Singlet
# TTGGGATAGTGGGAAA-1     MOMA17                         Singlet                         Singlet
# [1] "4 - MOMA302"
# orig.ident DF.classifications_0.25_0.09_222 DF.classifications_0.25_0.09_176
# AAGCATCCAGCTCTGG-1    MOMA302                          Singlet                          Singlet
# ACACAGTCAGCAGGAT-1    MOMA302                          Singlet                          Singlet
# AGAAATGAGATTAGAC-1    MOMA302                          Singlet                          Singlet
# AGGGAGTAGCGAGGAG-1    MOMA302                          Singlet                          Singlet
# AGTTAGCCACTAGGCC-1    MOMA302                          Singlet                          Singlet
# ATAGACCGTCAGACTT-1    MOMA302                          Singlet                          Singlet
# ATCACAGGTACGGATG-1    MOMA302                          Singlet                          Singlet
# ATGGGTTTCTCGTTTA-1    MOMA302                          Singlet                          Singlet
# CAAGCTATCGTCCATC-1    MOMA302                          Singlet                          Singlet
# CATAGACGTCAGTCGC-1    MOMA302                          Singlet                          Singlet
# GATGATCGTTGGGTTT-1    MOMA302                          Singlet                          Singlet
# GCAGGCTAGCAACTTC-1    MOMA302                          Singlet                          Singlet
# GGCGTCAGTCGAGCAA-1    MOMA302                          Singlet                          Singlet
# GGTAGAGAGTAAACGT-1    MOMA302                          Singlet                          Singlet
# GTATTGGTCCACAAGT-1    MOMA302                          Singlet                          Singlet
# GTTTACTGTCGAGCAA-1    MOMA302                          Singlet                          Singlet
# TACAACGTCTGCTGAA-1    MOMA302                          Singlet                          Singlet
# TACGGTAGTTGAAGTA-1    MOMA302                          Singlet                          Singlet
# TACTTCAAGTAGTCAA-1    MOMA302                          Singlet                          Singlet
# TAGTGCAAGCAAGTCG-1    MOMA302                          Singlet                          Singlet
# TCTATCAAGGTGCGAT-1    MOMA302                          Singlet                          Singlet
# TGAACGTCACGGCTAC-1    MOMA302                          Singlet                          Singlet
# TGAATCGCAGGCATGA-1    MOMA302                          Singlet                          Singlet
# [1] "5 - MOMA67"
# orig.ident DF.classifications_0.25_0.09_373 DF.classifications_0.25_0.09_310
# AAGGTAACATAATCCG-1     MOMA67                          Singlet                          Singlet
# AGATCGTCAGGTGACA-1     MOMA67                          Singlet                          Singlet
# AGGTAGGTCTACCCAC-1     MOMA67                          Singlet                          Singlet
# ATACTTCGTCTCTCCA-1     MOMA67                          Singlet                          Singlet
# ATGAGTCCAAGTGGGT-1     MOMA67                          Singlet                          Singlet
# ATGCGATGTATGAGCG-1     MOMA67                          Singlet                          Singlet
# CACCAAATCTAGTTCT-1     MOMA67                          Singlet                          Singlet
# CAGGGCTCAACCGGAA-1     MOMA67                          Doublet                          Singlet
# CGCCAGACACAATCTG-1     MOMA67                          Singlet                          Singlet
# CGTGATACACCTAAAC-1     MOMA67                          Singlet                          Singlet
# GAATCACAGTGTTCAC-1     MOMA67                          Singlet                          Singlet
# GAGTGTTCAATCTGCA-1     MOMA67                          Singlet                          Singlet
# GGGAGTACAAGTTCCA-1     MOMA67                          Singlet                          Singlet
# GTATTTCGTTTAGTCG-1     MOMA67                          Singlet                          Singlet
# GTGATGTAGCGCAATG-1     MOMA67                          Singlet                          Singlet
# TACTTGTCAGGTAGTG-1     MOMA67                          Singlet                          Singlet
# TCATCCGCAGACAAGC-1     MOMA67                          Singlet                          Singlet
# TCCCATGTCCCGATCT-1     MOMA67                          Singlet                          Singlet
# TCCGGGATCATTACTC-1     MOMA67                          Singlet                          Singlet
# TCGACGGTCACCGGGT-1     MOMA67                          Singlet                          Singlet
# TGAGCATTCCAAATGC-1     MOMA67                          Singlet                          Singlet
# TGAGGAGGTCACGTGC-1     MOMA67                          Singlet                          Singlet
# TGGATGTTCGGTGAAG-1     MOMA67                          Singlet                          Singlet
# TGGTTAGCAGGTGTGA-1     MOMA67                          Singlet                          Singlet
# TTAGGGTCACGACAGA-1     MOMA67                          Singlet                          Singlet
# TTTGATCAGGGCGAAG-1     MOMA67                          Doublet                          Singlet
# [1] "6 - MOMA68"
# orig.ident DF.classifications_0.25_0.09_206 DF.classifications_0.25_0.09_172
# AAAGGGCCAAATCGGG-1     MOMA68                          Singlet                          Singlet
# AAAGGTATCAGAGCGA-1     MOMA68                          Singlet                          Singlet
# AACGGGAGTCGAGCTC-1     MOMA68                          Singlet                          Singlet
# ACATCGAAGCTGAAGC-1     MOMA68                          Singlet                          Singlet
# ACTACGAAGACCGCCT-1     MOMA68                          Singlet                          Singlet
# ACTGTGATCCACCTCA-1     MOMA68                          Singlet                          Singlet
# ACTTATCCAGGGACTA-1     MOMA68                          Doublet                          Singlet
# ACTTTGTGTGTCTCCT-1     MOMA68                          Singlet                          Singlet
# AGAGAGCTCTAGTGAC-1     MOMA68                          Singlet                          Singlet
# AGCGCTGGTGTGTCGC-1     MOMA68                          Singlet                          Singlet
# AGGCCACCATTGCTGA-1     MOMA68                          Singlet                          Singlet
# AGGTTGTCAAGAAATC-1     MOMA68                          Singlet                          Singlet
# ATCCACCGTACTCGTA-1     MOMA68                          Singlet                          Singlet
# ATCCTATAGCGCTGAA-1     MOMA68                          Singlet                          Singlet
# ATGGGAGTCAAGTTGC-1     MOMA68                          Singlet                          Singlet
# ATTACCTCAACCGACC-1     MOMA68                          Singlet                          Singlet
# ATTCGTTGTCCACAGC-1     MOMA68                          Singlet                          Singlet
# ATTCTTGTCCACGGAC-1     MOMA68                          Singlet                          Singlet
# CAACCTCTCGCCGATG-1     MOMA68                          Singlet                          Singlet
# CAATTTCGTTTCGTTT-1     MOMA68                          Singlet                          Singlet
# CACATGATCGTTAGTG-1     MOMA68                          Singlet                          Singlet
# CACGTGGCACCAGGTC-1     MOMA68                          Singlet                          Singlet
# CAGATTGGTGAGGCAT-1     MOMA68                          Singlet                          Singlet
# CAGCAATAGTCCCGGT-1     MOMA68                          Singlet                          Singlet
# CAGTTCCTCCTTGGAA-1     MOMA68                          Singlet                          Singlet
# CCACCATTCCGTGTGG-1     MOMA68                          Singlet                          Singlet
# CCATCACTCTTTGGAG-1     MOMA68                          Singlet                          Singlet
# CCGAACGTCGCCCAGA-1     MOMA68                          Singlet                          Singlet
# CCTCCAAGTGAGTAAT-1     MOMA68                          Singlet                          Singlet
# CCTTGTGAGCGGATCA-1     MOMA68                          Singlet                          Singlet
# CCTTTGGTCCCATTTA-1     MOMA68                          Singlet                          Singlet
# CGATGCGGTCATGGCC-1     MOMA68                          Singlet                          Singlet
# CGGGTCATCAACGCTA-1     MOMA68                          Singlet                          Singlet
# CGTGCTTTCACTTGTT-1     MOMA68                          Singlet                          Singlet
# CGTTCTGTCGCACGGT-1     MOMA68                          Singlet                          Singlet
# CTACCCAAGGCGCTTC-1     MOMA68                          Singlet                          Singlet
# CTATCCGCACGGTCTG-1     MOMA68                          Singlet                          Singlet
# CTCAAGAAGATAACAC-1     MOMA68                          Singlet                          Singlet
# CTCAAGACATGGCCAC-1     MOMA68                          Singlet                          Singlet
# CTCAATTCAGCTTTGA-1     MOMA68                          Singlet                          Singlet
# CTCAATTGTCGTGCCA-1     MOMA68                          Singlet                          Singlet
# CTCCTCCCATGCCGCA-1     MOMA68                          Singlet                          Singlet
# CTGAATGTCAGTGATC-1     MOMA68                          Singlet                          Singlet
# CTGAGCGAGTCTGGTT-1     MOMA68                          Singlet                          Singlet
# CTGTACCCAGAGCTAG-1     MOMA68                          Singlet                          Singlet
# CTGTATTAGATCCCGC-1     MOMA68                          Singlet                          Singlet
# CTGTGAAAGCATGTTC-1     MOMA68                          Singlet                          Singlet
# GAAGTAAGTGTCATGT-1     MOMA68                          Singlet                          Singlet
# GAGTCATTCTACTGCC-1     MOMA68                          Singlet                          Singlet
# GAGTTGTTCACTACTT-1     MOMA68                          Singlet                          Singlet
# GATAGCTTCATTGCGA-1     MOMA68                          Singlet                          Singlet
# GATGCTACATAGCTGT-1     MOMA68                          Singlet                          Singlet
# GATGGAGCACGCTTAA-1     MOMA68                          Singlet                          Singlet
# GATGTTGAGGGTGAGG-1     MOMA68                          Singlet                          Singlet
# GATGTTGCAACTGCTA-1     MOMA68                          Doublet                          Singlet
# GCAGCCAGTTGTGCCG-1     MOMA68                          Singlet                          Singlet
# GCATCTCGTGCAGTGA-1     MOMA68                          Singlet                          Singlet
# GCCGATGAGGTAGTCG-1     MOMA68                          Singlet                          Singlet
# GGAAGTGAGCGCACAA-1     MOMA68                          Singlet                          Singlet
# GGAGAACAGCATCGAG-1     MOMA68                          Singlet                          Singlet
# GGGCTACCATGTTCGA-1     MOMA68                          Singlet                          Singlet
# GGGTCTGCATAACTCG-1     MOMA68                          Singlet                          Singlet
# GGTCACGGTGTGTACT-1     MOMA68                          Singlet                          Singlet
# GGTGTCGGTCATATGC-1     MOMA68                          Singlet                          Singlet
# GGTTGTATCACAAGAA-1     MOMA68                          Singlet                          Singlet
# GTGCAGCCAGGAGACT-1     MOMA68                          Singlet                          Singlet
# GTGGAAGCAAGATGGC-1     MOMA68                          Singlet                          Singlet
# GTTCATTGTGTCTCCT-1     MOMA68                          Singlet                          Singlet
# TACAACGAGACCTTTG-1     MOMA68                          Singlet                          Singlet
# TATATCCTCATTCACT-1     MOMA68                          Singlet                          Singlet
# TATCTGTTCTACTGCC-1     MOMA68                          Singlet                          Singlet
# TCATTTGCACTGCACG-1     MOMA68                          Singlet                          Singlet
# TCCTTTCGTTTAGTCG-1     MOMA68                          Singlet                          Singlet
# TCGAACACAGCAGTAG-1     MOMA68                          Singlet                          Singlet
# TCGACGGCACTCCACT-1     MOMA68                          Singlet                          Singlet
# TCGCTTGTCCAACCGG-1     MOMA68                          Singlet                          Singlet
# TCGGATAGTATTCCTT-1     MOMA68                          Singlet                          Singlet
# TCTATCAAGTGGCCTC-1     MOMA68                          Singlet                          Singlet
# TCTTGCGGTTTGGAGG-1     MOMA68                          Singlet                          Singlet
# TGAGACTCACAGCGCT-1     MOMA68                          Singlet                          Singlet
# TGCCGAGGTGGTTCTA-1     MOMA68                          Singlet                          Singlet
# TGCTGAAAGAAACCAT-1     MOMA68                          Singlet                          Singlet
# TGCTGAAGTCCTTTGC-1     MOMA68                          Singlet                          Singlet
# TGCTTCGCAAATAAGC-1     MOMA68                          Singlet                          Singlet
# TGTGGCGGTAATACCC-1     MOMA68                          Singlet                          Singlet
# TTCCAATAGGCTCTCG-1     MOMA68                          Singlet                          Singlet
# TTCTTGAGTGAATTAG-1     MOMA68                          Singlet                          Singlet
# TTGTTCATCGATACTG-1     MOMA68                          Singlet                          Singlet
# TTGTTTGCACAAATGA-1     MOMA68                          Singlet                          Singlet
# TTTATGCCATATTCGG-1     MOMA68                          Singlet                          Singlet
# TTTCACAAGTCTCGTA-1     MOMA68                          Singlet                          Singlet
# TTTCAGTGTCGAACAG-1     MOMA68                          Singlet                          Singlet
# [1] "7 - MOMA72"
# orig.ident DF.classifications_0.25_0.09_99 DF.classifications_0.25_0.09_75
# AATAGAGGTAATGATG-1     MOMA72                         Singlet                         Singlet
# ACGGTTAAGCTAGAGC-1     MOMA72                         Doublet                         Singlet
# ATCACAGAGCCAGAGT-1     MOMA72                         Singlet                         Singlet
# ATGGAGGGTCCTGTCT-1     MOMA72                         Singlet                         Singlet
# ATGGATCGTCTTTATC-1     MOMA72                         Singlet                         Singlet
# CCGGTGAGTACTAAGA-1     MOMA72                         Doublet                         Singlet
# CGACAGCAGAGTCAAT-1     MOMA72                         Singlet                         Singlet
# GACAGCCAGGCCTAGA-1     MOMA72                         Singlet                         Singlet
# GGATGTTCAGGTGGAT-1     MOMA72                         Singlet                         Singlet
# TCACTCGAGTACTCGT-1     MOMA72                         Singlet                         Singlet
# TCAGGGCTCGGCCTTT-1     MOMA72                         Singlet                         Singlet
# TCGCTTGAGCACCCAC-1     MOMA72                         Doublet                         Singlet
# TCGGGTGCACGCCACA-1     MOMA72                         Singlet                         Singlet
# TGGATCATCCGATAAC-1     MOMA72                         Singlet                         Singlet
# TTTCAGTGTTCTAAGC-1     MOMA72                         Doublet                         Singlet
# TTTGGAGGTATTTCCT-1     MOMA72                         Doublet                         Singlet

## NOTE
# Something goes wrong/ has gone wrong as the DF columns in the seurat.afterReclustering do not have the same name
# as the DF columns in the individual samples.
# That means that I can not use the names from the DF columns in the individual samples to extract the table for the Singlet/doublet
# assignment from the object after reclustering...
# It seems I have recalculated the DoubletFinding after I made the integrated sample ... but only for the first two samples, i.e. MOMA52 and MOMA57
# Probably after removal of TB cells???

sampleList <- c("MOMA52","MOMA57","MOMA17" ,"MOMA302","MOMA67" ,"MOMA68" ,"MOMA72")
doubletDir <- "3_doublet_detection"

DF_colReclustered <- colnames(seurat.afterReclustering@meta.data[grep("DF", colnames(seurat.afterReclustering@meta.data))])

DF_colSamples <- c()
for (i in 1:length(sampleList)){
  sc <- sampleList[i]
  print(paste0(i, " - ", sc))
  seurat <- readRDS(paste0(doubletDir,"/", sc, ".Mito10.doubletsRemoved.rds"))
  
  for (ix in grep("DF.classification", colnames(seurat@meta.data))){
    DF_colSamples <- c(DF_colSamples,colnames(seurat@meta.data)[ix])
  }
}

all(DF_colReclustered == DF_colSamples)

for (i in 1:length(sampleList)){
  sc <- sampleList[i]
  print(paste0(i, " - ", sc))
  seurat <- readRDS(paste0(doubletDir,"/", sc, ".Mito10.doubletsRemoved.rds"))
  
  for (ix in grep("DF.classification", colnames(seurat@meta.data))){
    colToUse <- colnames(seurat@meta.data)[ix]
    # Dirty hack, but only the first 4 columns deviate ...
    if (colnames(seurat@meta.data)[ix] == "DF.classifications_0.25_0.09_458"){ colToUse <- "DF.classifications_0.25_0.09_454"}
    if (colnames(seurat@meta.data)[ix] == "DF.classifications_0.25_0.09_401"){ colToUse <- "DF.classifications_0.25_0.09_402"}
    if (colnames(seurat@meta.data)[ix] == "DF.classifications_0.25_0.09_688"){ colToUse <- "DF.classifications_0.25_0.09_684"}
    if (colnames(seurat@meta.data)[ix] == "DF.classifications_0.25_0.09_609"){ colToUse <- "DF.classifications_0.25_0.09_605"}
    
    if ("Doublet" %in% unique(seurat.afterReclustering@meta.data[,colToUse])){
      print(colnames(seurat@meta.data)[ix])
      print(table(seurat.afterReclustering@meta.data$integrated_snn_res.1, seurat.afterReclustering@meta.data[,colToUse]))
    }
    a <- as.data.frame(table(seurat.afterReclustering@meta.data[,colToUse], seurat.afterReclustering@meta.data$integrated_snn_res.1))
  }
}

# Just plot the number of counts per cluster in a Violin plot
Idents(seurat.afterReclustering) <- "integrated_snn_res.1"

pdf(paste(figDir,"violinPlot_RNA_clusters_res1_nCount_RNA.pdf", sep="/"))

  p <- VlnPlot(object = seurat.afterReclustering, features=c("nCount_RNA"), pt.size = 1) + NoLegend() +
    theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 18), axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14))
  print(p)
  
  # And save as TIFF
  ggsave(paste(figDir, "violinPlot_RNA_clusters_res1_nCount_RNA.tiff", sep="/"), p, device='tiff',scale=2, width=5, height=5, dpi=360)

dev.off()

pdf(paste(figDir,"violinPlot_RNA_clusters_res1_nFeature_RNA.pdf", sep="/"))

  p <- VlnPlot(object = seurat.afterReclustering, features=c("nFeature_RNA"), pt.size = 1) + NoLegend() +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 18), axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14))
  print(p)

  # And save as TIFF
  ggsave(paste(figDir, "violinPlot_RNA_clusters_res1_nFeature_RNA.tiff", sep="/"), p, device='tiff', scale=2, width=5, height=5, dpi=360)

dev.off()

  
#### Myeloid cells and T cells ####
# -	Vergelijken clustering myeloid-enriched cells waar t cellen inzitten met subanalyse van alleen myeloide cellen.  
#   Door ons geannoteerde myeloide cellen pakken en die opnieuw clusteren -> laten zien dat annotatie niet af hangt van T cellen.
#   Bestaande integrated assay, zonder nieuwe annotatie, overlap van clusters?
dataDir <- "10_Reclustering_MP_and_MC"
seurat <- readRDS(paste(dataDir, "selectedClusters.subjectIntegrated.reclustered.rds", sep="/"))
DefaultAssay(seurat) <-"SCT"

# Make a new Ident
seurat@meta.data$CellTypes_res1 <- seurat@meta.data$integrated_snn_res.1
Idents(seurat) <- "CellTypes_res1"
levels(seurat)

# And rename
seurat <- RenameIdents(seurat, "0"="0: Non-classical MDMs", "1"="1: Classical MDMs", "2"="2: T cells",
                       "3"="3: cDCs","4"="4: Res-like C1q Mac","5"="5: T cells","6"="6: T & NK cells", "7"="7: cDCs & T cells",
                       "8"="8: SPP1 LAMs","9"="9: Overlap Mac", "10"="10: T cells","11"="11: Stress response non-classical MDMs",
                       "12"="12: cDCs")

# Take only myeloid cells and recluster?
DefaultAssay(seurat) <- "RNA"
selClusters.seurat <- subset(x = seurat, subset = integrated_snn_res.1 %in% c("0","1","4","8","9","11"))
dim(selClusters.seurat)
# [1] 22769  2356

table(selClusters.seurat$orig.ident)
#MOMA17 MOMA302  MOMA52  MOMA57  MOMA67  MOMA68  MOMA72 
#    76      96     197     870     167     766     184

originalClustering <- seurat@meta.data[, c(colnames(seurat@meta.data)[grep("integrated_snn", colnames(seurat@meta.data))], "CellTypes_res1")]

# Find neighbors
DefaultAssay(selClusters.seurat) <-"integrated"
seurat.clustered <- FindNeighbors(selClusters.seurat, reduction = "pca", dims = 1:30)

# Cluster at different resolution 
seurat.clustered <- FindClusters(seurat.clustered, resolution = 0.3)
# Number of communities: 5
seurat.clustered <- FindClusters(seurat.clustered, resolution = 0.6)
# Number of communities: 6
seurat.clustered <- FindClusters(seurat.clustered, resolution = 1.0)
# Number of communities: 8
seurat.clustered <- FindClusters(seurat.clustered, resolution = 1.4)
# Number of communities: 12

newClustering <- seurat.clustered@meta.data[, c(colnames(seurat.clustered@meta.data)[grep("integrated_snn", colnames(seurat.clustered@meta.data))])]

clustResults <- cbind(integrated_snn_res.1_old = originalClustering[rownames(newClustering),"CellTypes_res1"], newClustering)

clustResults$Freq = 1
clustResults2D = aggregate(Freq ~ integrated_snn_res.1_old + integrated_snn_res.0.3 + integrated_snn_res.0.6 + integrated_snn_res.1 + integrated_snn_res.1.4, data = clustResults, sum)

cols = rev(rainbow(nrow(clustResults2D), start = 0.1, end = 0.9))
pdf(paste0(figDir,"/Seurat_comparison_flowdiagram.pdf"), height=14, width=28)
alluvial(
  select(clustResults2D, integrated_snn_res.1_old, integrated_snn_res.0.3, integrated_snn_res.0.6, integrated_snn_res.1, integrated_snn_res.1.4),
  freq = clustResults2D$Freq,
  col = cols,
  alpha = 0.8,
  gap.width = 0.3,
  cw = 0.05,
  blocks = FALSE,
  axis_labels = sub("integrated_snn_res.","",colnames(clustResults2D)[1:5]),
  cex = 1.5,
  cex.axis = 1.5)  
dev.off()

# OBSERVATIONS:
# - Hardly any cells are clustered differently ...

## Clustering tree ##
# Clustree can not handle "06_old"
seurat.clustered@meta.data$integrated_snn_res.0.29 <- clustResults$integrated_snn_res.1_old

# Put the old clustering in front of the new ones
# And only map the original clustering (with 6 clusters) on the new clustering level with 6 clusters
indx1 <- which(grepl("integrated_snn_res.0.29", colnames(seurat.clustered@meta.data)))
indx2 <- which(grepl("integrated_snn_res.0.6", colnames(seurat.clustered@meta.data)))
indx <- c(indx1, indx2)
seurat.clustered@meta.data[,indx] <- seurat.clustered@meta.data[,c(indx[length(indx)],indx[1:(length(indx)-1)])]
colnames(seurat.clustered@meta.data)[indx] <- colnames(seurat.clustered@meta.data)[c(indx[length(indx)],indx[1:(length(indx)-1)])]

seurat.clustered@meta.data <- seurat.clustered@meta.data[indx]

pdf(paste0(figDir, "/Seurat_comparison_clustering_tree.pdf"), height=14, width=14)
p <- clustree(seurat.clustered, prefix = "integrated_snn_res.", node_text_size = 6, return = "plot",
              edge_width = 2.5, node_size_range = c(10,24))
p <- p + theme(legend.text = element_text(size=16), legend.key.size = unit(0.6, "cm"),
               legend.title = element_text(size=20)) +
  guides(color = guide_legend(override.aes = list(size=5))) +
  scale_color_hue(labels = c("1_old", "0.6")) # Get the right labels
print(p)
dev.off()

# And make a cross table
table(seurat.clustered$integrated_snn_res.0.29, seurat.clustered$integrated_snn_res.0.6)

# write.table(table(seurat.clustered$integrated_snn_res.0.29, seurat.clustered$integrated_snn_res.0.6), 
#             file = paste0(figDir, "/reclusteringWithoutNonMyeloids.tab"), sep="\t", col.names = TRUE, row.names = TRUE, quote=FALSE)

write.table(table(factor(seurat.clustered$integrated_snn_res.0.29, levels = c(1,0,4,8,9,11)), seurat.clustered$integrated_snn_res.0.6), 
                  file = paste0(figDir, "/reclusteringWithoutNonMyeloids.tab"), sep="\t", col.names = TRUE, row.names = TRUE, quote=FALSE)

tab <- table(factor(seurat.clustered$integrated_snn_res.0.29, levels = c(1,0,4,8,9,11)), seurat.clustered$integrated_snn_res.0.6)
diag(tab)/rowSums(tab)
#1         0         4         8         9        11 
#0.9923372 0.9712140 0.9780220 0.9839572 0.9691358 0.9508197 
100*diag(tab)/rowSums(tab)

## OBSERVATION
# - One-to-one correspondence!

#### Markers reclustering without non-myeloid cells ####
# Gleaned the below code from D:/Dropbox/Support/Yosta_Vegting/scRNASeq/Processing/20221220_AllSamples_v4/Code/10_reclustering_MP_and_MC.r
outDir <- paste0(figDir,"/markers_20240627_SCT")
dir.create(outDir)

# 20230221 - changed the directory as I made changes to the code (see mail Perry 20230206)
#          - redid the FindAllMarkers, now using 'latent.vars' and 'test.use = 'LR' '
outDir <- paste0(figDir,"/markers_20240627_LR")
dir.create(outDir)

# Make sure I do not have to change the code ....
seurat <- seurat.clustered

for (outDir in c(paste0(figDir,"/markers_20240627_SCT"),paste0(figDir,"/markers_20240627_LR") )){
  print("")
  print(outDir)
  for (clust in colnames(seurat@meta.data)[grep("integrated_snn", colnames(seurat@meta.data))]){
    res <- gsub("integrated_snn_","", clust)
    print(paste0("Now working on resolution ", res))
    
    # find markers for every cluster compared to all remaining cells, both negative and positive (only.pos = FALSE)
    Idents(seurat) <- clust
    
    if ( outDir == paste0(figDir,"/markers_integrated_20240627")){
      DefaultAssay(seurat) <- "integrated"
    } else {
      # Make sure to use 'SCT' and not 'integrated' as DefaultAssay, see mail Perry 06022023, but see above??
      DefaultAssay(seurat) <- "SCT"
    }
    
    print(paste0("... The following setting is used for DefaultAssay: ", DefaultAssay(seurat)))
    # if you set the DefaultAssay(seurat) to "SCT", you need to run PrepSCTFindMarkers or else you get an error message:
    #   Calculating cluster 0
    #   Calculating cluster 1
    #   Calculating cluster 2
    #   Calculating cluster 3
    #   Calculating cluster 4
    #   Warning: No DE genes identified
    #   Warning: The following tests were not performed: 
    #   Warning: When testing 0 versus all:
    #   Object contains multiple models with unequal library sizes. Run `PrepSCTFindMarkers()` before running `FindMarkers()`.
    
    if (DefaultAssay(seurat) == "SCT"){
      print("... PrepSCTFindMarkers is run")
      seurat <- PrepSCTFindMarkers(seurat) # Found 7 SCT models. Recorrecting SCT counts using minimum median counts: 3354
      # PrepSCTFindMarkers skips the re-correction the second time..
    }
    
    seurat@meta.data[,clust] <- factor(seurat@meta.data[,clust], levels=as.character(0:(length(unique(seurat@meta.data[,clust]))-1)))
    if (!file.exists(paste(outDir,paste0("clusterMarkers_", res,".txt"), sep="/"))){
      if ( outDir == paste0(figDir,"/markers_20240627_LR") ){
        print("... Using FindAllMarkers with latent.vars=\"orig.ident\" and test.use=\"LR\"")
        seurat.markers <- FindAllMarkers(seurat, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25,
                                         latent.vars = "orig.ident", test.use = 'LR')
      } else {
        print("... Using default FindAllMarkers")
        seurat.markers <- FindAllMarkers(seurat, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)
      }
    } else {
      seurat.markers <- read.table(paste(outDir,paste0("clusterMarkers_", res,".txt"), sep="/"), header = TRUE)
    }
    
    
    # Give top 5 per comparison based on p_adj_val (need top_n(n=-5) to get the right 5 genes!!)
    #top5.markers <- seurat.markers %>% group_by(cluster) %>% top_n(n = -5, wt = p_val_adj)  # wt defines ordering, avg_log2FC, or p_val_adj?
    # Hmm, returns more than 5 rows per clusters as there are ties in the 'p_val_adj' ...
    # So, also sort on avg_log2FC ...
    top5.markers <- seurat.markers %>% group_by(cluster) %>% arrange(p_val_adj, desc(avg_log2FC), .by_group=TRUE) %>% top_n(n = 5, wt = avg_log2FC)
    
    # Now add CD163, MRC1, MARCO and TREM2 (if they are not already in there)
    for (markerGene in c("CD163", "MRC1", "MARCO", "TREM2")){
      if ( !(markerGene %in% top5.markers$gene)){
        print(paste0("... Now adding ", markerGene, " to the list of marker genes"))
        top5.markers <- rbind(top5.markers, seurat.markers[seurat.markers$gene == markerGene, ])
      }
    }
    
    # Sort the top5 list on the cluster column
    top5.markers <- dplyr::arrange(top5.markers, cluster)
    
    # # For the DoHeatmap we need the scale.data slot ... do it here, instead of in one go for all genes, to be faster!! :-)
    if (DefaultAssay(seurat) != "SCT"){ # Added 20230208 - AJ: see https://github.com/satijalab/seurat/discussions/4259
      seurat <- ScaleData(seurat, features = top5.markers$gene, 
                          verbose = FALSE, vars.to.regress = c("nCount_RNA", "percent.mito")) 
    }
    
    # And plot in a heatmap and dotplot
    pdf(paste(outDir, paste0("markerPlots_Top5_", res,".pdf"), sep="/"), width = 28, height = 28)
    # Perhaps downsample to 100,200 cells??
    p <- DoHeatmap(seurat, features = top5.markers$gene, group.by = clust) + 
      theme(text = element_text(size=28), axis.text.y = element_text(size=14)) +
      NoLegend()  
    # 1: The following features were omitted as they were not found in the scale.data slot for the SCT assay: FCER1A, LYPD2
    # 2: The following features were omitted as they were not found in the scale.data slot for the SCT assay: FCER1A, LYPD2
    # 3: The following features were omitted as they were not found in the scale.data slot for the SCT assay: CCL17, IGLC2, FN1, LYPD2
    # 4: The following features were omitted as they were not found in the scale.data slot for the SCT assay: CCL17, IGLC2, FN1, FCER1A, LYPD2
    print(p)
    p <-  DotPlot(seurat, features = unique(top5.markers$gene), dot.scale = 6) + 
      coord_flip() +  
      theme(text = element_text(size=28), axis.text.y = element_text(size=14))
    print(p)
    dev.off()
    
    # Write file
    if (!file.exists(paste(outDir,paste0("clusterMarkers_", res,".txt"), sep="/"))){
      write.table(seurat.markers, file = paste(outDir,paste0("clusterMarkers_", res,".txt"), sep="/"), sep="\t", 
                  col.names = TRUE, row.names = FALSE, quote = FALSE)
    }
    
    # In case of resolution 1, only get clusters 0,1,4,8,9,11
    if ( clust == "integrated_snn_res.0.29"){
      selClusters.seurat <- subset(x = seurat, subset = integrated_snn_res.1 %in% c("0","1","4","8","9","11"))
      selClusters.seurat$integrated_snn_res.1 <- factor(selClusters.seurat$integrated_snn_res.1, levels=c("0","1","4","8","9","11"))
      
      pdf(paste(outDir, paste0("markerPlots_Top5_", res,".selectedClusters.pdf"), sep="/"), width = 28, height = 28)
      # Perhaps downsample to 100,200 cells??
      p <- DoHeatmap(selClusters.seurat, features = top5.markers$gene, group.by = clust) + 
        theme(text = element_text(size=28), axis.text.y = element_text(size=14)) +
        NoLegend()  
      print(p)
      p <-  DotPlot(selClusters.seurat, features = unique(top5.markers$gene), dot.scale = 6) + 
        coord_flip() +  
        theme(text = element_text(size=28), axis.text.y = element_text(size=14))
      print(p)
      dev.off()
      
      # And get a Heatmap of the ADTs
      # - We need first to Normalize and Scale the ADT data
      #   see https://github.com/satijalab/seurat/issues/3890, https://github.com/satijalab/seurat/issues/5089
      # - Do we take the ADT separate from the RNA data??
      # - We only have ADTs for two samples
      
      # The ADT has been normalized (see '2_filtering_normalization_and_scaling')
      adt.list <- SplitObject(seurat, split.by = "orig.ident")
      # Only get the two samples with actual ADT data
      adt.list <- adt.list[c("MOMA52", "MOMA57")]
      print("... Finding Varaiable Features in the ADT samples")
      adt.list <- lapply(X = adt.list, FUN = function(x) {
        #x <- NormalizeData(x, assay = "ADT", normalization.method = "CLR", verbose = FALSE)
        x <- FindVariableFeatures(x, assay = "ADT", verbose = FALSE)
      })
      print("... Selecting integration features")
      for ( j in length(adt.list)){
        print(length(VariableFeatures(adt.list[[j]])))
      }
      if ( length(VariableFeatures(adt.list[[1]])) == 0 ){
        features <- rownames(adt.list[[1]]@assays$ADT@counts)
      } else {
        features <- SelectIntegrationFeatures(object.list = adt.list, assay = c("ADT", "ADT"))
      }
      print("... Scaling and running PCA")
      adt.list <- lapply(X = adt.list, FUN = function(x) {
        x <- ScaleData(x, features = features, assay = "ADT", verbose = FALSE)
        x <- RunPCA(x, features = features, assay = "ADT", verbose = FALSE)
      })
      print("... Finding integration anchors (using RCPA) and integrating ")
      anchors <- FindIntegrationAnchors(object.list = adt.list, assay = c("ADT", "ADT"), reference = c(1, 2), 
                                        reduction = "rpca", dims = 1:30)
      adt.integrated <- IntegrateData(anchorset = anchors, dims = 1:30)
      print("... Scaling, running PCA and obtaining a UMAP of the integrated ADT data")
      adt.integrated <- ScaleData(adt.integrated, verbose = FALSE)
      adt.integrated <- RunPCA(adt.integrated, verbose = FALSE)
      adt.integrated <- RunUMAP(adt.integrated, dims = 1:30)
      
      # And add to the RNA object, but now you get an error as you only
      # have ADTs for two samples:
      #  Error: Cannot add a different number of cells than already present
      #seurat[["IADT"]] <- adt.integrated[["integrated"]]
      #seurat[["pca.adt"]] <- adt.integrated[["pca"]]
      #seurat[["umap.adt"]] <- adt.integrated[["umap"]]
      
      if ( outDir == paste0(figDir,"/markers_20240627_LR") ){
        print("... Integrated ADT data: Using FindAllMarkers with latent.vars=\"orig.ident\" and test.use=\"LR\"")
        adt.markers <- FindAllMarkers(adt.integrated, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25,
                                      latent.vars = "orig.ident", test.use = 'LR', assay = "ADT")
      } else {
        print("... Integrated ADT data: Using default FindAllMarkers")
        adt.markers <- FindAllMarkers(adt.integrated, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25 , assay = "ADT")
      }
      top5.adt_markers <- adt.markers %>% group_by(cluster) %>% arrange(p_val_adj, desc(avg_log2FC), .by_group=TRUE) %>% top_n(n = 5, wt = avg_log2FC)
      top5.adt_markers <- dplyr::arrange(top5.adt_markers, cluster)
      
      selClusters.adt <- subset(x = adt.integrated, subset = integrated_snn_res.0.29 %in% c("0","1","4","8","9","11"))
      selClusters.adt$integrated_snn_res.1 <- factor(selClusters.seurat$integrated_snn_res.0.29, levels=c("0","1","4","8","9","11"))
      
      pdf(paste(outDir, paste0("ADT_markerPlots_Top5_", res,".selectedClusters.pdf"), sep="/"), width = 28, height = 28)
      p <- DoHeatmap(selClusters.adt, features = top5.adt_markers$gene, assay = "integrated",group.by = clust) + 
        theme(text = element_text(size=28), axis.text.y = element_text(size=14)) +
        NoLegend()  
      print(p)
      p <-  DotPlot(selClusters.seurat, features = unique(top5.adt_markers$gene), dot.scale = 6, assay = "ADT") + 
        coord_flip() +  
        theme(text = element_text(size=28), axis.text.y = element_text(size=14))
      print(p)
      dev.off()
      
      
    }
    
  }
  
  # ViolinPlots of clusters
  for (clust in colnames(seurat@meta.data)[grep("integrated_snn", colnames(seurat@meta.data))]){
    res <- gsub("integrated_snn_r","R", clust)
    dir.create(paste(outDir,paste0("VlnPlots_", res), sep="/"))
    dir.create(paste(outDir,paste0("FeaturePlots_", res), sep="/"))
    
    Idents(seurat) <- clust
    maxY.RNA <- NULL # see mail Perry 06022023: original code maxY.RNA <- max(seurat@assays$RNA@data)
    # 21022023: Set the axis back in this case !!!
    
    if ( outDir == paste0(figDir,"/markers_20240627_LR") ){
      seurat.markers <- read.table(paste(figDir,"markers_20240627_LR",paste0("clusterMarkers_", res,".txt"), sep="/"), header = TRUE)
    } 
    if ( outDir == paste0(figDir,"/markers_20240627_SCT") ){
      seurat.markers <- read.table(paste(figDir,"markers_20240627_SCT",paste0("clusterMarkers_", res,".txt"), sep="/"), header = TRUE)
    } 
    # if ( outDir == paste0(dataDir,"/markers_integrated") ){
    #   seurat.markers <- read.table(paste("10_Reclustering_MP_and_MC/markers_integrated",paste0("clusterMarkers_", res,".txt"), sep="/"), header = TRUE)
    # } 
    
    #top5.markers <- seurat.markers %>% group_by(cluster) %>% top_n(n = -5, wt = p_val_adj)  # wt defines ordering, avg_log2FC, or p_val_adj?
    # Hmm, returns more than 5 rows per clusters as there are ties in the 'p_val_adj' ...
    # So, also sort on avg_log2FC ...
    top5.markers <- seurat.markers %>% group_by(cluster) %>% arrange(p_val_adj, desc(avg_log2FC), .by_group=TRUE) %>% top_n(n = 5, wt = avg_log2FC)
    
    # Sort the top5 list on the cluster column
    top5.markers <- dplyr::arrange(top5.markers, cluster)
    
    if ( outDir == paste0(figDir,"/markers_integrated")){
      DefaultAssay(seurat) <- "integrated"
    } else {
      # Make sure to use 'SCT' and not 'integrated' as DefaultAssay, see mail Perry 06022023
      DefaultAssay(seurat) <- "SCT"
    }
    
    print(paste0("... VlnPlots: The following setting is used for DefaultAssay: ", DefaultAssay(seurat)))
    
    for ( myClusters in unique(top5.markers$cluster) ){
      df <- top5.markers[top5.markers$cluster == myClusters,]
      pViolin <- list()
      pFeature <- list()
      for ( i in 1:nrow(df)){
        if ( i %in% c(2,4) & !(is.null(maxY.RNA)) ){
          pViolin[[i]] <- VlnPlot(object = seurat, features = df$gene[i], pt.size = 0.1, y.max = maxY.RNA, same.y.lims = TRUE) + NoLegend() +
            theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), 
                  axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) 
          pFeature[[i]] <- FeaturePlot(object = seurat, features = df$gene[i], pt.size = 0.1, reduction = "umap" ) + # see mail Perry 06022023 NoLegend() +
            theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), 
                  axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) 
        } else {
          pViolin[[i]] <- VlnPlot(object = seurat, features = df$gene[i], pt.size = 0.1, y.max = maxY.RNA, same.y.lims = TRUE) + NoLegend() +
            theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) 
          pFeature[[i]] <- FeaturePlot(object = seurat, features = df$gene[i], pt.size = 0.1, reduction = "umap" ) + # see mail Perry 06022023 NoLegend() +
            theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) 
          
        }
      }
      p1 <- pViolin[[1]] | pViolin[[2]] 
      p2 <- pViolin[[3]] | pViolin[[4]]
      p3 <- pViolin[[5]] | plot_spacer()
      p <- p1 / p2 / p3 + plot_layout(heights = c(1, 1, 1))
      pdf(paste(outDir, paste(paste0("VlnPlots_", res),paste0("Cluster_", myClusters,".pdf"), sep="/"), sep="/"))
      print(p)
      dev.off()
      
      p1 <- pFeature[[1]] | pFeature[[2]] 
      p2 <- pFeature[[3]] | pFeature[[4]]
      p3 <- pFeature[[5]] | plot_spacer()
      p <- p1 / p2 / p3 + plot_layout(heights = c(1, 1, 1))
      pdf(paste(outDir, paste(paste0("FeaturePlots_", res),paste0("Cluster_", myClusters,".pdf"), sep="/"), sep="/"))
      print(p)
      dev.off()
      
    }
  }
}

#### Violin plots CD163, S100A9, SPP1, TREM2 (Fig. 1D) ####
# - Violin plots CD163, S100A9, SPP1, TREM2 (e.g. Fig 1d) alleen cluster 0,1,4,8,9,11 (macrofaag clusters) and label 
#   clusters volgens onze naamgeving (kan ik ook in adobe doen).
dataDir <- "10_Reclustering_MP_and_MC/"
seurat <- readRDS(paste(dataDir, "selectedClusters.subjectIntegrated.reclustered.rds", sep="/"))

markerList <- c("CD14", "FCGR3A","CD163", "C1QC", "C1QA", "TREM2", "SPP1", "S100A9")

DefaultAssay(seurat) <- "SCT"

#for (clust in colnames(seurat@meta.data)[grep("integrated_snn", colnames(seurat@meta.data))]){
clust <- "integrated_snn_res.1"

res <- gsub("integrated_snn_","", clust)
seurat@meta.data[,clust] <- factor(seurat@meta.data[,clust], levels=as.character(0:(length(unique(seurat@meta.data[,clust]))-1)))

maxY.RNA <- NULL

# 
pdf(paste(figDir,paste0("violinPlot_RNA_selectedMarkers_", res,".pdf"), sep="/"))
Idents(seurat) <- clust

for ( marker in markerList){
  p <- VlnPlot(object = seurat, features = marker, idents = c(0,1,4,8,9,11), pt.size = 1, y.max = maxY.RNA, same.y.lims = TRUE) + NoLegend() +
    theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 18), axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14))
  print(p)
  
  # And save as TIFF
  ggsave(paste(figDir, paste0("/tiff/violinPlot_RNA_",marker,"_", res,".tiff"), sep="/"), p, scale=2, width=5, height=5, dpi=360,compression="lzw")
  
}

dev.off()

#### Check filtering on housekeeping genes ####
## Read the data ##
sampleList <- c("MOMA17", "MOMA52", "MOMA57", "MOMA67", "MOMA68", "MOMA72", "MOMA302")

tab <- c()

for ( sc in sampleList){
  seurat <- readRDS(paste0("../Data/scRNASeq/Processed/", sc, ".raw.annot.rds"))
  
  # Size before filtering
  origDim <- dim(seurat)
  print(paste0(sc, " - original number of cells: ", origDim[2]))
  # Perform some filtering based on nrGenes, nrExpHouseholdGenes, percMito
  print(paste0(sc, " - range of nFeatures_RNA: ", paste0(range(seurat@meta.data$nFeature_RNA), collapse=" - ")))
  # 
  
  ## Filtering ##
  for ( pctMito in c(10)){
    seuratFeatures <- subset(seurat, subset = nFeature_RNA >= 200 & nFeature_RNA < 6000)
    seuratMito <- subset(seurat, subset = nFeature_RNA >= 200 & nFeature_RNA < 6000 & percent.mito < pctMito)
    seuratMitoHK <- subset(seurat, subset = nFeature_RNA >= 200 & nFeature_RNA < 6000 & percent.mito < pctMito & n.exp.hkgenes > 55 )
    print(paste0("    Filtered on #Features                :", ncol(seuratFeatures), " cells remain"))
    print(paste0("    Filtered on #Features and pctMito    :", ncol(seuratMito), " cells remain"))
    print(paste0("    Filtered on #Features, pctMito and HK:", ncol(seuratMitoHK), " cells remain"))
    print(paste0("  Filtered an additional ",ncol(seuratMito)-ncol(seuratMitoHK), " cells by selecting for more than 55 housekeeping genes being present"))
    pctHK_Filtering <- round(((ncol(seuratMito)-ncol(seuratMitoHK))/ncol(seuratMito)) * 100, 2)
    print(paste0("  which corresponds to ", ((ncol(seuratMito)-ncol(seuratMitoHK))/ncol(seuratMito)) * 100, " %"))
    
    tab <- cbind(tab, c(ncol(seuratFeatures), ncol(seuratMito), ncol(seuratMitoHK), pctHK_Filtering ))
    
    pdf(paste0(figDir, "/", sc, "_nCountRNA_vs_nrHKgenes_beforeFiltering.pdf"))
      plot(seurat@meta.data$n.exp.hkgenes, seurat@meta.data$nCount_RNA, xlab="# housekeeping genes", ylab="nCount_RNA", main=sc, col="dodgerblue3")
      abline(v=55, col="red", lty=2, lwd=2) # Filtering on HK
      #plot(NULL, xlab="# housekeeping genes", ylab="nFeature_RNA", main=sc, xlim=c(0,max(seuratMito@meta.data$n.exp.hkgenes)+5), ylim=c(0, max(seuratMito@meta.data$nFeature_RNA))+25)
      #rect(0,200, max(seuratMito@meta.data$n.exp.hkgenes)+5, 6000, density=NA, col="slategray1")
      #points(seuratMito@meta.data$n.exp.hkgenes, seuratMito@meta.data$nFeature_RNA, col="dodgerblue3")
      plot(seurat@meta.data$n.exp.hkgenes, seurat@meta.data$nFeature_RNA, xlab="# housekeeping genes", ylab="nFeature_RNA", main=sc, col="dodgerblue3")
      abline(v=55, col="red", lty=2, lwd=2)
      #abline(h=200, col="blue", lty=2, lwd=2) # Filtering on nFeature (lower)
      #abline(h=6000, col="blue", lty=2, lwd=2) # Filtering on nFeature (upper)
      plot(seurat@meta.data$n.exp.hkgenes, seurat@meta.data$percent.mito, xlab="# housekeeping genes", ylab="pctMito_RNA", main=sc, col="dodgerblue3")
      abline(v=55, col="red", lty=2, lwd=2) # Filtering on HK
    dev.off()
    
    pdf(paste0(figDir, "/", sc, "_nCountRNA_vs_nrHKgenes_afterFiltering.pdf"))
      plot(seuratMito@meta.data$n.exp.hkgenes, seuratMito@meta.data$nCount_RNA, xlab="# housekeeping genes", ylab="nCount_RNA", main=sc, col="dodgerblue3")
      abline(v=55, col="red", lty=2, lwd=2) # Filtering on HK
      #plot(NULL, xlab="# housekeeping genes", ylab="nFeature_RNA", main=sc, xlim=c(0,max(seuratMito@meta.data$n.exp.hkgenes)+5), ylim=c(0, max(seuratMito@meta.data$nFeature_RNA))+25)
      #rect(0,200, max(seuratMito@meta.data$n.exp.hkgenes)+5, 6000, density=NA, col="slategray1")
      #points(seuratMito@meta.data$n.exp.hkgenes, seuratMito@meta.data$nFeature_RNA, col="dodgerblue3")
      plot(seuratMito@meta.data$n.exp.hkgenes, seuratMito@meta.data$nFeature_RNA, xlab="# housekeeping genes", ylab="nFeature_RNA", main=sc, col="dodgerblue3")
      abline(v=55, col="red", lty=2, lwd=2)
      #abline(h=200, col="blue", lty=2, lwd=2) # Filtering on nFeature (lower)
      #abline(h=6000, col="blue", lty=2, lwd=2) # Filtering on nFeature (upper)
      plot(seuratMito@meta.data$n.exp.hkgenes, seuratMito@meta.data$percent.mito, xlab="# housekeeping genes", ylab="pctMito_RNA", main=sc, col="dodgerblue3")
      abline(v=55, col="red", lty=2, lwd=2) # Filtering on HK
    dev.off()
  }
}

tab <- cbind(c("#Features", "#Features + pctMito", "#Features + pctMito + #HK", "pctHK_Filtered"), tab)
colnames(tab) <- c("Filtered on", sampleList)

rownames(tab) <- NULL

write.table(tab, file = paste0(figDir, "/FilteringSteps.tab"), sep="\t", col.names = TRUE, row.names = FALSE, quote=FALSE)


# [1] "MOMA17 - range of nFeatures_RNA: 24 - 8336"
# [1] "    Filtered on #Features and pctMito:1762 cells remain"
# [1] "    Filtered on #Features, pctMito and HK:1306 cells remain"
# [1] "  Filtered an additional 456 cells by selecting for more than 55 housekeeping genes being present"
# [1] " which corresponds to 25.8796821793417 %"
# [1] "MOMA52 - range of nFeatures_RNA: 28 - 6245"
# [1] "    Filtered on #Features and pctMito:6382 cells remain"
# [1] "    Filtered on #Features, pctMito and HK:6110 cells remain"
# [1] "  Filtered an additional 272 cells by selecting for more than 55 housekeeping genes being present"
# [1] " which corresponds to 4.26198683798182 %"
# [1] "MOMA57 - range of nFeatures_RNA: 23 - 5398"
# [1] "    Filtered on #Features and pctMito:9345 cells remain"
# [1] "    Filtered on #Features, pctMito and HK:9174 cells remain"
# [1] "  Filtered an additional 171 cells by selecting for more than 55 housekeeping genes being present"
# [1] " which corresponds to 1.82985553772071 %"
# [1] "MOMA67 - range of nFeatures_RNA: 28 - 8144"
# [1] "    Filtered on #Features and pctMito:5189 cells remain"
# [1] "    Filtered on #Features, pctMito and HK:4969 cells remain"
# [1] "  Filtered an additional 220 cells by selecting for more than 55 housekeeping genes being present"
# [1] " which corresponds to 4.2397379071112 %"
# [1] "MOMA68 - range of nFeatures_RNA: 47 - 7103"
# [1] "    Filtered on #Features and pctMito:3033 cells remain"
# [1] "    Filtered on #Features, pctMito and HK:2740 cells remain"
# [1] "  Filtered an additional 293 cells by selecting for more than 55 housekeeping genes being present"
# [1] " which corresponds to 9.66040224200462 %"
# [1] "MOMA72 - range of nFeatures_RNA: 36 - 7904"
# [1] "    Filtered on #Features and pctMito:1434 cells remain"
# [1] "    Filtered on #Features, pctMito and HK:1314 cells remain"
# [1] "  Filtered an additional 120 cells by selecting for more than 55 housekeeping genes being present"
# [1] " which corresponds to 8.36820083682008 %"
# [1] "MOMA302 - range of nFeatures_RNA: 46 - 7836"
# [1] "    Filtered on #Features and pctMito:3180 cells remain"
# [1] "    Filtered on #Features, pctMito and HK:2961 cells remain"
# [1] "  Filtered an additional 219 cells by selecting for more than 55 housekeeping genes being present"
# [1] " which corresponds to 6.88679245283019 %"


## OBSERVATIONS
# - Many more cells are filtered because of their mitochondrial content
# - Are the housekeeping genes that are now filtered, anyway cells that have low expression of genes??