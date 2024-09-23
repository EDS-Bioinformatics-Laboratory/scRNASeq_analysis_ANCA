# AJ - 20230104

# Reclustering of DC/Macrophages/Monocytes clusters

#### IMPORTANT NOTE ####
# 20230302 - Perry saw that used FindVariableFeatures after SCTransform, which
# should not have been done. Decided to leave it as is, since this would have
# changed all downstream analysis and the effect is extremely likely to be 
# minimal.

#### Set default for a good start ####
source("initAnalysis.r")
init()

# Create a directory to hold results and plots
outDir <- "10_Reclustering_MP_and_MC"
dir.create(outDir)

#### Library ####
library(Seurat)
library(cowplot)
library(patchwork)
library(ggplot2)
library(alluvial)
library(clustree)
library(dplyr)
library(corrplot)
library(tidyr)

#### Read the data ####
seurat <- readRDS("9_kidney_reference_mapping/Mito10.integrated.clustered.annot.kidney_immune_reference.rds")
DefaultAssay(seurat) <- "RNA"

# The monocytes and macrophages are in clusters 5, 7, 9 and 11 at resolution 0.6
# These clusters are assigned after integration, make sure you have selected the right assay (see line 31)
selClusters.seurat <- subset(x = seurat, subset = integrated_snn_res.0.6 %in% c("5","7","9","11"))
dim(selClusters.seurat)
# [1] 22769  4194

# Drop all superfluous factor levels
selClusters.seurat@meta.data <- droplevels(selClusters.seurat@meta.data)

# Clean the object
seurat <- DietSeurat(selClusters.seurat)

# And check
seurat
# An object of class Seurat 
# 46974 features across 4194 samples within 5 assays 
# Active assay: RNA (22769 features, 0 variable features)
# 4 other assays present: SCT, ADT, integrated, prediction.score.celltype

DefaultAssay(seurat) <- "RNA"

# Make a nicer sample annotation
seurat@meta.data$Condition <- ifelse(seurat@meta.data$orig.ident %in% c("MOMA17", "MOMA52", "MOMA67", "MOMA302", "MOMA72"), "ANCA",
                                     ifelse(seurat@meta.data$orig.ident == "MOMA57", "Control", "SLE"))

seurat@meta.data$subType <- ifelse(seurat@meta.data$orig.ident %in% c("MOMA17", "MOMA52", "MOMA302"), "PR3",
                                   ifelse(seurat@meta.data$orig.ident %in% c("MOMA57", "MOMA68"), "NONE", "MPO"))

seurat@meta.data$sampleAnnot <- paste0(gsub("MOMA", "s", seurat@meta.data$orig.ident),"_",seurat@meta.data$Condition,
                                       "_", seurat@meta.data$subType)

# Split the combined object
allSets.list <- SplitObject(seurat, split.by = "orig.ident")

# Perform standard preprocessing (SCTransform for RNA, CLR for ADT) and identify variable features per set
# We have enough features per cell, so take nfeatures=2000
keep <- c()
for (i in 1:length(allSets.list)) {
  # 
  cat(names(allSets.list)[i],dim(allSets.list[[i]]),"\n" )
  if (dim(allSets.list[[i]])[2] > 100 & !(names(allSets.list)[i] %in% c("Doublet","Negative")) ){
    allSets.list[[i]] <- SCTransform(allSets.list[[i]], vars.to.regress = "percent.mito", verbose = FALSE)
    allSets.list[[i]] <- FindVariableFeatures(allSets.list[[i]], selection.method = "vst", 
                                              nfeatures = 2000, verbose = FALSE)
    # And renormalize ADT???
    allSets.list[[i]] <- NormalizeData(allSets.list[[i]], normalization.method = "CLR", margin = 2, assay = "ADT")
    
    keep <- c(keep,i)
  } 
}

# MOMA17, MOMA302 have less than 250 cells, so put minimum at 100 cells ....

# MOMA52 22769 276 
# MOMA57 22769 1517 
# MOMA17 22769 107 
# MOMA302 22769 201 
# MOMA67 22769 292 
# MOMA68 22769 1538 
# MOMA72 22769 263

# Integrate 
features <- SelectIntegrationFeatures(object.list = allSets.list, nfeatures = 3000)
allSets.list <- PrepSCTIntegration(object.list = allSets.list, anchor.features = features)
anca.anchors <- FindIntegrationAnchors(object.list = allSets.list, normalization.method = "SCT",
                                       anchor.features = features)
anca.combined.sct <- IntegrateData(anchorset = anca.anchors, normalization.method = "SCT")

# Use 30 dimensions ...
anca.combined.sct <- RunPCA(anca.combined.sct, verbose = FALSE, npcs = 30)
anca.combined.sct <- RunUMAP(anca.combined.sct, reduction = "pca", dims = 1:30)

# And save
saveRDS(anca.combined.sct, paste(outDir, "selectedClusters.subjectIntegrated.rds", sep="/"))

# And clean
rm(anca.anchors, anca.combined.sct, features, seurat)

#### Plots ####
# Switch to integrated assay. 

# NOTE: Does the plotting take the raw counts or the scaled ones?
# Sometimes the y-axis looks a little odd ...

allSets.integrated <- readRDS(paste(outDir, "selectedClusters.subjectIntegrated.rds", sep="/"))
# allSets.integrated@assays$SCT@var.features
# logical(0)

# The variable features of this assay are automatically set during IntegrateData
DefaultAssay(allSets.integrated) <- "integrated"

# Create a directory to hold TIFFs
dir.create(paste0(outDir,"/tiff"))

pdf(paste(outDir, "dimPlot_afterIntegration.pdf", sep="/"))
Idents(allSets.integrated) <- "sampleAnnot"
p <- DimPlot(allSets.integrated, reduction = "umap", label = TRUE,  repel = TRUE) +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 3) ) )
  p[[1]]$layers[[1]]$aes_params$alpha = 0.6
  p[[1]]$layers[[1]]$mapping$alpha = 0.4
  print(p)
dev.off()

# And as TIFF - do something with the legend .... ??
p <- DimPlot(allSets.integrated, reduction = "umap", label = FALSE,  repel = TRUE) +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 3) ) )
  p[[1]]$layers[[1]]$aes_params$alpha = 0.6
  p[[1]]$layers[[1]]$mapping$alpha = 0.4

ggsave(paste(outDir, "/tiff/dimPlot_afterIntegration.tiff", sep="/"),p, scale=1, width=5, height=5, dpi=360,compression="lzw")

pdf(paste(outDir, "dimPlot_afterIntegration_groupBy_OldClusterRes0.6.pdf", sep="/"), width = 14, height = 7)
  Idents(allSets.integrated) <- "sampleAnnot"
  p <- DimPlot(allSets.integrated, reduction = "umap", group.by = "integrated_snn_res.0.6", label = TRUE,  repel = TRUE) +
    guides(color = guide_legend(override.aes = list(alpha = 1, size = 3) ) )
  p[[1]]$layers[[1]]$aes_params$alpha = 0.6
  p[[1]]$layers[[1]]$mapping$alpha = 0.4
  print(p)
  # Alternative
  p2 <- DimPlot(allSets.integrated, reduction = "umap", label = FALSE, pt.size = 1) +
    guides(color = guide_legend(override.aes = list(alpha = 1, size = 3) ) )
  p2[[1]]$layers[[1]]$aes_params$alpha = 0.6
  p2[[1]]$layers[[1]]$mapping$alpha = 0.4
  print(p2)
dev.off()

# And as TIFF
ggsave(paste(outDir, "/tiff/dimPlot_afterIntegration_groupBy_OldClusterRes0.6.p1.tiff", sep="/"),p, scale=1, width=5, height=5, dpi=360,compression="lzw")
ggsave(paste(outDir, "/tiff/dimPlot_afterIntegration_groupBy_OldClusterRes0.6.p2.tiff", sep="/"),p2, scale=1, width=5, height=5, dpi=360,compression="lzw")


## OBSERVATION:
# - Looks reassuringly like the UMAP with all data IMHO!!
# - Can clearly see the original 3 clusters as they appeared in 8_clustering/Mito10_custers/dimPlot_Clustered_res.0.6.pdf ?

#### Recluster ####
seurat <- readRDS(paste(outDir, "selectedClusters.subjectIntegrated.rds", sep="/"))
DefaultAssay(seurat)
# [1] "integrated"

#### Clustering ####
# Find neighbors
seurat.clustered <- FindNeighbors(seurat, reduction = "pca", dims = 1:30)

# Cluster at different resolution 
seurat.clustered <- FindClusters(seurat.clustered, resolution = 0.3)
# Number of communities: 6
seurat.clustered <- FindClusters(seurat.clustered, resolution = 0.6)
# Number of communities: 7
seurat.clustered <- FindClusters(seurat.clustered, resolution = 1.0)
# Number of communities: 12
seurat.clustered <- FindClusters(seurat.clustered, resolution = 1.4)
# Number of communities: 14

#And save
saveRDS(seurat.clustered, paste(outDir, "selectedClusters.subjectIntegrated.reclustered.rds", sep="/"))

# And clean
rm(seurat.clustered)

#### Plot clusters ####
seurat <- readRDS(paste(outDir, "selectedClusters.subjectIntegrated.reclustered.rds", sep="/"))

for (clust in colnames(seurat@meta.data)[grep("integrated_snn", colnames(seurat@meta.data))]){
  res <- gsub("integrated_snn_","", clust)
  
  pdf(paste(outDir, paste0("dimPlot_Clustered_", res,".pdf"), sep="/"))
    Idents(seurat) <- clust
    p <- DimPlot(seurat, reduction = "umap", label = TRUE, label.size = 6, repel = TRUE, pt.size = 1.5) +
      guides(colour=guide_legend(override.aes = list(size=2), ncol=1))
    print(p)
  dev.off()
  
  # And as TIFF
  ggsave(paste(outDir, paste0("/tiff/dimPlot_Clustered_", res,".tiff"), sep="/"),p, scale=2, width=5, height=5, dpi=360,compression="lzw")
  
}

# Per clustering resolution and per cluster get a table of how many cells from which subject are present
for (clust in colnames(seurat@meta.data)[grep("integrated_snn", colnames(seurat@meta.data))]){
  res <- gsub("integrated_snn_","", clust)
  
  t <- table(seurat@meta.data[,clust], seurat@meta.data$orig.ident)
  t.df <- as.data.frame.matrix(t)
  t.df$Cluster = rownames(t.df)
  t.df <- t.df[,c(5,1:4)]
  pct.df <- as.data.frame.matrix(round(proportions(as.matrix(t.df[,2:5]),1) * 100,2))
  colnames(pct.df) <- paste0("pct.",colnames(pct.df))
  t.df <- cbind(t.df, pct.df)
  
  # And normalize over columns
  norm.pct.df <- (scale(t.df[,2:5], center=FALSE, scale=colSums(as.data.frame(t.df[,2:5])))) * 100
  colnames(norm.pct.df) <- paste0("norm.",colnames(pct.df))
  t.df <- cbind(t.df, norm.pct.df)
  t.df[nrow(t.df)+1,] <- c("totals", colSums(as.data.frame(t.df[,2:5])), rep("-",4), colSums(norm.pct.df[,1:4]))
  
  # And save
  write.table(t.df, file=paste(outDir, paste0("dimPlot_Clustered_", res,".txt"), sep="/"), quote = FALSE,
              row.names = FALSE, col.names = TRUE, sep="\t")
  
  # And clean
  rm(t,t.df, pct.df)
}

# Violin plot per clustering resolution
# See: https://patchwork.data-imaginist.com/articles/guides/layout.html, https://patchwork.data-imaginist.com/articles/guides/assembly.html
# and https://www.rdocumentation.org/packages/grid/versions/3.6.2/topics/grid.text

# Common macrophage markers: CD14, CD16, CD64, CD68, CD71, CCR5
# See also https://www.frontiersin.org/articles/10.3389/fimmu.2019.01084/full, 
#          https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3448669/
# Non-classically activated MP: CD23

# For more markers see mail Yosta 13-12-2022:
# - Here only macrophages and interesting markers
#   macrophages: CD163, CD206 (= MRC1), MARCO, C1q, CD74 (Tissue-residency),
#                CD81 (Tissue-residency), CD14, CD16, HLA-DR, CD68
#   Interesting: C5aR1, C5aR2, CX3CR1, CCR2, TGF-b,	Procollagen1
#                IL1b, TNFa, MCP-1 (= CCL2), IL1R2 (erg upgereguleerd bij MPO monocyten actief), IL18R1

# See mail 20230414 - Add TREM2
## Macrophage markers
c("CD163", "CD206", "MRC1", "MARCO", "C1q", "CD74","CD81", "CD14", "CD16", "HLA-DR","CD68", "TREM2") %in% rownames(seurat@assays$SCT@data)
# [1]  TRUE FALSE  TRUE  TRUE FALSE  TRUE  TRUE  TRUE FALSE FALSE  TRUE TRUE

# CD206 == MRC1
# CD16 == FCGR3A, FCGR3B (both present!)
# CD64 == FCGR1A (present)
# CD68 == LAMP4, SCARD1, GP110 (NOT present)
# CD71 == TFRC (present)
# CD23 == FCER2 (present)
# HLA-DR == HLA-DRA
# C1Q == C1QA

c("CD163", "MRC1", "MARCO", "C1QA", "CD74","CD81", "CD14", "FCGR3A", "FCGR3B", "HLA-DRA","CD68", "TREM2") %in% rownames(seurat@assays$SCT@data)
#[1] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
# corresponding ADTs:
# c("anti-human CD163","anti-human CD81 (TAPA-1)","anti-human CD14", "anti-human CD16", "anti-human HLA-DR" )

## Interesting markers
# Procollagen1 == COL1A1
# TGF-b == TGFB1
# TNFa == TNF
c("C5AR1", "C5AR2", "CX3CR1", "CCR2", "TGFB1", "COL1A1","IL1B", "TNF", "CCL2", "IL1R2","IL18R1") %in% rownames(seurat@assays$SCT@data)
# [1]  TRUE  TRUE  TRUE  TRUE  TRUE FALSE  TRUE  TRUE  TRUE  TRUE  TRUE
# corresponding ADTs:
# c("anti-human CD88 (C5aR)","anti-human CX3CR1")

# T cel
c("CD3","CD4","CD8") %in% rownames(seurat@assays$SCT@data)
# [1] FALSE  TRUE FALSE
# CD3 == CD3E, CD3D, CD3G
# CD8 == CD8A, CD8B
c("CD3D", "CD3E", "CD3G","CD4","CD8A", "CD8B") %in% rownames(seurat@assays$SCT@data)
# [1] TRUE TRUE TRUE TRUE TRUE TRUE
# corresponding ADTs:
# c("anti-human CD3","anti-human CD4","anti-human CD8")

# B cel
c("CD19","CD20") %in% rownames(seurat@assays$SCT@data)
# [1]  TRUE FALSE
# CD20 == MS4A1
c("CD19","MS4A1") %in% rownames(seurat@assays$SCT@data)
# [1] TRUE TRUE
# corresponding ADTs:
# c("anti-human CD19","anti-human CD20")

# Plasma cell
c("CD138", "BCMA") %in% rownames(seurat@assays$SCT@data)
# [1] FALSE FALSE
# CD138 == SDC1, SYND1, SDC
# BCMA == TNFRSF17
c("SYND1", "TNFRSF17") %in% rownames(seurat@assays$SCT@data)
# [1] FALSE TRUE
# corresponding ADTs:
# -

# Neutrophil
c("CD16",  "CD15", "SSEA1", "CD66B") %in% rownames(seurat@assays$SCT@data)
# [1] FALSE FALSE FALSE FALSE
# CD16 == FCGR3A, FCGR3B (both present!)
# CD15 == FUT4 (SSEA-1)
# CD66b == CEACAM8 (CD67, CGM6)
c("FCGR3A", "FCGR3B",  "FUT4", "CEACAM8") %in% rownames(seurat@assays$SCT@data)
#[1]  TRUE  TRUE  TRUE FALSE
# corresponding ADTs:
# c("anti-human CD16")

# DC
c("CD11c","HLA-DR","CD1c","CLEG9a") %in% rownames(seurat@assays$SCT@data)
# [1] FALSE FALSE FALSE FALSE
# CD11c == CD11C, ITGAX
# HLA-DR == HLA-DRA
# CD1c == CD1C
# CLEG9a == CLEG9A, CD370, DNGR1
c("ITGAX","HLA-DRA","CD1C","DNGR1") %in% rownames(seurat@assays$SCT@data)
# [1]  TRUE  TRUE  TRUE FALSE
# corresponding ADTs:
# c("anti-human CD11c", "anti-human HLA-DR","anti-human CD1c")

# NK(T) cell
c("CD56") %in% rownames(seurat@assays$SCT@data)
# [1] FALSE
# CD56 == NCAM1
c("NCAM1") %in% rownames(seurat@assays$SCT@data)
# [1] TRUE
# corresponding ADTs:
# c("anti-human CD56")

markerList <- list(MP = c("CD163", "MRC1", "MARCO", "C1QA", "CD74","CD81", "CD14", "FCGR3A", "FCGR3B", "HLA-DRA","CD68", "TREM2"),
                   Interesting = c("C5AR1", "C5AR2", "CX3CR1", "CCR2", "TGFB1", "COL1A1","IL1B", "TNF", "CCL2", "IL1R2","IL18R1"),
                   T_cell = c("CD3D", "CD3E", "CD3G","CD4","CD8A", "CD8B"),
                   B_cell = c("CD19","MS4A1"),
                   Plasma_cell = c("SYND1", "TNFRSF17"),
                   Neutrophil = c("FCGR3A", "FCGR3B",  "FUT4", "CEACAM8"),
                   DC = c("ITGAX","HLA-DRA","CD1C","DNGR1"),
                   NK = c("NCAM1")
                  )

adtList <- list(MP = c("anti-human CD163","anti-human CD81 (TAPA-1)","anti-human CD14", "anti-human CD16", "anti-human HLA-DR" ),
                Interesting = c("anti-human CD88 (C5aR)","anti-human CX3CR1"),
                T_Cell = c("anti-human CD3","anti-human CD4","anti-human CD8"),
                B_cell = c("anti-human CD19","anti-human CD20"),
                Neutrophil = c("anti-human CD16"),
                DC = c("anti-human CD11c", "anti-human HLA-DR","anti-human CD1c"),
                NK = c("anti-human CD56")
               )

DefaultAssay(seurat) <- "SCT"

for (clust in colnames(seurat@meta.data)[grep("integrated_snn", colnames(seurat@meta.data))]){
  res <- gsub("integrated_snn_","", clust)
  seurat@meta.data[,clust] <- factor(seurat@meta.data[,clust], levels=as.character(0:(length(unique(seurat@meta.data[,clust]))-1)))
  
  # maxY.RNA <- max(seurat@assays$RNA@data[c("CD163", "MRC1", "MARCO", "C1QA", "CD74","CD81", "CD14", "FCGR3A", "FCGR3B", "HLA-DRA","CD68"),])
  maxY.RNA <- NULL
  
  # 
  pdf(paste(outDir,paste0("violinPlot_RNA_of_MP_markers_", res,".pdf"), sep="/"))
  Idents(seurat) <- clust
  p1 <- VlnPlot(object = seurat, features = c("CD163"), pt.size = 0.1, y.max = maxY.RNA, same.y.lims = TRUE) + NoLegend() +
    theme(axis.title.x = element_blank())
  p2 <- VlnPlot(object = seurat, features = c("MRC1"), pt.size = 0.1, y.max = maxY.RNA, same.y.lims = TRUE) + NoLegend() +
    theme(axis.title.x = element_blank())
  p3 <- VlnPlot(object = seurat, features = c("MARCO"), pt.size = 0.1, y.max = maxY.RNA, same.y.lims = TRUE) + NoLegend() +
    theme(axis.title.x = element_blank())
  p4 <- VlnPlot(object = seurat, features = c("C1QA"), pt.size = 0.1, y.max = maxY.RNA, same.y.lims = TRUE) + NoLegend() +
    theme(axis.title.x = element_blank())
  p5 <- VlnPlot(object = seurat, features = c("CD74"), pt.size = 0.1, y.max = maxY.RNA, same.y.lims = TRUE) + NoLegend() +
    theme(axis.title.x = element_blank())
  p6 <- VlnPlot(object = seurat, features = c("CD81"), pt.size = 0.1, y.max = maxY.RNA, same.y.lims = TRUE) + NoLegend() +
    theme(axis.title.x = element_blank())
  p7 <- VlnPlot(object = seurat, features = c("CD14"), pt.size = 0.1, y.max = maxY.RNA, same.y.lims = TRUE) + NoLegend() +
    theme(axis.title.x = element_blank())
  p8 <- VlnPlot(object = seurat, features = c("FCGR3A"), pt.size = 0.1, y.max = maxY.RNA, same.y.lims = TRUE) + NoLegend() +
    theme(axis.title.x = element_blank())
  p9 <- VlnPlot(object = seurat, features = c("FCGR3B"), pt.size = 0.1, y.max = maxY.RNA, same.y.lims = TRUE)  + NoLegend() +
    theme(axis.title.x = element_blank())
  p10 <- VlnPlot(object = seurat, features = c("HLA-DRA"), pt.size = 0.1, y.max = maxY.RNA, same.y.lims = TRUE)  + NoLegend() +
    theme(axis.title.x = element_blank())
  p11 <- VlnPlot(object = seurat, features = c("CD68"), pt.size = 0.1, y.max = maxY.RNA, same.y.lims = TRUE)  + NoLegend() +
    theme(axis.title.x = element_blank())
  p12 <- VlnPlot(object = seurat, features = c("TREM2"), pt.size = 0.1, y.max = maxY.RNA, same.y.lims = TRUE)  + NoLegend() +
    theme(axis.title.x = element_blank())
  # And print on three pages
  p <- p1 / p2 / p3 / p4  + plot_layout(heights = c(1, 1, 1, 1))
  print(p)
  
  p <- p5 / p6 / p7 / p8  + plot_layout(heights = c(1, 1, 1, 1))
  print(p)

  p <- p9 / p10 / p11 / p12  + plot_layout(heights = c(1, 1, 1, 1, 1))
  print(p)
  dev.off()
  
  # Plot some QC; fraction expression mitochondrial genes, housekeeping genes, number of features per cluster
  pdf(paste(outDir,paste0("violinPlot_qcStats_", res,".pdf"), sep="/"))
  p1 <- VlnPlot(seurat,features="percent.mito", pt.size = 0.1) + NoLegend() +
    ggtitle("fraction.mito") +
    theme(axis.title.x = element_blank())
  p2 <- VlnPlot(seurat,features="n.exp.hkgenes", pt.size = 0.1) + NoLegend() +
    theme(axis.title.x = element_blank())
  p3 <- VlnPlot(seurat,features="nFeature_RNA", pt.size = 0.1) + NoLegend() +
    theme(axis.title.x = element_blank())
  p <- p1 / p2 / p3  + plot_layout(heights = c(1, 1, 1))
  print(p)
  dev.off()
  
}

## OBSERVATIONS
# - Is the y-axis odd? i.e. don't you expect expression values, smaller numbers?

# Make these plots per marker list
# Make a function from the above plotting routine
# Or perhaps make a function just for the plotting in panels?
# - i.e. make the list with ggplots first, so you have more control?

plotPanels <- function(p, nrPanels, myTitle, nrow = 1){
  # Make sure the name of the list items is the desired name to be used in the plot
  pdf(paste(outDir,paste0(myTitle,".pdf"), sep="/"))
  # And print X panels per page
  if (length(p) >= nrPanels){
    for (i in seq(1, nrPanels*floor(length(p)/nrPanels), by=nrPanels)){
      listEnd <- (i+nrPanels)-1
      pp <- patchwork::wrap_plots(p[i:listEnd], nrow = nrow)  + plot_layout(heights = rep(1,nrPanels))
      print(pp)
    }
    i <- i + nrPanels
  } else {
    i <- 1
  }
  
  # And plot the remaining figures
  if ( length(p)%%nrPanels >= 1 ){
    p <- p[i:length(p)]
    for (j in (length(p)+1):nrPanels){
      p[[j]] <- plot_spacer()
    }
    pp <- patchwork::wrap_plots(p, nrow = nrow)  + plot_layout(heights = rep(1,nrPanels))
    print(pp)
  }
  dev.off()
}

plotMarkers_and_QC <- function(seurat, markers, my_ident, my_plotName){
  seurat@meta.data[,my_ident] <- factor(seurat@meta.data[,my_ident], levels=as.character(0:(length(unique(seurat@meta.data[,my_ident]))-1)))
  
  maxY.RNA <- NULL
  
  # Make list with plots
  Idents(seurat) <- my_ident
  p <- list()
  
  # Check whether all markers can be found
  markers <- markers[markers %in% rownames(seurat@assays$SCT@data)]
  
  for (i in 1:length(markers)){
    p[[i]] <- VlnPlot(object = seurat, features = markers[i], pt.size = 0.1, y.max = maxY.RNA, same.y.lims = TRUE) + NoLegend() +
      theme(axis.title.x = element_blank())
  }
  
  pdf(paste(outDir,paste0("violinPlot_RNA_of_", my_plotName,".pdf"), sep="/"))
    # And print 4 panels per page
    if (length(markers) >= 4){
      for (i in seq(1, 4*floor(length(markers)/4), by=4)){
        pp <- p[[i]] / p[[i+1]] / p[[i+2]] / p[[i+3]]  + plot_layout(heights = c(1, 1, 1, 1))
        print(pp)
      }
      i <- i + 4
    } else {
      i <- 1
    }
  
    if ( length(markers)%%4 == 1 ){
      pp <- p[[i]] / plot_spacer() / plot_spacer() / plot_spacer()  + plot_layout(heights = c(1, 1, 1, 1, 1))
      print(pp)
    }
    if ( length(markers)%%4 == 2 ){
      pp <- p[[i]] / p[[i+1]] / plot_spacer() / plot_spacer()  + plot_layout(heights = c(1, 1, 1, 1, 1))
      print(pp)
    }
    if ( length(markers)%%4 == 3 ){
      pp <- p[[i]] / p[[i+1]] / p[[i+2]] / plot_spacer()  + plot_layout(heights = c(1, 1, 1, 1, 1))
      print(pp)
    }
  dev.off()

  # And clean
  rm(pp,p)
}

# # And now for every resolution and every marker set
# for (clust in colnames(seurat@meta.data)[grep("integrated_snn", colnames(seurat@meta.data))]){
#   res <- gsub("integrated_snn_","", clust)
#   
#   for (markers in names(markerList)){
#     plotName <- paste0(markers,"_markers_",res)
#     plotMarkers_and_QC(seurat, markerList[[markers]], clust, plotName )
#   }
#   
#   # Plot some QC; fraction expression mitochondrial genes, housekeeping genes, number of features per cluster
#   pdf(paste(outDir,paste0("violinPlot_qcStats_", res,".pdf"), sep="/"))
#   p1 <- VlnPlot(seurat,features="percent.mito", pt.size = 0.1) + NoLegend() +
#     ggtitle("fraction.mito") +
#     theme(axis.title.x = element_blank())
#   p2 <- VlnPlot(seurat,features="n.exp.hkgenes", pt.size = 0.1) + NoLegend() +
#     theme(axis.title.x = element_blank())
#   p3 <- VlnPlot(seurat,features="nFeature_RNA", pt.size = 0.1) + NoLegend() +
#     theme(axis.title.x = element_blank())
#   p <- p1 / p2 / p3  + plot_layout(heights = c(1, 1, 1))
#   print(p)
#   dev.off()
#   
# }

# Alternatively, more flexible for what kind of plot you want ...
for (clust in colnames(seurat@meta.data)[grep("integrated_snn", colnames(seurat@meta.data))]){
  res <- gsub("integrated_snn_","", clust)
  seurat@meta.data[,clust] <- factor(seurat@meta.data[,clust], levels=as.character(0:(length(unique(seurat@meta.data[,clust]))-1)))

  maxY.RNA <- NULL

  # Make list with plots
  Idents(seurat) <- clust
  
  for (markerSet in names(markerList)){
    markers <- markerList[[markerSet]]
    # Check whether all markers can be found
    markers <- markers[markers %in% rownames(seurat@assays$SCT@data)]

    p <- list()
    
    for (i in 1:length(markers)){
      p[[i]] <- VlnPlot(object = seurat, features = markers[i], pt.size = 0.1, y.max = maxY.RNA, same.y.lims = TRUE) + NoLegend() +
                theme(axis.title.x = element_blank())
      
    }
    myTitle <- paste0("violinPlot_RNA_of_", markerSet,"_markers_",res)
    plotPanels(p, 4, myTitle, nrow=4)
  }
  
  # Plot some QC; fraction expression mitochondrial genes, housekeeping genes, number of features per cluster
  pdf(paste(outDir,paste0("violinPlot_qcStats_", res,".pdf"), sep="/"))
    p1 <- VlnPlot(seurat,features="percent.mito", pt.size = 0.1) + NoLegend() +
      ggtitle("fraction.mito") +
      theme(axis.title.x = element_blank())
    p2 <- VlnPlot(seurat,features="n.exp.hkgenes", pt.size = 0.1) + NoLegend() +
      theme(axis.title.x = element_blank())
    p3 <- VlnPlot(seurat,features="nFeature_RNA", pt.size = 0.1) + NoLegend() +
      theme(axis.title.x = element_blank())
    p <- p1 / p2 / p3  + plot_layout(heights = c(1, 1, 1))
    print(p)
  dev.off()
  
}

### And similarly for the selected ADT markers
DefaultAssay(seurat) <- "ADT"
for (clust in colnames(seurat@meta.data)[grep("integrated_snn", colnames(seurat@meta.data))]){
  res <- gsub("integrated_snn_","", clust)
  seurat@meta.data[,clust] <- factor(seurat@meta.data[,clust], levels=as.character(0:(length(unique(seurat@meta.data[,clust]))-1)))
  
  maxY.RNA <- NULL
  
  # Make list with plots
  Idents(seurat) <- clust
  
  for (markerSet in names(adtList)){
    markers <- adtList[[markerSet]]
    # Check whether all ADT markers can be found
    markers <- markers[markers %in% rownames(seurat@assays$ADT@data)]
    
    p <- list()
    
    for (i in 1:length(markers)){
      p[[i]] <- VlnPlot(object = seurat, features = markers[i], pt.size = 0.1, y.max = maxY.RNA, same.y.lims = TRUE) + NoLegend() +
        theme(axis.title.x = element_blank())
      
    }
    myTitle <- paste0("violinPlot_RNA_of_", markerSet,"_ADT_markers_",res)
    plotPanels(p, 4, myTitle, nrow=4)
  }
  
}

#### Reference mapping ####
# see also 9_kidney_reference_mapping
# Mapping is against Stewart et al 2019

# Coloring function
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

## Import the Kidney Immune reference ##
# Copied from ../../20220623_AllSamples_InitialAnalysis/Code/8_kidney_reference_mapping
if(file.exists(paste(outDir,"kidney_immune_reference.rds", sep="/"))){
  reference <- readRDS(paste(outDir,"kidney_immune_reference.rds",sep="/"))
} else {
  reference <- readRDS(paste("../../20220623_AllSamples_InitialAnalysis/Code/8_kidney_reference_mapping","kidney_immune_reference.rds",sep="/"))
}

# Read the data
query1 <- readRDS(paste(outDir, "selectedClusters.subjectIntegrated.reclustered.rds", sep="/"))
if(file.exists(paste(outDir,"kidney_immune_anchors.rds", sep="/"))){
  kidney.anchors <- readRDS(paste(outDir,"kidney_immune_anchors.rds", sep="/"))
} else {
  kidney.anchors <- FindTransferAnchors(reference = reference, query = query1, dims = 1:30, reference.reduction = "pca")
  saveRDS(kidney.anchors, paste(outDir,"kidney_immune_anchors.rds", sep="/"))
}

predictions <- TransferData(anchorset = kidney.anchors, refdata = reference$celltype, dims = 1:30)
query1 <- AddMetaData(query1, metadata = predictions)

write.table(table(predictions$predicted.id, query1$integrated_snn_res.0.3),paste(outDir,"kidney_immune_query_cluster_res_0.3.txt", sep="/"), sep = "\t", quote = FALSE, col.names = NA)
write.table(table(predictions$predicted.id, query1$integrated_snn_res.0.6),paste(outDir,"kidney_immune_query_cluster_res_0.6.txt", sep="/"), sep = "\t", quote = FALSE, col.names = NA)
write.table(table(predictions$predicted.id, query1$integrated_snn_res.1),paste(outDir,"kidney_immune_query_cluster_res_1.txt", sep="/"), sep = "\t", quote = FALSE, col.names = NA)
write.table(table(predictions$predicted.id, query1$integrated_snn_res.1.4),paste(outDir,"kidney_immune_query_cluster_res_1.4.txt", sep="/"), sep = "\t", quote = FALSE, col.names = NA)

# Check whether there already is a 'umap' slot
if(!("^umap" %in% names(reference@reductions))){
  reference <- RunUMAP(reference, dims = 1:30, reduction = "pca", return.model = TRUE)
}
query1 <- MapQuery(anchorset = kidney.anchors, reference = reference, query = query1,
                   refdata = list(celltype = "celltype"), reference.reduction = "pca", reduction.model = "umap")

# And save
saveRDS(query1, paste(outDir,"selectedClusters.subjectIntegrated.reclustered.kidney_immune_reference.rds", sep="/"))
saveRDS(reference, paste(outDir,"kidney_immune_reference.rds", sep="/"))

# Read the data
reference <- readRDS(paste(outDir,"kidney_immune_reference.rds", sep="/"))
query1 <- readRDS(paste(outDir,"selectedClusters.subjectIntegrated.reclustered.kidney_immune_reference.rds", sep="/"))

# Get common colors
refColors <- gg_color_hue(length(levels(factor(reference@meta.data$celltype))))

# Make the levels order the same
query1@meta.data$predicted.celltype <- factor(query1@meta.data$predicted.celltype, levels = levels(factor(reference@meta.data$celltype)))

pdf(paste(outDir,"umap_reference_query.pdf", sep="/"), width=14)
p1 <- DimPlot(reference, reduction = "umap", group.by = "celltype", label = TRUE, label.size = 3,
              repel = TRUE, cols = refColors) + NoLegend() + ggtitle("Reference annotations")
p2 <- DimPlot(query1, reduction = "ref.umap", group.by = "predicted.celltype", label = TRUE,
              label.size = 3, repel = TRUE, cols = refColors) + NoLegend() + ggtitle("Query transferred labels")
p <- p1 + p2
print(p)
dev.off()

pdf(paste(outDir,"umap_query_kidney_immune.pdf", sep="/"))
p <- DimPlot(query1, reduction = "ref.umap", group.by = "predicted.celltype", label = TRUE, pt.size = 1.7,
        label.size = 5, repel = TRUE, cols = refColors) + NoLegend() + ggtitle("Query transferred labels: immune cells") +
  xlim(-10,16) + ylim(-17,10) 
print(p)
dev.off()


pdf(paste(outDir,"umap_query_kidney_immune_cluster_res_0.3.pdf", sep="/"), width=14)
p1 <- DimPlot(query1, reduction = "ref.umap", group.by = "predicted.celltype", label = TRUE,
              label.size = 3, repel = TRUE, cols = refColors) + NoLegend() + ggtitle("Query transferred labels")
p2 <- DimPlot(object = query1, reduction = "ref.umap", group.by = "integrated_snn_res.0.3", label = TRUE,
              label.size = 5, repel = TRUE) + NoLegend() + ggtitle("Cluster labels (resolution = 0.3)")
p <- p1 + p2
print(p)
dev.off()

pdf(paste(outDir,"umap_query_kidney_immune_cluster_res_0.6.pdf", sep="/"), width=14)
p1 <- DimPlot(query1, reduction = "ref.umap", group.by = "predicted.celltype", label = TRUE,
              label.size = 3, repel = TRUE, cols = refColors) + NoLegend() + ggtitle("Query transferred labels")
p2 <- DimPlot(object = query1, reduction = "ref.umap", group.by = "integrated_snn_res.0.6", label = TRUE,
              label.size = 5, repel = TRUE) + NoLegend() + ggtitle("Cluster labels (resolution = 0.6)")
p <- p1 + p2
print(p)
dev.off()

pdf(paste(outDir,"umap_query_kidney_immune_cluster_res_1.pdf", sep="/"), width=14)
p1 <- DimPlot(query1, reduction = "ref.umap", group.by = "predicted.celltype", label = TRUE,
              label.size = 3, repel = TRUE, cols = refColors) + NoLegend() + ggtitle("Query transferred labels")
p2 <- DimPlot(object = query1, reduction = "ref.umap", group.by = "integrated_snn_res.1", label = TRUE,
              label.size = 5, repel = TRUE) + NoLegend() + ggtitle("Cluster labels (resolution = 1)")
p <- p1 + p2
print(p)
dev.off()

pdf(paste(outDir,"umap_query_kidney_immune_cluster_res_1.4.pdf", sep="/"), width=14)
p1 <- DimPlot(query1, reduction = "ref.umap", group.by = "predicted.celltype", label = TRUE,
              label.size = 3, repel = TRUE, cols = refColors) + NoLegend() + ggtitle("Query transferred labels")
p2 <- DimPlot(object = query1, reduction = "ref.umap", group.by = "integrated_snn_res.1.4", label = TRUE,
              label.size = 5, repel = TRUE) + NoLegend() + ggtitle("Cluster labels (resolution = 1.4)")
p <- p1 + p2
print(p)
dev.off()

# flow diagram
clustResults <- query1@meta.data[ , c("predicted.celltype", "integrated_snn_res.0.3", "integrated_snn_res.1.4")]

# Column 'Freq' to store the total count of all combinations
clustResults$Freq <- 1
clustResults2D <- aggregate(Freq ~ ., data = clustResults, sum)

cols <- rev(rainbow(nrow(clustResults2D), start = 0.1, end = 0.9))
pdf(paste(outDir,"kidney_immune_celltype_clustering_flowdiagram.pdf", sep="/"), height=14, width=14)
alluvial(
  select(clustResults2D, predicted.celltype, integrated_snn_res.0.3, integrated_snn_res.1.4),
  freq = clustResults2D$Freq,
  col = cols,
  alpha = 0.8,
  gap.width = 0.5,
  cw = 0.2,
  blocks = FALSE,
  axis_labels = sub("integrated_snn_res.","",colnames(clustResults2D)[1:3]),
  cex = 1.5,
  cex.axis = 1.5)  
dev.off()

# Predict compartment: {lymphoid,myeloid}
predictions.compartment <- TransferData(anchorset = kidney.anchors, refdata = reference$compartment, dims = 1:30)
table(predictions.compartment$predicted.id, query1$integrated_snn_res.0.3)
#              0   1   2   3   4
#   lymphoid 591  44  79 149 205
#   myeloid  349 840 796 603 538

table(predictions.compartment$predicted.id, query1$integrated_snn_res.0.6)
#             0   1   2   3   4   5   6   7   8
#  lymphoid  38  72 280 103  74 319 129  46   7
#  myeloid  803 755 350 403 332  33 207 196  47

# Clean as you get the dreaded 'Error: cannot allocate vector of size .... '
rm(query1, reference, p1, p2, predictions, predictions.compartment, clustResults, clustResults2D, kidney.anchors)
gc()

#### CellTypist with majority-voting ####
# Using majority-voting within CellTypist just on this subset of data will have a different effect
# as using it on the whole data set.
# So, extract the raw counts from the subset, run Celltypist and add the resulting data back to the subset
# for analysis and visualisation

# Read the data
seurat <- readRDS(paste(outDir, "selectedClusters.subjectIntegrated.reclustered.rds", sep="/"))

# Get the raw counts #
write.table(as.data.frame(seurat@assays$RNA@counts), paste(outDir, "selectedClusters.final.rawcounts.csv", sep="/"),
            sep = ',', row.names = T, col.names = T, quote = F)

# Run CellTypist via NoteBook: ../NoteBooks/runCellTypist_selectedClusters.ipynb

# Get the resulting data
# - The 'probability_matrix.csv' contains the probabilities per cell of to which celltype it would fit.
annot.df <- read.csv("./10_Reclustering_MP_and_MC/selectedClusters/majorityVoting_predicted_labels.csv")
probMatrix.df <- read.csv("./10_Reclustering_MP_and_MC/selectedClusters/majorityVoting_probability_matrix.csv", check.names = FALSE)

# Make sure to keep the rownames !!
seurat@meta.data$barcodes.orig <- rownames(seurat@meta.data)

# now, merge by rowname and celltypistName
colnames(annot.df) <- c("cellName","predicted.labels.celltypist.reclustered", "over_clustering.celltypist.reclustered", 
                        "majority_voting.celltypist.reclustered")
#seurat@meta.data <- merge(seurat@meta.data, annot.df, by.x="celltypistName",by.y="cellName", all.x = TRUE, sort = FALSE )
seurat@meta.data <- merge(seurat@meta.data, annot.df, by.x="barcodes.orig",by.y="cellName", all.x = TRUE, sort = FALSE )

# Put the rownames back
rownames(seurat@meta.data) <- seurat@meta.data$barcodes.orig

# And add the probability scores ... becoming a very large object :-)
colnames(probMatrix.df) <- paste0(colnames(probMatrix.df),".celltypist.reclustered")
seurat@meta.data <- merge(seurat@meta.data, probMatrix.df, by.x="barcodes.orig",by.y=".celltypist.reclustered", all.x = TRUE, sort = FALSE )

# Put the rownames back
rownames(seurat@meta.data) <- seurat@meta.data$barcodes.orig

# Create a celltypist.score column holding the probability of the assignment (the highest one in case two or more
# celltypes have been assigned. Then the celltypes are separated by a "|")
seurat@meta.data$celltypist.max_probability <- 0
for (i in 1:nrow(seurat@meta.data)){
  if (seurat$majority_voting.celltypist.reclustered[i] != "Unassigned"){
    celltypes <- strsplit(seurat@meta.data$majority_voting.celltypist.reclustered[i],"[|]")[[1]]
    seurat@meta.data$celltypist.max_probability[i] <- max(seurat@meta.data[i,paste0(celltypes,".celltypist.reclustered")])
  } 
}

## And plot ##
# Use the majority votings as labels
pdf(paste(outDir, "umap_celltypist_annotation_res0.3.pdf", sep="/"), width = 14, height = 7)
  p1 <- DimPlot(object = seurat, reduction = "umap", group.by = "majority_voting.celltypist.reclustered", label = TRUE) + NoLegend()
  p2 <- DimPlot(object = seurat, reduction = "umap", group.by = "integrated_snn_res.0.3", label = TRUE) + NoLegend()
  plot_grid(p1, p2)
  # Get the legends on a separate page using the cowplot package
  p1 <- DimPlot(object = seurat, reduction = "umap", group.by = "majority_voting.celltypist.reclustered", label = TRUE)
  legend_p1 <- cowplot::get_legend(p1)
  p2 <- DimPlot(object = seurat, reduction = "umap", group.by = "integrated_snn_res.0.3", label = TRUE)
  legend_p2 <- cowplot::get_legend(p2)
  plot_grid(legend_p1, NULL, legend_p2,ncol = 3)
dev.off()

pdf(paste(outDir, "umap_celltypist_annotation_res0.6.pdf", sep="/"), width = 14, height = 7)
p1 <- DimPlot(object = seurat, reduction = "umap", group.by = "majority_voting.celltypist.reclustered", label = TRUE) + NoLegend()
p2 <- DimPlot(object = seurat, reduction = "umap", group.by = "integrated_snn_res.0.6", label = TRUE) + NoLegend()
plot_grid(p1, p2)
# Get the legends on a separate page using the cowplot package
p1 <- DimPlot(object = seurat, reduction = "umap", group.by = "majority_voting.celltypist.reclustered", label = TRUE)
legend_p1 <- cowplot::get_legend(p1)
p2 <- DimPlot(object = seurat, reduction = "umap", group.by = "integrated_snn_res.0.6", label = TRUE)
legend_p2 <- cowplot::get_legend(p2)
plot_grid(legend_p1, NULL, legend_p2,ncol = 3)
dev.off()

pdf(paste(outDir, "umap_celltypist_annotation_res1.pdf", sep="/"), width = 14, height = 7)
p1 <- DimPlot(object = seurat, reduction = "umap", group.by = "majority_voting.celltypist.reclustered", label = TRUE) + NoLegend()
p2 <- DimPlot(object = seurat, reduction = "umap", group.by = "integrated_snn_res.1", label = TRUE) + NoLegend()
plot_grid(p1, p2)
# Get the legends on a separate page using the cowplot package
p1 <- DimPlot(object = seurat, reduction = "umap", group.by = "majority_voting.celltypist.reclustered", label = TRUE)
legend_p1 <- cowplot::get_legend(p1)
p2 <- DimPlot(object = seurat, reduction = "umap", group.by = "integrated_snn_res.1", label = TRUE)
legend_p2 <- cowplot::get_legend(p2)
plot_grid(legend_p1, NULL, legend_p2,ncol = 3)
dev.off()

pdf(paste(outDir, "umap_celltypist_probabilities.pdf", sep="/"), width = 14, height = 7)
  p1 <- FeaturePlot(object = seurat, features = "celltypist.max_probability", reduction = "umap")
  plot_grid(p1)
dev.off()

pdf(paste(outDir, "violin_celltypist_probabilities.pdf", sep="/"), width = 14, height = 14)
  p1 <- VlnPlot(object = seurat, features = "celltypist.max_probability", group.by = "majority_voting.celltypist.reclustered") + NoLegend()
  plot_grid(p1)
dev.off()

pdf(paste(outDir, "celltype_celltypist_annotation.pdf", sep="/"), width = 14, height =14)
  tab <- table(seurat$integrated_snn_res.0.3, seurat$majority_voting.celltypist.reclustered)
  corrplot(tab/rowSums(tab), is.corr = FALSE, col = 
           colorRampPalette(c("#FFFFFF", "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC", "#053061"))(200))
dev.off()

pdf(paste(outDir, "celltype_celltypist_annotation_res1.pdf", sep="/"), width = 14, height =14)
tab <- table(seurat$integrated_snn_res.1, seurat$majority_voting.celltypist.reclustered)
corrplot(tab/rowSums(tab), is.corr = FALSE, col = 
           colorRampPalette(c("#FFFFFF", "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC", "#053061"))(200))
dev.off()

# And get some tables 
# - of the distribution of the different celltypes over the different clusters
# - of the distribution of the different celltypes over the different subjects
# - of the distribution of the different celltypes over the different subjects and clusters

# A lot of 'Unassigned'??

# We do not use the CellTypist annotation any more ... but still here
tab <- table(seurat$integrated_snn_res.0.3, seurat$majority_voting.celltypist.reclustered)
write.table(tab, paste(outDir,"celltypes_per_cluster_res0_3.txt", sep="/"), quote=FALSE, sep="\t", row.names = TRUE, col.names = TRUE)

tab <- table(seurat$integrated_snn_res.0.6, seurat$majority_voting.celltypist.reclustered)
write.table(tab, paste(outDir,"celltypes_per_cluster_res0_6.txt", sep="/"), quote=FALSE, sep="\t", row.names = TRUE, col.names = TRUE)

tab <- table(seurat$integrated_snn_res.1, seurat$majority_voting.celltypist.reclustered)
write.table(tab, paste(outDir,"celltypes_per_cluster_res1.txt", sep="/"), quote=FALSE, sep="\t", row.names = TRUE, col.names = TRUE)

tab <- table(seurat$sampleAnnot, seurat$majority_voting.celltypist.reclustered)
write.table(tab, paste(outDir,"celltypes_per_subject.txt", sep="/"), quote=FALSE, sep="\t", row.names = TRUE, col.names = TRUE)

tab_0_3 <- table(seurat$integrated_snn_res.0.3, seurat$sampleAnnot, seurat$majority_voting.celltypist.reclustered)
write.table(tab_0_3, paste(outDir,"celltypes_per_subject_and_cluster_res0_3.txt", sep="/"), quote=FALSE, sep="\t", row.names = TRUE, col.names = TRUE)

tab_0_6 <- table(seurat$integrated_snn_res.0.6, seurat$sampleAnnot, seurat$majority_voting.celltypist.reclustered)
write.table(tab_0_6, paste(outDir,"celltypes_per_subject_and_cluster_res0_6.txt", sep="/"), quote=FALSE, sep="\t", row.names = TRUE, col.names = TRUE)

tab_1 <- table(seurat$integrated_snn_res.1, seurat$sampleAnnot, seurat$majority_voting.celltypist.reclustered)
write.table(tab_1, paste(outDir,"celltypes_per_subject_and_cluster_res1.txt", sep="/"), quote=FALSE, sep="\t", row.names = TRUE, col.names = TRUE)

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

for (res in c("0_3","0_6", "1")){
  tab <- read.table(paste(outDir,paste0("celltypes_per_subject_and_cluster_res",res,".txt"), sep="/"), sep="\t", header = TRUE, row.names = 1)
  df <- as.data.frame(tab)
  colnames(df) <- c("cluster", "subject", "celltype", "frequency")

  pdf(paste(outDir, paste0("distribution_per_subject_and_cluster_res", res,".pdf"), sep="/"), width = 14, height =14)
    p <- ggplot(df, aes(x=cluster, y=frequency))+
      geom_bar(stat='identity', fill="forest green")+
      facet_wrap(~subject,  ncol=1, strip.position = "left")
    print(p)
  dev.off()

  df$celltype <- factor(df$celltype)
  #nrColors <- length(levels(factor(df$celltype)))
  #myColors = rainbow(nrColors)
  myColors <- gg_color_hue(length(unique(df$celltype)))
  df$colors <- myColors[factor(df$celltype)]

  pdf(paste(outDir, paste0("distribution_per_subject_and_cluster_stacked_res", res,".pdf"), sep="/"), width = 14, height =14)
    p <- ggplot(df, aes(x=cluster, y=frequency, fill=celltype)) +
      geom_bar(stat='identity', position="stack")+
      facet_wrap(~subject,  ncol=1, strip.position = "left") +
      scale_fill_manual(values = myColors) +
      guides(fill=guide_legend(ncol =1))
    print(p)
    legend_p <- cowplot::get_legend(p)
    p <- ggplot(df, aes(x=cluster, y=frequency, fill=celltype)) +
      geom_bar(stat='identity', position="stack")+
      facet_wrap(~subject,  ncol=1, strip.position = "left") +
      scale_fill_manual(values = myColors) + 
      NoLegend()
    print(p)
    plot_grid(legend_p)
  dev.off()

  # percentages per subject
  df2 <- group_by(df, c(subject)) %>% mutate(percent = frequency/sum(frequency))
  df2$celltype <- factor(df2$celltype)
  #nrColors <- length(levels(factor(df2$celltype)))
  #myColors = rainbow(nrColors)
  myColors <- gg_color_hue(length(unique(df2$celltype)))
  df2$colors <- myColors[factor(df2$celltype)]

  pdf(paste(outDir, paste0("distribution_per_subject_and_cluster_stacked_percentages_res", res,".pdf"), sep="/"), width = 14, height =14)
  p <- ggplot(df2, aes(x=cluster, y=percent, fill=celltype)) +
    geom_bar(stat='identity', position="stack")+
    facet_wrap(~subject,  ncol=1, strip.position = "left") +
    scale_fill_manual(values = myColors) +
    guides(fill=guide_legend(ncol =1))
  print(p)
  legend_p <- cowplot::get_legend(p)
  p <- ggplot(df2, aes(x=cluster, y=percent, fill=celltype)) +
    geom_bar(stat='identity', position="stack")+
    facet_wrap(~subject,  ncol=1, strip.position = "left") +
    scale_fill_manual(values = myColors) + 
    NoLegend()
  print(p)
  plot_grid(legend_p)
  dev.off()
  
  # percentages per subject
  df2 <- group_by(df, c(subject)) %>% mutate(percent = frequency/sum(frequency))
  sum(df2[grep("s17_ANCA_PR3", df2$subject), "percent"])
  #[1] 1
  df2$celltype <- factor(df2$celltype)
  #nrColors <- length(levels(df2$celltype))
  #rainbowColors = rainbow(nrColors)
  myColors <- gg_color_hue(length(unique(df2$celltype)))
  df2$colors <- myColors[factor(df2$celltype)]
  
  df2$subject <- factor(df2$subject, 
                        levels = c("s17_ANCA_PR3", "s302_ANCA_PR3","s52_ANCA_PR3", 
                                   "s72_ANCA_MPO", "s67_ANCA_MPO", "s68_SLE_NONE", "s57_Control_NONE"),
                        labels = c("s1_ANCA_PR3", "s2_ANCA_PR3","s3_ANCA_PR3", 
                                   "s4_ANCA_MPO", "s5_ANCA_MPO", "SLE", "HC"))
  
  levels(df2$celltype)
  #[1] "Classical monocytes"  "DC1"  "DC2"  "Intermediate macrophages" "Non-classical monocytes" 
  #[5] "Unassigned" 
  
  # Celltypes (CellTypist annotation - not used anymore - 20230302)
  pdf(paste(outDir, paste0("distribution_per_subject_stacked_percentages_res", res,".pdf"), sep="/"), width = 14, height =14)
  p <- ggplot(df2, aes(x=subject, y=percent, fill=celltype)) +
    geom_bar(stat='identity', position="stack")+
    scale_fill_manual(values = myColors) +
    guides(fill=guide_legend(ncol =1)) +
    theme_cowplot() +
    theme(axis.text.x=element_text(angle=90,hjust=1))
  print(p)
  legend_p <- cowplot::get_legend(p)
  p2 <- ggplot(df2, aes(x=subject, y=percent, fill=celltype)) +
    geom_bar(stat='identity', position="stack")+
    scale_fill_manual(values = myColors) +
    guides(fill=guide_legend(ncol =1)) +
    theme_cowplot() +
    theme(axis.text.x=element_text(angle=90,hjust=1)) +
    NoLegend()
  print(p2)
  # And print the legend
  plot_grid(legend_p)
  dev.off()
  
  # And as TIFF
  ggsave(paste(outDir, paste0("/tiff/distribution_per_subject_stacked_percentages_res",res,".tiff"), sep="/"), p, scale=2, width=5, height=5, dpi=360,compression="lzw")
  ggsave(paste(outDir, paste0("/tiff/distribution_per_subject_stacked_percentages_res", res,".NoLegend.tiff"), sep="/"), p2, scale=2, width=5, height=5, dpi=360,compression="lzw")
  ggsave(paste(outDir, paste0("/tiff/distribution_per_subject_stacked_percentages_res", res,".Legend.tiff"), sep="/"), legend_p, scale=0.5, width=5, height=5, dpi=360,compression="lzw")
  
  
  # Cluster percentages
  tab <- table(seurat@meta.data[,paste0("integrated_snn_res.",gsub("_",".",res))], seurat$sampleAnnot)
  write.table(tab, paste(outDir,paste0("clusters_per_subject_freq_res", res,".txt"), sep="/"), quote=FALSE, sep="\t", row.names = TRUE, col.names = TRUE)
  
  df <- as.data.frame(tab)
  colnames(df) <- c("cluster", "subject", "frequency")
  
  if( res == "1"){
    myColors <- gg_color_hue(length(unique(df$cluster)))
    
    df_res1 <- df[df$cluster %in% c("0","1","4","8"), ]
    df2 <- group_by(df_res1, c(subject)) %>% mutate(percent = frequency/sum(frequency))
    sum(df2[grep("s17_ANCA_PR3", df2$subject), "percent"])
    #[1] 1
    df2$cluster <- factor(df2$cluster, levels=c("0","1","4","8"))
    
    myColors = myColors[c(0,1,4,8)+1]
    df2$colors <- myColors[factor(df2$cluster)]
    
    df2$subject <- factor(df2$subject, 
                          levels = c("s17_ANCA_PR3", "s302_ANCA_PR3","s52_ANCA_PR3", 
                                     "s72_ANCA_MPO", "s67_ANCA_MPO", "s68_SLE_NONE", "s57_Control_NONE"),
                          labels = c("s1_ANCA_PR3", "s2_ANCA_PR3","s3_ANCA_PR3", 
                                     "s4_ANCA_MPO", "s5_ANCA_MPO", "SLE", "HC"))
    
    # And write with percentages
    df3 <- df2[, c("cluster","subject", "percent")]
    df3 <- df3 %>% pivot_wider(names_from = subject, values_from = percent)
    write.table(as.data.frame(df3), paste(outDir,paste0("clusters_per_subject_percentages_res", res,".selectedClusters.txt"), sep="/"), quote=FALSE, sep="\t", row.names = FALSE, col.names = TRUE)
    
    pdf(paste(outDir, paste0("distribution_per_subject_stacked_percentages_cluster_res", res,".selectedClusters.pdf"), sep="/"), width = 14, height =14)
    p <- ggplot(df2, aes(x=subject, y=percent, fill=cluster)) +
      geom_bar(stat='identity', position="stack")+
      scale_fill_manual(values = myColors) +
      guides(fill=guide_legend(ncol =1)) +
      theme_cowplot() +
      theme(axis.text.x=element_text(angle=90,hjust=1))
    print(p)
    legend_p <- cowplot::get_legend(p)
    p2 <- ggplot(df2, aes(x=subject, y=percent, fill=cluster)) +
      geom_bar(stat='identity', position="stack")+
      scale_fill_manual(values = myColors) +
      guides(fill=guide_legend(ncol =1)) +
      theme_cowplot() +
      theme(axis.text.x=element_text(angle=90,hjust=1)) +
      NoLegend()
    print(p2)
    # And print the legend
    plot_grid(legend_p)
    dev.off()
    
    # And as TIFF
    ggsave(paste(outDir, paste0("/tiff/distribution_per_subject_stacked_percentages_clusters_res",res,".selectedClusters.tiff"), sep="/"), p, scale=2, width=5, height=5, dpi=360,compression="lzw")
    ggsave(paste(outDir, paste0("/tiff/distribution_per_subject_stacked_percentages_clusters_res", res,".selectedClusters.NoLegend.tiff"), sep="/"), p2, scale=2, width=5, height=5, dpi=360,compression="lzw")
    ggsave(paste(outDir, paste0("/tiff/distribution_per_subject_stacked_percentages_clusters_res", res,".selectedClusters.Legend.tiff"), sep="/"), legend_p, scale=0.5, width=5, height=5, dpi=360,compression="lzw")
    
  }
  
  df2 <- group_by(df, c(subject)) %>% mutate(percent = frequency/sum(frequency))
  sum(df2[grep("s17_ANCA_PR3", df2$subject), "percent"])
  #[1] 1
  df2$cluster <- factor(df2$cluster)
  #nrColors <- length(levels(df2$cluster))
  #myColors = rainbow(nrColors)
  myColors <- gg_color_hue(length(unique(df2$cluster)))
  df2$colors <- myColors[factor(df2$cluster)]
  
  df2$subject <- factor(df2$subject, 
                        levels = c("s17_ANCA_PR3", "s302_ANCA_PR3","s52_ANCA_PR3", 
                                   "s72_ANCA_MPO", "s67_ANCA_MPO", "s68_SLE_NONE", "s57_Control_NONE"),
                        labels = c("s1_ANCA_PR3", "s2_ANCA_PR3","s3_ANCA_PR3", 
                                   "s4_ANCA_MPO", "s5_ANCA_MPO", "SLE", "HC"))
  
  # And write with percentages
  df3 <- df2[, c("cluster","subject", "percent")]
  df3 <- df3 %>% pivot_wider(names_from = subject, values_from = percent)
  write.table(as.data.frame(df3), paste(outDir,paste0("clusters_per_subject_percentages_res", res,".txt"), sep="/"), quote=FALSE, sep="\t", row.names = FALSE, col.names = TRUE)
  
  pdf(paste(outDir, paste0("distribution_per_subject_stacked_percentages_cluster_res", res,".pdf"), sep="/"), width = 14, height =14)
  p <- ggplot(df2, aes(x=subject, y=percent, fill=cluster)) +
    geom_bar(stat='identity', position="stack")+
    scale_fill_manual(values = myColors) +
    guides(fill=guide_legend(ncol =1)) +
    theme_cowplot() +
    theme(axis.text.x=element_text(angle=90,hjust=1))
  print(p)
  legend_p <- cowplot::get_legend(p)
  p2 <- ggplot(df2, aes(x=subject, y=percent, fill=cluster)) +
    geom_bar(stat='identity', position="stack")+
    scale_fill_manual(values = myColors) +
    guides(fill=guide_legend(ncol =1)) +
    theme_cowplot() +
    theme(axis.text.x=element_text(angle=90,hjust=1)) +
    NoLegend()
  print(p2)
  # And print the legend
  plot_grid(legend_p)
  dev.off()
  
  # And as TIFF
  ggsave(paste(outDir, paste0("/tiff/distribution_per_subject_stacked_percentages_clusters_res",res,".tiff"), sep="/"), p, scale=2, width=5, height=5, dpi=360,compression="lzw")
  ggsave(paste(outDir, paste0("/tiff/distribution_per_subject_stacked_percentages_clusters_res", res,".NoLegend.tiff"), sep="/"), p2, scale=2, width=5, height=5, dpi=360,compression="lzw")
  ggsave(paste(outDir, paste0("/tiff/distribution_per_subject_stacked_percentages_clusters_res", res,".Legend.tiff"), sep="/"), legend_p, scale=0.5, width=5, height=5, dpi=360,compression="lzw")
  
}

#### Flow diagram using alluvial ####
# Extract the data frame with the clusterings at 5 different resolutions
clustResults <- seurat@meta.data[ , grep("integrated_snn",colnames(seurat@meta.data))]

# Column 'Freq' to store the total count of all cluster combinations
clustResults$Freq = 1
clustResults2D = aggregate(Freq ~ integrated_snn_res.0.3 + integrated_snn_res.0.6 + integrated_snn_res.1 + 
                             integrated_snn_res.1.4, data = clustResults, sum)

cols = rev(rainbow(nrow(clustResults2D), start = 0.1, end = 0.9))
pdf(paste0(outDir, "/Seurat_clustering_flowdiagram.pdf"), height=14, width=28)
alluvial(
  select(clustResults2D, integrated_snn_res.0.3, integrated_snn_res.0.6, integrated_snn_res.1, integrated_snn_res.1.4),
  freq = clustResults2D$Freq,
  col = cols,
  alpha = 0.8,
  gap.width = 0.3,
  cw = 0.05,
  blocks = FALSE,
  axis_labels = sub("integrated_snn_res.","",colnames(clustResults2D)[1:4]),
  cex = 1.5,
  cex.axis = 1.5)  
dev.off()

#### Clustering tree ####

pdf(paste0(outDir, "/Seurat_clustering_tree.pdf"), height=14, width=14)
p <- clustree(seurat, prefix = "integrated_snn_res.", node_text_size = 6, return = "plot",
              edge_width = 2.5, node_size_range = c(10,24))
p <- p + theme(legend.text = element_text(size=16), legend.key.size = unit(0.6, "cm"),
               legend.title = element_text(size=20)) +
  guides(color = guide_legend(override.aes = list(size=5)))
print(p)

# Colour by gene expression
p <- clustree(seurat, prefix = "integrated_snn_res.", node_colour = "CD14", node_colour_aggr = "median",
              node_text_size = 6, return = "plot", node_text_colour =  "white",
              edge_width = 2.5, node_size_range = c(10,24))
p <- p + theme(legend.text = element_text(size=16), legend.key.size = unit(0.6, "cm"),
               legend.title = element_text(size=20)) +
  guides(color = guide_legend(override.aes = list(size=5)))
print(p)

# Overlay clustering tree on UMAP
p <- clustree_overlay(seurat, prefix = "integrated_snn_res.", red_dim = "umap", x_value = "umap1", y_value = "umap2",
                      plot_sides = FALSE, node_size_range = c(6, 18))
p <- p + theme(legend.text = element_text(size=16), legend.key.size = unit(0.6, "cm"),
               legend.title = element_text(size=20),
               axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 16), 
               axis.title.x = element_text(size = 24), axis.title.y = element_text(size = 24)  ) +
  guides(color = guide_legend(override.aes = list(size=3)))
print(p)

dev.off()

#### Find markers per cluster ####
# Take care on which clustering you now take!!
# Could make a function out of this (specifying nr. markers to plot etc.)

# Read the data
dataDir <- "10_Reclustering_MP_and_MC/"
seurat <- readRDS(paste(dataDir, "selectedClusters.subjectIntegrated.reclustered.rds", sep="/"))

# Original (BUT WRONG!!) results, i.e. DefaultAssay was set to 'integrated'
outDir <- paste0(dataDir,"/markers_integrated")
dir.create(outDir)

# 20230207 - changed the directory as I made changes to the code (see mail Perry 20230206)
# - i.e. DefaultAssay set to 'SCT'
outDir <- paste0(dataDir,"/markers_20230207_SCT")
dir.create(outDir)

# 20230221 - changed the directory as I made changes to the code (see mail Perry 20230206)
#          - redid the FindAllMarkers, now using 'latent.vars' and 'test.use = 'LR' '
outDir <- paste0(dataDir,"/markers_20230221_LR")
dir.create(outDir)

for (outDir in c(paste0(dataDir,"/markers_integrated"), paste0(dataDir,"/markers_20230207_SCT"),paste0(dataDir,"/markers_20230221_LR") )){
  for (clust in colnames(seurat@meta.data)[grep("integrated_snn", colnames(seurat@meta.data))]){
    res <- gsub("integrated_snn_","", clust)
    print(outDir)
    
    # find markers for every cluster compared to all remaining cells, both negative and positive (only.pos = FALSE)
    Idents(seurat) <- clust
    
    if ( outDir == paste0(dataDir,"/markers_integrated")){
      DefaultAssay(seurat) <- "integrated"
    } else {
      # Make sure to use 'SCT' and not 'integrated' as DefaultAssay, see mail Perry 06022023, but see above??
      DefaultAssay(seurat) <- "SCT"
    }
    
    print(paste0("The following setting is used for DefaultAssay: ", DefaultAssay(seurat)))
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
      print("PrepSCTFindMarkers is run ...")
      seurat <- PrepSCTFindMarkers(seurat) # Found 7 SCT models. Recorrecting SCT counts using minimum median counts: 3354
      # PrepSCTFindMarkers skips the re-correction the second time..
    }
    
    seurat@meta.data[,clust] <- factor(seurat@meta.data[,clust], levels=as.character(0:(length(unique(seurat@meta.data[,clust]))-1)))
    if (!file.exists(paste(outDir,paste0("clusterMarkers_", res,".txt"), sep="/"))){
      if ( outDir == paste0(dataDir,"/markers_20230221_LR") ){
        print("Using FindAllMarkers with latent.vars=\"orig.ident\" and test.use=\"LR\"")
        seurat.markers <- FindAllMarkers(seurat, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25,
                                         latent.vars = "orig.ident", test.use = 'LR')
      } else {
        print("Using default FindAllMarkers")
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
  }
  
  # ViolinPlots of clusters
  for (clust in colnames(seurat@meta.data)[grep("integrated_snn", colnames(seurat@meta.data))]){
    res <- gsub("integrated_snn_r","R", clust)
    dir.create(paste(outDir,paste0("VlnPlots_", res), sep="/"))
    dir.create(paste(outDir,paste0("FeaturePlots_", res), sep="/"))
    
    Idents(seurat) <- clust
    maxY.RNA <- NULL # see mail Perry 06022023: original code maxY.RNA <- max(seurat@assays$RNA@data)
    # 21022023: Set the axis back in this case !!!
    
    if ( outDir == paste0(dataDir,"/markers_20230221_LR") ){
      seurat.markers <- read.table(paste("10_Reclustering_MP_and_MC/markers_20230221_LR",paste0("clusterMarkers_", res,".txt"), sep="/"), header = TRUE)
    } 
    if ( outDir == paste0(dataDir,"/markers_20230207") ){
      seurat.markers <- read.table(paste("10_Reclustering_MP_and_MC/markers_20230207_SCT",paste0("clusterMarkers_", res,".txt"), sep="/"), header = TRUE)
    } 
    if ( outDir == paste0(dataDir,"/markers_integrated") ){
      seurat.markers <- read.table(paste("10_Reclustering_MP_and_MC/markers_integrated",paste0("clusterMarkers_", res,".txt"), sep="/"), header = TRUE)
    } 
    
    #top5.markers <- seurat.markers %>% group_by(cluster) %>% top_n(n = -5, wt = p_val_adj)  # wt defines ordering, avg_log2FC, or p_val_adj?
    # Hmm, returns more than 5 rows per clusters as there are ties in the 'p_val_adj' ...
    # So, also sort on avg_log2FC ...
    top5.markers <- seurat.markers %>% group_by(cluster) %>% arrange(p_val_adj, desc(avg_log2FC), .by_group=TRUE) %>% top_n(n = 5, wt = avg_log2FC)
    
    # Sort the top5 list on the cluster column
    top5.markers <- dplyr::arrange(top5.markers, cluster)
    
    for ( myClusters in unique(top5.markers$cluster) ){
      if ( outDir == paste0(dataDir,"/markers_integrated")){
        DefaultAssay(seurat) <- "integrated"
      } else {
        # Make sure to use 'SCT' and not 'integrated' as DefaultAssay, see mail Perry 06022023
        DefaultAssay(seurat) <- "SCT"
      }

      print(paste0("VlnPlots: The following setting is used for DefaultAssay: ", DefaultAssay(seurat)))

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
# Clean
rm(seurat, p, p1, p2, p3, pFeature,pViolin,top5.markers, seurat.markers, res, maxY.RNA)

#### Find markers per cluster - UPDATE 20230418 ####
# - Make markerplot only for 0,1,4,8,9,11
# - Take top5 + CD163, MRC1, MARCO and TREM2

# Also the top 5 of the surface protein levels for clusters 0,1,4 and 8

# Take care on which clustering you now take!!
# Could make a function out of this (specifying nr. markers to plot etc.)

# Read the data
dataDir <- "10_Reclustering_MP_and_MC/"
seurat <- readRDS(paste(dataDir, "selectedClusters.subjectIntegrated.reclustered.rds", sep="/"))

seurat@assays$ADT
# Assay data with 137 features for 4194 cells
# First 10 features:
#   anti-human CD86, anti-human CD274 (B7-H1, PD-L1), anti-human CD270 (HVEM, TR2), anti-human CD155 (PVR), anti-human CD112
# (Nectin-2), anti-human CD47, anti-human CD48, anti-human CD40, anti-human CD154, anti-human CD52 

# Original (BUT WRONG!!) results, i.e. DefaultAssay was set to 'integrated'
outDir <- paste0(dataDir,"/markers_integrated_20230418")
dir.create(outDir)

# 20230207 - changed the directory as I made changes to the code (see mail Perry 20230206)
# - i.e. DefaultAssay set to 'SCT'
outDir <- paste0(dataDir,"/markers_20230418_SCT")
dir.create(outDir)

# 20230221 - changed the directory as I made changes to the code (see mail Perry 20230206)
#          - redid the FindAllMarkers, now using 'latent.vars' and 'test.use = 'LR' '
outDir <- paste0(dataDir,"/markers_20230418_LR")
dir.create(outDir)

for (outDir in c(paste0(dataDir,"/markers_integrated_20230418"), paste0(dataDir,"/markers_20230418_SCT"),paste0(dataDir,"/markers_20230418_LR") )){
  print("")
  print(outDir)
  for (clust in colnames(seurat@meta.data)[grep("integrated_snn", colnames(seurat@meta.data))]){
    res <- gsub("integrated_snn_","", clust)
    print(paste0("Now working on resolution ", res))
    
    # find markers for every cluster compared to all remaining cells, both negative and positive (only.pos = FALSE)
    Idents(seurat) <- clust
    
    if ( outDir == paste0(dataDir,"/markers_integrated_20230418")){
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
      if ( outDir == paste0(dataDir,"/markers_20230418_LR") ){
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
    if ( clust == "integrated_snn_res.1"){
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
      
      if ( outDir == paste0(dataDir,"/markers_20230418_LR") ){
        print("... Integrated ADT data: Using FindAllMarkers with latent.vars=\"orig.ident\" and test.use=\"LR\"")
        adt.markers <- FindAllMarkers(adt.integrated, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25,
                                         latent.vars = "orig.ident", test.use = 'LR', assay = "ADT")
      } else {
        print("... Integrated ADT data: Using default FindAllMarkers")
        adt.markers <- FindAllMarkers(adt.integrated, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25 , assay = "ADT")
      }
      top5.adt_markers <- adt.markers %>% group_by(cluster) %>% arrange(p_val_adj, desc(avg_log2FC), .by_group=TRUE) %>% top_n(n = 5, wt = avg_log2FC)
      top5.adt_markers <- dplyr::arrange(top5.adt_markers, cluster)
      
      selClusters.adt <- subset(x = adt.integrated, subset = integrated_snn_res.1 %in% c("0","1","4","8","9","11"))
      selClusters.adt$integrated_snn_res.1 <- factor(selClusters.seurat$integrated_snn_res.1, levels=c("0","1","4","8","9","11"))
      
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
    
    if ( outDir == paste0(dataDir,"/markers_20230221_LR") ){
      seurat.markers <- read.table(paste("10_Reclustering_MP_and_MC/markers_20230221_LR",paste0("clusterMarkers_", res,".txt"), sep="/"), header = TRUE)
    } 
    if ( outDir == paste0(dataDir,"/markers_20230207") ){
      seurat.markers <- read.table(paste("10_Reclustering_MP_and_MC/markers_20230207_SCT",paste0("clusterMarkers_", res,".txt"), sep="/"), header = TRUE)
    } 
    if ( outDir == paste0(dataDir,"/markers_integrated") ){
      seurat.markers <- read.table(paste("10_Reclustering_MP_and_MC/markers_integrated",paste0("clusterMarkers_", res,".txt"), sep="/"), header = TRUE)
    } 
    
    #top5.markers <- seurat.markers %>% group_by(cluster) %>% top_n(n = -5, wt = p_val_adj)  # wt defines ordering, avg_log2FC, or p_val_adj?
    # Hmm, returns more than 5 rows per clusters as there are ties in the 'p_val_adj' ...
    # So, also sort on avg_log2FC ...
    top5.markers <- seurat.markers %>% group_by(cluster) %>% arrange(p_val_adj, desc(avg_log2FC), .by_group=TRUE) %>% top_n(n = 5, wt = avg_log2FC)
    
    # Sort the top5 list on the cluster column
    top5.markers <- dplyr::arrange(top5.markers, cluster)
    
    if ( outDir == paste0(dataDir,"/markers_integrated")){
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
# Clean
rm(seurat, p, p1, p2, p3, pFeature,pViolin,top5.markers, seurat.markers, res, maxY.RNA)

#### Comparison with previous clustering ####
#Code by Perry Moerland - 20210331
library(alluvial)
library(dplyr)
library(Seurat)

outDir <- "10_Reclustering_MP_and_MC"
dir.create(outDir)

seurat <- readRDS("9_kidney_reference_mapping/Mito10.integrated.clustered.annot.kidney_immune_reference.rds")
Idents(seurat) <- "integrated_snn_res.0.6"
seurat <- subset(x = seurat, idents = c(5,7,9,11), slot = 'counts')
clustResults.old <- seurat@meta.data[ , grep("integrated_snn",colnames(seurat@meta.data))]

seurat.new <- readRDS(paste(outDir, "selectedClusters.subjectIntegrated.reclustered.rds", sep="/"))
clustResults.new <- seurat.new@meta.data[ , grep("integrated_snn",colnames(seurat.new@meta.data))]

clustResults <- cbind(integrated_snn_res.0.6_old = clustResults.old$integrated_snn_res.0.6, clustResults.new)

#### Flow diagram using alluvial ####

# Column 'Freq' to store the total count of all cluster combinations
clustResults$Freq = 1
clustResults2D = aggregate(Freq ~ integrated_snn_res.0.6_old + integrated_snn_res.0.3 + integrated_snn_res.0.6 + integrated_snn_res.1 + integrated_snn_res.1.4, data = clustResults, sum)

cols = rev(rainbow(nrow(clustResults2D), start = 0.1, end = 0.9))
pdf(paste0(outDir,"/Seurat_comparison_flowdiagram.pdf"), height=14, width=28)
alluvial(
  select(clustResults2D, integrated_snn_res.0.6_old, integrated_snn_res.0.3, integrated_snn_res.0.6, integrated_snn_res.1, integrated_snn_res.1.4),
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
# Warning messages:
# 1: In min(x) : no non-missing arguments to min; returning Inf
# 2: In max(x) : no non-missing arguments to max; returning -Inf
# ... 32 in total


#### Clustering tree ####
# Clustree can not handle "06_old"
seurat.new@meta.data$integrated_snn_res.0.29 <- clustResults.old$integrated_snn_res.0.6

# Put the old clustering in front of the new ones
indx <- which(grepl("integrated", colnames(seurat.new@meta.data)))
seurat.new@meta.data[,indx] <- seurat.new@meta.data[,c(indx[length(indx)],indx[1:(length(indx)-1)])]
colnames(seurat.new@meta.data)[indx] <- colnames(seurat.new@meta.data)[c(indx[length(indx)],indx[1:(length(indx)-1)])]

pdf(paste0(outDir, "/Seurat_comparison_clustering_tree.pdf"), height=14, width=14)
p <- clustree(seurat.new, prefix = "integrated_snn_res.", node_text_size = 6, return = "plot",
              edge_width = 2.5, node_size_range = c(10,24))
p <- p + theme(legend.text = element_text(size=16), legend.key.size = unit(0.6, "cm"),
               legend.title = element_text(size=20)) +
  guides(color = guide_legend(override.aes = list(size=5))) +
  scale_color_hue(labels = c("0.6_old", "0.3","0.6",'1.0',"1.4")) # Get the right labels
print(p)

# Colour by gene expression
p <- clustree(seurat.new, prefix = "integrated_snn_res.", node_colour = "CD14", node_colour_aggr = "median",
              node_text_size = 6, return = "plot", node_text_colour =  "white",
              edge_width = 2.5, node_size_range = c(10,24))
p <- p + theme(legend.text = element_text(size=16), legend.key.size = unit(0.6, "cm"),
               legend.title = element_text(size=20)) +
  guides(color = guide_legend(override.aes = list(size=5))) 
print(p)

# Overlay clustering tree on UMAP
p <- clustree_overlay(seurat.new, prefix = "integrated_snn_res.", red_dim = "umap", x_value = "umap1", y_value = "umap2",
                      plot_sides = FALSE, node_size_range = c(6, 18))
p <- p + theme(legend.text = element_text(size=16), legend.key.size = unit(0.6, "cm"),
               legend.title = element_text(size=20),
               axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 16), 
               axis.title.x = element_text(size = 24), axis.title.y = element_text(size = 24)  ) +
  guides(color = guide_legend(override.aes = list(size=3))) +
  scale_fill_discrete(labels = c("0.6_old", "0.3","0.6",'1.0',"1.4")) +  # Get the right labels
  scale_color_discrete(labels= c("0.6_old", "0.3","0.6",'1.0',"1.4"))
print(p)

dev.off()

#### Per clustering resolution and per cluster get a table of how many cells from which sample are present ####

# Read the data
seurat <- readRDS(paste(outDir, "selectedClusters.subjectIntegrated.reclustered.rds", sep="/"))

for (clust in colnames(seurat@meta.data)[grep("integrated_snn", colnames(seurat@meta.data))]){
  res <- gsub("integrated_snn_","", clust)
  
  pdf(paste(outDir, paste0("dimPlot_Clustered_", res,"_MOMA.pdf"), sep="/"))
  Idents(seurat) <- clust
  p <- DimPlot(seurat, reduction = "umap", label = FALSE,  group.by = "sampleAnnot", repel = TRUE) +
    guides(colour=guide_legend(override.aes = list(size=2), ncol=1))
  print(p)

  dev.off()
}

for (clust in colnames(seurat@meta.data)[grep("integrated_snn", colnames(seurat@meta.data))]){
  res <- gsub("integrated_snn_","", clust)
  
  t <- table(seurat@meta.data[,clust], seurat@meta.data$sampleAnnot)
  t.df <- as.data.frame.matrix(t)
  t.df$Cluster = rownames(t.df)
  t.df <- t.df[,c(3,1:2)]
  pct.df <- as.data.frame.matrix(round(proportions(as.matrix(t.df[,2:3]),1) * 100,2))
  colnames(pct.df) <- paste0("pct.",colnames(pct.df))
  t.df <- cbind(t.df, pct.df)
  
  # And normalize over columns
  norm.pct.df <- (scale(t.df[,2:3], center=FALSE, scale=colSums(as.data.frame(t.df[,2:3])))) * 100
  colnames(norm.pct.df) <- paste0("norm.",colnames(pct.df))
  t.df <- cbind(t.df, norm.pct.df)
  t.df[nrow(t.df)+1,] <- c("totals", colSums(as.data.frame(t.df[,2:3])), rep("-",2), colSums(norm.pct.df[,1:2]))
  
  # And save
  write.table(t.df, file=paste(outDir, paste0("dimPlot_Clustered_", res,"_MOMA.txt"), sep="/"), quote = FALSE,
              row.names = FALSE, col.names = TRUE, sep="\t")
  
  # And clean
  rm(t,t.df, pct.df)
}

#### Conserved markers ####
# see: https://satijalab.org/seurat/articles/integration_introduction.html
dataDir <- "10_Reclustering_MP_and_MC"
#dir.create(outDir)

# For performing differential expression after integration, we switch back to the original data
# Should it be the SCT, i.e. normalized data???
seurat <- readRDS(paste(dataDir, "selectedClusters.subjectIntegrated.reclustered.rds", sep="/"))
DefaultAssay(seurat) <- "RNA"

# Create a directory to hold the results
dir.create(paste0(dataDir,"/conserved_markers"))
outDir <- paste0(dataDir,"/conserved_markers")

# Create new directory as I changed some code (see mail Perry 20230206)
# - including 'latent.vars' and 'test.use = "LR"' does not work as for some conditions, you only have one factor
# - We use 'grouping.var = "orig.ident" to look at the conserved markers per sample (in a way trying to handle batch effects)
dir.create(paste0(dataDir,"/conserved_markers_20230221"))
outDir <- paste0(dataDir,"/conserved_markers_20230221")

# Only do 0.3, 0.6, 1.0, 1.4??
for (res in c(0.3, 0.6, 1, 1.4)){
  #res <- gsub("integrated_snn_","", clustRes)
  clustRes <- paste0("integrated_snn_res.", res)
  
  dir.create(paste(outDir,paste0("resolution_", gsub("\\.", "_", res)), sep="/"))
  
  Idents(seurat) <- clustRes
  
  for (cluster in levels(seurat@meta.data[,clustRes])){
    cons.markers <- FindConservedMarkers(seurat, ident.1 = cluster, 
                                          grouping.var = "Condition", verbose = FALSE)
       
    # And save
    write.table(cons.markers, file=paste(paste(outDir,paste0("resolution_", gsub("\\.", "_", res)), sep="/"), paste0("conservedMarkers_cluster", cluster,".Condition.txt"), sep="/"), quote = FALSE,
              row.names = TRUE, col.names = TRUE, sep="\t")
    
    seurat@meta.data$Condition_subType <- factor(paste(seurat@meta.data$Condition,seurat@meta.data$subType, sep="_"))
    
    cons.markers <- FindConservedMarkers(seurat, ident.1 = cluster, 
                                           grouping.var = "Condition_subType", verbose = FALSE)
  
    # And save
    write.table(cons.markers, file=paste(paste(outDir,paste0("resolution_", gsub("\\.", "_", res)), sep="/"), paste0("conservedMarkers_cluster", cluster,".Condition_SubType.txt"), sep="/"), quote = FALSE,
              row.names = TRUE, col.names = TRUE, sep="\t")
    
    cons.markers <- FindConservedMarkers(seurat, ident.1 = cluster, 
                                         grouping.var = "orig.ident", verbose = FALSE)
    
    # And save
    write.table(cons.markers, file=paste(paste(outDir,paste0("resolution_", gsub("\\.", "_", res)), sep="/"), paste0("conservedMarkers_cluster", cluster,".orig_ident.txt"), sep="/"), quote = FALSE,
                row.names = TRUE, col.names = TRUE, sep="\t")
    
  }
  
  # Feature plots?
  # How to select genes that stick out?
}  

# Warning: Identity: 8 not present in group ANCA. Skipping ANCA
# Warning: Identity: 8 not present in group ANCA_MPO. Skipping ANCA_MPO
# Warning: Identity: 8 not present in group ANCA_PR3. Skipping ANCA_PR3
# Warning: ANCA has fewer than 3 cells in Identity: 11. Skipping ANCA
# Warning: ANCA_MPO has fewer than 3 cells in Identity: 11. Skipping ANCA_MPO
# Warning: Identity: 11 not present in group ANCA_PR3. Skipping ANCA_PR3
# Warning: ANCA has fewer than 3 cells in Identity: 12. Skipping ANCA
# Warning: ANCA_MPO has fewer than 3 cells in Identity: 12. Skipping ANCA_MPO
# Warning: Identity: 12 not present in group ANCA_PR3. Skipping ANCA_PR3
# Warning: Identity: 12 not present in group ANCA_PR3. Skipping ANCA_PR3
# Warning: ANCA has fewer than 3 cells in Identity: 14. Skipping ANCA
# Warning: ANCA_MPO has fewer than 3 cells in Identity: 14. Skipping ANCA_MPO
# Warning: Identity: 14 not present in group ANCA_PR3. Skipping ANCA_PR3
# Warning: Identity: 15 not present in group ANCA. Skipping ANCA
# Warning: Identity: 15 not present in group ANCA_MPO. Skipping ANCA_MPO
# Warning: Identity: 15 not present in group ANCA_PR3. Skipping ANCA_PR3
# ....

# Featureplot for resolution 0.3
# Highlight some of the specific markers - these were manual selected (TBD)
Idents(seurat) <- "integrated_snn_res.0.3"

pdf(paste(outDir, "featurePlot_Clustered_0_3.pdf", sep="/")) # see mail Perry 06022023
  DefaultAssay(seurat) <- 'SCT' # see mail Perry 06022023
  p <- FeaturePlot(seurat, features = c("GNLY", "FCER1A", "FCGR3A", "VCAN", "CD3E", "CD2"), min.cutoff = "q9")
  print(p)
  
  p <- FeaturePlot(seurat, features = c("S100A8", "CD36", "EEF1A1", "FCN1","GZMA","LST1"), min.cutoff = "q9")
  print(p)
  
  p <- FeaturePlot(seurat, features = c("SMIM25", "LILRB2", "CD14", "C1QC", "CD9", "HLA-DRA"), min.cutoff = "q9")
  print(p)
  
  p <- FeaturePlot(seurat, features = c("CD69", "PLD4", "CCL5", "IL18","RPL32","RNASE6"), min.cutoff = "q9")
  print(p)
  
dev.off()

# Kidney-resident macrophages
# Stewart et al, 2019, RUNX1 and CD206 (== MRC1, CD206 not found in data...): see mail Yosta 29082022 with marker lists
pdf(paste(outDir, "featurePlot_Clustered_0_3_KidneyResM.pdf", sep="/"))
  p <- FeaturePlot(seurat, features = c("RUNX1","MRC1","SEPP1","CSF2RA","CD74","OR2A25","LINC01644","TRBV10-2"), min.cutoff = "q9")
  print(p)
dev.off()
# Warning messages:
#   1: In FetchData.Seurat(object = object, vars = c(dims, "ident", features),  :
#               The following requested variables were not found: SEPP1, OR2A25
#   2: In FeaturePlot(seurat, features = c("SEPP1", "CSF2RA", "CD74", "OR2A25",  :
#               All cells have the same value (1) of LINC01644.
#   3: In FeaturePlot(seurat, features = c("SEPP1", "CSF2RA", "CD74", "OR2A25",  :
#               All cells have the same value (13) of TRBV10-2.                                   

