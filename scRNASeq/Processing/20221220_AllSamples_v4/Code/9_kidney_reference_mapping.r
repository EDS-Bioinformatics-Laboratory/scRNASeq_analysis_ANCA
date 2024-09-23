# Use Kidney reference sets from Stewart et al

#### Set default for a good start ####
source("initAnalysis.r")
init()

# Create a directory to hold results and plots
outDir <- "9_kidney_reference_mapping"
dir.create(outDir)

#### Library ####
library(alluvial)
library(dplyr)
library(ggplot2)
library(Seurat)
library(SeuratDisk) 
#library(zellkonverter)
library(tidyverse)

dataDir <- "8_clustering//"

# Coloring function
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

#### Data ####
# See ../../20220623_AllSamples_InitialAnalysis/Code/8_kidney_reference_mapping.r
# for the conversion of the data from Stewart et al. (https://science.sciencemag.org/content/365/6460/1461) 
# into something useable 

#### Kidney Mature Immune ####
# Copied from ../../20220623_AllSamples_InitialAnalysis/Code/8_kidney_reference_mapping
if(file.exists(paste(outDir,"kidney_immune_reference.rds", sep="/"))){
  reference <- readRDS(paste(outDir,"kidney_immune_reference.rds",sep="/"))
} else {
  reference <- readRDS(paste("../../20220623_AllSamples_InitialAnalysis/Code/8_kidney_reference_mapping","kidney_immune_reference.rds",sep="/"))
  saveRDS(reference, paste(outDir,"kidney_immune_reference.rds",sep="/"))
}

# Make a plot of the reference map
pdf(paste(outDir, "umap_kidney_immune_reference.pdf", sep="/"), width=14)
p <- DimPlot(reference, reduction = "umap", group.by = "celltype", label = TRUE, label.size = 3,
        repel = TRUE) + ggtitle("Reference annotations kidney immune set")
print(p)
dev.off()

# Import Yosta's data from 8_clustering
query1 <- readRDS(paste(dataDir,"Mito10.integrated.clustered.rds", sep="/"))
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

# Check whether there already is a 'umap' slot
if(!("^umap" %in% names(reference@reductions))){
  reference <- RunUMAP(reference, dims = 1:30, reduction = "pca", return.model = TRUE)
}
query1 <- MapQuery(anchorset = kidney.anchors, reference = reference, query = query1,
                   refdata = list(celltype = "celltype"), reference.reduction = "pca", reduction.model = "umap")
# NOTE:
# Predicting cell labels
#Error in x$.self$finalize() : attempt to apply non-function
#|                                                  | 0 % ~calculating  


# And save
saveRDS(query1, paste(outDir,"Mito10.integrated.clustered.annot.kidney_immune_reference.rds", sep="/"))
saveRDS(reference, paste(outDir,"kidney_immune_reference.rds", sep="/"))

# Read the data
reference <- readRDS(paste(outDir,"kidney_immune_reference.rds", sep="/"))
query1 <- readRDS(paste(outDir,"Mito10.integrated.clustered.annot.kidney_immune_reference.rds", sep="/"))

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


# flow diagram
clustResults <- query1@meta.data[ , c("predicted.celltype", "integrated_snn_res.0.3", "integrated_snn_res.0.6")]

# Column 'Freq' to store the total count of all combinations
clustResults$Freq <- 1
clustResults2D <- aggregate(Freq ~ ., data = clustResults, sum)

cols <- rev(rainbow(nrow(clustResults2D), start = 0.1, end = 0.9))
pdf(paste(outDir,"kidney_immune_celltype_clustering_flowdiagram.pdf", sep="/"), height=14, width=14)
alluvial(
  select(clustResults2D, predicted.celltype, integrated_snn_res.0.3, integrated_snn_res.0.6),
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
print(table(predictions.compartment$predicted.id, query1$integrated_snn_res.0.3))
#             0    1    2    3    4    5    6    7    8    9   10   11
# lymphoid 5754 5053 3878 3049   87 1275  131 1109  613  726   84    0
# myeloid     2    4    7   34 2141    4 1031    0  426   14    0   63

print(table(predictions.compartment$predicted.id, query1$integrated_snn_res.0.6))
#             0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16   17   18   19
# lymphoid 3880 3165 3044 3036 2375  106 1170  105 1115    4  971  413  640  646  639  166  114  109    0   61
# myeloid     7    1   34    1    1 1175    0 1016    0  980    0  395   14    2    4   33    0    0   63    0

# Clean as you get the dreaded 'Error: cannot allocate vector of size .... '
rm(query1, reference, matureImmuneKidney.seurat, p1, p2, predictions, predictions.compartment, sce, counts, clustResults, clustResults2D, kidney.anchors)
gc()

#### Kidney Mature Full ####
# 
if(file.exists(paste(outDir,"kidney_full_reference.rds", sep="/"))){
  reference <- readRDS(paste(outDir,"kidney_full_reference.rds",sep="/"))
} else {
  reference <- readRDS(paste("../../20220623_AllSamples_InitialAnalysis/Code/8_kidney_reference_mapping","kidney_full_reference.rds",sep="/"))
  saveRDS(reference,paste(outDir,"kidney_full_reference.rds",sep="/") )
}

# Get common colors
refColors <- gg_color_hue(length(levels(factor(reference@meta.data$label.main))))

pdf(paste(outDir,"umap_kidney_full_reference.pdf", sep="/"), width=14)
p <- DimPlot(reference, reduction = "umap", group.by = "label.main", label = TRUE, label.size = 3,
        repel = TRUE, cols = refColors) + ggtitle("Reference annotations kidney full set")
print(p)
dev.off()

# Import Yosta's data from 8_clustering
query1 <- readRDS(paste(dataDir,"Mito10.integrated.clustered.rds", sep="/"))
if(file.exists(paste(outDir,"kidney_full_anchors.rds", sep="/"))){
  kidney.anchors <- readRDS(paste(outDir,"kidney_full_anchors.rds", sep="/"))
} else {
  kidney.anchors <- FindTransferAnchors(reference = reference, query = query1, dims = 1:30, reference.reduction = "pca")
  saveRDS(kidney.anchors, paste(outDir,"kidney_full_anchors.rds", sep="/"))
}

predictions <- TransferData(anchorset = kidney.anchors, refdata = reference$label.main, dims = 1:30)
query1 <- AddMetaData(query1, metadata = predictions)

write.table(table(predictions$predicted.id, query1$integrated_snn_res.0.3),paste(outDir,"kidney_full_query_cluster_res_0.3.txt", sep="/"), sep = "\t", quote = FALSE, col.names = NA)
write.table(table(predictions$predicted.id, query1$integrated_snn_res.0.6),paste(outDir,"kidney_full_query_cluster_res_0.6.txt", sep="/"), sep = "\t", quote = FALSE, col.names = NA)
write.table(table(predictions$predicted.id, query1$integrated_snn_res.1),paste(outDir,"kidney_full_query_cluster_res_1.txt", sep="/"), sep = "\t", quote = FALSE, col.names = NA)

if(!("^umap" %in% names(reference@reductions))){
  reference <- RunUMAP(reference, dims = 1:30, reduction = "pca", return.model = TRUE)
}
query1 <- MapQuery(anchorset = kidney.anchors, reference = reference, query = query1,
                   refdata = list(celltype = "label.main"), reference.reduction = "pca", reduction.model = "umap")

# And save
saveRDS(query1, paste(outDir,"Mito10.integrated.clustered.annot.kidney_full_reference.rds", sep="/"))
saveRDS(reference, paste(outDir,"kidney_full_reference.rds", sep="/"))

# Read the data
reference <- readRDS(paste(outDir,"kidney_full_reference.rds", sep="/"))
query1 <- readRDS(paste(outDir,"Mito10.integrated.clustered.annot.kidney_full_reference.rds", sep="/"))

# Make the levels order the same
query1@meta.data$predicted.celltype <- factor(query1@meta.data$predicted.celltype)

# Get common colors
refColors <- gg_color_hue(length(levels(factor(reference@meta.data$label.main))))
names(refColors) <- levels(factor(reference@meta.data$label.main))

queryColors <- refColors[names(refColors) %in% levels(query1@meta.data$predicted.celltype)]

pdf(paste(outDir,"umap_full_reference_query.pdf", sep="/"), width=14)
p1 <- DimPlot(reference, reduction = "umap", group.by = "label.main", label = TRUE, label.size = 3,
              repel = TRUE, cols = refColors) + NoLegend() + ggtitle("Reference annotations")
p2 <- DimPlot(query1, reduction = "ref.umap", group.by = "predicted.celltype", label = TRUE,
              label.size = 3, repel = TRUE, cols = queryColors) + NoLegend() + ggtitle("Query transferred labels")
p <- p1 + p2
print(p)
dev.off()

pdf(paste(outDir,"umap_query_kidney_full.pdf", sep="/"))
p <- DimPlot(query1, reduction = "ref.umap", group.by = "predicted.celltype", label = TRUE, pt.size = 1.7,
        label.size = 5, repel = TRUE, cols = queryColors) + NoLegend() + ggtitle("Query transferred labels: full kidney set") #+
  # xlim(-10,16) + ylim(-17,10) 
print(p)
dev.off()


pdf(paste(outDir,"umap_query_kidney_full_cluster_res_0.3.pdf", sep="/"), width=14)
p1 <- DimPlot(query1, reduction = "ref.umap", group.by = "predicted.celltype", label = TRUE,
              label.size = 3, repel = TRUE, cols = queryColors) + NoLegend() + ggtitle("Query transferred labels")
p2 <- DimPlot(object = query1, reduction = "ref.umap", group.by = "integrated_snn_res.0.3", label = TRUE,
              label.size = 5, repel = TRUE) + NoLegend() + ggtitle("Cluster labels (resolution = 0.3)")
p <- p1 + p2
print(p)
dev.off()

pdf(paste(outDir,"umap_query_kidney_full_cluster_res_0.6.pdf", sep="/"), width=14)
p1 <- DimPlot(query1, reduction = "ref.umap", group.by = "predicted.celltype", label = TRUE,
              label.size = 3, repel = TRUE, cols = queryColors) + NoLegend() + ggtitle("Query transferred labels")
p2 <- DimPlot(object = query1, reduction = "ref.umap", group.by = "integrated_snn_res.0.6", label = TRUE,
              label.size = 5, repel = TRUE) + NoLegend() + ggtitle("Cluster labels (resolution = 0.6)")
p <- p1 + p2
print(p)
dev.off()

pdf(paste(outDir,"umap_query_kidney_full_cluster_res_1.pdf", sep="/"), width=14)
p1 <- DimPlot(query1, reduction = "ref.umap", group.by = "predicted.celltype", label = TRUE,
              label.size = 3, repel = TRUE, cols = queryColors) + NoLegend() + ggtitle("Query transferred labels")
p2 <- DimPlot(object = query1, reduction = "ref.umap", group.by = "integrated_snn_res.1", label = TRUE,
              label.size = 5, repel = TRUE) + NoLegend() + ggtitle("Cluster labels (resolution = 1)")
p <- p1 + p2
print(p)
dev.off()

# flow diagram
clustResults <- query1@meta.data[ , c("predicted.celltype", "integrated_snn_res.0.3", "integrated_snn_res.0.6")]

# Column 'Freq' to store the total count of all combinations
clustResults$Freq <- 1
clustResults2D <- aggregate(Freq ~ ., data = clustResults, sum)

cols <- rev(rainbow(nrow(clustResults2D), start = 0.1, end = 0.9))
pdf(paste(outDir,"kidney_full_celltype_clustering_flowdiagram.pdf", sep="/"), height=14, width=14)
alluvial(
  select(clustResults2D, predicted.celltype, integrated_snn_res.0.3, integrated_snn_res.0.6),
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

# Predict compartment: {lymphoid,myeloid, non_PT, PT }
predictions.compartment <- TransferData(anchorset = kidney.anchors, refdata = reference$label.compartment, dims = 1:30)
print(table(predictions.compartment$predicted.id, query1$integrated_snn_res.0.3))
#               0    1    2    3    4    5    6    7    8    9   10   11
#   lymphoid 5740 5052 3866 3059   62 1273   26 1108  459  725   74    8
#   myeloid     7    4   19   24 2164    6 1126    0  507   15   10   54
#   non_PT      0    0    0    0    0    0    0    0    7    0    0    0
#   PT          9    1    0    0    2    0   10    1   66    0    0    1

print(table(predictions.compartment$predicted.id, query1$integrated_snn_res.0.6))
#             0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16   17   18   19
# lymphoid 3867 3155 3054 3034 2372   67 1169   19 1114    0  971  348  643  644  640  157   86  104    8    0
# myeloid    20    6   24    1    2 1212    0 1092    0  984    0  458   11    4    3   36   23    5   54    1
# non_PT      0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    7
# PT          0    5    0    2    2    2    1   10    1    0    0    2    0    0    0    6    5    0    1   53

# And clean 
rm(query1, reference, p1, p2, predictions, predictions.compartment, clustResults, clustResults2D, kidney.anchors)
gc()
