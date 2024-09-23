#### Set default for a good start ####
source("initAnalysis.r")
init()

# Create a directory to hold results and plots
outDir <- "finalFigures"
dir.create(outDir)

dir.create(paste0(outDir,"/20230606_Wensenlijst"))
figDir <- paste0(outDir,"/20230606_Wensenlijst")

# And create a directory for holding the TIFFs?
dir.create(paste0(figDir,"/tiff"))

#### Markerplot @res1.0 ####
# See '10_reclustering_MP_and_MC.r' - line 1139
#
# Perhaps interesting?
# - https://github.com/satijalab/seurat/issues/2201
library(Seurat)
library(dplyr)
library(tidyr)
library(cowplot)
library(ggplot2)

# Color function
# Also get the correct colors
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

# Read the data
dataDir <- "10_Reclustering_MP_and_MC/"
seurat <- readRDS(paste(dataDir, "selectedClusters.subjectIntegrated.reclustered.rds", sep="/"))

# Set the default assay slot correct
DefaultAssay(seurat) <- "SCT"

# We are going to show the marker plot for resolution 1 and only cluster 0, 1, 4, 8, 9 and 11
# and thus also only the top5 marker for these clusters and add some other ones
# See mails 20 November 2023: Perry noted that some genes (C1QB etc.) were actually missing from the heatmap
# It turned out that there were two 'clusterMarkers'. file for resolution 1 (check in 10_Reclustering_MP_and_MC.r' why that would be so and which is now the correct/latest one)!

# Below command leads to pdf(paste(figDir, "markerPlots_Top5_res1_LR.selectedClusters.pdf", sep="/"), width = 28, height = 28)
#seurat.markers <- read.table(paste(dataDir,paste0("markers_20230418_LR/clusterMarkers_1.txt"), sep="/"), header = TRUE)

# This commands leads pdf(paste(figDir, "markerPlots_Top5_res1_LR.selectedClusters.20231121.pdf", sep="/"), width = 28, height = 28)
seurat.markers <- read.table(paste(dataDir,paste0("markers_20230418_LR/clusterMarkers_res.1.txt"), sep="/"), header = TRUE)

# Give top 5 per comparison based on p_adj_val (need top_n(n=-5) to get the right 5 genes!!)
#top5.markers <- seurat.markers %>% group_by(cluster) %>% top_n(n = -5, wt = p_val_adj)  # wt defines ordering, avg_log2FC, or p_val_adj?
# Hmm, returns more than 5 rows per clusters as there are ties in the 'p_val_adj' ...
# So, also sort on avg_log2FC ...
top5.markers <- seurat.markers %>% group_by(cluster) %>% dplyr::arrange(p_val_adj, desc(avg_log2FC), .by_group=TRUE) %>% top_n(n = 5, wt = avg_log2FC)

# Select only the genes in the cluster 0,1,4,8,9,11
top5.markers <- top5.markers[top5.markers$cluster %in% c(0,1,4,8,9,11),]

# Now get the data for CD163, MRC1, MARCO, TREM2, CD14, PLIN2, C1QC and SPP1 (if they are not already in there)
# AJ 20230606 - Not C1QA as previously??
# 20231121 - save them first in a seprate dataframe, so you can select the entry with the highest avg_logFC
for (markerGene in c("CD163", "MRC1", "MARCO", "TREM2", "PLIN2", "CD14", "C1QC", "SPP1")){
  if ( !(markerGene %in% top5.markers$gene)){
    print(paste0("... Now adding ", markerGene, " to the list of marker genes"))
    if (exists("selectedGenes.df")){
      selectedGenes.df <- rbind(selectedGenes.df, seurat.markers[seurat.markers$gene == markerGene, ])
    } else {
      selectedGenes.df <- seurat.markers[seurat.markers$gene == markerGene, ]
    }
  }
}

# Note that MARCO is not in the final plot because it is not in the list of marker genes

# Now, select only those entries that are in the desired clusters
selectedGenes.df <- selectedGenes.df[selectedGenes.df$cluster %in% c(0,1,4,8,9,11),]

# And get the entry with the highest avg_log2FC
selectedGenes.df <- selectedGenes.df %>% group_by(gene) %>% dplyr::arrange(desc(avg_log2FC), .by_group=TRUE) %>% top_n(n = 1, wt = avg_log2FC)

# And now add these to the top5 list
# First remove them, if they are present, to avoid duplicates?
top5.markers <- top5.markers[!(top5.markers$gene %in% selectedGenes.df$gene), ]
top5.markers <- rbind(top5.markers, selectedGenes.df)
rm(selectedGenes.df)

# Remove duplicate entries of genes based on the avg_log2FC?
top5.markers <- top5.markers %>% group_by(gene) %>% dplyr::arrange(top5.markers, avg_log2FC, .by_group = TRUE) %>% top_n(n = 1, wt = avg_log2FC)

# Sort the top5 list on the cluster column
top5.markers <- dplyr::arrange(top5.markers, cluster) 

# Plot a heatmap and dotplot
# resolution 1 -> only get clusters 0,1,4,8,9,11
clust = "integrated_snn_res.1"
selClusters.seurat <- subset(x = seurat, subset = integrated_snn_res.1 %in% c("0","1","4","8","9","11"))
selClusters.seurat$integrated_snn_res.1 <- factor(selClusters.seurat$integrated_snn_res.1, levels=c("0","1","4","8","9","11"))

myColors <- gg_color_hue(length(unique(seurat$integrated_snn_res.1)))
myColors = myColors[c(0,1,4,8,9,11)+1]

#pdf(paste(figDir, "markerPlots_Top5_res1_LR.selectedClusters.pdf", sep="/"), width = 28, height = 28)
pdf(paste(figDir, "markerPlots_Top5_res1_LR.selectedClusters.20231121.pdf", sep="/"), width = 28, height = 28)
  # Perhaps downsample to 100,200 cells??
  p1 <- DoHeatmap(selClusters.seurat, assay = "SCT", features = top5.markers$gene, group.by = clust, 
               size = 28, angle=0, hjust=0.5, group.colors = myColors) + 
    theme(axis.text.y = element_text(size=28, face = "bold.italic")) +
    theme(legend.key.size = unit(6, 'cm'), 
          legend.key.height= unit(4, 'cm'),
          legend.key.width= unit(4, 'cm'),
          legend.title = element_text(size=40),legend.text = element_text(size=32) ) +
    guides(colour = guide_legend(override.aes = list(size=30, alpha =1))) 
  # Warning message:
  #   In DoHeatmap(selClusters.seurat, features = top5.markers$gene, group.by = clust,  :
  #                  The following features were omitted as they were not found in the scale.data slot for the SCT assay: FN1, IGLC2, LYPD2
  # Get the legend
  legend_p <- cowplot::get_legend(p1)
  print(p1)
  # Now plot it without legend
  p2 <- DoHeatmap(selClusters.seurat, features = top5.markers$gene, group.by = clust, 
               size = 28, angle=0, hjust=0.5, group.colors = myColors) + 
    theme(axis.text.y = element_text(size=28, face = "bold.italic")) +
    NoLegend()
  print(p2)
  # And plot the legend
  plot_grid(legend_p)

dev.off()

# And as TIFF
# ggsave(paste(figDir, "/tiff/markerPlots_Top5_res1_LR.selectedClusters.tiff", sep="/"), p1, scale=4, width=5, height=5, dpi=360,compression="lzw")
# ggsave(paste(figDir, "/tiff/markerPlots_Top5_res1_LR.selectedClusters.NoLegend.tiff", sep="/"), p2, scale=4, width=5, height=5, dpi=360,compression="lzw")
# ggsave(paste(figDir, "/tiff/markerPlots_Top5_res1_LR.selectedClusters.LegendOnly.tiff", sep="/"), legend_p, scale=4, width=5, height=5, dpi=360,compression="lzw")
ggsave(paste(figDir, "/tiff/markerPlots_Top5_res1_LR.selectedClusters.20231121.tiff", sep="/"), p1, scale=4, width=5, height=5, dpi=360,compression="lzw")
ggsave(paste(figDir, "/tiff/markerPlots_Top5_res1_LR.selectedClusters.NoLegend.20231121.tiff", sep="/"), p2, scale=4, width=5, height=5, dpi=360,compression="lzw")
ggsave(paste(figDir, "/tiff/markerPlots_Top5_res1_LR.selectedClusters.LegendOnly.20231121.tiff", sep="/"), legend_p, scale=4, width=5, height=5, dpi=360,compression="lzw")

# Why does FTL end up so low?
# - It is a marker for cluster 9 ...
#   Together with IGLC2 (but where is that one ???)
#.  And also LYPD2 and FN1 are lacking ...
"IGLC2" %in% rownames(selClusters.seurat@assays$SCT$counts)
#[1] TRUE
"IGLC2" %in% rownames(selClusters.seurat@assays$SCT$data)
#[1] TRUE
"IGLC2" %in% rownames(selClusters.seurat@assays$SCT$scale.data)
#[1] FALSE
top5.markers$gene %in% rownames(selClusters.seurat@assays$SCT$scale.data)
#[1]  TRUE  TRUE  TRUE  TRUE FALSE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
#[26]  TRUE  TRUE FALSE  TRUE FALSE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
top5.markers$gene
#[1] "LST1"     "FCGR3A"   "SMIM25"   "SERPINA1" "LYPD2"    "S100A9"   "S100A8"   "VCAN"     "S100A12"  "CD14"     "CD14"     "EREG"     "C1QC"     "C1QA"    
#[15] "C1QB"     "APOE"     "TREM2"    "TREM2"    "RNASE1"   "CD163"    "CD163"    "MRC1"     "MRC1"     "FABP5"    "C15orf48" "SPP1"     "APOC1"    "FN1"     
#[29] "FTL"      "IGLC2"    "HSP90AA1" "DNAJB1"   "HSPA1A"   "HSPA1B"   "PLIN2"    "PLIN2"    "HSPB1"   
# All genes are present in the 'data' slot

# DoHeatmap takes the 'scale.data'. slot of the Seurat object

# NOTE
# Compare paste(dataDir,paste0("markers_20230418_LR/clusterMarkers_res.1.txt"), sep="/") and
# paste(dataDir,paste0("markers_20230221_LR/clusterMarkers_res.1.txt"), sep="/"), just to make sure
library(tools)
md5sum(paste(dataDir,paste0("markers_20230418_LR/clusterMarkers_res.1.txt"), sep="/")) == md5sum(paste(dataDir,paste0("markers_20230221_LR/clusterMarkers_res.1.txt"), sep="/"))
# 10_Reclustering_MP_and_MC//markers_20230418_LR/clusterMarkers_res.1.txt 
#                                                                    TRUE

# NOTE
# We can rename the clusters to the cell types using 'levels()' if wanted
# And then perhaps rotate the labels??
# But they might become large etc. so perhaps do this manual later??

# Just try:
p1 <- DoHeatmap(selClusters.seurat, features = top5.markers$gene, group.by = clust,
                size = 28, label=TRUE, group.colors = myColors, group.bar.height = 0.02) + 
  theme(axis.text.y = element_text(size=28, face = "bold.italic")) +
  theme(legend.key.size = unit(6, 'cm'), 
        legend.key.height= unit(4, 'cm'),
        legend.key.width= unit(4, 'cm'),
        legend.title = element_text(size=40),legend.text = element_text(size=32) ) +
  guides(colour = guide_legend(override.aes = list(size=30, alpha =1)))

pbuild <- ggplot_build(plot = p1)

clusterLabels.df <- pbuild$data[[4]]
# Replace the cluster numbers by their celltypes
require(plyr)
clusterLabels.df$label <- mapvalues(clusterLabels.df$label, 
                               from=c("0","1","4","8","9","11"), 
                               to=c("Non-classical MDMs","Classical MDMs", "Res-like C1q Mac",
                                    "SPP1 LAMs","Overlap Mac", "Stress response \n non-classical MDMs"))

p1 <- DoHeatmap(selClusters.seurat, features = top5.markers$gene, group.by = clust,
                size = 28, label=FALSE, group.colors = myColors, group.bar.height = 0.02) +
  theme(axis.text.y = element_text(size=28, face = "bold.italic")) +
  theme(legend.key.size = unit(6, 'cm'), 
        legend.key.height= unit(4, 'cm'),
        legend.key.width= unit(4, 'cm'),
        legend.title = element_text(size=40),legend.text = element_text(size=32) ) +
  guides(colour = guide_legend(override.aes = list(size=30, alpha =1)))

# But now there is no space reserved for the labels ...

group.bar.height = 0.02

y.range <- diff(x = pbuild$layout$panel_params[[1]]$y.range)
y.pos <- max(pbuild$layout$panel_params[[1]]$y.range) + y.range * 0.015
y.max <- y.pos + group.bar.height * y.range

p1 <- p1 + 
  theme(plot.margin = margin(12,1,1,1, "cm")) +
  geom_text(stat = "identity", data=clusterLabels.df, 
            aes_string(label="label", x="x"), 
            y = y.max - 5 , angle=60, size=14, 
            hjust=0, vjust=1, lineheight = 0.6, color="black") 

pdf(paste(figDir, "markerPlots_Top5_res1_LR.selectedClusters.Alternative.pdf", sep="/"), width = 28, height = 28)
  print(p1)
dev.off()

# The celltypes are cumbersome ... "Stress response non-classical MDMs" is too long
# Putting in a line break does help in combination with the 'lineheight' specification in 'geom_text'
# But you still have a lot of whitespace, so it might be better to use one legend with the color coding?
# You also need to increase the 'margin' using 'plot.margin'

#### UMAP CD45+ cells ####
# Use celltypes instead of numbers ... T cells/ B cells/ myeloid cells/ granulocytes etc.
# The figure corresponds to 9_kidney_reference_mapping/umap_query_kidney_immune_cluster_res_0.6.pdf

# Read the data
dataDir <- "9_kidney_reference_mapping"
reference <- readRDS(paste(dataDir,"kidney_immune_reference.rds", sep="/"))
query1 <- readRDS(paste(dataDir,"Mito10.integrated.clustered.annot.kidney_immune_reference.rds", sep="/"))

# Get common colors
refColors <- gg_color_hue(length(levels(factor(reference@meta.data$celltype))))

# Make the levels order the same
query1@meta.data$predicted.celltype <- factor(query1@meta.data$predicted.celltype, levels = levels(factor(reference@meta.data$celltype)))

pdf(paste(figDir,"umap_query_kidney_immune_cluster_res_0.6.pdf", sep="/"), width=14)
# Original code 9_kidney_reference_mapping.r
# p1 <- DimPlot(query1, reduction = "ref.umap", group.by = "predicted.celltype", label = TRUE,
#               label.size = 3, repel = TRUE, cols = queryColors) + NoLegend() + ggtitle("Query transferred labels")
# p2 <- DimPlot(object = query1, reduction = "ref.umap", group.by = "integrated_snn_res.0.6", label = TRUE,
#               label.size = 5, repel = TRUE) + NoLegend() + ggtitle("Cluster labels (resolution = 0.6)")
#p <- p1 + p2
  p1 <- DimPlot(object = query1, reduction = "ref.umap", group.by = "integrated_snn_res.0.6", label = FALSE,
                label.size = 8, repel = TRUE) + NoLegend() + ggtitle("Cluster labels (resolution = 0.6)")
  p1 <- LabelClusters(p1, id = "integrated_snn_res.0.6",  fontface = "bold")
  p2 <- DimPlot(object = query1, reduction = "ref.umap", group.by = "predicted.celltype", label = FALSE,
                label.size = 6, repel = TRUE) + NoLegend() + ggtitle("Cluster labels\n(resolution = 0.6, query transferred labels)") 
  p2 <- LabelClusters(p2, id = "predicted.celltype",  fontface = "bold")
  p <- p1 + p2
  print(p)
  print(p1)
  print(p2)
  
  # And make them without labels and get the legends
  p1b <- DimPlot(object = query1, reduction = "ref.umap", group.by = "integrated_snn_res.0.6", label = FALSE,
                label.size = 8, repel = TRUE) + ggtitle("Cluster labels (resolution = 0.6)")
  print(p1b)
  # Get the legend
  legend_p1 <- cowplot::get_legend(p1b)
  p1c <- DimPlot(object = query1, reduction = "ref.umap", group.by = "integrated_snn_res.0.6", label = FALSE,
                label.size = 8, repel = TRUE) +  NoLegend() + ggtitle("Cluster labels (resolution = 0.6)")
  print(p1c)
  plot_grid(legend_p1)
  #
  p2b <- DimPlot(object = query1, reduction = "ref.umap", group.by = "predicted.celltype", label = FALSE,
                label.size = 6, repel = TRUE) + ggtitle("Cluster labels\n(resolution = 0.6, query transferred labels)") 
  print(p2b)
  legend_p2 <- cowplot::get_legend(p2b)
  p2c <- DimPlot(object = query1, reduction = "ref.umap", group.by = "predicted.celltype", label = FALSE,
                label.size = 6, repel = TRUE) + NoLegend() + ggtitle("Cluster labels\n(resolution = 0.6, query transferred labels)") 
  print(p2c)
  plot_grid(legend_p2)
  
dev.off()

# And as TIFF
ggsave(paste(figDir, "/tiff/umap_query_kidney_immune_cluster_res_0.6.all.tiff", sep="/"), p, scale=2, width=5, height=5, dpi=360,compression="lzw")
ggsave(paste(figDir, "/tiff/umap_query_kidney_immune_cluster_res_0.6.clusters.tiff", sep="/"), p1, scale=1.5, width=5, height=5, dpi=360,compression="lzw")
ggsave(paste(figDir, "/tiff/umap_query_kidney_immune_cluster_res_0.6.transferredLabels.tiff", sep="/"), p2, scale=1.5, width=5, height=5, dpi=360,compression="lzw")
ggsave(paste(figDir, "/tiff/umap_query_kidney_immune_cluster_res_0.6.clusters_NoLabel.tiff", sep="/"), p1b, scale=1.5, width=5, height=5, dpi=360,compression="lzw")
ggsave(paste(figDir, "/tiff/umap_query_kidney_immune_cluster_res_0.6.clusters_NoLabel_NoLegend.tiff", sep="/"), p1c, scale=1.5, width=5, height=5, dpi=360,compression="lzw")
ggsave(paste(figDir, "/tiff/umap_query_kidney_immune_cluster_res_0.6.clusters_Legend.tiff", sep="/"), legend_p1, scale=1.25, width=5, height=5, dpi=360,compression="lzw")
ggsave(paste(figDir, "/tiff/umap_query_kidney_immune_cluster_res_0.6.transferredLabels_NoLabel.tiff", sep="/"), p2b, scale=1.5, width=5, height=5, dpi=360,compression="lzw")
ggsave(paste(figDir, "/tiff/umap_query_kidney_immune_cluster_res_0.6.transferredLabels_NoLabel_NoLegend.tiff", sep="/"), p2c, scale=1.5, width=5, height=5, dpi=360,compression="lzw")
ggsave(paste(figDir, "/tiff/umap_query_kidney_immune_cluster_res_0.6.transferredLabels_Legend.tiff", sep="/"), legend_p2, scale=0.75, width=5, height=5, dpi=360,compression="lzw")

pdf(paste(figDir,"umap_query_kidney_immune_cluster_res_0.6_v2.pdf", sep="/"))
# Original code 9_kidney_reference_mapping.r
# p1 <- DimPlot(query1, reduction = "ref.umap", group.by = "predicted.celltype", label = TRUE,
#               label.size = 3, repel = TRUE, cols = queryColors) + NoLegend() + ggtitle("Query transferred labels")
# p2 <- DimPlot(object = query1, reduction = "ref.umap", group.by = "integrated_snn_res.0.6", label = TRUE,
#               label.size = 5, repel = TRUE) + NoLegend() + ggtitle("Cluster labels (resolution = 0.6)")
#p <- p1 + p2
p1 <- DimPlot(object = query1, reduction = "ref.umap", group.by = "integrated_snn_res.0.6", label = FALSE,
              label.size = 8, repel = TRUE) + NoLegend() + ggtitle("Cluster labels (resolution = 0.6)")
p1 <- LabelClusters(p1, id = "integrated_snn_res.0.6",  fontface = "bold")
p2 <- DimPlot(object = query1, reduction = "ref.umap", group.by = "predicted.celltype", label = FALSE,
              label.size = 6, repel = TRUE) + NoLegend() + ggtitle("Cluster labels\n(resolution = 0.6, query transferred labels)") 
p2 <- LabelClusters(p2, id = "predicted.celltype",  fontface = "bold")
#p <- p1 + p2
#print(p)
print(p1)
print(p2)

# And make them without labels and get the legends
p1b <- DimPlot(object = query1, reduction = "ref.umap", group.by = "integrated_snn_res.0.6", label = FALSE,
               label.size = 8, repel = TRUE) + ggtitle("Cluster labels (resolution = 0.6)")
print(p1b)
# Get the legend
legend_p1 <- cowplot::get_legend(p1b)
p1c <- DimPlot(object = query1, reduction = "ref.umap", group.by = "integrated_snn_res.0.6", label = FALSE,
               label.size = 8, repel = TRUE) +  NoLegend() + ggtitle("Cluster labels (resolution = 0.6)")
print(p1c)
plot_grid(legend_p1)
#
p2b <- DimPlot(object = query1, reduction = "ref.umap", group.by = "predicted.celltype", label = FALSE,
               label.size = 6, repel = TRUE) + ggtitle("Cluster labels\n(resolution = 0.6, query transferred labels)") 
print(p2b)
legend_p2 <- cowplot::get_legend(p2b)
p2c <- DimPlot(object = query1, reduction = "ref.umap", group.by = "predicted.celltype", label = FALSE,
               label.size = 6, repel = TRUE) + NoLegend() + ggtitle("Cluster labels\n(resolution = 0.6, query transferred labels)") 
print(p2c)
plot_grid(legend_p2)

dev.off()

#### UMAP @res1.0 after reclustering ####
# Use celltypes instead of numbers ... 
# The figure corresponds to 10_Reclustering_MP_and_MC/dimPlot_Clustered_res.1.pdf
dataDir <- "10_Reclustering_MP_and_MC"
seurat <- readRDS(paste(dataDir, "selectedClusters.subjectIntegrated.reclustered.rds", sep="/"))

# Make a new Ident
seurat@meta.data$CellTypes_res1 <- seurat@meta.data$integrated_snn_res.1
Idents(seurat) <- "CellTypes_res1"
levels(seurat)

# And rename
seurat <- RenameIdents(seurat, "0"="0: Non-classical MDMs", "1"="1: Classical MDMs", "2"="2: T cells",
                       "3"="3: cDCs","4"="4: Res-like C1q Mac","5"="5: T cells","6"="6: T & NK cells", "7"="7: cDCs & T cells",
                       "8"="8: SPP1 LAMs","9"="9: Overlap Mac", "10"="10: T cells","11"="11: Stress response non-classical MDMs",
                       "12"="12: cDCs")

clust <- "CellTypes_res1"
res <- gsub("CellTypes_res","", clust)
  
pdf(paste(figDir, paste0("dimPlot_Clustered_", res,".pdf"), sep="/"))
  p <- DimPlot(seurat, reduction = "umap", label = FALSE, label.size =14, repel = TRUE, pt.size = 1.5) +
    guides(colour=guide_legend(override.aes = list(size=2), ncol=1))
  # Get the legend
  legend_p <- cowplot::get_legend(p)
  
  Idents(seurat) <- "integrated_snn_res.1"
  p <- DimPlot(seurat, reduction = "umap", label = FALSE, label.size =14, repel = TRUE, pt.size = 1.5) +
    guides(colour=guide_legend(override.aes = list(size=2), ncol=1))
  p <- LabelClusters(p, id = "ident",  fontface = "bold", size=10)
  print(p)
  # Now plot it without legend
  # And use the cluster numbers
  p2 <- p + NoLegend()
  print(p2)
  # And plot the legend
  plot_grid(legend_p)
  
dev.off()

# And save as TIFF
ggsave(paste(figDir, "/tiff/dimPlot_Clustered_1.tiff", sep="/"), p, scale=2, width=5, height=5, dpi=360,compression="lzw")
ggsave(paste(figDir, "/tiff/dimPlot_Clustered_1.NoLegend.tiff", sep="/"), p2, scale=2, width=5, height=5, dpi=360,compression="lzw")
ggsave(paste(figDir, "/tiff/dimPlot_Clustered_1.LegendOnly.tiff", sep="/"), legend_p, scale=2, width=5, height=5, dpi=360,compression="lzw")


#### Distribution PR3 vs MPO ####
# Using the same UMAP
dataDir <- "10_Reclustering_MP_and_MC"
seurat <- readRDS(paste(dataDir, "selectedClusters.subjectIntegrated.reclustered.rds", sep="/"))

# Use the 'subType' column to subset
seurat <- subset(seurat, subset = subType %in% c("PR3", "MPO"))
Idents(seurat) <- "subType"

# And visualize
pdf(paste(figDir, "dimPlot_PR3_vs_MPO.pdf", sep="/"))
  p <- DimPlot(seurat, reduction = "umap", label = FALSE, label.size = 6, repel = TRUE, pt.size = 1.5) +
    guides(colour=guide_legend(override.aes = list(size=2), ncol=1))
  print(p)
  # Get the legend
  legend_p <- cowplot::get_legend(p)
  # Now plot it without legend
  p2 <- DimPlot(seurat, reduction = "umap", label = FALSE, label.size = 6, repel = TRUE, pt.size = 1.5) +
    guides(colour=guide_legend(override.aes = list(size=2), ncol=1)) + NoLegend()
  print(p2)
  # And plot the legend
  plot_grid(legend_p)
dev.off()

# And save as TIFF
ggsave(paste(figDir, "/tiff/dimPlot_PR3_vs_MPO.tiff", sep="/"), p, scale=2, width=5, height=5, dpi=360,compression="lzw")
ggsave(paste(figDir, "/tiff/dimPlot_PR3_vs_MPO.NoLegend.tiff", sep="/"), p2, scale=2, width=5, height=5, dpi=360,compression="lzw")
ggsave(paste(figDir, "/tiff/dimPlot_PR3_vs_MPO.LegendOnly.tiff", sep="/"), legend_p, scale=1, width=5, height=5, dpi=360,compression="lzw")


#### UMAP @res0.3 after reclustering ####
# Use celltypes instead of numbers ... 
# The figure corresponds to 10_Reclustering_MP_and_MC/dimPlot_Clustered_res.0.3.pdf
dataDir <- "10_Reclustering_MP_and_MC"
seurat <- readRDS(paste(dataDir, "selectedClusters.subjectIntegrated.reclustered.rds", sep="/"))

# Make a new Ident
seurat@meta.data$CellTypes_res0.3 <- seurat@meta.data$integrated_snn_res.0.3
Idents(seurat) <- "CellTypes_res0.3"
levels(seurat)

# And rename
seurat <- RenameIdents(seurat, "0"="T cells", "1"="Classical MDMs", "2"="Non-classical MDMs",
                       "3"="Res-like Mac","4"="DCs")

clust <- "CellTypes_res0.3"
res <- gsub("CellTypes_res","", clust)

pdf(paste(figDir, paste0("dimPlot_Clustered_", res,".pdf"), sep="/"))
  p <- DimPlot(seurat, reduction = "umap", label = FALSE, label.size = 6, repel = TRUE, pt.size = 1.5) +
    guides(colour=guide_legend(override.aes = list(size=2), ncol=1))
  # Get the legend
  legend_p <- cowplot::get_legend(p)
  
  # Set the ident to integrated_snn_res.0.3
  Idents(seurat) <- "integrated_snn_res.0.3"
  p <- DimPlot(seurat, reduction = "umap", label = FALSE, label.size =14, repel = TRUE, pt.size = 1.5) +
    guides(colour=guide_legend(override.aes = list(size=2), ncol=1))
  p <- LabelClusters(p, id = "ident",  fontface = "bold", size=10)
  print(p)
  # Now plot it without legend
  # And use the cluster numbers
  p2 <- p + NoLegend()
  print(p2)
  # And plot the legend
  plot_grid(legend_p)
  
dev.off()

# And save as TIFF
ggsave(paste(figDir, "/tiff/dimPlot_Clustered_0.3.tiff", sep="/"), p, scale=2, width=5, height=5, dpi=360,compression="lzw")
ggsave(paste(figDir, "/tiff/dimPlot_Clustered_0.3.NoLegend.tiff", sep="/"), p2, scale=2, width=5, height=5, dpi=360,compression="lzw")
ggsave(paste(figDir, "/tiff/dimPlot_Clustered_0.3.LegendOnly.tiff", sep="/"), legend_p, scale=0.5, width=5, height=5, dpi=360,compression="lzw")

#### Distribution clusters at res 0.3 ####
# see 10_Reclustering_MP_and_MC/distribution_per_subject_stacked_percentages_cluster_res0_3.pdf
dataDir <- "10_Reclustering_MP_and_MC"

df3 <- read.table(paste(dataDir,"clusters_per_subject_percentages_res0_3.txt", sep="/"), sep="\t", header=TRUE)

# df3 <- df2[, c("cluster","subject", "percent")]
# df3 <- df3 %>% pivot_wider(names_from = subject, values_from = percent)
# 
# Now, do the reverse, i.e. pivot_longer
df2 <- df3 %>% pivot_longer(!cluster, names_to = "subject", values_to = "percent")

# Get the right colors per cluster
myColors <- gg_color_hue(length(unique(df2$cluster)))
df2$colors <- myColors[factor(df2$cluster)]

# Rename: AGN1_PR3, AGN2_PR3, AGN3_PR3, AGN4_MPO, AGN5_MPO, LN, NC
df2$subject <- factor(df2$subject, 
                      levels = c("s1_ANCA_PR3", "s2_ANCA_PR3","s3_ANCA_PR3", 
                                 "s4_ANCA_MPO", "s5_ANCA_MPO", "SLE", "HC"),
                      labels = c("AGN1_PR3", "AGN2_PR3","AGN3_PR3", 
                                 "AGN4_MPO", "AGN5_MPO", "LN", "NC"))

# And make a factor of the colors
# TO DO: We can also rename the levels if wanted ... see above
#    "0"="T cells", "1"="Classical MDMs", "2"="Non-classical MDMs","3"="Res-like Mac","4"="DCs"
df2$cluster <- factor(df2$cluster,
                      levels = c(0,1,2,3,4),
                      labels = c("T cells", "Classical MDMs", "Non-classical MDMs",
                                 "Res-like Mac","DCs"))

pdf(paste(figDir, paste0("distribution_per_subject_stacked_percentages_cluster_res0_3.pdf"), sep="/"), width = 14, height =14)
  p <- ggplot(df2, aes(x=subject, y=percent, fill=cluster)) +
    geom_bar(stat='identity', position="stack") +
    xlab("") + ylab("Proportion") +
    scale_fill_manual("Cell type", values = myColors) +
    guides(fill=guide_legend(ncol=1)) +
    theme_cowplot() +
    theme(axis.text.x=element_text(angle=90,hjust=1, size=24)) +
    theme(axis.text.y=element_text(size=24)) +
    theme(axis.title.y = element_text(size=24)) +
    theme(legend.title = element_text(size=28), #change legend title font size
          legend.text = element_text(size=24)) #change legend text font size
  print(p)
  # Get the legend as a separate figure
  legend_p <- cowplot::get_legend(p)
  p2 <- ggplot(df2, aes(x=subject, y=percent, fill=cluster)) +
    geom_bar(stat='identity', position="stack") +
    xlab("") + ylab("Proportion") +
    scale_fill_manual("Cell type",values = myColors) +
    guides(fill=guide_legend(ncol=1)) +
    theme_cowplot() +
    theme(axis.text.x=element_text(angle=90,hjust=1, size=18)) +
    theme(axis.text.y=element_text(size=16)) +
    theme(axis.title.y = element_text(size=18)) +
    NoLegend()
  print(p2)
  # And print the legend
  plot_grid(legend_p)
dev.off()

# And as TIFF
ggsave(paste(figDir, "/tiff/distribution_per_subject_stacked_percentages_clusters_res0_3.tiff", sep="/"), p, scale=2, width=5, height=5, dpi=360,compression="lzw")
ggsave(paste(figDir, "/tiff/distribution_per_subject_stacked_percentages_clusters_res0_3.NoLegend.tiff", sep="/"), p2, scale=2, width=5, height=5, dpi=360,compression="lzw")
ggsave(paste(figDir, "/tiff/distribution_per_subject_stacked_percentages_clusters_res0_3.LegendOnly.tiff", sep="/"), legend_p, scale=0.75, width=5, height=5, dpi=360,compression="lzw")

#### Markerplot @res0.3 ####
dataDir <- "10_Reclustering_MP_and_MC/"
seurat <- readRDS(paste(dataDir, "selectedClusters.subjectIntegrated.reclustered.rds", sep="/"))

# Set the default assay slot correct
DefaultAssay(seurat) <- "SCT"

# We are going to show the marker plot for resolution 1 and only cluster 0, 1, 4, 8, 9 and 11
# and thus also only the top5 marker for these clusters and add some other ones
seurat.markers <- read.table(paste(dataDir,paste0("markers_20230418_LR/clusterMarkers_res.0.3.txt"), sep="/"), header = TRUE)

# Give top 5 per comparison based on p_adj_val (need top_n(n=-5) to get the right 5 genes!!)
#top5.markers <- seurat.markers %>% group_by(cluster) %>% top_n(n = -5, wt = p_val_adj)  # wt defines ordering, avg_log2FC, or p_val_adj?
# Hmm, returns more than 5 rows per clusters as there are ties in the 'p_val_adj' ...
# So, also sort on avg_log2FC ...
top5.markers <- seurat.markers %>% group_by(cluster) %>% dplyr::arrange(p_val_adj, desc(avg_log2FC), .by_group=TRUE) %>% top_n(n = 5, wt = avg_log2FC)

# # Select only the genes in the cluster 0,1,4,8,9,11
# top5.markers <- top5.markers[top5.markers$cluster %in% c(0,1,4,8,9,11),]

# NOTE: Should we add these
# Now add CD163, MRC1, MARCO, TREM2, CD14, PLIN2, C1QC and SPP1 (if they are not already in there)
# AJ 20230606 - Not C1QA as previously??
# for (markerGene in c("CD163", "MRC1", "MARCO", "TREM2", "PLIN2", "CD14", "C1QC", "SPP1")){
#   if ( !(markerGene %in% top5.markers$gene)){
#     print(paste0("... Now adding ", markerGene, " to the list of marker genes"))
#     top5.markers <- rbind(top5.markers, seurat.markers[seurat.markers$gene == markerGene, ])
#   }
# }

# Sort the top5 list on the cluster column
top5.markers <- dplyr::arrange(top5.markers, cluster)

# Plot a heatmap and dotplot
# resolution 1 -> only get clusters 0,1,4,8,9,11
clust = "integrated_snn_res.0.3"

# Also get the correct colors
myColors <- gg_color_hue(length(unique(seurat.markers$cluster)))

pdf(paste(figDir, "markerPlots_Top5_res0_3_LR.pdf", sep="/"), width = 28, height = 28)
  # Perhaps downsample to 100,200 cells??
  p1 <- DoHeatmap(seurat, features = top5.markers$gene, group.by = clust, 
                size = 28, angle=0, hjust=0.5, group.colors = myColors) + 
    theme(axis.text.y = element_text(size=36, face = "bold.italic")) +
    theme(legend.key.size = unit(2, 'cm'), 
        legend.key.height= unit(1, 'cm'),
        legend.key.width= unit(1, 'cm'),
        legend.title = element_text(size=28),legend.text = element_text(size=24) ) +
    guides(colour = guide_legend(override.aes = list(size=28, alpha =1))) 
  # Get the legend
  legend_p <- cowplot::get_legend(p1)
  print(p1)
  # Now plot it without legend
  p2 <- DoHeatmap(seurat, features = top5.markers$gene, group.by = clust, 
                size = 28, angle=0, hjust=0.5, group.colors = myColors) + 
    theme(axis.text.y = element_text(size=36, face = "bold.italic")) +
    NoLegend()
  print(p2)
  # And plot the legend
  plot_grid(legend_p)

dev.off()

# And as TIFF
ggsave(paste(figDir, "/tiff/markerPlots_Top5_res0_3_LR.tiff", sep="/"), p1, scale=3, width=5, height=5, dpi=360,compression="lzw")
ggsave(paste(figDir, "/tiff/markerPlots_Top5_res0_3_LR.NoLegend.tiff", sep="/"), p2, scale=3, width=5, height=5, dpi=360,compression="lzw")
ggsave(paste(figDir, "/tiff/markerPlots_Top5_res0_3_LR.LegendOnly.tiff", sep="/"), legend_p, scale=3, width=5, height=5, dpi=360,compression="lzw")

#### Markerplot Surface proteins (ADT) ####
# See 10_Reclustering_MP_and_MC/ADT_markerPlots_Top5_res.1.selectedClusters.pdf
# i.e. 10_Reclustering_MP_and_MC #'Find markers per cluster - UPDATE 20230418'
# Increase font names and add a legend
dataDir <- "10_Reclustering_MP_and_MC/"
seurat <- readRDS(paste(dataDir, "selectedClusters.subjectIntegrated.reclustered.rds", sep="/"))

clust <- "integrated_snn_res.1"

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

#if ( outDir == paste0(dataDir,"/markers_20230418_LR") ){
  print("... Integrated ADT data: Using FindAllMarkers with latent.vars=\"orig.ident\" and test.use=\"LR\"")
  Idents(adt.integrated) <- clust
  adt.markers <- FindAllMarkers(adt.integrated, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25,
                                latent.vars = "orig.ident", test.use = 'LR', assay = "ADT")
  ## IMPORTANT ##
  # We get 16 clusters out ... so, there is no match any more with the 12 clusters that we get for resolution 1.0
  # on the gene expression level. But we restore that later on (line 509)
  
  ## SOLUTION - 20230721
  # set Idents(adt.integrated) <- clust, i.e. integrated_snn_res.1
  # Probably that will give you another figure as well ...
  # Ran this code again at 20230721 and added the date in the name
  
#} else {
#  print("... Integrated ADT data: Using default FindAllMarkers")
#  adt.markers <- FindAllMarkers(adt.integrated, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25 , assay = "ADT")
#}
top5.adt_markers <- adt.markers %>% group_by(cluster) %>% dplyr::arrange(p_val_adj, desc(avg_log2FC), .by_group=TRUE) %>% top_n(n = 5, wt = avg_log2FC)
top5.adt_markers <- dplyr::arrange(top5.adt_markers, cluster)

selClusters.adt <- subset(x = adt.integrated, subset = integrated_snn_res.1 %in% c("0","1","4","8","9","11"))
selClusters.adt$integrated_snn_res.1 <- factor(selClusters.adt$integrated_snn_res.1, levels=c("0","1","4","8","9","11"))

# And also make the selection in the top ADT markers
top5.adt_markers <- top5.adt_markers[top5.adt_markers$cluster %in% c("0","1","4","8","9","11"),]
top5.adt_markers$cluster <- factor(top5.adt_markers$cluster, levels=c("0","1","4","8","9","11"))
# and order by cluster
top5.adt_markers <- dplyr::arrange(top5.adt_markers, cluster)

# Also get the correct colors
myColors <- gg_color_hue(length(unique(seurat$integrated_snn_res.1)))
myColors <- myColors[c(0,1,4,8,9,11)+1]

#pdf(paste(figDir, paste0("ADT_markerPlots_Top5_res.1.selectedClusters.pdf"), sep="/"), width = 28, height = 28)
pdf(paste(figDir, paste0("ADT_markerPlots_Top5_res.1.selectedClusters.20230721.pdf"), sep="/"), width = 28, height = 28)
  p1 <- DoHeatmap(selClusters.adt, features = top5.adt_markers$gene, assay = "integrated", group.by = clust,
                 size=12, angle=0, hjust=0.5, group.colors = myColors) + 
    scale_y_discrete(labels=function(x)gsub("anti-human ", "", x)) + 
    theme(axis.text.y = element_text(size=36, face = "bold.italic")) +
    theme(legend.key.size = unit(4, 'cm'), 
          legend.key.height= unit(3, 'cm'),
          legend.key.width= unit(2, 'cm'),
          legend.title = element_text(size=24),legend.text = element_text(size=24) ) +
    guides(colour = guide_legend(override.aes = list(size=18, alpha=1)))
  # Get the legend
  legend_p <- cowplot::get_legend(p1)
  print(p1)
  # Now plot it without legend
  p2 <- DoHeatmap(selClusters.adt, features = top5.adt_markers$gene, assay = "integrated",group.by = clust,
                 size=12, angle=0, hjust=0.5, group.colors = myColors) + 
    scale_y_discrete(labels=function(x)gsub("anti-human ", "", x)) + 
    theme(axis.text.y = element_text(size=36, face = "bold.italic")) +
    NoLegend()  
  print(p2)
  # And plot the legend
  plot_grid(legend_p)
  
dev.off()

# And as TIFF
ggsave(paste(figDir, "/tiff/ADT_markerPlots_Top5_res.1.selectedClusters.tiff", sep="/"), p1, scale=4, width=5, height=5, dpi=360,compression="lzw")
ggsave(paste(figDir, "/tiff/ADT_markerPlots_Top5_res.1.selectedClusters.NoLegend.tiff", sep="/"), p2, scale=4, width=5, height=5, dpi=360,compression="lzw")
ggsave(paste(figDir, "/tiff/ADT_markerPlots_Top5_res.1.selectedClusters.LegendOnly.tiff", sep="/"), legend_p, scale=4, width=5, height=5, dpi=360,compression="lzw")

#### Markerplot Surface proteins (ADT) - macrophage markers  ####
# Mail Yosta 20072023
# - She wants a list of 59 markers not just the top 5
#   See Yosta_Vegting\scRNASeq\Data\Dataset_1\Meta\TotalSeq_C_Human_Universal_Cocktail_v1_137_Antibodies_399905_Barcodes.xlsx - the markers in yellow

# See 10_Reclustering_MP_and_MC/ADT_markerPlots_Top5_res.1.selectedClusters.pdf
# i.e. 10_Reclustering_MP_and_MC #'Find markers per cluster - UPDATE 20230418'
# Increase font names and add a legend
dataDir <- "10_Reclustering_MP_and_MC/"
seurat <- readRDS(paste(dataDir, "selectedClusters.subjectIntegrated.reclustered.rds", sep="/"))

clust <- "integrated_snn_res.1"

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

# Now, get the 59 markers that Yosta wants
library(xlsx)

wb     <- loadWorkbook("../../../Data/Dataset_1/Meta/TotalSeq_C_Human_Universal_Cocktail_v1_137_Antibodies_399905_Barcodes.xlsx") 
sheet1 <- getSheets(wb)[[1]]
rows   <- getRows(sheet1)
cells  <- getCells(rows)
styles <- sapply(cells, getCellStyle)
cellColor <- function(style) 
{
  fg  <- style$getFillForegroundXSSFColor()
  rgb <- tryCatch(fg$getRgb(), error = function(e) NULL)
  rgb <- paste(rgb, collapse = "")
  return(rgb)
}
mycolor <- (yellow = "ffff00")
m     <- match(sapply(styles, cellColor), mycolor)

idx <- which(m==1)
length(idx)
#[1] 73

# And get the values
adt.markers <- unlist(lapply(cells[idx], getCellValue))
names(adt.markers) <- adt.markers
# These are the gene names, unfortunately not the rownames of the ADT's...
# Have to get these form the XLSX
df <- read.xlsx("../../../Data/Dataset_1/Meta/TotalSeq_C_Human_Universal_Cocktail_v1_137_Antibodies_399905_Barcodes.xlsx", 1,header=TRUE,colClasses=c("character"),encoding="UTF-8")
adt.markers <- df[df$`Gene.name` %in% adt.markers, "Description"]
length(adt.markers)
# [1] 74     -- Ah, funny, 1 more ! Why is this?

# Change some of the 'special characters'
# https://stackoverflow.com/questions/20916070/importing-an-excel-file-with-greek-characters-into-r-in-the-correct-encoding
tolatin=function(compname) {
  n=as.character(compname,encoding="UTF-8")
  n=gsub("\u03B1","alpha",n)
  n=gsub("\u03B2","beta",n)
  n=gsub("\u03B3","gamma",n)
  n=gsub("\u03B4","delta",n)
  n=gsub("\u03B5","epsilon",n)
  n
}

adt.markers <- tolatin(adt.markers)

# NOTE:
# Something might have gone wrong earlier in naming the ADT's
# FCERIA has been denoted as 'Fc eta RI beta' ...
adt.markers <- gsub("FcepsilonRIalpha", "Fc eta RI beta", 
                    gsub("Integrin alphaE", "Integrin E", adt.markers))

# Subset the clusters
selClusters.adt <- subset(x = adt.integrated, subset = integrated_snn_res.1 %in% c("0","1","4","8","9","11"))
selClusters.adt$integrated_snn_res.1 <- factor(selClusters.adt$integrated_snn_res.1, levels=c("0","1","4","8","9","11"))

# Also get the correct colors
myColors <- gg_color_hue(length(unique(seurat$integrated_snn_res.1)))
myColors <- myColors[c(0,1,4,8,9,11)+1]

# But now there would be no clustering ....
# https://www.biostars.org/p/9564526/
# https://bioinformatics.stackexchange.com/questions/8973/doheatmap-hierarchical-clustering-seurat
#if ( outDir == paste0(dataDir,"/markers_20230418_LR") ){
print("... Integrated ADT data: Using FindAllMarkers with latent.vars=\"orig.ident\" and test.use=\"LR\"")
Idents(adt.integrated) <- clust     # Essential!!
all.markers <- FindAllMarkers(adt.integrated, only.pos = FALSE, min.pct = 0.0, logfc.threshold = 0.0,
                              latent.vars = "orig.ident", test.use = 'LR', assay = "ADT", return.thresh = 1.0)

# Do not use any cutoff, just to get all back ...

## IMPORTANT ##
# You have to correctly set the ident!!
# Otherwise you will get 16 clusters out ... and there is no match any more with the 12 clusters that we get for resolution 1.0
# on the gene expression level. 

#} else {
#  print("... Integrated ADT data: Using default FindAllMarkers")
#  adt.markers <- FindAllMarkers(adt.integrated, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25 , assay = "ADT")
#}
all.markers <- all.markers %>% group_by(cluster) %>% dplyr::arrange(desc(avg_log2FC), .by_group=TRUE)

all.markers <- dplyr::arrange(all.markers, cluster)

# And also make the selection in the top ADT markers
all.markers <- all.markers[all.markers$cluster %in% c("0","1","4","8","9","11"),]
all.markers$cluster <- factor(all.markers$cluster, levels=c("0","1","4","8","9","11"))
# and order by cluster
all.markers <- dplyr::arrange(all.markers, cluster)

# Is every selected ADT marker present?
match(adt.markers, all.markers$gene)
# [1]   1   7   3 132 131   4   6 135   5 136 123   2 110   8 137  12 115 126 109 114  14 111  10  16 118 122 100 120 105 117 128 106  17
#[34] 129 108  25 104 127  15 103 124  88  95  82 101  19  32  84  30  24  43  38  97  67  22  18  21  27  40  49  51  55  58  60  63  68
#[67]  72  73  75  76  85  86  91  99
# All found ... :-)

# And now put the selected markers in this order and put the NA's last??
adt.markers <- adt.markers[order(match(adt.markers, all.markers$gene))]
adt.markers[1:5]

pdf(paste(figDir, paste0("ADT_markerPlots_MacrophageMarkers_res.1.selectedClusters.pdf"), sep="/"), width = 28, height = 28)
p1 <- DoHeatmap(selClusters.adt, features = adt.markers, assay = "integrated", group.by = clust,
                size=12, angle=0, hjust=0.5, group.colors = myColors) + 
  scale_y_discrete(labels=function(x)gsub("anti-human ", "", x)) + 
  theme(axis.text.y = element_text(size=28, face = "bold.italic")) +
  theme(legend.key.size = unit(4, 'cm'), 
        legend.key.height= unit(3, 'cm'),
        legend.key.width= unit(2, 'cm'),
        legend.title = element_text(size=24),legend.text = element_text(size=24) ) +
  guides(colour = guide_legend(override.aes = list(size=18, alpha=1)))
# In DoHeatmap(selClusters.adt, features = adt.markers, assay = "integrated",  :
#    The following features were omitted as they were not found in the scale.data slot for the integrated assay: 
#    anti-human CD127 (IL-7Ra), anti-human CD124 (IL-4Ra), anti-human FceRIa, anti-human CD119 (IFN-<U+03B3> R a chain), 
#    anti-human/mouse integrin ?7, anti-human CD103 (Integrin aE)
# Has to do with special characters ... see code above to remedy this!!!

# Get the legend
legend_p <- cowplot::get_legend(p1)
print(p1)
# Now plot it without legend
p2 <- DoHeatmap(selClusters.adt, features = adt.markers, assay = "integrated",group.by = clust,
                size=12, angle=0, hjust=0.5, group.colors = myColors) + 
  scale_y_discrete(labels=function(x)gsub("anti-human ", "", x)) + 
  theme(axis.text.y = element_text(size=28, face = "bold.italic")) +
  NoLegend()  
print(p2)
# And plot the legend
plot_grid(legend_p)

dev.off()

# And as TIFF
ggsave(paste(figDir, "/tiff/ADT_markerPlots_MacrophageMarkers_res.1.selectedClusters.tiff", sep="/"), p1, scale=4, width=5, height=5, dpi=360,compression="lzw")
ggsave(paste(figDir, "/tiff/ADT_markerPlots_MacrophageMarkers_res.1.selectedClusters.NoLegend.tiff", sep="/"), p2, scale=4, width=5, height=5, dpi=360,compression="lzw")
ggsave(paste(figDir, "/tiff/ADT_markerPlots_MacrophageMarkers_res.1.selectedClusters.LegendOnly.tiff", sep="/"), legend_p, scale=4, width=5, height=5, dpi=360,compression="lzw")


#### Distribution clusters at res 1 ####
#Extra desired figures - 11-06-2023 

# see 10_Reclustering_MP_and_MC/distribution_per_subject_stacked_percentages_cluster_res0_3.pdf
dataDir <- "10_Reclustering_MP_and_MC"

df3 <- read.table(paste(dataDir,"clusters_per_subject_percentages_res1.txt", sep="/"), sep="\t", header=TRUE)

# df3 <- df2[, c("cluster","subject", "percent")]
# df3 <- df3 %>% pivot_wider(names_from = subject, values_from = percent)
# 
# Now, do the reverse, i.e. pivot_longer
df2 <- df3 %>% pivot_longer(!cluster, names_to = "subject", values_to = "percent")

# Get the right colors per cluster
myColors <- gg_color_hue(length(unique(df2$cluster)))
df2$colors <- myColors[factor(df2$cluster)]

# Rename: AGN1_PR3, AGN2_PR3, AGN3_PR3, AGN4_MPO, AGN5_MPO, LN, NC
df2$subject <- factor(df2$subject, 
                      levels = c("s1_ANCA_PR3", "s2_ANCA_PR3","s3_ANCA_PR3", 
                                 "s4_ANCA_MPO", "s5_ANCA_MPO", "SLE", "HC"),
                      labels = c("AGN1_PR3", "AGN2_PR3","AGN3_PR3", 
                                 "AGN4_MPO", "AGN5_MPO", "LN", "NC"))

# And make a factor of the colors
# TO DO: We can also rename the levels if wanted ... see above
#    "0"="0: Non-classical MDMs", "1"="1: Classical MDMs", "2"="2: T cells",
#    "3"="3: cDCs","4"="4: Res-like C1q Mac","5"="5: T cells","6"="6: T & NK cells", "7"="7: cDCs & T cells",
#    "8"="8: SPP1 LAMs","9"="9: Overlap Mac", "10"="10: T cells","11"="11: Stress response non-classical MDMs",
#    "12"="12: cDCs"

df2$cluster <- factor(df2$cluster,
                      levels = 0:12,
                      labels = c("0: Non-classical MDMs", "1: Classical MDMs", "2: T cells",
                        "3: cDCs","4: Res-like C1q Mac","5: T cells","6: T & NK cells", "7: cDCs & T cells",
                        "8: SPP1 LAMs","9: Overlap Mac", "10: T cells","11: Stress response non-classical MDMs",
                        "12: cDCs"))

pdf(paste(figDir, paste0("distribution_per_subject_stacked_percentages_cluster_res1.pdf"), sep="/"), width = 14, height =14)
p <- ggplot(df2, aes(x=subject, y=percent, fill=cluster)) +
  geom_bar(stat='identity', position="stack") +
  xlab("") + ylab("Proportion") +
  scale_fill_manual("Cell type", values = myColors) +
  guides(fill=guide_legend(ncol=1)) +
  theme_cowplot() +
  theme(axis.text.x=element_text(angle=90,hjust=1, size=24)) +
  theme(axis.text.y=element_text(size=24)) +
  theme(axis.title.y = element_text(size=24)) +
  theme(legend.title = element_text(size=28), #change legend title font size
        legend.text = element_text(size=24)) #change legend text font size
print(p)
# Get the legend as a separate figure
legend_p <- cowplot::get_legend(p)
p2 <- ggplot(df2, aes(x=subject, y=percent, fill=cluster)) +
  geom_bar(stat='identity', position="stack") +
  xlab("") + ylab("Proportion") +
  scale_fill_manual("Cell type",values = myColors) +
  guides(fill=guide_legend(ncol=1)) +
  theme_cowplot() +
  theme(axis.text.x=element_text(angle=90,hjust=1, size=18)) +
  theme(axis.text.y=element_text(size=16)) +
  theme(axis.title.y = element_text(size=18)) +
  NoLegend()
print(p2)
# And print the legend
plot_grid(legend_p)
dev.off()

# And as TIFF
ggsave(paste(figDir, "/tiff/distribution_per_subject_stacked_percentages_clusters_res1.tiff", sep="/"), p, scale=2, width=5, height=5, dpi=360,compression="lzw")
ggsave(paste(figDir, "/tiff/distribution_per_subject_stacked_percentages_clusters_res1.NoLegend.tiff", sep="/"), p2, scale=2, width=5, height=5, dpi=360,compression="lzw")
ggsave(paste(figDir, "/tiff/distribution_per_subject_stacked_percentages_clusters_res1.LegendOnly.tiff", sep="/"), legend_p, scale=0.75, width=5, height=5, dpi=360,compression="lzw")

#### ViolinPlots at res1 ####
dataDir <- "10_Reclustering_MP_and_MC/"
seurat <- readRDS(paste(dataDir, "selectedClusters.subjectIntegrated.reclustered.rds", sep="/"))

markerList <- c("CD14", "FCGR3A","CD163", "C1QC", "C1QA", "TREM2", "SPP1", "S100A9")

DefaultAssay(seurat) <- "SCT"

#for (clust in colnames(seurat@meta.data)[grep("integrated_snn", colnames(seurat@meta.data))]){
clust <- "integrated_snn_res.1"

res <- gsub("integrated_snn_","", clust)
seurat@meta.data[,clust] <- factor(seurat@meta.data[,clust], levels=as.character(0:(length(unique(seurat@meta.data[,clust]))-1)))

# maxY.RNA <- max(seurat@assays$RNA@data[c("CD163", "MRC1", "MARCO", "C1QA", "CD74","CD81", "CD14", "FCGR3A", "FCGR3B", "HLA-DRA","CD68"),])
maxY.RNA <- NULL

# 
pdf(paste(figDir,paste0("violinPlot_RNA_selectedMarkers_", res,".pdf"), sep="/"))
Idents(seurat) <- clust

for ( marker in markerList){
  p <- VlnPlot(object = seurat, features = marker, pt.size = 1, y.max = maxY.RNA, same.y.lims = TRUE) + NoLegend() +
    theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 18), axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14))
  print(p)
  
  # And save as TIFF
  ggsave(paste(figDir, paste0("/tiff/violinPlot_RNA_",marker,"_", res,".tiff"), sep="/"), p, scale=2, width=5, height=5, dpi=360,compression="lzw")
  
}

dev.off()

#### Plots before reclustering ####
seurat <- readRDS("9_kidney_reference_mapping/Mito10.integrated.clustered.annot.kidney_immune_reference.rds")
DefaultAssay(seurat) <- "RNA"

# Create a directory to store these results
baseOutDir <- outDir
outDir <- paste0(outDir,"/beforeReclustering")
dir.create(outDir)
dir.create(paste0(outDir,"/tiff"))

## 20230926 - AJ: changed seurat$predicted.labels.celltypist.after to seurat$predicted.celltype and changed outDir to 'finalFigures/beforeReclustering/predicted_celltypes

# Make a nicer sample annotation
seurat@meta.data$Condition <- ifelse(seurat@meta.data$orig.ident %in% c("MOMA17", "MOMA52", "MOMA67", "MOMA302", "MOMA72"), "ANCA",
                                     ifelse(seurat@meta.data$orig.ident == "MOMA57", "Control", "SLE"))

seurat@meta.data$subType <- ifelse(seurat@meta.data$orig.ident %in% c("MOMA17", "MOMA52", "MOMA302"), "PR3",
                                   ifelse(seurat@meta.data$orig.ident %in% c("MOMA57", "MOMA68"), "NONE", "MPO"))

seurat@meta.data$sampleAnnot <- paste0(gsub("MOMA", "s", seurat@meta.data$orig.ident),"_",seurat@meta.data$Condition,
                                       "_", seurat@meta.data$subType)


tab_0_3 <- table(seurat$integrated_snn_res.0.3, seurat$sampleAnnot, seurat$predicted.celltype)
write.table(tab_0_3, paste(outDir,"celltypes_per_subject_and_cluster_res0_3.txt", sep="/"), quote=FALSE, sep="\t", row.names = TRUE, col.names = TRUE)

tab_0_6 <- table(seurat$integrated_snn_res.0.6, seurat$sampleAnnot, seurat$predicted.celltype)
write.table(tab_0_6, paste(outDir,"celltypes_per_subject_and_cluster_res0_6.txt", sep="/"), quote=FALSE, sep="\t", row.names = TRUE, col.names = TRUE)

tab_1 <- table(seurat$integrated_snn_res.1, seurat$sampleAnnot, seurat$predicted.celltype)
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
  
  # levels(df2$celltype)
  # #[1] "Classical monocytes"  "DC1"  "DC2"  "Intermediate macrophages" "Non-classical monocytes" 
  # #[5] "Unassigned" 
  
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

## UMAP - see 8_clustering.r
for (res in c("0_3","0_6", "1")){
  pdf(paste(outDir, paste0("dimPlot_Clustered_", res,".pdf"), sep="/"))
    p <- DimPlot(seurat, reduction = "umap", label = FALSE, label.size =14, repel = TRUE, pt.size = NULL) +
         guides(colour=guide_legend(override.aes = list(size=2), ncol=1))
    # Get the legend
    legend_p <- cowplot::get_legend(p)

    Idents(seurat) <- paste0("integrated_snn_res.",gsub("_",".", res))
    p <- DimPlot(seurat, reduction = "umap", label = FALSE, label.size =14, repel = TRUE, pt.size = NULL) +
         guides(colour=guide_legend(override.aes = list(size=2), ncol=1))
    p <- LabelClusters(p, id = "ident",  fontface = "bold", size=6)
    print(p)
    # Now plot it without legend
    # And use the cluster numbers
    p2 <- p + NoLegend()
    print(p2)
    # And plot the legend
    plot_grid(legend_p)

  dev.off()

  # And save as TIFF
  ggsave(paste(outDir, paste0("/tiff/dimPlot_Clustered_", res, ".tiff"), sep="/"), p, scale=2, width=5, height=5, dpi=360,compression="lzw")
  ggsave(paste(outDir, paste0("/tiff/dimPlot_Clustered_", res, ".NoLegend.tiff"), sep="/"), p2, scale=2, width=5, height=5, dpi=360,compression="lzw")
  ggsave(paste(outDir, paste0("/tiff/dimPlot_Clustered_", res, ".LegendOnly.tiff"), sep="/"), legend_p, scale=2, width=5, height=5, dpi=360,compression="lzw")
}

## Top 5 markerplot
# See 10_reclustering_MP_and_MC.r
# - Only use LR version
dataDir <- outDir
outDir <- paste0(dataDir,"/markers_20230901_LR")
dir.create(outDir)

#for (outDir in c(paste0(dataDir,"/markers_integrated_20230418"), paste0(dataDir,"/markers_20230418_SCT"),paste0(dataDir,"/markers_20230418_LR") )){
for (outDir in c(paste0(dataDir,"/markers_20230901_LR") )){
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
      if ( outDir == paste0(dataDir,"/markers_20230901_LR") ){
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
    
  #   # In case of resolution 1, only get clusters 0,1,4,8,9,11
  #   if ( clust == "integrated_snn_res.1"){
  #     selClusters.seurat <- subset(x = seurat, subset = integrated_snn_res.1 %in% c("0","1","4","8","9","11"))
  #     selClusters.seurat$integrated_snn_res.1 <- factor(selClusters.seurat$integrated_snn_res.1, levels=c("0","1","4","8","9","11"))
  #     
  #     pdf(paste(outDir, paste0("markerPlots_Top5_", res,".selectedClusters.pdf"), sep="/"), width = 28, height = 28)
  #     # Perhaps downsample to 100,200 cells??
  #     p <- DoHeatmap(selClusters.seurat, features = top5.markers$gene, group.by = clust) + 
  #       theme(text = element_text(size=28), axis.text.y = element_text(size=14)) +
  #       NoLegend()  
  #     print(p)
  #     p <-  DotPlot(selClusters.seurat, features = unique(top5.markers$gene), dot.scale = 6) + 
  #       coord_flip() +  
  #       theme(text = element_text(size=28), axis.text.y = element_text(size=14))
  #     print(p)
  #     dev.off()
  #     
  #     # And get a Heatmap of the ADTs
  #     # - We need first to Normalize and Scale the ADT data
  #     #   see https://github.com/satijalab/seurat/issues/3890, https://github.com/satijalab/seurat/issues/5089
  #     # - Do we take the ADT separate from the RNA data??
  #     # - We only have ADTs for two samples
  #     
  #     # The ADT has been normalized (see '2_filtering_normalization_and_scaling')
  #     adt.list <- SplitObject(seurat, split.by = "orig.ident")
  #     # Only get the two samples with actual ADT data
  #     adt.list <- adt.list[c("MOMA52", "MOMA57")]
  #     print("... Finding Varaiable Features in the ADT samples")
  #     adt.list <- lapply(X = adt.list, FUN = function(x) {
  #       #x <- NormalizeData(x, assay = "ADT", normalization.method = "CLR", verbose = FALSE)
  #       x <- FindVariableFeatures(x, assay = "ADT", verbose = FALSE)
  #     })
  #     print("... Selecting integration features")
  #     for ( j in length(adt.list)){
  #       print(length(VariableFeatures(adt.list[[j]])))
  #     }
  #     if ( length(VariableFeatures(adt.list[[1]])) == 0 ){
  #       features <- rownames(adt.list[[1]]@assays$ADT@counts)
  #     } else {
  #       features <- SelectIntegrationFeatures(object.list = adt.list, assay = c("ADT", "ADT"))
  #     }
  #     print("... Scaling and running PCA")
  #     adt.list <- lapply(X = adt.list, FUN = function(x) {
  #       x <- ScaleData(x, features = features, assay = "ADT", verbose = FALSE)
  #       x <- RunPCA(x, features = features, assay = "ADT", verbose = FALSE)
  #     })
  #     print("... Finding integration anchors (using RCPA) and integrating ")
  #     anchors <- FindIntegrationAnchors(object.list = adt.list, assay = c("ADT", "ADT"), reference = c(1, 2), 
  #                                       reduction = "rpca", dims = 1:30)
  #     adt.integrated <- IntegrateData(anchorset = anchors, dims = 1:30)
  #     print("... Scaling, running PCA and obtaining a UMAP of the integrated ADT data")
  #     adt.integrated <- ScaleData(adt.integrated, verbose = FALSE)
  #     adt.integrated <- RunPCA(adt.integrated, verbose = FALSE)
  #     adt.integrated <- RunUMAP(adt.integrated, dims = 1:30)
  #     
  #     # And add to the RNA object, but now you get an error as you only
  #     # have ADTs for two samples:
  #     #  Error: Cannot add a different number of cells than already present
  #     #seurat[["IADT"]] <- adt.integrated[["integrated"]]
  #     #seurat[["pca.adt"]] <- adt.integrated[["pca"]]
  #     #seurat[["umap.adt"]] <- adt.integrated[["umap"]]
  #     
  #     if ( outDir == paste0(dataDir,"/markers_20230418_LR") ){
  #       print("... Integrated ADT data: Using FindAllMarkers with latent.vars=\"orig.ident\" and test.use=\"LR\"")
  #       adt.markers <- FindAllMarkers(adt.integrated, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25,
  #                                     latent.vars = "orig.ident", test.use = 'LR', assay = "ADT")
  #     } else {
  #       print("... Integrated ADT data: Using default FindAllMarkers")
  #       adt.markers <- FindAllMarkers(adt.integrated, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25 , assay = "ADT")
  #     }
  #     top5.adt_markers <- adt.markers %>% group_by(cluster) %>% arrange(p_val_adj, desc(avg_log2FC), .by_group=TRUE) %>% top_n(n = 5, wt = avg_log2FC)
  #     top5.adt_markers <- dplyr::arrange(top5.adt_markers, cluster)
  #     
  #     selClusters.adt <- subset(x = adt.integrated, subset = integrated_snn_res.1 %in% c("0","1","4","8","9","11"))
  #     selClusters.adt$integrated_snn_res.1 <- factor(selClusters.seurat$integrated_snn_res.1, levels=c("0","1","4","8","9","11"))
  #     
  #     pdf(paste(outDir, paste0("ADT_markerPlots_Top5_", res,".selectedClusters.pdf"), sep="/"), width = 28, height = 28)
  #     p <- DoHeatmap(selClusters.adt, features = top5.adt_markers$gene, assay = "integrated",group.by = clust) + 
  #       theme(text = element_text(size=28), axis.text.y = element_text(size=14)) +
  #       NoLegend()  
  #     print(p)
  #     p <-  DotPlot(selClusters.seurat, features = unique(top5.adt_markers$gene), dot.scale = 6, assay = "ADT") + 
  #       coord_flip() +  
  #       theme(text = element_text(size=28), axis.text.y = element_text(size=14))
  #     print(p)
  #     dev.off()
  #     
  #     
  #   }
  #   
    }
  # 
  # # ViolinPlots of clusters
  # for (clust in colnames(seurat@meta.data)[grep("integrated_snn", colnames(seurat@meta.data))]){
  #   res <- gsub("integrated_snn_r","R", clust)
  #   dir.create(paste(outDir,paste0("VlnPlots_", res), sep="/"))
  #   dir.create(paste(outDir,paste0("FeaturePlots_", res), sep="/"))
  #   
  #   Idents(seurat) <- clust
  #   maxY.RNA <- NULL # see mail Perry 06022023: original code maxY.RNA <- max(seurat@assays$RNA@data)
  #   # 21022023: Set the axis back in this case !!!
  #   
  #   if ( outDir == paste0(dataDir,"/markers_20230221_LR") ){
  #     seurat.markers <- read.table(paste("10_Reclustering_MP_and_MC/markers_20230221_LR",paste0("clusterMarkers_", res,".txt"), sep="/"), header = TRUE)
  #   } 
  #   if ( outDir == paste0(dataDir,"/markers_20230207") ){
  #     seurat.markers <- read.table(paste("10_Reclustering_MP_and_MC/markers_20230207_SCT",paste0("clusterMarkers_", res,".txt"), sep="/"), header = TRUE)
  #   } 
  #   if ( outDir == paste0(dataDir,"/markers_integrated") ){
  #     seurat.markers <- read.table(paste("10_Reclustering_MP_and_MC/markers_integrated",paste0("clusterMarkers_", res,".txt"), sep="/"), header = TRUE)
  #   } 
  #   
  #   #top5.markers <- seurat.markers %>% group_by(cluster) %>% top_n(n = -5, wt = p_val_adj)  # wt defines ordering, avg_log2FC, or p_val_adj?
  #   # Hmm, returns more than 5 rows per clusters as there are ties in the 'p_val_adj' ...
  #   # So, also sort on avg_log2FC ...
  #   top5.markers <- seurat.markers %>% group_by(cluster) %>% arrange(p_val_adj, desc(avg_log2FC), .by_group=TRUE) %>% top_n(n = 5, wt = avg_log2FC)
  #   
  #   # Sort the top5 list on the cluster column
  #   top5.markers <- dplyr::arrange(top5.markers, cluster)
  #   
  #   if ( outDir == paste0(dataDir,"/markers_integrated")){
  #     DefaultAssay(seurat) <- "integrated"
  #   } else {
  #     # Make sure to use 'SCT' and not 'integrated' as DefaultAssay, see mail Perry 06022023
  #     DefaultAssay(seurat) <- "SCT"
  #   }
  #   
  #   print(paste0("... VlnPlots: The following setting is used for DefaultAssay: ", DefaultAssay(seurat)))
  #   
  #   for ( myClusters in unique(top5.markers$cluster) ){
  #     df <- top5.markers[top5.markers$cluster == myClusters,]
  #     pViolin <- list()
  #     pFeature <- list()
  #     for ( i in 1:nrow(df)){
  #       if ( i %in% c(2,4) & !(is.null(maxY.RNA)) ){
  #         pViolin[[i]] <- VlnPlot(object = seurat, features = df$gene[i], pt.size = 0.1, y.max = maxY.RNA, same.y.lims = TRUE) + NoLegend() +
  #           theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), 
  #                 axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) 
  #         pFeature[[i]] <- FeaturePlot(object = seurat, features = df$gene[i], pt.size = 0.1, reduction = "umap" ) + # see mail Perry 06022023 NoLegend() +
  #           theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), 
  #                 axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) 
  #       } else {
  #         pViolin[[i]] <- VlnPlot(object = seurat, features = df$gene[i], pt.size = 0.1, y.max = maxY.RNA, same.y.lims = TRUE) + NoLegend() +
  #           theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) 
  #         pFeature[[i]] <- FeaturePlot(object = seurat, features = df$gene[i], pt.size = 0.1, reduction = "umap" ) + # see mail Perry 06022023 NoLegend() +
  #           theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) 
  #         
  #       }
  #     }
  #     p1 <- pViolin[[1]] | pViolin[[2]] 
  #     p2 <- pViolin[[3]] | pViolin[[4]]
  #     p3 <- pViolin[[5]] | plot_spacer()
  #     p <- p1 / p2 / p3 + plot_layout(heights = c(1, 1, 1))
  #     pdf(paste(outDir, paste(paste0("VlnPlots_", res),paste0("Cluster_", myClusters,".pdf"), sep="/"), sep="/"))
  #     print(p)
  #     dev.off()
  #     
  #     p1 <- pFeature[[1]] | pFeature[[2]] 
  #     p2 <- pFeature[[3]] | pFeature[[4]]
  #     p3 <- pFeature[[5]] | plot_spacer()
  #     p <- p1 / p2 / p3 + plot_layout(heights = c(1, 1, 1))
  #     pdf(paste(outDir, paste(paste0("FeaturePlots_", res),paste0("Cluster_", myClusters,".pdf"), sep="/"), sep="/"))
  #     print(p)
  #     dev.off()
  #     
  #   }
  # }
}


# And put the variable 'outDir' back to its original value
outDir <- baseOutDir


#### Wensenlijst 3.0 - 20231016 ####

outDir <- "finalFigures"
dir.create(outDir)

dir.create(paste0(outDir,"/20231016_Wensenlijst"))
figDir <- paste0(outDir,"/20231016_Wensenlijst")

# And create a directory for holding the TIFFs?
dir.create(paste0(figDir,"/tiff"))

#### Markerplot @res1.0 ####
# See '10_reclustering_MP_and_MC.r' - line 1139
#
# Perhaps interesting?
# - https://github.com/satijalab/seurat/issues/2201
library(Seurat)
library(dplyr)
library(tidyr)
library(cowplot)
library(ggplot2)

# Color function
# Also get the correct colors
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

# Read the data
dataDir <- "10_Reclustering_MP_and_MC/"
seurat <- readRDS(paste(dataDir, "selectedClusters.subjectIntegrated.reclustered.rds", sep="/"))

# Set the default assay slot correct
DefaultAssay(seurat) <- "SCT"

# We are going to show the marker plot for resolution 1 and only cluster 0, 1, 4, 8, 9 and 11
# and thus also only the top5 marker for these clusters and add some other ones
seurat.markers <- read.table(paste(dataDir,paste0("markers_20230418_LR/clusterMarkers_1.txt"), sep="/"), header = TRUE)

# Give top 5 per comparison based on p_adj_val (need top_n(n=-5) to get the right 5 genes!!)
#top5.markers <- seurat.markers %>% group_by(cluster) %>% top_n(n = -5, wt = p_val_adj)  # wt defines ordering, avg_log2FC, or p_val_adj?
# Hmm, returns more than 5 rows per clusters as there are ties in the 'p_val_adj' ...
# So, also sort on avg_log2FC ...
top5.markers <- seurat.markers %>% group_by(cluster) %>% dplyr::arrange(p_val_adj, desc(avg_log2FC), .by_group=TRUE) %>% top_n(n = 5, wt = avg_log2FC)

# Select only the genes in the cluster 0,1,4,8,9,11
top5.markers <- top5.markers[top5.markers$cluster %in% c(0,1,4,8,9,11),]

# Now add CD163, MRC1, MARCO, TREM2, CD14, PLIN2, C1QC and SPP1 (if they are not already in there)
# AJ 20230606 - Not C1QA as previously??
# Wenselijst 3.0.pptx (mail 20231016) - add LYPD2, CD31 (PECAM1), CCR2, CD81, APOE,FOLR2, LPL, CD36, FN1, HSP90AA1
for (markerGene in c("CD163", "MRC1", "MARCO", "TREM2", "PLIN2", "CD14", "C1QC", "SPP1",
                     "LYPD2", "CD31", "CCR2", "CD81", "APOE","FOLR2", "LPL", "CD36", "FN1", "HSP90AA1")){
  if ( !(markerGene %in% top5.markers$gene)){
    print(paste0("... Now adding ", markerGene, " to the list of marker genes"))
    top5.markers <- rbind(top5.markers, seurat.markers[seurat.markers$gene == markerGene, ])
  }
}
# Note that MARCO is not in he final plot because it is not in the list of marker genes
# LYPD2 is already there?? ...

# Sort the top5 list on the cluster column
top5.markers <- dplyr::arrange(top5.markers, cluster)

# Plot a heatmap and dotplot
# resolution 1 -> only get clusters 0,1,4,8,9,11
clust = "integrated_snn_res.1"
selClusters.seurat <- subset(x = seurat, subset = integrated_snn_res.1 %in% c("0","1","4","8","9","11"))
selClusters.seurat$integrated_snn_res.1 <- factor(selClusters.seurat$integrated_snn_res.1, levels=c("0","1","4","8","9","11"))

myColors <- gg_color_hue(length(unique(seurat$integrated_snn_res.1)))
myColors = myColors[c(0,1,4,8,9,11)+1]

pdf(paste(figDir, "markerPlots_Top5_res1_LR.selectedClusters.pdf", sep="/"), width = 28, height = 28)
# Perhaps downsample to 100,200 cells??
p1 <- DoHeatmap(selClusters.seurat, features = top5.markers$gene, group.by = clust, 
                size = 28, angle=0, hjust=0.5, group.colors = myColors) + 
  theme(axis.text.y = element_text(size=28, face = "bold.italic")) +
  theme(legend.key.size = unit(6, 'cm'), 
        legend.key.height= unit(4, 'cm'),
        legend.key.width= unit(4, 'cm'),
        legend.title = element_text(size=40),legend.text = element_text(size=32) ) +
  guides(colour = guide_legend(override.aes = list(size=30, alpha =1))) 
# Warning message:
#In DoHeatmap(selClusters.seurat, features = top5.markers$gene, group.by = clust,  :
#               The following features were omitted as they were not found in the scale.data slot for the SCT assay: 
#    LPL, FN1, FCER1A, LYPD2

# These genes were not selected as variable genes and hence not used in the 'ScaleData' function
# (or similar as used in the SCTransform workflow)
# Can I add these later ????

# Get the legend
legend_p <- cowplot::get_legend(p1)
print(p1)
# Now plot it without legend
p2 <- DoHeatmap(selClusters.seurat, features = top5.markers$gene, group.by = clust, 
                size = 28, angle=0, hjust=0.5, group.colors = myColors) + 
  theme(axis.text.y = element_text(size=28, face = "bold.italic")) +
  NoLegend()
print(p2)
# And plot the legend
plot_grid(legend_p)

dev.off()

# And as TIFF
ggsave(paste(figDir, "/tiff/markerPlots_Top5_res1_LR.selectedClusters.tiff", sep="/"), p1, scale=4, width=5, height=5, dpi=360,compression="lzw")
ggsave(paste(figDir, "/tiff/markerPlots_Top5_res1_LR.selectedClusters.NoLegend.tiff", sep="/"), p2, scale=4, width=5, height=5, dpi=360,compression="lzw")
ggsave(paste(figDir, "/tiff/markerPlots_Top5_res1_LR.selectedClusters.LegendOnly.tiff", sep="/"), legend_p, scale=4, width=5, height=5, dpi=360,compression="lzw")

#### Fig 1F - Markerplot Surface proteins (ADT) with extra markers ####
# See 10_Reclustering_MP_and_MC/ADT_markerPlots_Top5_res.1.selectedClusters.pdf
# i.e. 10_Reclustering_MP_and_MC #'Find markers per cluster - UPDATE 20230418'
# Increase font names and add a legend
dataDir <- "10_Reclustering_MP_and_MC/"
seurat <- readRDS(paste(dataDir, "selectedClusters.subjectIntegrated.reclustered.rds", sep="/"))

clust <- "integrated_snn_res.1"

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

#if ( outDir == paste0(dataDir,"/markers_20230418_LR") ){
print("... Integrated ADT data: Using FindAllMarkers with latent.vars=\"orig.ident\" and test.use=\"LR\"")
Idents(adt.integrated) <- clust
adt.markers <- FindAllMarkers(adt.integrated, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25,
                              latent.vars = "orig.ident", test.use = 'LR', assay = "ADT")

# 20231027 - I now get 12 clusters .... see remark below
# Warning message:
#   glm.fit: fitted probabilities numerically 0 or 1 occurred 

## IMPORTANT ##
# We get 16 clusters out ... so, there is no match any more with the 12 clusters that we get for resolution 1.0
# on the gene expression level. But we restore that later on (line 509)

## SOLUTION - 20230721
# set Idents(adt.integrated) <- clust, i.e. integrated_snn_res.1
# Probably that will give you another figure as well ...
# Ran this code again at 20230721 and added the date in the name

#} else {
#  print("... Integrated ADT data: Using default FindAllMarkers")
#  adt.markers <- FindAllMarkers(adt.integrated, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25 , assay = "ADT")
#}
top5.adt_markers <- adt.markers %>% group_by(cluster) %>% dplyr::arrange(p_val_adj, desc(avg_log2FC), .by_group=TRUE) %>% top_n(n = 5, wt = avg_log2FC)
top5.adt_markers <- dplyr::arrange(top5.adt_markers, cluster)

# 20231027 - Add CD14, CD163 and CX3CR1 (see mail Yosta 20231026)
adt.markers.formatted <- adt.markers %>% group_by(cluster) %>% dplyr::arrange(p_val_adj, desc(avg_log2FC), .by_group=TRUE)
top5.adt_markers <- rbind(top5.adt_markers, adt.markers.formatted[grepl("CD14$|CD163|CX3CR1", adt.markers.formatted$gene),])

selClusters.adt <- subset(x = adt.integrated, subset = integrated_snn_res.1 %in% c("0","1","4","8","9","11"))
selClusters.adt$integrated_snn_res.1 <- factor(selClusters.adt$integrated_snn_res.1, levels=c("0","1","4","8","9","11"))

# And also make the selection in the top ADT markers
top5.adt_markers <- top5.adt_markers[top5.adt_markers$cluster %in% c("0","1","4","8","9","11"),]
top5.adt_markers$cluster <- factor(top5.adt_markers$cluster, levels=c("0","1","4","8","9","11"))
# and order by cluster
top5.adt_markers <- dplyr::arrange(top5.adt_markers, cluster)

# Also get the correct colors
myColors <- gg_color_hue(length(unique(seurat$integrated_snn_res.1)))
myColors <- myColors[c(0,1,4,8,9,11)+1]

#pdf(paste(figDir, paste0("ADT_markerPlots_Top5_res.1.selectedClusters.pdf"), sep="/"), width = 28, height = 28)
#pdf(paste(figDir, paste0("ADT_markerPlots_Top5_res.1.selectedClusters.20230721.pdf"), sep="/"), width = 28, height = 28)
pdf(paste(figDir, paste0("Fig1F_ADT_markerPlots_Top5_res.1.selectedClusters.20231027.pdf"), sep="/"), width = 28, height = 28)
p1 <- DoHeatmap(selClusters.adt, features = top5.adt_markers$gene, assay = "integrated", group.by = clust,
                size=12, angle=0, hjust=0.5, group.colors = myColors) + 
  scale_y_discrete(labels=function(x)gsub("anti-human ", "", x)) + 
  theme(axis.text.y = element_text(size=36, face = "bold.italic")) +
  theme(legend.key.size = unit(4, 'cm'), 
        legend.key.height= unit(3, 'cm'),
        legend.key.width= unit(2, 'cm'),
        legend.title = element_text(size=24),legend.text = element_text(size=24) ) +
  guides(colour = guide_legend(override.aes = list(size=18, alpha=1)))
# Get the legend
legend_p <- cowplot::get_legend(p1)
print(p1)
# Now plot it without legend
p2 <- DoHeatmap(selClusters.adt, features = top5.adt_markers$gene, assay = "integrated",group.by = clust,
                size=12, angle=0, hjust=0.5, group.colors = myColors) + 
  scale_y_discrete(labels=function(x)gsub("anti-human ", "", x)) + 
  theme(axis.text.y = element_text(size=36, face = "bold.italic")) +
  NoLegend()  
print(p2)
# And plot the legend
plot_grid(legend_p)

dev.off()

# And as TIFF
ggsave(paste(figDir, "/tiff/Fig_1F_ADT_markerPlots_Top5_res.1.selectedClusters.20231027.tiff", sep="/"), p1, scale=4, width=5, height=5, dpi=360,compression="lzw")
ggsave(paste(figDir, "/tiff/Fig_1F_ADT_markerPlots_Top5_res.1.selectedClusters.NoLegend.20231027.tiff", sep="/"), p2, scale=4, width=5, height=5, dpi=360,compression="lzw")
ggsave(paste(figDir, "/tiff/Fig_1F_ADT_markerPlots_Top5_res.1.selectedClusters.LegendOnly.20231027.tiff", sep="/"), legend_p, scale=4, width=5, height=5, dpi=360,compression="lzw")

#### Markerplot Surface proteins (ADT) - macrophage markers  ####
# Mail Yosta 20072023
# - She wants a list of 59 markers not just the top 5
#   See Yosta_Vegting\scRNASeq\Data\Dataset_1\Meta\TotalSeq_C_Human_Universal_Cocktail_v1_137_Antibodies_399905_Barcodes.xlsx - the markers in yellow

# See 10_Reclustering_MP_and_MC/ADT_markerPlots_Top5_res.1.selectedClusters.pdf
# i.e. 10_Reclustering_MP_and_MC #'Find markers per cluster - UPDATE 20230418'
# Increase font names and add a legend
dataDir <- "10_Reclustering_MP_and_MC/"
seurat <- readRDS(paste(dataDir, "selectedClusters.subjectIntegrated.reclustered.rds", sep="/"))

clust <- "integrated_snn_res.1"

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

# Now, get the 59 markers that Yosta wants
# Does not work on Mac, could not get 'rJava' compiled ...
# Used code to generate the list and save that on Windows
# Or perhaps use 'tidyxl' package??
# See https://higgi13425.github.io/medical_r/posts/2021-01-13-extracting-highlighting-as-data-from-excel/
library(tidyxl)
wb <- xlsx_cells("../../../Data/Dataset_1/Meta/TotalSeq_C_Human_Universal_Cocktail_v1_137_Antibodies_399905_Barcodes.xlsx")
wb.format <- xlsx_formats("../../../Data/Dataset_1/Meta/TotalSeq_C_Human_Universal_Cocktail_v1_137_Antibodies_399905_Barcodes.xlsx")

location <- wb %>% 
  filter(local_format_id %in% 
           which(wb.format$local$fill$patternFill$fgColor$rgb == "FFFFFF00")) %>%
  select(address, data_type)

location
# This provides the cells

# Now get the value of these cells
genes <- wb %>% 
  filter(address %in% location$address) %>% 
  select(address, character)

# And get the name of the ADT in the rows these cells are in
# Actually, column C, so replace every 'G' in the location table with a 'C' ?
location$address <- gsub('G','C', location$address)
adt.markers <- as.vector(wb %>% 
  filter(address %in% location$address) %>% 
  select(character))

adt.markers <- unlist(adt.markers)
names(adt.markers) <- adt.markers

## On Windows ##
if (Sys.info()['sysname'] == "Windows"){
library(xlsx)

wb     <- loadWorkbook("../../../Data/Dataset_1/Meta/TotalSeq_C_Human_Universal_Cocktail_v1_137_Antibodies_399905_Barcodes.xlsx") 
sheet1 <- getSheets(wb)[[1]]
rows   <- getRows(sheet1)
cells  <- getCells(rows)
styles <- sapply(cells, getCellStyle)
cellColor <- function(style) 
{
  fg  <- style$getFillForegroundXSSFColor()
  rgb <- tryCatch(fg$getRgb(), error = function(e) NULL)
  rgb <- paste(rgb, collapse = "")
  return(rgb)
}
mycolor <- (yellow = "ffff00")
m     <- match(sapply(styles, cellColor), mycolor)

idx <- which(m==1)
length(idx)
#[1] 73

# And get the values
adt.markers <- unlist(lapply(cells[idx], getCellValue))
names(adt.markers) <- adt.markers
# These are the gene names, unfortunately not the rownames of the ADT's...
# Have to get these from the XLSX
df <- read.xlsx("../../../Data/Dataset_1/Meta/TotalSeq_C_Human_Universal_Cocktail_v1_137_Antibodies_399905_Barcodes.xlsx", 1,header=TRUE,colClasses=c("character"),encoding="UTF-8")
adt.markers <- df[df$`Gene.name` %in% adt.markers, "Description"]
length(adt.markers)
# [1] 74     -- Ah, funny, 1 more ! Why is this?
}


# Change some of the 'special characters'
# https://stackoverflow.com/questions/20916070/importing-an-excel-file-with-greek-characters-into-r-in-the-correct-encoding
tolatin=function(compname) {
  n=as.character(compname,encoding="UTF-8")
  n=gsub("\u03B1","alpha",n)
  n=gsub("\u03B2","beta",n)
  n=gsub("\u03B3","gamma",n)
  n=gsub("\u03B4","delta",n)
  n=gsub("\u03B5","epsilon",n)
  n
}

adt.markers <- tolatin(adt.markers)

# NOTE:
# Something might have gone wrong earlier in naming the ADT's
# FCERIA has been denoted as 'Fc eta RI beta' ...
adt.markers <- gsub("FcepsilonRIalpha", "Fc eta RI beta", 
                    gsub("Integrin alphaE", "Integrin E", adt.markers))

write.table(adt.markers, paste0(figDir,"/adt_markers.txt"), quote=FALSE)

# Subset the clusters
selClusters.adt <- subset(x = adt.integrated, subset = integrated_snn_res.1 %in% c("0","1","4","8","9","11"))
selClusters.adt$integrated_snn_res.1 <- factor(selClusters.adt$integrated_snn_res.1, levels=c("0","1","4","8","9","11"))

# Also get the correct colors
myColors <- gg_color_hue(length(unique(seurat$integrated_snn_res.1)))
myColors <- myColors[c(0,1,4,8,9,11)+1]

# But now there would be no clustering ....
# https://www.biostars.org/p/9564526/
# https://bioinformatics.stackexchange.com/questions/8973/doheatmap-hierarchical-clustering-seurat
#if ( outDir == paste0(dataDir,"/markers_20230418_LR") ){
print("... Integrated ADT data: Using FindAllMarkers with latent.vars=\"orig.ident\" and test.use=\"LR\"")
Idents(adt.integrated) <- clust     # Essential!!
all.markers <- FindAllMarkers(adt.integrated, only.pos = FALSE, min.pct = 0.0, logfc.threshold = 0.0,
                              latent.vars = "orig.ident", test.use = 'LR', assay = "ADT", return.thresh = 1.0)

# Do not use any cutoff, just to get all back ...

## IMPORTANT ##
# You have to correctly set the ident!!
# Otherwise you will get 16 clusters out ... and there is no match any more with the 12 clusters that we get for resolution 1.0
# on the gene expression level. 

#} else {
#  print("... Integrated ADT data: Using default FindAllMarkers")
#  adt.markers <- FindAllMarkers(adt.integrated, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25 , assay = "ADT")
#}
all.markers <- all.markers %>% group_by(cluster) %>% dplyr::arrange(desc(avg_log2FC), .by_group=TRUE)

all.markers <- dplyr::arrange(all.markers, cluster)

# And also make the selection in the top ADT markers
all.markers <- all.markers[all.markers$cluster %in% c("0","1","4","8","9","11"),]
all.markers$cluster <- factor(all.markers$cluster, levels=c("0","1","4","8","9","11"))
# and order by cluster
all.markers <- dplyr::arrange(all.markers, cluster)

# Is every selected ADT marker present?
match(adt.markers, all.markers$gene)
# [1]   1   7   3 132 131   4   6 135   5 136 123   2 110   8 137  12 115 126 109 114  14 111  10  16 118 122 100 120 105 117 128 106  17
#[34] 129 108  25 104 127  15 103 124  88  95  82 101  19  32  84  30  24  43  38  97  67  22  18  21  27  40  49  51  55  58  60  63  68
#[67]  72  73  75  76  85  86  91  99
# All found ... :-)

# And now put the selected markers in this order and put the NA's last??
adt.markers <- adt.markers[order(match(adt.markers, all.markers$gene))]
adt.markers[1:5]

pdf(paste(figDir, paste0("ADT_markerPlots_MacrophageMarkers_res.1.selectedClusters.pdf"), sep="/"), width = 28, height = 28)
p1 <- DoHeatmap(selClusters.adt, features = adt.markers, assay = "integrated", group.by = clust,
                size=12, angle=0, hjust=0.5, group.colors = myColors) + 
  scale_y_discrete(labels=function(x)gsub("anti-human ", "", x)) + 
  theme(axis.text.y = element_text(size=28, face = "bold.italic")) +
  theme(legend.key.size = unit(4, 'cm'), 
        legend.key.height= unit(3, 'cm'),
        legend.key.width= unit(2, 'cm'),
        legend.title = element_text(size=24),legend.text = element_text(size=24) ) +
  guides(colour = guide_legend(override.aes = list(size=18, alpha=1)))
# In DoHeatmap(selClusters.adt, features = adt.markers, assay = "integrated",  :
#    The following features were omitted as they were not found in the scale.data slot for the integrated assay: 
#    anti-human CD127 (IL-7Ra), anti-human CD124 (IL-4Ra), anti-human FceRIa, anti-human CD119 (IFN-<U+03B3> R a chain), 
#    anti-human/mouse integrin ?7, anti-human CD103 (Integrin aE)
# Has to do with special characters ... see code above to remedy this!!!

# Get the legend
legend_p <- cowplot::get_legend(p1)
print(p1)
# Now plot it without legend
p2 <- DoHeatmap(selClusters.adt, features = adt.markers, assay = "integrated",group.by = clust,
                size=12, angle=0, hjust=0.5, group.colors = myColors) + 
  scale_y_discrete(labels=function(x)gsub("anti-human ", "", x)) + 
  theme(axis.text.y = element_text(size=28, face = "bold.italic")) +
  NoLegend()  
print(p2)
# And plot the legend
plot_grid(legend_p)

dev.off()

# And as TIFF
ggsave(paste(figDir, "/tiff/ADT_markerPlots_MacrophageMarkers_res.1.selectedClusters.tiff", sep="/"), p1, scale=4, width=5, height=5, dpi=360,compression="lzw")
ggsave(paste(figDir, "/tiff/ADT_markerPlots_MacrophageMarkers_res.1.selectedClusters.NoLegend.tiff", sep="/"), p2, scale=4, width=5, height=5, dpi=360,compression="lzw")
ggsave(paste(figDir, "/tiff/ADT_markerPlots_MacrophageMarkers_res.1.selectedClusters.LegendOnly.tiff", sep="/"), legend_p, scale=4, width=5, height=5, dpi=360,compression="lzw")

#### Stacked bar plot immune cells acc. Stewart @0.6 ####
# MNP/Neutrophil, DC, Mast cell as 'myeloid cells' 
# see UMAP CD45+ cells #

# Use celltypes instead of numbers ... T cells/ B cells/ myeloid cells/ granulocytes etc.
# The figure corresponds to 9_kidney_reference_mapping/umap_query_kidney_immune_cluster_res_0.6.pdf

# Read the data
dataDir <- "9_kidney_reference_mapping"
reference <- readRDS(paste(dataDir,"kidney_immune_reference.rds", sep="/"))
query1 <- readRDS(paste(dataDir,"Mito10.integrated.clustered.annot.kidney_immune_reference.rds", sep="/"))

# Make the levels order the same
query1@meta.data$predicted.celltype <- factor(query1@meta.data$predicted.celltype, 
                                              levels = levels(factor(reference@meta.data$celltype)))

# Collapse MNP, Neutrophils, DC and Mast cells into 'myeloid cells'
# And remove the corresponding colors from refColors (myeloid cells taking the first color)
testVec <- as.vector(reference@meta.data$celltype)
reference@meta.data$celltype <- ifelse(grepl('^MNP|Neutrophil|dendritic|Mast',testVec), 'Myeloid cell', testVec)
rm(testVec)

reference@meta.data$celltype <- factor(reference@meta.data$celltype,
                                       levels = levels(factor(reference@meta.data$celltype)))

# Get the colors; am I using those???
refColors <- gg_color_hue(length(levels(factor(reference@meta.data$celltype))))

# Make the levels order the same
testVec <- as.vector(query1@meta.data$predicted.celltype)
query1@meta.data$predicted.celltype <- ifelse(grepl('^MNP|Neutrophil|dendritic|Mast',testVec), 'Myeloid cell', testVec)
rm(testVec)

query1@meta.data$predicted.celltype <- factor(query1@meta.data$predicted.celltype, 
                                              levels = levels(factor(reference@meta.data$celltype)))

pdf(paste(figDir,"umap_query_kidney_immune_cluster_res_0.6.pdf", sep="/"), width=14)
p1 <- DimPlot(object = query1, reduction = "ref.umap", group.by = "integrated_snn_res.0.6", label = FALSE,
              label.size = 8, repel = TRUE) + NoLegend() + ggtitle("Cluster labels (resolution = 0.6)")
p1 <- LabelClusters(p1, id = "integrated_snn_res.0.6",  fontface = "bold")
p2 <- DimPlot(object = query1, reduction = "ref.umap", group.by = "predicted.celltype", label = FALSE,
              label.size = 6, repel = TRUE, cols = refColors) + NoLegend() + ggtitle("Cluster labels\n(resolution = 0.6, query transferred labels)") 
p2 <- LabelClusters(p2, id = "predicted.celltype",  fontface = "bold")
p <- p1 + p2
print(p)
print(p1)
print(p2)

# And make them without labels and get the legends
p1b <- DimPlot(object = query1, reduction = "ref.umap", group.by = "integrated_snn_res.0.6", label = FALSE,
               label.size = 8, repel = TRUE) + ggtitle("Cluster labels (resolution = 0.6)")
print(p1b)
# Get the legend
legend_p1 <- cowplot::get_legend(p1b)
p1c <- DimPlot(object = query1, reduction = "ref.umap", group.by = "integrated_snn_res.0.6", label = FALSE,
               label.size = 8, repel = TRUE) +  NoLegend() + ggtitle("Cluster labels (resolution = 0.6)")
print(p1c)
plot_grid(legend_p1)
#
p2b <- DimPlot(object = query1, reduction = "ref.umap", group.by = "predicted.celltype", label = FALSE,
               label.size = 6, repel = TRUE, cols = refColors) + ggtitle("Cluster labels\n(resolution = 0.6, query transferred labels)") 
print(p2b)
legend_p2 <- cowplot::get_legend(p2b)
p2c <- DimPlot(object = query1, reduction = "ref.umap", group.by = "predicted.celltype", label = FALSE,
               label.size = 6, repel = TRUE, cols = refColors) + NoLegend() + ggtitle("Cluster labels\n(resolution = 0.6, query transferred labels)") 
print(p2c)
plot_grid(legend_p2)

dev.off()

# And as TIFF
ggsave(paste(figDir, "/tiff/umap_query_kidney_immune_cluster_res_0.6.all.tiff", sep="/"), p, scale=2, width=5, height=5, dpi=360,compression="lzw")
ggsave(paste(figDir, "/tiff/umap_query_kidney_immune_cluster_res_0.6.clusters.tiff", sep="/"), p1, scale=1.5, width=5, height=5, dpi=360,compression="lzw")
ggsave(paste(figDir, "/tiff/umap_query_kidney_immune_cluster_res_0.6.transferredLabels.tiff", sep="/"), p2, scale=1.5, width=5, height=5, dpi=360,compression="lzw")
ggsave(paste(figDir, "/tiff/umap_query_kidney_immune_cluster_res_0.6.clusters_NoLabel.tiff", sep="/"), p1b, scale=1.5, width=5, height=5, dpi=360,compression="lzw")
ggsave(paste(figDir, "/tiff/umap_query_kidney_immune_cluster_res_0.6.clusters_NoLabel_NoLegend.tiff", sep="/"), p1c, scale=1.5, width=5, height=5, dpi=360,compression="lzw")
ggsave(paste(figDir, "/tiff/umap_query_kidney_immune_cluster_res_0.6.clusters_Legend.tiff", sep="/"), legend_p1, scale=1.25, width=5, height=5, dpi=360,compression="lzw")
ggsave(paste(figDir, "/tiff/umap_query_kidney_immune_cluster_res_0.6.transferredLabels_NoLabel.tiff", sep="/"), p2b, scale=1.5, width=5, height=5, dpi=360,compression="lzw")
ggsave(paste(figDir, "/tiff/umap_query_kidney_immune_cluster_res_0.6.transferredLabels_NoLabel_NoLegend.tiff", sep="/"), p2c, scale=1.5, width=5, height=5, dpi=360,compression="lzw")
ggsave(paste(figDir, "/tiff/umap_query_kidney_immune_cluster_res_0.6.transferredLabels_Legend.tiff", sep="/"), legend_p2, scale=0.75, width=5, height=5, dpi=360,compression="lzw")

pdf(paste(figDir,"umap_query_kidney_immune_cluster_res_0.6_v2.pdf", sep="/"))
# Original code 9_kidney_reference_mapping.r
# p1 <- DimPlot(query1, reduction = "ref.umap", group.by = "predicted.celltype", label = TRUE,
#               label.size = 3, repel = TRUE, cols = queryColors) + NoLegend() + ggtitle("Query transferred labels")
# p2 <- DimPlot(object = query1, reduction = "ref.umap", group.by = "integrated_snn_res.0.6", label = TRUE,
#               label.size = 5, repel = TRUE) + NoLegend() + ggtitle("Cluster labels (resolution = 0.6)")
#p <- p1 + p2
p1 <- DimPlot(object = query1, reduction = "ref.umap", group.by = "integrated_snn_res.0.6", label = FALSE,
              label.size = 8, repel = TRUE) + NoLegend() + ggtitle("Cluster labels (resolution = 0.6)")
p1 <- LabelClusters(p1, id = "integrated_snn_res.0.6",  fontface = "bold")
p2 <- DimPlot(object = query1, reduction = "ref.umap", group.by = "predicted.celltype", label = FALSE,
              label.size = 6, repel = TRUE, cols = refColors) + NoLegend() + ggtitle("Cluster labels\n(resolution = 0.6, query transferred labels)") 
p2 <- LabelClusters(p2, id = "predicted.celltype",  fontface = "bold")
#p <- p1 + p2
#print(p)
print(p1)
print(p2)

# And make them without labels and get the legends
p1b <- DimPlot(object = query1, reduction = "ref.umap", group.by = "integrated_snn_res.0.6", label = FALSE,
               label.size = 8, repel = TRUE) + ggtitle("Cluster labels (resolution = 0.6)")
print(p1b)
# Get the legend
legend_p1 <- cowplot::get_legend(p1b)
p1c <- DimPlot(object = query1, reduction = "ref.umap", group.by = "integrated_snn_res.0.6", label = FALSE,
               label.size = 8, repel = TRUE) +  NoLegend() + ggtitle("Cluster labels (resolution = 0.6)")
print(p1c)
plot_grid(legend_p1)
#
p2b <- DimPlot(object = query1, reduction = "ref.umap", group.by = "predicted.celltype", label = FALSE,
               label.size = 6, repel = TRUE, cols = refColors) + ggtitle("Cluster labels\n(resolution = 0.6, query transferred labels)") 
print(p2b)
legend_p2 <- cowplot::get_legend(p2b)
p2c <- DimPlot(object = query1, reduction = "ref.umap", group.by = "predicted.celltype", label = FALSE,
               label.size = 6, repel = TRUE, cols = refColors) + NoLegend() + ggtitle("Cluster labels\n(resolution = 0.6, query transferred labels)") 
print(p2c)
plot_grid(legend_p2)

dev.off()

#### ... and the corresponding distribution ####
query1@meta.data$Condition <- ifelse(query1@meta.data$orig.ident %in% c("MOMA17", "MOMA52", "MOMA67", "MOMA302", "MOMA72"), "ANCA",
                                     ifelse(query1@meta.data$orig.ident == "MOMA57", "Control", "SLE"))

query1@meta.data$subType <- ifelse(query1@meta.data$orig.ident %in% c("MOMA17", "MOMA52", "MOMA302"), "PR3",
                                   ifelse(query1@meta.data$orig.ident %in% c("MOMA57", "MOMA68"), "NONE", "MPO"))

query1@meta.data$sampleAnnot <- paste0(gsub("MOMA", "s", query1@meta.data$orig.ident),"_",query1@meta.data$Condition,
                                       "_", query1@meta.data$subType)

t <- table(query1@meta.data[,"predicted.celltype"], query1$sampleAnnot)
write.table(t, paste(figDir,"clusters_per_subject_freq_predictedCelltype.txt", sep="/"), quote=FALSE, sep="\t", row.names = TRUE, col.names = TRUE)

t.df <- as.data.frame.matrix(t)
t.df$Cluster = rownames(t.df)
t.df <- t.df[,c(8,1:7)]
pct.df <- as.data.frame.matrix(round(proportions(as.matrix(t.df[,2:8]),1) * 100,2))
colnames(pct.df) <- paste0("pct.",colnames(pct.df))
t.df <- cbind(t.df, pct.df)

# And normalize over columns
norm.pct.df <- (scale(t.df[,9:15], center=FALSE, scale=colSums(as.data.frame(t.df[,9:15])))) * 100
colnames(norm.pct.df) <- paste0("norm.",colnames(pct.df))
t.df <- cbind(t.df, norm.pct.df)
t.df[nrow(t.df)+1,] <- c("totals", colSums(as.data.frame(t.df[,2:8])), rep("-",7), colSums(norm.pct.df[,1:7]))
rownames(t.df)[7] <- "totals"

# And save
write.table(t.df, file=paste(figDir, "dimPlot_pctPerCelltype.txt", sep="/"), quote = FALSE,
            row.names = FALSE, col.names = TRUE, sep="\t")

t.df$cluster = rownames(t.df)
t.df.longer <- pivot_longer(data=t.df[1:6,], cols = starts_with("norm"), names_to = "subject", values_to = "percent")
df2 <- t.df.longer[,c('cluster', 'subject', 'percent')]

# Get the right colors per cluster
#myColors <- gg_color_hue(length(unique(df2$cluster)))
df2$colors <- refColors[factor(df2$cluster)]
df2$percent <- as.numeric(df2$percent)

# Rename: AGN1_PR3, AGN2_PR3, AGN3_PR3, AGN4_MPO, AGN5_MPO, LN, NC
# s17_ANCA_PR3     == s1_ANCA_PR3 == AGN1_PR3
# s302_ANCA_PR3    == s2_ANCA_PR3 == AGN2_PR3
# s52_ANCA_PR3     == s3_ANCA_PR3 == AGN3_PR3
# s72_ANCA_MPO     == s4_ANCA_MPO == AGN4_MPO
# s67_ANCA_MPO     == s5_ANCA_MPO == AGN5_MPO
# s68_SLE_NONE     == SLE         == LN
# s57_Control_NONE == HC          == NC

df2$subject <- gsub("norm.pct\\.", "", df2$subject)
df2$subject <- factor(df2$subject, 
                      levels = c("s17_ANCA_PR3", "s302_ANCA_PR3","s52_ANCA_PR3", 
                                 "s72_ANCA_MPO", "s67_ANCA_MPO", "s68_SLE_NONE", "s57_Control_NONE"),
                      labels = c("AGN1_PR3", "AGN2_PR3","AGN3_PR3", 
                                 "AGN4_MPO", "AGN5_MPO", "LN", "NC"))

pdf(paste(figDir, paste0("distribution_per_subject_stacked_percentages_cluster_res0_6.pdf"), sep="/"), width = 14, height =14)
p <- ggplot(df2, aes(x=subject, y=percent, fill=cluster)) +
  geom_bar(stat='identity', position="stack") +
  xlab("") + ylab("Proportion") +
  scale_fill_manual("Cell type", values = refColors) +
  guides(fill=guide_legend(ncol=1)) +
  theme_cowplot() +
  theme(axis.text.x=element_text(angle=90,hjust=1, size=24)) +
  theme(axis.text.y=element_text(size=24)) +
  theme(axis.title.y = element_text(size=24)) +
  theme(legend.title = element_text(size=28), #change legend title font size
        legend.text = element_text(size=24)) #change legend text font size
print(p)
# Get the legend as a separate figure
legend_p <- cowplot::get_legend(p)
p2 <- ggplot(df2, aes(x=subject, y=percent, fill=cluster)) +
  geom_bar(stat='identity', position="stack") +
  xlab("") + ylab("Proportion") +
  scale_fill_manual("Cell type",values = refColors) +
  guides(fill=guide_legend(ncol=1)) +
  theme_cowplot() +
  theme(axis.text.x=element_text(angle=90,hjust=1, size=18)) +
  theme(axis.text.y=element_text(size=16)) +
  theme(axis.title.y = element_text(size=18)) +
  NoLegend()
print(p2)
# And print the legend
plot_grid(legend_p)
dev.off()

# And as TIFF
ggsave(paste(figDir, "/tiff/distribution_per_subject_stacked_percentages_clusters_res0_6.tiff", sep="/"), p, scale=2, width=5, height=5, dpi=360,compression="lzw")
ggsave(paste(figDir, "/tiff/distribution_per_subject_stacked_percentages_clusters_res0_6.NoLegend.tiff", sep="/"), p2, scale=2, width=5, height=5, dpi=360,compression="lzw")
ggsave(paste(figDir, "/tiff/distribution_per_subject_stacked_percentages_clusters_res0_6.LegendOnly.tiff", sep="/"), legend_p, scale=0.75, width=5, height=5, dpi=360,compression="lzw")

#### dimPlots after reclustering with patient info ####
# AJ - 22012024: Reviewers want to have UMAPs per patient to compare
# Look at the code after integration of the complete data set, i.e 6_after_integration.r and
# after reclustering, 10_reclustering_MP_and_MC.r
dir.create(paste0(outDir,"/20240122"))
figDir <- paste(outDir,"20240122", sep="/")

# Get the data from this directory
dataDir <- "10_Reclustering_MP_and_MC"

# Libraries
library(Seurat)
library(dplyr)
library(tidyr)
library(cowplot)
library(ggplot2)

# Color function
# Also get the correct colors
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

allSets.integrated <- readRDS(paste(dataDir, "selectedClusters.subjectIntegrated.reclustered.rds", sep="/"))
# allSets.integrated@assays$SCT@var.features
# logical(0)

# The variable features of this assay are automatically set during IntegrateData
DefaultAssay(allSets.integrated) <- "integrated"

# Create a directory to hold TIFFs
#dir.create(paste0(outDir,"/tiff"))

pdf(paste(figDir, "dimPlot_byPatient_afterReclustering_res1.pdf", sep="/"))
Idents(allSets.integrated) <- "sampleAnnot"
# Get common colors
refColors <- gg_color_hue(length(levels(factor(allSets.integrated@meta.data$integrated_snn_res.1))))

p <- DimPlot(allSets.integrated, reduction = "umap", label = TRUE,  repel = TRUE,
             group.by = 'integrated_snn_res.1', cols = refColors) +
  ggtitle("All samples") + 
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 3) ) )
p[[1]]$layers[[1]]$aes_params$alpha = 0.6
p[[1]]$layers[[1]]$mapping$alpha = 0.4
print(p)

p <- DimPlot(allSets.integrated, reduction = "umap", label = TRUE,  repel = TRUE,
             group.by = 'sampleAnnot') +
  ggtitle("All samples") + 
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 3) ) )
p[[1]]$layers[[1]]$aes_params$alpha = 0.6
p[[1]]$layers[[1]]$mapping$alpha = 0.4
print(p)

Idents(allSets.integrated) <- "Condition"
seurat.subset <- subset(allSets.integrated, idents = "ANCA")
p <- DimPlot(seurat.subset, reduction = "umap", label = TRUE,  repel = TRUE,
             group.by = 'integrated_snn_res.1', cols = refColors) +
  ggtitle("ANCA samples") + 
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 3) ) )
p[[1]]$layers[[1]]$aes_params$alpha = 0.6
p[[1]]$layers[[1]]$mapping$alpha = 0.4
print(p)

seurat.subset <- subset(allSets.integrated, idents = "SLE")
p <- DimPlot(seurat.subset, reduction = "umap", label = TRUE,  repel = TRUE,
             group.by = 'integrated_snn_res.1', cols = refColors) +
  ggtitle("LN sample") + 
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 3) ) )
p[[1]]$layers[[1]]$aes_params$alpha = 0.6
p[[1]]$layers[[1]]$mapping$alpha = 0.4
print(p)

seurat.subset <- subset(allSets.integrated, idents = "Control")
p <- DimPlot(seurat.subset, reduction = "umap", label = TRUE,  repel = TRUE,
             group.by = 'integrated_snn_res.1', cols = refColors) +
  ggtitle("HC samples") + 
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 3) ) )
p[[1]]$layers[[1]]$aes_params$alpha = 0.6
p[[1]]$layers[[1]]$mapping$alpha = 0.4
print(p)

Idents(allSets.integrated) <- "sampleAnnot"
for (sample in unique(allSets.integrated@meta.data$sampleAnnot)){
  seurat.subset <- subset(allSets.integrated, idents = sample)
  p <- DimPlot(seurat.subset, reduction = "umap", label = TRUE,  repel = TRUE,
             group.by = 'integrated_snn_res.1', cols = refColors) +
    ggtitle(sample) + 
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 3) ) )
  p[[1]]$layers[[1]]$aes_params$alpha = 0.6
  p[[1]]$layers[[1]]$mapping$alpha = 0.4
  print(p)
}
dev.off()


#### Heatmap plots comparisons ####
# see '15_pathway_visualization_20230419.r @Heatmap like plot of all comparisons line 1000

#### SCTransform with all features in the output - 20231023 ####
# As I missed some features in the 'scale.data' slot I decided to redo the SCTransform part of the pipeline.
# I copied the code from '10_reclustering_MP_and_MC.r' and added the setting 'return.only.var.genes = FALSE' to
# see whether that woud give me the genes that I now missed (LPL, FN1, FCER1A, LYPD2, see line 1576) in the output

## Read the data ##
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
    # set 'return.only.var.genes = FALSE' to get all genes in the 'scaled.data' slot back ...
    allSets.list[[i]] <- SCTransform(allSets.list[[i]], vars.to.regress = "percent.mito", verbose = FALSE,
                                     return.only.var.genes = FALSE)
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
# You need Matrix v1.5.4
# See https://github.com/satijalab/seurat/issues/6746 as I got this error using Matrix v1.6.1 ...
remove.packages('Matrix')
library(remotes)
install_version("Matrix", "1.5.4")

features <- SelectIntegrationFeatures(object.list = allSets.list, nfeatures = 3000)
allSets.list.prep <- PrepSCTIntegration(object.list = allSets.list, anchor.features = features)
anca.anchors <- FindIntegrationAnchors(object.list = allSets.list.prep, normalization.method = "SCT",
                                       anchor.features = features)
anca.combined.sct <- IntegrateData(anchorset = anca.anchors, normalization.method = "SCT")

# Use 30 dimensions ...
anca.combined.sct <- RunPCA(anca.combined.sct, verbose = FALSE, npcs = 30)
anca.combined.sct <- RunUMAP(anca.combined.sct, reduction = "pca", dims = 1:30)

# Check whether we still have all features in the 'scaled.data' slot
dim(anca.combined.sct@assays$SCT@scale.data)
# [1] 3000 4194
# Hmmm, I might need all features!!
dim(anca.combined.sct@assays$RNA@counts)
# [1] 22769  4194

# Are LPL, FN1, FCER1A, LYPD2 in the 'scale.data' slot ?
c('LPL', 'FN1', 'FCER1A', 'LYPD2') %in% rownames(anca.combined.sct@assays$SCT@scale.data)
# Nope ...
# [1] FALSE FALSE FALSE FALSE

features <- SelectIntegrationFeatures(object.list = allSets.list, nfeatures = 3000)
features <- c(features, c('LPL', 'FN1', 'FCER1A', 'LYPD2'))

allSets.list.prep <- PrepSCTIntegration(object.list = allSets.list, anchor.features = features)
# chor.features, ] : subscript out of bounds
# In addition: Warning message:
#   The following requested features are not present in any models: FCER1A, LYPD2 

# Hmm, it seems that I loose a couple of genes alreay early on ... so, why dit Yosta select these?
# On what basis; literature or our data?

# anca.anchors <- FindIntegrationAnchors(object.list = allSets.list.prep, normalization.method = "SCT",
#                                        anchor.features = features)
# anca.combined.sct <- IntegrateData(anchorset = anca.anchors, normalization.method = "SCT")
# 
# # Use 30 dimensions ...
# anca.combined.sct <- RunPCA(anca.combined.sct, verbose = FALSE, npcs = 30)
# anca.combined.sct <- RunUMAP(anca.combined.sct, reduction = "pca", dims = 1:30)
# 
# # And save
# saveRDS(anca.combined.sct, paste(outDir, "selectedClusters.subjectIntegrated.v20231023.rds", sep="/"))
# 
# # And clean
# rm(anca.anchors, anca.combined.sct, features, seurat)
