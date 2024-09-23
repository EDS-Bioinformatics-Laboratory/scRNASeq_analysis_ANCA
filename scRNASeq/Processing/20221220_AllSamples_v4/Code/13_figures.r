# AJ - 20230206

#### Set default for a good start ####
source("initAnalysis.r")
init()

# Create a directory to hold results and plots
outDir <- "13_figures"
dir.create(outDir)

#### Library ####
library(Seurat)
library(cowplot)
library(patchwork)
library(ggplot2)
library(dplyr)

#### Helper functions ####
plotPanels <- function(p, nrPanels, myTitle, nrow = 1, out = outDir){
  # Make sure the name of the list items is the desired name to be used in the plot
  pdf(paste(out,paste0(myTitle,".pdf"), sep="/"))
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

#### Marker genes per celltype ####
# see also https://satijalab.org/seurat/articles/visualization_vignette.html
#
# see mail 06022023:
# "Perry en ik hebben vorige week uitgebreid naar de sc-data gekeken en dachten dat het mooi zou zijn om
# een markerplot te maken met de belangrijkste identificatiegenen per celtype over de clusters heen, 
# dus heb een lijstje met belangrijkste genen opgesteld. Ik denk dat we op resolutie 0.6 of 1.0 gaan uitkomen.
#
# Classical monocyte-derived Macrofagen:
#   - HLA-DRB5, CD68, CD36, S100A9, LYZ, VCAN, CD14
# 
# Non-classical monocyte-derived macrofagen
# - HLA-DRB5, CD68, FCGR3A, CX3CR1, low CCR2
# 
# Tissue resident macrofaag
# - HLA-DRB5, CD68, CD74, CD81, C1Qa
# 
# cDC
# - HLA-DRB5, CD1C, FLT3, CLEC10A
# - CD14 negative 
# 
# pDC
# - HLA-DRB5, IL3RA, CLEC4C, GZMB
# 
# T cells
# - CD3E, CD4, CD8A, GZMB, GZMK, CD27
# 
# NK cells
# - NCAM1,  CD16,  CD3E negative, HLA-DRB5 pos and negative"

# and Perry 20230206:
# "Ik bedoelde dat het mogelijk interessant is om deze als signatures aan VISION toe te voegen of via 
# Seurat::AddModuleScore een gecombineerde score uit te rekenen (en daarna te visualiseren met FeaturePlot of VlnPlot)."


markerList <- list(Classical = c("HLA-DRB5", "CD68", "CD14", "CD36","S100A9","LYZ","VCAN"),
                   NonClassical = c("HLA-DRB5", "CD68", "FCGR3A", "CX3CR1", "CCR2"),
                   TissueRes = c("HLA-DRB5", "CD68", "CD74", "CD81", "C1QA"),
                   T_cell = c("CD3E", "CD4","CD8A", "GZMB", "GZMK", "CD27"),
                   cDC = c("HLA-DRB5","CD1C","FLT3","CLEC10A", "CD14"),
                   pDC = c("HLA-DRB5","IL3RA","CLEC4C","GZMB"),
                   NK = c("NCAM1", "FCGR3A", "HLA-DRB5", "CD3E")
)

#### Read the data ####
dir.create("afterReclustering")

seurat <- readRDS("10_Reclustering_MP_and_MC/selectedClusters.subjectIntegrated.reclustered.rds")
DefaultAssay(seurat) <- "integrated"

setwd("afterReclustering")

# 
for (clust in colnames(seurat@meta.data)[grep("integrated_snn", colnames(seurat@meta.data))]){
  res <- gsub("integrated_snn_","", clust)
  seurat@meta.data[,clust] <- factor(seurat@meta.data[,clust], levels=as.character(0:(length(unique(seurat@meta.data[,clust]))-1)))
  
  markerPlotDir <- paste(outDir, paste0("markerPlots_",res), sep="/")
  dir.create(markerPlotDir)
  
  maxY.RNA <- NULL
  
  # Make list with plots
  Idents(seurat) <- clust
  
  for (cellTypes in names(markerList)){
    print(cellTypes)
    
    DefaultAssay(seurat) <- "SCT"
    
    markers <- markerList[[cellTypes]]
    # Check whether all markers can be found
    markers <- markers[markers %in% rownames(seurat@assays$SCT@data)]
    
    p <- list()
    
    p[[1]] <- DimPlot(seurat, reduction = "umap", label = TRUE, label.size = 6, repel = TRUE, pt.size = 0.5) +
      guides(colour=guide_legend(override.aes = list(size=2), ncol=1))
    
    for (i in 1:length(markers)){
      #p[[i]] <- VlnPlot(object = seurat, features = markers[i], pt.size = 0.1, y.max = maxY.RNA, same.y.lims = TRUE) + NoLegend() +
      #  theme(axis.title.x = element_blank())
      p[[(i+1)]] <- FeaturePlot(seurat, features = markers[i], min.cutoff = "q05", max.cutoff = "q95", order = TRUE,
                              cols = c("lightgrey", "blue"), pt.size = 0.5)
      
    }
    myTitle <- paste0("featurePlot_RNA_of_", cellTypes,"_markers")
    plotPanels(p, 4, myTitle, nrow=2, out = markerPlotDir)
  }
}

# And return to the base directory
setwd("..")

#### All Data ####
dir.create("allData")

seurat <- readRDS("10_Reclustering_MP_and_MC/selectedClusters.subjectIntegrated.reclustered.rds")
DefaultAssay(seurat) <- "integrated"

setwd("afterReclustering")

