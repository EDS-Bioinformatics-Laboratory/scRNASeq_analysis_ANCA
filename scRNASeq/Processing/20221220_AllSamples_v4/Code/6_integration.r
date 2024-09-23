# AJ - 20230103
#

# Integration of TB cells removed, doubles removed, tubulus cells removed and SCTransformed normalized datasets - per mito content
# see: https://satijalab.org/seurat/articles/integration_introduction.html#performing-integration-on-datasets-normalized-with-sctransform-1

#### Set default for a good start ####
source("initAnalysis.r")
init()

# Create a directory to hold results and plots
outDir <- "6_integration"
dir.create(outDir)

# And create to subdirectories to hold plots from before integration
outDir_beforeIntegration <- paste(outDir,"beforeIntegration", sep ="/")
dir.create(outDir_beforeIntegration)

#### Library ####
library(Seurat)
library(ggplot2)
library(cowplot)
library(patchwork)
library(grid)

#### Helper function ####
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

# Integration of datasets per mito content
pcs <- list(Mito10=30)

#### Some plots before integration ####
for ( i in 1:length(pcs)){
  
  file_names <- dir("./5_removal_tubuluscells/", pattern = "*.rds")
  file_names <- file_names[grepl(names(pcs)[i],file_names)]
  # For the two samples with VDJ info get the object from which the TB cells were removed
  file_names_vdj <- file_names[grepl("noTB_cells",file_names)]
  file_names <- c(file_names_vdj,file_names[grepl("MOMA17|MOMA67|MOMA68|MOMA72|MOMA302",file_names)])
  
  # Get all objects in a list - 7 in total
  p <- list()
  for (sc in file_names){
    seurat <- readRDS(paste0("./5_removal_tubuluscells//",sc))
    Idents(seurat) <- "SCT_snn_res.0.3"
    
    p[[sc]] <- DimPlot(seurat, reduction = "umap", label = TRUE,  repel = TRUE) + NoLegend()
    p[[sc]] <- p[[sc]] + ggtitle(seurat@meta.data$orig.ident)
    
  }
  
  # And plot two rows of 4?
  pdf(paste(outDir_beforeIntegration,paste0("dimPlot_beforeIntegration.",names(pcs)[i],".pdf"), sep="/"), width = 21, height = 21)
    p1 <- p[[1]] + p[[2]] + p[[3]] + p[[4]]
    p2 <- p[[5]] + p[[6]] + p[[7]] + plot_spacer()
    pp <- p1 | p2
    print(pp)
  dev.off()
}


#### Batch effects w.r.t. ADTs ####
# Just check our two samples for the most important ADTs (see mail Yosta 13-12-2022, below 'Plots - ADT' and 10_reclustering_MP_and_MC.r)
# MOMA52, MOMA57

adtList <- list(MP = c("anti-human CD163","anti-human CD81 (TAPA-1)","anti-human CD14", "anti-human CD16", "anti-human HLA-DR" ),
                Interesting = c("anti-human CD88 (C5aR)","anti-human CX3CR1"),
                T_cell = c("anti-human CD3","anti-human CD4","anti-human CD8"),
                B_cell = c("anti-human CD19","anti-human CD20"),
                Neutrophil = c("anti-human CD16"),
                DC = c("anti-human CD11c", "anti-human HLA-DR","anti-human CD1c"),
                NK = c("anti-human CD56")
)

MOMA52.sc <- readRDS("./5_removal_tubuluscells/MOMA52.Mito10.noTB_cells.doubletsRemoved.annot.rds")
MOMA57.sc <- readRDS("./5_removal_tubuluscells/MOMA57.Mito10.noTB_cells.doubletsRemoved.annot.rds")

p <- list()
j <- 1

DefaultAssay(MOMA52.sc) <- "ADT"
DefaultAssay(MOMA57.sc) <- "ADT"

for (cellTypes in names(adtList)){
  for (i in 1:length(adtList[[cellTypes]])){
    p[[j]] <- FeaturePlot(MOMA52.sc, features = adtList[[cellTypes]][i], min.cutoff = "q05", max.cutoff = "q95", 
                  cols = c("lightgrey", "blue"), pt.size = 0.5) + ggtitle(label = adtList[[cellTypes]][i], subtitle = "MOMA52" )
    p[[(j+1)]] <- FeaturePlot(MOMA57.sc, features = adtList[[cellTypes]][i], min.cutoff = "q05", max.cutoff = "q95", 
                          cols = c("lightgrey", "blue"), pt.size = 0.5) + ggtitle(label = adtList[[cellTypes]][i], subtitle = "MOMA57" )
    j <- j + 2
  }
  
  myTitle <- paste0("featurePlot_ADT_beforeIntegration_", cellTypes,"_markers")
  plotPanels(p, 4, myTitle, nrow=2, out = outDir_beforeIntegration)
  
}

# And clean
rm(MOMA52.sc, MOMA57.sc, p, j, myTitle)

# @ Memory issues?? see https://github.com/satijalab/seurat/issues/1720
# But for now, just restarting R worked against the 'Calloc' error w.r.t. the integration

for ( i in 1:length(pcs)){
  
  file_names <- dir("./5_removal_tubuluscells/", pattern = "*.rds")
  file_names <- file_names[grepl(names(pcs)[i],file_names)]
  # For the two samples with VDJ info get the object from which the TB cells were removed
  file_names_vdj <- file_names[grepl("noTB_cells",file_names)]
  file_names <- c(file_names_vdj,file_names[grepl("MOMA17|MOMA67|MOMA68|MOMA72|MOMA302",file_names)])
  
  # Get all objects in a list
  scSampleList <- list()
  for (sc in file_names){
    scSampleList[[sc]] <- readRDS(paste0("./5_removal_tubuluscells//",sc))
  }
  
  features <- SelectIntegrationFeatures(object.list = scSampleList, nfeatures = 3000)
  scSampleList <- PrepSCTIntegration(object.list = scSampleList, anchor.features = features)
  # anca.anchors <- FindIntegrationAnchors(object.list = scSampleList, normalization.method = "SCT",
  #                                        anchor.features = features)  # Perhaps use 'reduction = 'rpca' '??
  #
  # As I got a memory allocation error with our normal code I used
  # https://satijalab.org/seurat/articles/integration_large_datasets.html
  # to tackle this ...
  # Use reduction="rpca" and use MOMA52 as reference -> way faster!!
  #
  # But need to compare this with results from original method (use cluster??)??
  #
  anca.anchors <- FindIntegrationAnchors(object.list = scSampleList, reference = c(2), normalization.method = "SCT",
                                         anchor.features = features, reduction = "rpca")
  anca.combined.sct <- IntegrateData(anchorset = anca.anchors, normalization.method = "SCT")
  rm(anca.anchors, features, scSampleList)
  anca.combined.sct <- RunPCA(anca.combined.sct, verbose = FALSE, npcs = pcs[i][[1]])
  anca.combined.sct <- RunUMAP(anca.combined.sct, reduction = "pca", dims = 1:pcs[i][[1]])
  
  # And save the integrated object
  saveRDS(anca.combined.sct, paste0(outDir,"/", names(pcs)[i],".integrated.rds"))
  
  # And clean
  rm(anca.combined.sct)
} 

#### Plots ####
# List with number of PCs
pcs <- list(Mito10=30)
i = 1
seurat <- readRDS(paste0(outDir,"/",names(pcs)[i],".integrated.rds"))
seurat
# An object of class Seurat 
# 46962 features across 25485 samples within 4 assays 
# Active assay: integrated (3000 features, 3000 variable features)
# 3 other assays present: RNA, SCT, ADT
# 2 dimensional reductions calculated: pca, umap

# Should I use the SCT assay as default??

# Coloring function
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

colors <- gg_color_hue(length(unique(seurat@meta.data$orig.ident)))
names(colors) <- unique(seurat@meta.data$orig.ident)

seurat@meta.data$Condition[seurat@meta.data$orig.ident %in% c("MOMA17","MOMA52","MOMA67","MOMA72","MOMA302")] <- "ANCA"
seurat@meta.data$Condition[seurat@meta.data$orig.ident %in% c("MOMA68")] <- "SLE"
seurat@meta.data$Condition[seurat@meta.data$orig.ident %in% c("MOMA57")] <- "Control"

# See also mail with Yosta 20220919
seurat@meta.data$anca.type <- "NONE"
seurat@meta.data$anca.type[seurat@meta.data$orig.ident %in% c("MOMA17","MOMA52", "MOMA302")] <- "PR3"
seurat@meta.data$anca.type[seurat@meta.data$orig.ident %in% c("MOMA67", "MOMA72")] <- "MPO"

# To  look at each condition combine seurat@meta.data$Condition and subType
seurat@meta.data$sampleAnnot <- paste(seurat@meta.data$Condition,seurat@meta.data$anca.type, sep="_")

dir.create(paste0(outDir,"/tiff"))

for ( i in 1:length(pcs)){
  # Read the data
#  seurat <- readRDS(paste0(outDir,"/",names(pcs)[i],".integrated.rds"))
  
  pdf(paste(outDir, "dimPlot_afterIntegration.pdf", sep="/"))
    Idents(seurat) <- "orig.ident"
    p <- DimPlot(seurat, reduction = "umap", label = FALSE, pt.size = 1) +
      guides(color = guide_legend(override.aes = list(alpha = 1, size = 3) ) )
    p[[1]]$layers[[1]]$aes_params$alpha = 0.6
    p[[1]]$layers[[1]]$mapping$alpha = 0.4  #
    print(p)
  dev.off()

  # And as TIFF
  ggsave(paste(outDir, "/tiff/dimPlot_afterIntegration.tiff", sep="/"),p, scale=2, width=5, height=5, dpi=360,compression="lzw")
  
  pdf(paste(outDir, "dimPlot_afterIntegration_splitBy_Subject.pdf", sep="/"), width = 14, height = 7)
    Idents(seurat) <- "orig.ident"
    p <- DimPlot(seurat, reduction = "umap", split.by = "orig.ident",label = TRUE,  repel = TRUE, pt.size =1 ) +
      guides(color = guide_legend(override.aes = list(alpha = 1, size = 3) ) )
    p[[1]]$layers[[1]]$aes_params$alpha = 0.6
    p[[1]]$layers[[1]]$mapping$alpha = 0.4
    print(p)
  dev.off()
  
  # And as TIFF
  ggsave(paste(outDir, "/tiff/dimPlot_afterIntegration_splitBy_Subject.tiff", sep="/"),p, scale=2, width=5, height=5, dpi=360,compression="lzw")
  
  for ( moma in unique(seurat@meta.data$orig.ident)){
    seurat.subset <- subset(x = seurat, subset = orig.ident == moma )
    Idents(seurat.subset) <- "orig.ident"
    
    pdf(paste(outDir, paste0("dimPlot_afterIntegration_", moma,".pdf"), sep="/"), width = 14, height = 7)
      p <- DimPlot(seurat.subset, reduction = "umap", label = FALSE,  repel = TRUE, pt.size =1, cols = colors[moma] ) +
        guides(color = guide_legend(override.aes = list(alpha = 1, size = 3) ) )
      p[[1]]$layers[[1]]$aes_params$alpha = 0.6
      p[[1]]$layers[[1]]$mapping$alpha = 0.4
      print(p)
    dev.off()
    
    # And as TIFF
    ggsave(paste(outDir, paste0("/tiff/dimPlot_afterIntegration_", moma,".tiff"), sep="/"),p, scale=2, width=5, height=5, dpi=360,compression="lzw")
  }
  
  # And per condition:
  # ANCA_MPO, ANCA_PR3, SLE_NONE, 
  for ( cond in unique(seurat@meta.data$Condition)){
    seurat.subset <- subset(x = seurat, subset = Condition == cond )
    seurat.subset@meta.data$orig.ident <- factor(seurat.subset@meta.data$orig.ident)
    colorsNew <- colors[levels(seurat.subset@meta.data$orig.ident)]
    Idents(seurat.subset) <- "orig.ident"
    
    pdf(paste(outDir, paste0("dimPlot_afterIntegration_", cond,".pdf"), sep="/"), width = 14, height = 7)
      p <- DimPlot(seurat.subset, reduction = "umap", label = FALSE,  repel = TRUE, pt.size = 1, cols = colorsNew ) +
        guides(color = guide_legend(override.aes = list(alpha = 1, size = 3) ) )
      p[[1]]$layers[[1]]$aes_params$alpha = 0.6
      p[[1]]$layers[[1]]$mapping$alpha = 0.4
      print(p)
    dev.off()
    
    # And as TIFF
    ggsave(paste(outDir, paste0("/tiff/dimPlot_afterIntegration_", cond,".tiff"), sep="/"),p, scale=2, width=5, height=5, dpi=360,compression="lzw")
    
    if (cond == "ANCA"){
      seurat.subset@meta.data$anca.type <- factor(seurat.subset@meta.data$anca.type)
      Idents(seurat.subset) <- "anca.type"
      
      pdf(paste(outDir, paste0("dimPlot_afterIntegration_", cond,"_subTypes.pdf"), sep="/"), width = 14, height = 7)
        p <- DimPlot(seurat.subset, reduction = "umap", label = FALSE,  repel = TRUE, pt.size =1, cols = c("red","blue") ) +
          guides(color = guide_legend(override.aes = list(alpha = 1, size = 3) ) )
        p[[1]]$layers[[1]]$aes_params$alpha = 0.6
        p[[1]]$layers[[1]]$mapping$alpha = 0.4
        print(p)
      dev.off()
      
      # And as TIFF
      ggsave(paste(outDir, paste0("/tiff/dimPlot_afterIntegration_", cond,"_subTypes.tiff"), sep="/"),p, scale=2, width=5, height=5, dpi=360,compression="lzw")
    
    }
    
  }
  
  # And clean
  # rm(seurat)
}  

# Color by condition and subtype, i.e. sampleAnnot
Idents(seurat) <- "sampleAnnot"
g_ANCA_PR3 <- WhichCells(seurat, idents = "ANCA_PR3")
g_ANCA_MPO <- WhichCells(seurat, idents = "ANCA_MPO")
g_SLE <- WhichCells(seurat, idents = "SLE_NONE")
g_CTRL <- WhichCells(seurat, idents = "Control_NONE")

my_colors <- gg_color_hue(4)
p <- DimPlot(seurat, label=TRUE, label.box = TRUE, repel = TRUE, group.by="sampleAnnot", pt.size = 1.5)
        # cells.highlight= list(g_ANCA_PR3, g_ANCA_MPO, g_SLE, g_CTRL), 
        # cols.highlight = my_colors, cols= "grey")
p <- p + ggtitle("Sample origin")
  p[[1]]$layers[[1]]$aes_params$alpha = 0.6
  p[[1]]$layers[[1]]$mapping$alpha = 0.4

pdf(paste(outDir, "dimPlot_afterIntegration_sampleOrigin.pdf", sep="/"), width = 14, height = 7)
  print(p)
  p1 <- p + NoLegend()
  print(p1)
  p2 <- DimPlot(seurat, label=FALSE, group.by="sampleAnnot" , pt.size = 1.5)
  # cells.highlight= list(g_ANCA_PR3, g_ANCA_MPO, g_SLE, g_CTRL), 
  # cols.highlight = my_colors, cols= "grey")
  p2 <- p2 + ggtitle("Sample origin")
  p2[[1]]$layers[[1]]$aes_params$alpha = 0.6
  p2[[1]]$layers[[1]]$mapping$alpha = 0.4
  print(p2)
dev.off()

# And as TIFF
ggsave(paste(outDir, "/tiff/dimPlot_afterIntegration_sampleOrigin.v1.tiff", sep="/"), p, scale=2, width=5, height=5, dpi=360,compression="lzw")
ggsave(paste(outDir, "/tiff/dimPlot_afterIntegration_sampleOrigin.v2.tiff", sep="/"), p1, scale=2, width=5, height=5, dpi=360,compression="lzw")
ggsave(paste(outDir, "/tiff/dimPlot_afterIntegration_sampleOrigin.v3.tiff", sep="/"), p2, scale=2, width=5, height=5, dpi=360,compression="lzw")


#### Plots - ADT ####
# Only two samples contain ADT data and it seems as if the order of plotting is determined by the order of the factor you
# have the 'Idents(seurat)' in. Hence, the last sample is 'overlayed' on the previous ones..
# In FeaturePlot put order = TRUE ?
# 
dim(seurat[['ADT']])
# [1] 137 25485

# NOTE: MRC1 == CD206, but there is no corresponding 'anti-human' ADT, so I took CD45 (i.c.w. PTPRC)
pdf(paste(outDir, "featurePlot_afterIntegration_ADT_RNA.pdf", sep="/"), width = 14, height = 7)
  p <- FeaturePlot(seurat, features = c("anti-human CD14", "anti-human CD16", "anti-human CD163", "anti-human CD45",
                                    "CD14", "FCGR3B", "CD163", "PTPRC"), min.cutoff = "q05", max.cutoff = "q95", 
            ncol = 4, cols = c("lightgrey", "blue"), pt.size = 0.5)
  print(p)
dev.off()

# And as TIFF
ggsave(paste(outDir, "/tiff/featurePlot_afterIntegration_ADT_RNA.tiff", sep="/"), p, scale=2, width=5, height=5, dpi=360,compression="lzw")

## Gleaned froom 10_reclustering_MP_and_MC.r
# For more markers see mail Yosta 13-12-2022:
# - Here only macrophages and interesting markers
#   macrophages: CD163, CD206 (= MRC1), MARCO, C1q, CD74 (Tissue-residency),
#                CD81 (Tissue-residency), CD14, CD16, HLA-DR, CD68
#   Interesting: C5aR1, C5aR2, CX3CR1, CCR2, TGF-b,	Procollagen1
#                IL1b, TNFa, MCP-1 (= CCL2), IL1R2 (erg upgereguleerd bij MPO monocyten actief), IL18R1

## Macrophage markers
c("CD163", "CD206", "MRC1", "MARCO", "C1q", "CD74","CD81", "CD14", "CD16", "HLA-DR","CD68") %in% rownames(seurat@assays$SCT@data)
# [1]  TRUE FALSE  TRUE  TRUE FALSE  TRUE  TRUE  TRUE FALSE FALSE  TRUE

# CD206 == MRC1
# CD16 == FCGR3A, FCGR3B (both present!)
# CD64 == FCGR1A (present)
# CD68 == LAMP4, SCARD1, GP110 (NOT present)
# CD71 == TFRC (present)
# CD23 == FCER2 (present)
# HLA-DR == HLA-DRA
# C1Q == C1QA

c("CD163", "MRC1", "MARCO", "C1QA", "CD74","CD81", "CD14", "FCGR3A", "FCGR3B", "HLA-DRA","CD68") %in% rownames(seurat@assays$SCT@data)
#[1] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
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

# I have put the genes and ADTs in order per celltype .. easier with plotting later on!
markerList <- list(MP = c("CD163", "CD81", "CD14", "FCGR3A", "HLA-DRA", "FCGR3B", "MRC1", "MARCO", "C1QA", "CD74","CD68"),
                   Interesting = c("C5AR1", "C5AR2", "CX3CR1", "CCR2", "TGFB1", "COL1A1","IL1B", "TNF", "CCL2", "IL1R2","IL18R1"),
                   T_cell = c("CD3E", "CD4","CD8A", "CD3D", "CD3G","CD8B"),
                   B_cell = c("CD19","MS4A1"),
                   Plasma_cell = c("SYND1", "TNFRSF17"),
                   Neutrophil = c("FCGR3A", "FCGR3B",  "FUT4", "CEACAM8"),
                   DC = c("ITGAX","HLA-DRA","CD1C","DNGR1"),
                   NK = c("NCAM1")
)

adtList <- list(MP = c("anti-human CD163","anti-human CD81 (TAPA-1)","anti-human CD14", "anti-human CD16", "anti-human HLA-DR" ),
                Interesting = c("anti-human CD88 (C5aR)","anti-human CX3CR1"),
                T_cell = c("anti-human CD3","anti-human CD4","anti-human CD8"),
                B_cell = c("anti-human CD19","anti-human CD20"),
                Neutrophil = c("anti-human CD16"),
                DC = c("anti-human CD11c", "anti-human HLA-DR","anti-human CD1c"),
                NK = c("anti-human CD56")
)

# Create a separate directory for these plots
markerPlotDir <- paste(outDir, "markerPlots", sep="/")
dir.create(markerPlotDir)

# Plotting 4 figures per page; 2 ADTs in the upper row, 2 corresponding RNAs in the lower 
# And then fill with RNA as I have more of those
for (cellTypes in names(adtList)){
  print(cellTypes)
  
  DefaultAssay(seurat) <- "ADT"
  
  adt_markers <- adtList[[cellTypes]]
  # Check whether all markers can be found
  adt_markers <- adt_markers[adt_markers %in% rownames(seurat@assays$ADT@data)]
  
  p_adt <- list()
  
  # maxY.RNA <- max(seurat@assays$RNA@data[c("CD163", "MRC1", "MARCO", "C1QA", "CD74","CD81", "CD14", "FCGR3A", "FCGR3B", "HLA-DRA","CD68"),])
  maxY.RNA <- NULL
  
  # 
  for (i in 1:length(adt_markers)){
    p_adt[[i]] <- FeaturePlot(seurat, features = adt_markers[i], min.cutoff = "q05", max.cutoff = "q95", order = TRUE,
                          cols = c("lightgrey", "blue"), pt.size = 0.5)
    
    # p[[i]] <- VlnPlot(object = seurat, features = adt_markers[i], pt.size = 0.1, y.max = maxY.RNA, same.y.lims = TRUE) + NoLegend() +
    #   theme(axis.title.x = element_blank())
  }
  
  # Get the corresponding ones from the RNA
  DefaultAssay(seurat) <- "SCT"
  
  markers <- markerList[[cellTypes]]
  # Check whether all markers can be found
  markers <- markers[markers %in% rownames(seurat@assays$SCT@data)]
  
  p_markers <- list()
  
  for (j in 1:length(markers)){
    p_markers[[j]] <- FeaturePlot(seurat, features = markers[j], min.cutoff = "q05", max.cutoff = "q95", order = TRUE,
                              cols = c("lightgrey", "blue"), pt.size = 0.5)
    
  }
  
  # Rearrange as we would like to have the ADT in the first row and the corresponding RNA below
  pp <- list()
  
  x <- seq(1,length(adt_markers),2)
  if (length(x) > 1){
   x <- x[-length(x)]
  }
  print(x)
  
  j <- 1
  for (i in x){
    pp[[j]] <- p_adt[[i]]
    if ( length(adt_markers) > 1 ){
      pp[[(j+1)]] <- p_adt[[(i+1)]]
    } else {
      pp[[(j+1)]] <- plot_spacer()
    }
    pp[[(j+2)]] <- p_markers[[i]]
    if ( length(markers) > 1 ){
      pp[[(j+3)]] <- p_markers[[(i+1)]]
    } else {
      pp[[(j+3)]] <- plot_spacer()
    }
    j <- j + 4
  }
  
  j <- length(pp) + 1
  if ( length(adt_markers)%%2 == 1 & length(adt_markers) > 1){
    pp[[j]] <- p_adt[[(i+2)]]
    pp[[(j+1)]] <- plot_spacer()
    if ( length(markers) >= length(adt_markers)) {
      pp[[(j+2)]] <- p_markers[[(i+2)]]
      pp[[(j+3)]] <- plot_spacer()
    } else {
      pp[[(j+2)]] <- plot_spacer()
      pp[[(j+3)]] <- plot_spacer()
    }
  }
  
  j <- length(pp)
  if ( length(markers) > length(adt_markers) ){
    for (i in 1:(length(markers) - length(adt_markers))){
      pp[[(j+i)]] <- p_markers[[(length(adt_markers)+i)]]
    }
  }
  
  myTitle <- paste0("featurePlot_ADT_RNA_of_", cellTypes,"_markers")
  plotPanels(pp, 4, myTitle, nrow=2, out = markerPlotDir)
  
}

for (cellTypes in names(markerList)[!(names(markerList) %in% names(adtList))]){
  print(cellTypes)
  
  DefaultAssay(seurat) <- "SCT"
  
  markers <- markerList[[cellTypes]]
  # Check whether all markers can be found
  markers <- markers[markers %in% rownames(seurat@assays$SCT@data)]
  
  p <- list()
  print(length(markers))
  for (i in 1:length(markers)){
    p[[(i)]] <- FeaturePlot(seurat, features = markers[i], min.cutoff = "q05", max.cutoff = "q95", order = TRUE,
                              cols = c("lightgrey", "blue"), pt.size = 0.5)
    
  }
  
  myTitle <- paste0("featurePlot_RNA_of_", cellTypes,"_markers")
  plotPanels(p, 2, myTitle, nrow=2, out = markerPlotDir)
  
}

# ViolinPlot w.r.t. Condition ..?
for (cellTypes in names(adtList)){
  print(cellTypes)
  
  DefaultAssay(seurat) <- "ADT"
  Idents(seurat) <- "sampleAnnot"
  
  adt_markers <- adtList[[cellTypes]]
  # Check whether all markers can be found
  adt_markers <- adt_markers[adt_markers %in% rownames(seurat@assays$ADT@data)]
  
  p_adt <- list()
  
  # maxY.RNA <- max(seurat@assays$RNA@data[c("CD163", "MRC1", "MARCO", "C1QA", "CD74","CD81", "CD14", "FCGR3A", "FCGR3B", "HLA-DRA","CD68"),])
  maxY.RNA <- NULL
  
  # 
  for (i in 1:length(adt_markers)){
    p_adt[[i]] <- VlnPlot(object = seurat, features = adt_markers[i], pt.size = 0.1, y.max = maxY.RNA, same.y.lims = TRUE) + NoLegend() +
      theme(axis.title.x = element_blank())
  }
  
  # Get the corresponding ones from the RNA
  DefaultAssay(seurat) <- "SCT"
  
  markers <- markerList[[cellTypes]]
  # Check whether all markers can be found
  markers <- markers[markers %in% rownames(seurat@assays$SCT@data)]
  
  p_markers <- list()
  
  for (j in 1:length(markers)){
    p_markers[[j]] <- VlnPlot(object = seurat, features = markers[j], pt.size = 0.1, y.max = maxY.RNA, same.y.lims = TRUE) + NoLegend() +
      theme(axis.title.x = element_blank())
    
  }
  
  # Rearrange as we would like to have the ADT in the first row and the corresponding RNA below
  pp <- list()
  
  x <- seq(1,length(adt_markers),2)
  if (length(x) > 1){
    x <- x[-length(x)]
  }
  print(x)
  
  j <- 1
  for (i in x){
    pp[[j]] <- p_adt[[i]]
    if ( length(adt_markers) > 1 ){
      pp[[(j+1)]] <- p_adt[[(i+1)]]
    } else {
      pp[[(j+1)]] <- plot_spacer()
    }
    pp[[(j+2)]] <- p_markers[[i]]
    if ( length(markers) > 1 ){
      pp[[(j+3)]] <- p_markers[[(i+1)]]
    } else {
      pp[[(j+3)]] <- plot_spacer()
    }
    j <- j + 4
  }
  
  j <- length(pp) + 1
  if ( length(adt_markers)%%2 == 1 & length(adt_markers) > 1){
    pp[[j]] <- p_adt[[(i+2)]]
    pp[[(j+1)]] <- plot_spacer()
    if ( length(markers) >= length(adt_markers)) {
      pp[[(j+2)]] <- p_markers[[(i+2)]]
      pp[[(j+3)]] <- plot_spacer()
    } else {
      pp[[(j+2)]] <- plot_spacer()
      pp[[(j+3)]] <- plot_spacer()
    }
  }
  
  j <- length(pp)
  if ( length(markers) > length(adt_markers) ){
    for (i in 1:(length(markers) - length(adt_markers))){
      pp[[(j+i)]] <- p_markers[[(length(adt_markers)+i)]]
    }
  }
  
  myTitle <- paste0("violinPlot_ADT_RNA_of_", cellTypes,"_markers")
  plotPanels(pp, 4, myTitle, nrow=2, out = markerPlotDir)
  
}

for (cellTypes in names(markerList)[!(names(markerList) %in% names(adtList))]){
  print(cellTypes)
  
  DefaultAssay(seurat) <- "SCT"
  Idents(seurat) <- "sampleAnnot"
  
  markers <- markerList[[cellTypes]]
  # Check whether all markers can be found
  markers <- markers[markers %in% rownames(seurat@assays$SCT@data)]
  
  p <- list()
  
  for (i in 1:length(markers)){
    p[[(i)]] <- VlnPlot(object = seurat, features = markers[], pt.size = 0.1, y.max = maxY.RNA, same.y.lims = TRUE) + NoLegend() +
      theme(axis.title.x = element_blank())
    
  }
  
  myTitle <- paste0("violinPlot_RNA_of_", cellTypes,"_markers")
  plotPanels(p, 4, myTitle, nrow=2, out = markerPlotDir)
  
}

#### Plots - VDJ ####
# Just show all cells that have data, i.e. do not have <NA> in fx. umis_vdj_t or umis_vdj_b
Idents(seurat) <- "sampleAnnot"

idx <- which(!is.na(seurat@meta.data$umis_vdj_t))
selCells_T <- rownames(seurat@meta.data)[idx]
pdf(paste(outDir, "VDJ_T_cells_dimPlot_afterIntegration.pdf", sep="/"), width = 14, height = 7)
  p <- DimPlot(seurat, cells.highlight = list("T cells" = selCells_T), cols.highlight = "blue", cols = "lightgrey", pt.size = 0.5)
  print(p)
dev.off()

# And as TIFF
ggsave(paste(outDir, "tiff", "VDJ_T_cells_featurePlot_afterIntegration.tiff", sep="/"), p, scale=2, width=5, height=5, dpi=360,compression="lzw")

# And same for B cells
idx <- which(!is.na(seurat@meta.data$umis_vdj_b))
selCells_B <- rownames(seurat@meta.data)[idx]
pdf(paste(outDir,"VDJ_B_cells_dimPlot_afterIntegration.pdf", sep="/"), width = 14, height = 7)
  p <- DimPlot(seurat, cells.highlight = list("B cells"=selCells_B), cols.highlight = "blue", cols = "lightgrey", pt.size = 0.5)
  print(p)
dev.off()

# And as TIFF
ggsave(paste(outDir, "tiff", "VDJ_B_cells_featurePlot_afterIntegration.tiff", sep="/"), p, scale=2, width=5, height=5, dpi=360,compression="lzw")

# And combine?
T_and_B <- intersect(selCells_B, selCells_T)
# Remove the cells present in both T and B cell sets from the original ones
selCells_B <- setdiff(selCells_B, selCells_T)
selCells_B <- selCells_B[!(selCells_B %in% T_and_B)]
selCells_T <- setdiff(selCells_T, selCells_B)
selCells_T <- selCells_T[!(selCells_T %in% T_and_B)]

seurat@meta.data$TB <- ifelse(rownames(seurat@meta.data) %in% selCells_B, "B", 
                              ifelse(rownames(seurat@meta.data) %in% selCells_T, "T",
                                     ifelse(rownames(seurat@meta.data) %in% T_and_B, "T_and_B", "unknown")))

print(table(seurat@meta.data$TB))
#   B       T unknown 
# 790    8511   16184

write.table(as.data.frame(table(seurat@meta.data$TB)), file=paste(outDir, "distribution_T_and_B_cells.txt", sep="/"), sep="\t", quote=FALSE,
            col.names = NA)

# and plot, but see https://github.com/satijalab/seurat/issues/3750
pdf(paste(outDir, "VDJ_All_cells_dimPlot_afterIntegration.pdf", sep="/"), width = 14, height = 7)
# p <- DimPlot(seurat, cells.highlight = list("T cells only"=selCells_T,"B cells only"=selCells_B, "T and B cells"=T_and_B), cols.highlight = c("blue", "red", "purple"), cols = "lightgrey", pt.size = 0.5) 
# there are no TB cells any longer: 'T_and_B' = 'purple'
  p1 <- DimPlot(seurat, group.by = "TB", pt.size = 0.5, cols = c('B' = 'red', 'T' = 'blue', 'unknown' = 'grey'),
                order=c("uknown","B", "T"))
  p1 <- p1 + ggtitle("T and B cells")
  print(p1)
  p2 <- DimPlot(seurat, group.by = "TB", split.by = "TB",
                pt.size = 0.5, cols = c('B' = 'red', 'T' = 'blue', 'unknown' = 'grey'), order=c("uknown","B", "T"))
  p2 <- p2 + ggtitle("T and B cells")
  print(p2)
dev.off()

# And as TIFF
ggsave(paste(outDir, "tiff", "VDJ_All_cells_dimPlot_afterIntegration.tiff", sep="/"), p1, scale=2, width=5, height=5, dpi=360,compression="lzw")
ggsave(paste(outDir, "tiff", "VDJ_All_cells_dimPlot_afterIntegration_splitBy_Type.tiff", sep="/"), p2, scale=2, width=5, height=5, dpi=360,compression="lzw")

# And clean
rm(seurat, seurat.subset)

# MOMA17 is a little the odd one out ...
# The other samples seem to resemble each other quite well!
# I can not see a difference (directly) between MPO, PR3 or SLE and Control


#### Integrate data for all genes ####
# In order to get a count matrix with batch corrected values for all genes and not just the features used to integrate?
# I am always stumbling upon this ... and does it now make sense?
library(Seurat)
library(ggplot2)
library(cowplot)
library(patchwork)
library(grid)

# Integration of datasets per mito content
pcs <- list(Mito10=30)

# @ Memory issues?? see https://github.com/satijalab/seurat/issues/1720
# But for now, just restarting R worked against the 'Calloc' error w.r.t. the integration

for ( i in 1:length(pcs)){
  
  file_names <- dir("./5_removal_tubuluscells/", pattern = "*.rds")
  file_names <- file_names[grepl(names(pcs)[i],file_names)]
  # For the two samples with VDJ info get the object from which the TB cells were removed
  file_names_vdj <- file_names[grepl("noTB_cells",file_names)]
  file_names <- c(file_names_vdj,file_names[grepl("MOMA17|MOMA67|MOMA68|MOMA72|MOMA302",file_names)])
  
  # Get all objects in a list
  scSampleList <- list()
  for (sc in file_names){
    scSampleList[[sc]] <- readRDS(paste0("./5_removal_tubuluscells//",sc))
  }
  
  features <- rownames(scSampleList[[1]])
  scSampleList <- PrepSCTIntegration(object.list = scSampleList, anchor.features = features)
  # subscript out of bounds
  # In addition: Warning message:
  #   The following requested features are not present in any models: C1QTNF12, SCNN1D, TAS1R3, 
  # TTC34, AL365255.1, AL034417.3, CORT, AL139423.1, EXOSC10-AS1, LINC01772, HTR6, AL445471.1, SDC3, 
  # AC114488.3, ZMYM4-AS1, AL591845.1, BMP8B, AL604028.1, TAL1, AL050343.2, GLIS1, CYP2J2, GNG12, 
  # PKN2-AS1, LINC02609, BRDT, NTNG1, KCNA2, AP4B1-AS1, AL445231.1, HIST2H4A, C2CD4D-AS1, MUC1, 
  # AL590666.2, NHLH1, NOS1AP, SELP, LHX4, AL162431.1, AL096803.2, NR5A2, INAVA, ACBD3-AS1, MIXL1, 
  # AL731702.1, AL117350.1, AL591848.3, GCSAML, TPO, SDC1, TDRD15, AC018742.1, AC104699.1, OTOF, 
  # AC013403.2, AC104695.3, AL121652.1, LINC00211, AC074366.1, AC008280.3, LINC01800, MEIS1, 
  # ATOH8, AC015971.1, IGKV3OR2-268, AC104134.1, IGKV5-2, IGKV1-6, IGKV6-21, IGKV2D-30, IGKV6D-21, 
  # IGKV3D-20, IGKV1D-17, IGKV3D-15, IGKV3D-11, IGKV1D-8, CNNM3-DT, LINC01918, AC068491.3, RGPD8, 
  # AC079753.2, AC016745.1, IGKV1OR2-108, TMEM37, LYPD6B, AC012443.2, DAPL1, CERS6-AS1, MYO3B, 
  # ITGA6-AS1, AC010894.4, LINC01116, FRZB, AC009315.1, TFPI, AC013468.1,  [... truncated] 
  
  anca.anchors <- FindIntegrationAnchors(object.list = scSampleList, reference = c(2), normalization.method = "SCT",
                                         anchor.features = features, reduction = "rpca")
  anca.combined.sct <- IntegrateData(anchorset = anca.anchors, normalization.method = "SCT")
  rm(anca.anchors, features, scSampleList)
  anca.combined.sct <- RunPCA(anca.combined.sct, verbose = FALSE, npcs = pcs[i][[1]])
  anca.combined.sct <- RunUMAP(anca.combined.sct, reduction = "pca", dims = 1:pcs[i][[1]])
  
  # And save the integrated object
  saveRDS(anca.combined.sct, paste0(outDir,"/", names(pcs)[i],".integrated_allFeatures.rds"))
  
  # And clean
  rm(anca.combined.sct)
} 
