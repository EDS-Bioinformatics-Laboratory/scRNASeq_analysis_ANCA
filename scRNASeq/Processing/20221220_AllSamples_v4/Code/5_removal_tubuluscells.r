# Cell annotation with SingleR and CellTypist

#### Set default for a good start ####
source("initAnalysis.r")
init()

# Create a directory to hold results and plots
outDir <- "5_removal_tubuluscells"
dir.create(outDir)

#### Library ####
library(Seurat)
library(ggplot2)
library(patchwork)
library(grid)

dataDir <- "4_cell_annotation_before_integration"

# Sample list
sampleList <- c("MOMA17", "MOMA52", "MOMA57", "MOMA67", "MOMA68", "MOMA72", "MOMA302")

# List with number of PCs (only for Mito10)
pcs <- list(Mito10=30)

# First, we have 'trained' SingleR with the data from Stewart (see Yosta_Vegting\scRNASeq\Processing\20210630_Analysis_MitoContent\Code)
# to identify the Tubulus cells

# Now, we remove these cells and do some final QC ...

for ( i in 1:length(sampleList) ){
  print(sampleList[i])

  for ( mito in 1:length(pcs) ){
    print(names(pcs)[mito])
     
    if (sampleList[i] %in% c("MOMA52","MOMA57")){
      seurat <- readRDS(paste0(dataDir,"/",sampleList[i],".",names(pcs)[mito],".noTB_cells.doubletsRemoved.annot.rds"))
      
      # Make sure the default assay is the normalized one, i.e. 'SCT'
      DefaultAssay(seurat) <- "SCT"
      
      cat(paste0(sampleList[i],".",names(pcs)[mito],".noTB_cells.doubletsRemoved.annot.rds")," contains ",dim(seurat)[2]," cells before filtering for Proximal tubulus cells\n")
      
      seurat <- subset(seurat, subset = matureKidney.main.before  != "Proximal tubule")
      cat("... and ",dim(seurat)[2]," cells after filtering\n")
      
      # Save the filtered seurat object
      saveRDS(seurat, paste0(outDir,"/",sampleList[i],".",names(pcs)[mito],".noTB_cells.doubletsRemoved.annot.rds"))
      
      # # Let us do some initial clustering to see whether we have cleaned data, i.e. no clusters with suspect low numbers of housekeeping genes,
      # no clusters with low number of expressed genes etc.
      
      #### Clustering ####
      # Find neighbors
      seurat <- FindNeighbors(seurat, reduction = "pca", dims = 1:pcs[[mito]])
      
      # Cluster at different resolution
      seurat <- FindClusters(seurat, resolution = 0.3)
      seurat <- FindClusters(seurat, resolution = 0.6)
      
      #### Plot clusters ####
      clustDir <- paste0(outDir,"/", paste0(sampleList[i],"_", names(pcs)[mito], "_noTB_cells_clusters"))
      dir.create(clustDir)
      
      for (clust in colnames(seurat@meta.data)[grep("SCT_snn", colnames(seurat@meta.data))]){
        res <- gsub("SCT_snn_","", clust)
        
        pdf(paste(clustDir, paste0("dimPlot_Clustered_", res,".pdf"), sep="/"))
        Idents(seurat) <- clust
        p <- DimPlot(seurat, reduction = "umap", label = TRUE,  repel = TRUE) +
          guides(colour=guide_legend(override.aes = list(size=2), ncol=1))
        print(p)
        dev.off()
      }
      
      # Violin plot per clustering resolution
      # See: https://patchwork.data-imaginist.com/articles/guides/layout.html, https://patchwork.data-imaginist.com/articles/guides/assembly.html
      # and https://www.rdocumentation.org/packages/grid/versions/3.6.2/topics/grid.text
      for (clust in colnames(seurat@meta.data)[grep("SCT_snn", colnames(seurat@meta.data))]){
        res <- gsub("SCT_snn_","", clust)
        seurat@meta.data[,clust] <- factor(seurat@meta.data[,clust], levels=as.character(0:(length(unique(seurat@meta.data[,clust]))-1)))
        
        Idents(seurat) <- clust
        
        # Plot some QC; fraction expression mitochondrial genes, housekeeping genes, number of features per cluster
        pdf(paste(clustDir,paste0("violinPlot_qcStats_", res,".pdf"), sep="/"))
        p1 <- VlnPlot(seurat,features="percent.mito", pt.size = 0.1) + NoLegend() +
          ggtitle("percent.mito") +
          theme(axis.title.x = element_blank())
        p2 <- VlnPlot(seurat,features="n.exp.hkgenes", pt.size = 0.1) + NoLegend() +
          theme(axis.title.x = element_blank())
        p3 <- VlnPlot(seurat,features="nFeature_RNA", pt.size = 0.1) + NoLegend() +
          theme(axis.title.x = element_blank())
        p <- p1 / p2 / p3 + plot_layout(heights = c(1, 1, 1))
        print(p)
        dev.off()
        
      }
    }
    
    # Just the normal samples
    seurat <- readRDS(paste0(dataDir,"/",sampleList[i],".",names(pcs)[mito],".doubletsRemoved.annot.rds"))
    
    # Make sure the default assay is the normalized one, i.e. 'SCT'
    DefaultAssay(seurat) <- "SCT"
    
    cat(paste0(sampleList[i],".",names(pcs)[mito],".doubletsRemoved.annot.rds")," contains ",dim(seurat)[2]," cells before filtering for Proximal tubulus cells\n")

    seurat <- subset(seurat, subset = matureKidney.main.before  != "Proximal tubule")
    cat("... and ",dim(seurat)[2]," cells after filtering\n")

    # Save the filtered seurat object
    saveRDS(seurat, paste0(outDir,"/",sampleList[i],".",names(pcs)[mito],".doubletsRemoved.annot.rds"))
    
    # # Let us do some initial clustering to see whether we have cleaned data, i.e. no clusters with suspect low numbers of housekeeping genes,
    # no clusters with low number of expressed genes etc.
    
    #### Clustering ####
    # Find neighbors
    seurat <- FindNeighbors(seurat, reduction = "pca", dims = 1:pcs[[mito]])
    
    # Cluster at different resolution
    seurat <- FindClusters(seurat, resolution = 0.3)
    seurat <- FindClusters(seurat, resolution = 0.6)
    
    #### Plot clusters ####
    clustDir <- paste0(outDir,"/", paste0(sampleList[i],"_", names(pcs)[mito], "_clusters"))
    dir.create(clustDir)

    for (clust in colnames(seurat@meta.data)[grep("SCT_snn", colnames(seurat@meta.data))]){
      res <- gsub("SCT_snn_","", clust)
      
      pdf(paste(clustDir, paste0("dimPlot_Clustered_", res,".pdf"), sep="/"))
        Idents(seurat) <- clust
        p <- DimPlot(seurat, reduction = "umap", label = TRUE,  repel = TRUE) +
          guides(colour=guide_legend(override.aes = list(size=2), ncol=1))
        print(p)
      dev.off()
    }
    
    # Violin plot per clustering resolution
    # See: https://patchwork.data-imaginist.com/articles/guides/layout.html, https://patchwork.data-imaginist.com/articles/guides/assembly.html
    # and https://www.rdocumentation.org/packages/grid/versions/3.6.2/topics/grid.text
    for (clust in colnames(seurat@meta.data)[grep("SCT_snn", colnames(seurat@meta.data))]){
      res <- gsub("SCT_snn_","", clust)
      seurat@meta.data[,clust] <- factor(seurat@meta.data[,clust], levels=as.character(0:(length(unique(seurat@meta.data[,clust]))-1)))
      
      Idents(seurat) <- clust
      
      # Plot some QC; fraction expression mitochondrial genes, housekeeping genes, number of features per cluster
      pdf(paste(clustDir,paste0("violinPlot_qcStats_", res,".pdf"), sep="/"))
      p1 <- VlnPlot(seurat,features="percent.mito", pt.size = 0.1) + NoLegend() +
        ggtitle("percent.mito") +
        theme(axis.title.x = element_blank())
      p2 <- VlnPlot(seurat,features="n.exp.hkgenes", pt.size = 0.1) + NoLegend() +
        theme(axis.title.x = element_blank())
      p3 <- VlnPlot(seurat,features="nFeature_RNA", pt.size = 0.1) + NoLegend() +
        theme(axis.title.x = element_blank())
      p <- p1 / p2 / p3 + plot_layout(heights = c(1, 1, 1))
      print(p)
      dev.off()
      
    }
  }
}

# [1] "MOMA17"
# [1] "Mito10"
# MOMA17.Mito10.doubletsRemoved.annot.rds  contains  1225  cells before filtering for Proximal tubulus cells
# ... and  466  cells after filtering
# 
# [1] "MOMA52"
# [1] "Mito10"
# MOMA52.Mito10.noTB_cells.doubletsRemoved.annot.rds  contains  5656  cells before filtering for Proximal tubulus cells
# ... and  5553  cells after filtering
# 
# MOMA52.Mito10.doubletsRemoved.annot.rds  contains  5709  cells before filtering for Proximal tubulus cells
# ... and  5600  cells after filtering
# 
# [1] "MOMA57"
# [1] "Mito10"
# MOMA57.Mito10.noTB_cells.doubletsRemoved.annot.rds  contains  8517  cells before filtering for Proximal tubulus cells
# ... and  8359  cells after filtering
# 
# MOMA57.Mito10.doubletsRemoved.annot.rds  contains  8565  cells before filtering for Proximal tubulus cells
# ... and  8406  cells after filtering
# 
# [1] "MOMA67"
# [1] "Mito10"
# MOMA67.Mito10.doubletsRemoved.annot.rds  contains  4659  cells before filtering for Proximal tubulus cells
# ... and  4583  cells after filtering
# 
# [1] "MOMA68"
# [1] "Mito10"
# MOMA68.Mito10.doubletsRemoved.annot.rds  contains  2568  cells before filtering for Proximal tubulus cells
# ... and  2546  cells after filtering
# 
# [1] "MOMA72"
# [1] "Mito10"
# MOMA72.Mito10.doubletsRemoved.annot.rds  contains  1239  cells before filtering for Proximal tubulus cells
# ... and  1231  cells after filtering
# 
# [1] "MOMA302"
# [1] "Mito10"
# MOMA302.Mito10.doubletsRemoved.annot.rds  contains  2785  cells before filtering for Proximal tubulus cells
# ... and  2747  cells after filtering
