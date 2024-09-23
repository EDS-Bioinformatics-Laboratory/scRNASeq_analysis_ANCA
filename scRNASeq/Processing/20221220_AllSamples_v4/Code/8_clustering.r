# Clustering
#

#### Set default for a good start ####
source("initAnalysis.r")
init()

# Create a directory to hold results and plots
outDir <- "8_clustering"
dir.create(outDir)

#### Library ####
library(Seurat)
library(ggplot2)
library(patchwork)
library(grid)

dataDir <- "7_cellTypist_annotation_after_integration/"

# List with number of PCs
pcs <- list(Mito10=30)

for ( i in 1:length(pcs)){
  clustDir <- paste0(outDir,"/", paste0(names(pcs)[i], "_clusters"))
  dir.create(clustDir)
  
  # Read the data
  seurat <- readRDS(paste0(dataDir,"/",names(pcs)[i],".integrated.rds"))
  
  # Make sure to use the 'integrated' assay object to do the clustering upon
  # https://github.com/satijalab/seurat/issues/1717
  DefaultAssay(seurat) <- 'integrated'

  #### Clustering ####
  # Find neighbors
  seurat <- FindNeighbors(seurat, reduction = "pca", dims = 1:pcs[i][[1]])

  # Cluster at different resolution
  seurat <- FindClusters(seurat, resolution = 0.3)
  seurat <- FindClusters(seurat, resolution = 0.6)
  seurat <- FindClusters(seurat, resolution = 1)

  #And save
  saveRDS(seurat, paste0(outDir,"/",names(pcs)[i],".integrated.clustered.rds"))
  
  #### Plot clusters ####
  for (clust in colnames(seurat@meta.data)[grep("integrated_snn", colnames(seurat@meta.data))]){
    res <- gsub("integrated_snn_","", clust)
    
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
  for (clust in colnames(seurat@meta.data)[grep("integrated_snn", colnames(seurat@meta.data))]){
    res <- gsub("integrated_snn_","", clust)
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
    
    # And get some annotation in (see 7_celltypist_annotation_after_integration.r), also
    # for resolution 0.6
    # Make plots for bestMatchNoMajorityVoting and probMatchMajorityVoting
    for (scenario in c("bestMatchNoMajorityVoting", "probMatchMajorityVoting")){
      pdf(paste(clustDir, paste0("umap_celltypist_", scenario,"_annot.", names(pcs)[i],".",res,".pdf"), sep="/"), width = 14, height = 7)
        p1 <- DimPlot(object = seurat, reduction = "umap", group.by = paste0(scenario,".labels.celltypist.after"), label = TRUE) + NoLegend()
        p2 <- DimPlot(object = seurat, reduction = "umap", group.by = clust, label = TRUE) + NoLegend()
        print((plot_grid(p1, p2)))
        # Get the legends on a separate page using the cowplot package
        p1 <- DimPlot(object = seurat, reduction = "umap", group.by = paste0(scenario,".labels.celltypist.after"), label = TRUE) +
          guides(color = guide_legend(override.aes = list(size=3), nrow = 30) ) +
          theme(legend.text = element_text(size=9))
        legend_p1 <- cowplot::get_legend(p1)
        p2 <- DimPlot(object = seurat, reduction = "umap", group.by = clust, label = TRUE)
        legend_p2 <- cowplot::get_legend(p2)
        print(plot_grid(NULL, legend_p1, NULL, NULL, legend_p2,ncol = 5))
      dev.off()
      
      pdf(paste(clustDir, paste0(scenario,"_celltypist_annotation.", names(pcs)[i],".",res,".pdf"), sep="/"), width = 14, height =14)
        tab <- table(seurat@meta.data[,clust], seurat@meta.data[,paste0(scenario,".labels.celltypist.after")])
        corrplot(tab/rowSums(tab), is.corr = FALSE, col = 
                 colorRampPalette(c("#FFFFFF", "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC", "#053061"))(200))
      dev.off()
      
      # And get some tables 
      # - of the distribution of the different celltypes over the different clusters
      
      tab <- table(seurat@meta.data[,clust], seurat@meta.data[,paste0(scenario,".labels.celltypist.after")])
      write.table(tab, paste(clustDir,paste0(scenario,"_celltypes_per_cluster.", names(pcs)[i],".", res,".txt"), sep="/"), quote=FALSE, sep="\t", row.names = TRUE, col.names = TRUE)
    }
  }
}
# And clean
rm(seurat, p, p1, p2, p3, res)

