# Cell after integration with CellTypist

#### Set default for a good start ####
source("initAnalysis.r")
init()

# Create a directory to hold results and plots
outDir <- "7_cellTypist_annotation_after_integration"
dir.create(outDir)

#### Library ####
library(Seurat)
library(ggplot2)
library(cowplot)
library(corrplot)

#### CellTypist - after integration ####
dataDir <- "6_integration"

# List with number of PCs (only for Mito10)
pcs <- list(Mito10=30)

## Read the data ##
for (j in 1:length(pcs)){
  seurat <- readRDS(paste0(dataDir,"/",names(pcs)[j],".integrated.rds"))
  DefaultAssay(seurat) <- "integrated"

  #### Get the raw counts ####
  write.table(as.data.frame(seurat@assays$RNA@counts), paste0(outDir,"/",names(pcs)[j], "_countsAfterIntegration.csv"),
                sep = ',', row.names = T, col.names = T, quote = F)
}

#### CellTypist - Run ####
# see ../NoteBooks/runCellTypist_afterIntegration.ipynb
#
# I tried different scenario's
# 1. Just starting from the raw counts and without majority voting, but with 'best match'
#    - This should be the same before and after integration
# 2. Apply majority voting
#    - This has a dependency on neighbouring cells, so will be different before and after integration

#### CellTypist - Read predicted labels ####
# see ../NoteBooks/runCellTypist_afterIntegration.ipynb
# 
# see https://www.celltypist.org/tutorials/onlineguide
#
# I ran both scenario's on every subject
# - saved in their own subdirectory with the sample directory
#   i.e. integratedSamples/bestMatch_noMajorityvoting
#        integratedSamples/probMatch_majorityVoting

for ( i in 1:length(pcs) ){
  print(names(pcs[i]))
  sampleDir <- paste(outDir,"integratedSamples", sep="/")
  
  # Read data
  if ( file.exists(paste0(outDir,"/",names(pcs)[j],".integrated.rds"))) {
    seurat <- readRDS(paste0(outDir,"/",names(pcs)[j],".integrated.rds"))
  } else {
    seurat <- readRDS(paste0(dataDir,"/",names(pcs)[j],".integrated.rds"))
  }
    
  if (length(grep("celltypist.after", colnames(seurat@meta.data))) > 0){
    # Make a fresh start, remove all columns containing 'celltypist'in their name ...
    print("Removing all columns holding celltypist information from after integration ...")
    seurat@meta.data <- seurat@meta.data[,-c(grep("celltypist.after", colnames(seurat@meta.data)))]
  }
  
  # The 'probability_matrix.csv' contains the probabilities per cell of to which celltype it would fit.
  # NOTE: the 'probability_matrix' is of course the same for both the 'best match' and 'majority voting' scenario
  annot.bestMatch.df <- read.csv(paste0(sampleDir,"/bestMatch_noMajorityVoting/predicted_labels.csv"))
  annot.majorityVoting.df <- read.csv(paste0(sampleDir,"/probMatch_majorityVoting/predicted_labels.csv"))
  probMatrix.df <- read.csv(paste0(sampleDir,"/bestMatch_noMajorityVoting/probability_matrix.csv"), check.names = FALSE)
  
  # check whether cells match
  all(colnames(seurat) == annot.bestMatch.df$X)
  all(colnames(seurat) == annot.majorityVoting.df$X)
  
  # # Make sure to keep the rownames !!
  seurat@meta.data$barcodes.orig <- rownames(seurat@meta.data)
  
  # # now, merge by rowname and celltypistName
  colnames(annot.bestMatch.df) <- c("cellName","bestMatchNoMajorityVoting.labels.celltypist.after")
  if (ncol(annot.majorityVoting.df) > 2){
    colnames(annot.majorityVoting.df) <- c("cellName","predicted.labels.celltypist.after", "over_clustering.celltypist.after", "probMatchMajorityVoting.labels.celltypist.after")
  } else {
    colnames(annot.majorityVoting.df) <- c("cellName","probMatchMajorityVoting.labels.celltypist.after")
  } 
  seurat@meta.data <- merge(seurat@meta.data, annot.bestMatch.df, by.x="row.names",by.y="cellName", all.x = TRUE, sort = FALSE )
  seurat@meta.data <- merge(seurat@meta.data, annot.majorityVoting.df, by.x="Row.names",by.y="cellName", all.x = TRUE, sort = FALSE )
  # 
  
  # And add the probability scores ... becoming a very large object :-)
  colnames(probMatrix.df) <- paste0(colnames(probMatrix.df),".celltypist.after")
  celltypist.max_probability <- 0
  celltypist.max_probability <- apply(probMatrix.df[,-1], 1, function(x) max(x))
  
  seurat@meta.data <- merge(seurat@meta.data, probMatrix.df, by.x="Row.names",by.y=".celltypist.after", all.x = TRUE, sort = FALSE )
  
  # Put the rownames back
  rownames(seurat@meta.data) <- seurat@meta.data$barcodes.orig
  seurat@meta.data$Row.names <- NULL
  seurat@meta.data$barcodes.orig <- NULL
  # And add the column with the maximum probability
  seurat@meta.data$celltypist.max_probability <- celltypist.max_probability
  
  # And save
  saveRDS(seurat, paste0(outDir,"/",names(pcs)[i],".integrated.rds"))
}

#### CellTypist - Plots and tables ####
plotsDir <- paste(outDir, "plotsAndTables", sep="/" )
dir.create(plotsDir)

for ( i in 1:length(pcs) ){
  print(names(pcs[i]))
  
  # Read data
  if ( file.exists(paste0(outDir,"/",names(pcs)[i],".integrated.rds"))) {
    seurat <- readRDS(paste0(outDir,"/",names(pcs)[i],".integrated.rds"))
  } else {
    seurat <- readRDS(paste0(dataDir,"/",names(pcs)[i],".integrated.rds"))
  }
  
  # Make plots for bestMatchNoMajorityVoting and probMatchMajorityVoting
  for (scenario in c("bestMatchNoMajorityVoting", "probMatchMajorityVoting")){
    pdf(paste(plotsDir, paste0("umap_celltypist_", scenario,"_annot.", names(pcs)[i],".pdf"), sep="/"), width = 14, height = 7)
    p1 <- DimPlot(object = seurat, reduction = "umap", group.by = paste0(scenario,".labels.celltypist.after"), label = TRUE) + NoLegend()
    p2 <- DimPlot(object = seurat, reduction = "umap", group.by = "SCT_snn_res.0.3", label = TRUE) + NoLegend()
    print((plot_grid(p1, p2)))
    # Get the legends on a separate page using the cowplot package
    p1 <- DimPlot(object = seurat, reduction = "umap", group.by = paste0(scenario,".labels.celltypist.after"), label = TRUE) +
      guides(color = guide_legend(override.aes = list(size=3), nrow = 30) ) +
      theme(legend.text = element_text(size=9))
    legend_p1 <- cowplot::get_legend(p1)
    p2 <- DimPlot(object = seurat, reduction = "umap", group.by = "SCT_snn_res.0.3", label = TRUE)
    legend_p2 <- cowplot::get_legend(p2)
    print(plot_grid(NULL, legend_p1, NULL, NULL, legend_p2,ncol = 5))
    dev.off()
    
    pdf(paste(plotsDir, paste0("umap_celltypist_maxProbability.", names(pcs)[i],".pdf"), sep="/"), width = 14, height = 7)
    p1 <- FeaturePlot(object = seurat, features = "celltypist.max_probability", reduction = "umap")
    print(plot_grid(p1))
    dev.off()
    
    pdf(paste(plotsDir, paste0("violin_celltypist_probabilities.", names(pcs)[i],".pdf"), sep="/"), width = 14, height = 14)
    p1 <- VlnPlot(object = seurat, features = "celltypist.max_probability", group.by = paste0(scenario,".labels.celltypist.after")) + NoLegend()
    print(plot_grid(p1))
    dev.off()
    
    pdf(paste(plotsDir, paste0(scenario,"_celltypist_annotation.", names(pcs)[i],".pdf"), sep="/"), width = 14, height =14)
    tab <- table(seurat$SCT_snn_res.0.3, seurat@meta.data[,paste0(scenario,".labels.celltypist.after")])
    corrplot(tab/rowSums(tab), is.corr = FALSE, col = 
               colorRampPalette(c("#FFFFFF", "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC", "#053061"))(200))
    dev.off()
    
    # And get some tables 
    # - of the distribution of the different celltypes over the different clusters
    
    tab <- table(seurat$SCT_snn_res.0.3, seurat@meta.data[,paste0(scenario,".labels.celltypist.after")])
    write.table(tab, paste(plotsDir,paste0(scenario,"_celltypes_per_cluster_res03.", names(pcs)[i],".txt"), sep="/"), quote=FALSE, sep="\t", row.names = TRUE, col.names = TRUE)
    
    df <- as.data.frame(tab)
    colnames(df) <- c("cluster", "celltype", "frequency")
    
    df$celltype <- factor(df$celltype)
    nrColors <- length(levels(factor(df$celltype)))
    rainbowColors = rainbow(nrColors)
    df$colors <- rainbowColors[factor(df$celltype)]
    
    pdf(paste(plotsDir, paste0("distribution_per_cluster_stack.", scenario,".", names(pcs)[i],".pdf"), sep="/"), width = 14, height =14)
    p <- ggplot(df, aes(x=cluster, y=frequency, fill=celltype)) +
      geom_bar(stat='identity', position="stack")+
      scale_fill_manual(values = rainbowColors) +
      guides(fill=guide_legend(ncol =1))
    print(p)
    legend_p <- cowplot::get_legend(p)
    p <- ggplot(df, aes(x=cluster, y=frequency, fill=celltype)) +
      geom_bar(stat='identity', position="stack")+
      scale_fill_manual(values = rainbowColors) + 
      NoLegend()
    print(p)
    print(plot_grid(legend_p))
    dev.off()
  }
}
