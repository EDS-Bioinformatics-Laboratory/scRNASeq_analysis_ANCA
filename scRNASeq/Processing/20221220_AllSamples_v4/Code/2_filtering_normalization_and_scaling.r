#### Filter ####
#
# Previously for Patient1 data
# 20210630 - Based on paper by Schroeder et al (https://www.biorxiv.org/content/10.1101/2020.03.02.973925v1.full) allow more 
#            mitochondrial reads/counts
# 20210716 - Decide to go for 20, 50 and 80% on filtered CellRanger data
# 20221117 - We only used 10% cut-off in the 20220916_AllSamples_v2 analysis
#          - We include filtering on housekeeping genes (> 55) to get rid of the
#            cluster with cells that have a low number of housekeeping cells expressed
#            (see cluster 4 in 20220916_AllSamples_v2/Code/6_clustering/Mito10_clusters/violinPlot_qcStats_res.0.3.pdf)

#### Set default for a good start ####
source("initAnalysis.r")
init()

# Create a directory to hold results and plots
outDir <- "2_filtering_normalization_and_scaling"
dir.create(outDir)

#### Library ####
library(Seurat)
library(viridis)
library(ggplot2)
library(patchwork)

#### Read the data ####
sampleList <- c("MOMA17", "MOMA52", "MOMA57", "MOMA67", "MOMA68", "MOMA72", "MOMA302")

for ( sc in sampleList){
  seurat <- readRDS(paste0("../Data/scRNASeq/Processed/", sc, ".raw.annot.rds"))
  
  # Size before filtering
  origDim <- dim(seurat)
  
  # Perform some filtering based on nrGenes, nrExpHouseholdGenes, percMito
  print(paste0(sc, " - range of nFeatures_RNA: ", paste0(range(seurat@meta.data$nFeature_RNA), collapse=" - ")))
  # 
  
  #### Filtering ####
  # The Workshop of the Broad Institute filtered also on a maximum of the nr of genes...
  #seurat <- subset(seurat, subset = nFeature_RNA > 350 & nFeature_RNA < 5000 & percent.mito < 0.1 & n.expressed.hkgenes > 55 )
  # see 'filtering_scenarios.r' - CreateSeuratObject filers first on min.features and then on the resulting object on min.cells
  # (via CreateAssayObject)
  # Schroeder et al:
  # "we expanded the final threshold to keep cells with <80% mitochondrial content and >200 genes per cell in an effort to 
  #  maximize our cellular population. Additionally, cells with >6,000 genes per cell were excluded to eliminate the likelihood
  #  of bias introduced by the presence of doublets."
  for ( pctMito in c(10)){
    seuratMito <- subset(seurat, subset = nFeature_RNA >= 200 & nFeature_RNA < 6000 & percent.mito < pctMito & n.exp.hkgenes > 55 )
    
    print(seuratMito)  #AJ - 20221130; Check to see whether ADT slot is still there...
    # YEP, still there!
    
    # A gene should be present in at least 3 (i.e. more than 2) cells
    cells_per_gene <- Matrix::rowSums(seuratMito@assays$RNA@counts > 0)
    cells_per_gene.gt.2 <- names(which(cells_per_gene > 2))
    
    #seuratMito <- subset(seuratMito, features = c(cells_per_gene.gt.2))
    #print(seuratMito)  #AJ - 20221130; Check to see whether ADT slot is still there...
    # And here we loose the ADT slot .... we have to use:
    seuratMito@assays$RNA <- subset(seuratMito@assays$RNA, features = c(cells_per_gene.gt.2))
    
    # Size after filtering
    filteredDim <- dim(seuratMito)
    
    cat(paste0("... pctMito <", pctMito,"%: From the original ", origDim[1], " features, ", origDim[1]-filteredDim[1], " were filtered, leaving ", filteredDim[1], " features\n"))
    cat(paste0("... pctMito <", pctMito,"%: From the original ", origDim[2], " cells, ", origDim[2]-filteredDim[2], " were filtered, leaving ", filteredDim[2], " cells\n"))
    
    # OBSERVATION
    # 
    
    # And save
    saveRDS(seuratMito, paste0(outDir,"/", sc, ".filtMito",pctMito,".rds"))
    
    # And remove
    rm(seuratMito)
  }
}

## NOTE
# It seems that here we loose the ADT slot!!
# In the subsetting you really have to specify the Assay to work on!

# [1] "MOMA17 - range of nFeatures_RNA: 24 - 8336"
# An object of class Seurat 
# 36738 features across 1306 samples within 2 assays 
# Active assay: RNA (36601 features, 0 variable features)
# 1 other assay present: ADT
# ... pctMito <10%: From the original 36601 features, 18133 were filtered, leaving 18468 features
# ... pctMito <10%: From the original 7242 cells, 5936 were filtered, leaving 1306 cells
# [1] "MOMA52 - range of nFeatures_RNA: 28 - 6245"
# An object of class Seurat 
# 36738 features across 6110 samples within 2 assays 
# Active assay: RNA (36601 features, 0 variable features)
# 1 other assay present: ADT
# ... pctMito <10%: From the original 36601 features, 19228 were filtered, leaving 17373 features
# ... pctMito <10%: From the original 6695 cells, 585 were filtered, leaving 6110 cells
# [1] "MOMA57 - range of nFeatures_RNA: 23 - 5398"
# An object of class Seurat 
# 36738 features across 9174 samples within 2 assays 
# Active assay: RNA (36601 features, 0 variable features)
# 1 other assay present: ADT
# ... pctMito <10%: From the original 36601 features, 18391 were filtered, leaving 18210 features
# ... pctMito <10%: From the original 9605 cells, 431 were filtered, leaving 9174 cells
# [1] "MOMA67 - range of nFeatures_RNA: 28 - 8144"
# An object of class Seurat 
# 36738 features across 4969 samples within 2 assays 
# Active assay: RNA (36601 features, 0 variable features)
# 1 other assay present: ADT
# ... pctMito <10%: From the original 36601 features, 17259 were filtered, leaving 19342 features
# ... pctMito <10%: From the original 5752 cells, 783 were filtered, leaving 4969 cells
# [1] "MOMA68 - range of nFeatures_RNA: 47 - 7103"
# An object of class Seurat 
# 36738 features across 2740 samples within 2 assays 
# Active assay: RNA (36601 features, 0 variable features)
# 1 other assay present: ADT
# ... pctMito <10%: From the original 36601 features, 19486 were filtered, leaving 17115 features
# ... pctMito <10%: From the original 3128 cells, 388 were filtered, leaving 2740 cells
# [1] "MOMA72 - range of nFeatures_RNA: 36 - 7904"
# An object of class Seurat 
# 36738 features across 1314 samples within 2 assays 
# Active assay: RNA (36601 features, 0 variable features)
# 1 other assay present: ADT
# ... pctMito <10%: From the original 36601 features, 20185 were filtered, leaving 16416 features
# ... pctMito <10%: From the original 1663 cells, 349 were filtered, leaving 1314 cells
# [1] "MOMA302 - range of nFeatures_RNA: 46 - 7836"
# An object of class Seurat 
# 36738 features across 2961 samples within 2 assays 
# Active assay: RNA (36601 features, 0 variable features)
# 1 other assay present: ADT
# ... pctMito <10%: From the original 36601 features, 18990 were filtered, leaving 17611 features
# ... pctMito <10%: From the original 3549 cells, 588 were filtered, leaving 2961 cells

# What happened with the TB_cells?
# MOMA52
seurat <- readRDS(paste0(outDir,"/", "MOMA52.filtMito10.rds"))
idx <- which(!is.na(seurat@meta.data$umis_vdj_t))
selCells_T <- rownames(seurat@meta.data)[idx]
idx <- which(!is.na(seurat@meta.data$umis_vdj_b))
selCells_B <- rownames(seurat@meta.data)[idx]
T_and_B <- intersect(selCells_B, selCells_T)
# Remove the cells present in both T and B cell sets from the original ones
selCells_B <- setdiff(selCells_B, selCells_T)
selCells_B <- selCells_B[!(selCells_B %in% T_and_B)]
selCells_T <- setdiff(selCells_T, selCells_B)
selCells_T <- selCells_T[!(selCells_T %in% T_and_B)]
seurat@meta.data$TB <- ifelse(rownames(seurat@meta.data) %in% selCells_B, "B", 
                              ifelse(rownames(seurat@meta.data) %in% selCells_T, "T",
                                     ifelse(rownames(seurat@meta.data) %in% T_and_B, "T_and_B", "unknown")))
table(seurat@meta.data$TB)
#   B       T T_and_B unknown 
# 422    3926      52    1710  -> 6110

# Still 52 ...

# MOMA57
seurat <- readRDS(paste0(outDir,"/", "MOMA57.filtMito10.rds"))
idx <- which(!is.na(seurat@meta.data$umis_vdj_t))
selCells_T <- rownames(seurat@meta.data)[idx]
idx <- which(!is.na(seurat@meta.data$umis_vdj_b))
selCells_B <- rownames(seurat@meta.data)[idx]
T_and_B <- intersect(selCells_B, selCells_T)
# Remove the cells present in both T and B cell sets from the original ones
selCells_B <- setdiff(selCells_B, selCells_T)
selCells_B <- selCells_B[!(selCells_B %in% T_and_B)]
selCells_T <- setdiff(selCells_T, selCells_B)
selCells_T <- selCells_T[!(selCells_T %in% T_and_B)]
seurat@meta.data$TB <- ifelse(rownames(seurat@meta.data) %in% selCells_B, "B", 
                              ifelse(rownames(seurat@meta.data) %in% selCells_T, "T",
                                     ifelse(rownames(seurat@meta.data) %in% T_and_B, "T_and_B", "unknown")))
table(seurat@meta.data$TB)
#   B       T T_and_B unknown 
# 417    5019      52    3686  -> 9174

# Also , still 52

#### Plots - after filtering ####
file_names <- dir(outDir, pattern = "*.rds")
file_names <- file_names[grepl("Mito",file_names)]

for ( rds in file_names ){
  # Read data
  seurat <- readRDS(paste0(outDir,"/",rds))
  
  # File name
  newName <- gsub(".rds","", rds)
  
  # Create useful vectors
  counts_per_cell <- Matrix::colSums(seurat@assays$RNA@counts)
  counts_per_gene <- Matrix::rowSums(seurat@assays$RNA@counts)
  genes_per_cell <- Matrix::colSums(seurat@assays$RNA@counts > 0) # count gene only if it has non-zero reads mapped.
  cells_per_gene <- Matrix::rowSums(seurat@assays$RNA@counts > 0) # count number of cells a gene is present in, ie has at least 1 count
  
  pdf(paste(outDir,paste0(newName, "_statsAfterFiltering_Luecken_Figure2.pdf"), sep="/"))
  # A
  counts_per_cell.df <- as.data.frame(counts_per_cell)
  
  hist(counts_per_cell, main="", breaks=100, col='cornflowerblue', xlab = "Count depth")
  hist(counts_per_cell, main="", xlim=c(0,5000), breaks=1000, col='cornflowerblue', xlab = "Count depth")
  
  histLarge.df <- hist(counts_per_cell, main=sc, breaks=100, col='cornflowerblue', xlab = "Count depth", plot = FALSE)
  binWidthLarge <- histLarge.df$breaks[2]-histLarge.df$breaks[1]
  breaksLarge <- histLarge.df$breaks
  maxLarge <- max(histLarge.df$counts)
  
  histDetail.df <- hist(counts_per_cell, main=sc, xlim=c(0,5000), breaks=1000, col='cornflowerblue', xlab = "Count depth", plot=FALSE)
  binWidthDetail <- histDetail.df$breaks[2]-histDetail.df$breaks[1]
  breaksDetail <- histDetail.df$breaks
  maxDetail <- max(histDetail.df$counts)
  
  p.overview <- ggplot(counts_per_cell.df,aes(x=counts_per_cell)) +
    geom_histogram(binwidth = binWidthLarge, color='black', fill='cornflowerblue') +
    #    scale_x_continuous(breaks = breaksLarge) +
    labs(x="Count depth", y="Frequency") +
    theme_classic()
  p.detail <- ggplot(counts_per_cell.df,aes(x=counts_per_cell)) +
    geom_histogram(binwidth = binWidthDetail, color='black', fill='cornflowerblue') +
    xlim(0,5000) +
    #    scale_x_continuous(breaks = breaksDetail) +
    labs(x="Count depth", y="Frequency") +
    theme_classic()
  p.all <- p.overview + annotation_custom(ggplotGrob(p.detail),
                                          xmin = max(counts_per_cell.df)/2, ymin = maxLarge/2,
                                          xmax = max(counts_per_cell.df)) + ggtitle(sc)
  
  print(p.all)
  rm(counts_per_cell.df)
  
  # B
  hist(genes_per_cell, main="", breaks=100, col='cornflowerblue', xlab="Number of genes", xaxt='n')
  axis(side=1, at=c(seq(0,999,250),seq(1000,6000,500)), labels=c(seq(0,999,250),seq(1000,6000,500)), las=2)
  # Q: Where to filter? at 500, 750, 800
  abline(v=c(500,750,800), col="red", lwd=2)
  
  # C
  plot(sort(counts_per_cell, decreasing = TRUE), log='y', main='', xlab="Cell rank", ylab="Count depth",
       yaxt='n')
  axis(2, at=c(seq(0,900,100),seq(1000,9000, 1000),seq(10000,70000, 10000)), labels=FALSE)
  labels <- sapply(c(3,4), function(i) as.expression(bquote(10^ .(i))))
  axis(2, at=c(1000, 10000), labels=labels)
  abline(h=1500, col="red", lwd=2) # as in Luecken & Theis
  
  # D
  #Create a function to generate a continuous color palette
  #rbPal <- colorRampPalette(c('red','blue'))
  
  #This adds a column of color values based on the y values
  range(seurat@meta.data$percent.mito)
  # [1] 
  gradCol <- seurat@meta.data$percent.mito
  #myCol <- viridis(10)[1:5][as.numeric(cut(gradCol,breaks = 5))]
  myCol <- viridis(10)[as.numeric(cut(gradCol,breaks = seq(0,100,10)))]
  layout(matrix(1:2,ncol=2), width = c(3,1),height = c(1,1))
  par(mar=c(5.1,4.1,4.1,0.1))
  plot(counts_per_cell, genes_per_cell, log='xy', col=myCol, xlab="Count depth", ylab="Number of genes", las=2, pch=16, cex=0.6)
  legend_image <- as.raster(matrix(rev(viridis(10)), ncol=1))
  par(mar=c(5.1,1,4.1,2.1))
  plot(c(0,2),c(0,1), type = 'n', axes = F, xaxt='n', xlab = '', ylab='', main = '')
  axis(4, at=seq(0,1,0.1), labels=FALSE, pos=0.5)
  text(x=0.9, y = seq(0,1,0.1), labels = seq(0,100,10), adj=0, cex = 0.6)
  mtext("Percentage mitochondrial counts", side=4, line=-1)
  rasterImage(legend_image, -1, 0, 0.5, 1)
  
  # And reset plotting margins,
  par(mar=c(5.1,4.1,4.1,2.1))
  par(mfrow=c(1,1))
  
  # Extra
  plot(sort(genes_per_cell), log='y', main='genes per cell (ordered)', xlab="Cell rank", ylab="Number of genes", las=2)
  
  # Cell cycle markers
  FeatureScatter(object = seurat, feature1 = "S.Score", feature2 = "nFeature_RNA")
  
  dev.off()
  
  # And plot
  pdf(paste(outDir,paste0(newName,"_violinPlot_afterFiltering.pdf"), sep="/"))
  Idents(seurat) <- "hash.ID"
  p1 <- VlnPlot(object = seurat, features = c("nFeature_RNA"), pt.size = 0.1)
  p2 <- VlnPlot(object = seurat, features = c("nCount_RNA"), pt.size = 0.1)
  p <- p1|p2
  print(p)
  p3 <- VlnPlot(object = seurat, features = c("percent.mito"), pt.size = 0.1) + ggtitle("percentage mito")
  p4 <- VlnPlot(object = seurat, features = c("n.exp.hkgenes"), pt.size = 0.1)
  p <- p3|p4
  print(p)
  dev.off()
  
  pdf(paste(outDir,paste0(newName, "_nrGenes_vs_Counts_afterFiltering.pdf"), sep="/"))
  Idents(seurat) <- "hash.ID"
  p <- FeatureScatter(object = seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  print(p)
  #
  dev.off()
  
  # And clean
  rm(counts_per_cell, counts_per_gene, genes_per_cell, cells_per_gene)
  rm(p1,p2,p3,p4, p.overview, p.detail, p.all)
  
  # And clean
  rm(seurat)
}

#### Normalize RNA ####
# Normalize RNA data with SCTransform
file_names <- dir(outDir, pattern = "*.rds")
file_names <- file_names[grepl("Mito",file_names)]

for ( rds in file_names ){
  # Read data
  seurat <- readRDS(paste0(outDir,"/",rds))
  
  DefaultAssay(seurat) <- "RNA"
  
  # File name
  newName <- gsub(".rds","", rds)
  
  # Also regress out "S.score" etc. ??
  seurat <- SCTransform(seurat, vars.to.regress = "percent.mito")
  # ncells	              
  #    Number of subsampling cells used to build NB regression; default is 5000
  # residual.features	    
  #    Genes to calculate residual features for; default is NULL (all genes). If specified, will be set to VariableFeatures of the returned object.
  # variable.features.n	  
  #    Use this many features as variable features after ranking by residual variance; default is 3000. Only applied if residual.features is not set.
  ###  NOTE: AJ - I do not see this in the output ... it looks like it is still using 2000 features:
  #      Variance stabilizing transformation of count matrix of size 17034 by 9525
  #      Model formula is y ~ log_umi
  #      Get Negative Binomial regression parameters per gene
  #      Using 2000 genes, 5000 cells
  #    BUT: ?SCTransform indeed says variable.features.n = 3000
  #
  # vars.to.regress	
  #    Variables to regress out in a second non-regularized linear regression. For example, percent.mito. Default is NULL
  # do.scale	
  #    Whether to scale residuals to have unit variance; default is FALSE
  # do.center	
  #    Whether to center residuals to have mean zero; default is TRUE
  
  # https://satijalab.org/seurat/articles/sctransform_vignette.html
  # The latest version of sctransform also supports using glmGamPoi package which substantially improves the speed
  # of the learning procedure. It can be invoked by specifying method="glmGamPoi".

  # And save
  saveRDS(seurat, paste0(outDir,"/", newName, ".SCTransformed.rds"))
  
  print(seurat) # And yes, now ADT is still present!
}

## OBSERVATION
# - Generally, between 60-110 outliers are found, that will be ignored in the fitting/regularization step
# - MOMA72 has much more of a 'knee' in the 'counts vs. cells (ranked)' and 'Number of genes vs. cells' plots
#   than MOMA302 - probably as there are less cells sequenced??

#### And normalize ADT ####
# Normalize using CLR
file_names <- dir(outDir, pattern = "SCTransformed.rds")
file_names <- file_names[grepl("Mito",file_names)]

for ( rds in file_names ){
  # Read data
  seurat <- readRDS(paste0(outDir,"/",rds))
  
  DefaultAssay(seurat) <- "ADT"
  seurat <- NormalizeData(seurat, normalization.method = "CLR", margin = 2)

  # And save
  saveRDS(seurat, paste0(outDir,"/", rds))
}