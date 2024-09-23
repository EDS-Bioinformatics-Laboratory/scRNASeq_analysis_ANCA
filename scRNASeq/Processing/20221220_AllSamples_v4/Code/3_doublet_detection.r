# RNA analysis
#

#### IMPORTANT NOTE ####
# 20230515 - Perry saw that I did not adjust the 'pK' value in the 'doubletFinder_v3' function
#            and that I also used the 'unadjusted' (for homotypic doublets) classification when removing the 'Doublets'
#            see line 175
# See '3_doublet_detection_20230514_CHECK.r' for some checks..

#### Set default for a good start ####
source("initAnalysis.r")
init()

# Create a directory to hold results and plots
outDir <- "3_doublet_detection"
dir.create(outDir)

#### Library ####
library(Seurat)
library(ggplot2)
library(cowplot)
library(patchwork)
library(grid)
library(DoubletFinder)
library(qpdf)

# Sample list
sampleList <- c("MOMA17", "MOMA52", "MOMA57", "MOMA67", "MOMA68", "MOMA72", "MOMA302")

#### Extract singlets ####
# Read the data from the previous step
file_names <- dir("./2_filtering_normalization_and_scaling", pattern = "*.rds")
file_names <- file_names[grepl("SCTransformed",file_names)]

for ( rds in file_names ){
  # Read data
  seurat <- readRDS(paste0("2_filtering_normalization_and_scaling/",rds))
  
  # Set the DefaultAssay to SCT, i.e. the normalized data
  DefaultAssay(seurat) <- "SCT"
  
  # File name
  newName <- gsub(".rds","", rds)
  
  # Detect doublets
  # - https://github.com/EDePasquale/DoubletDecon
  #   Problems installing loads of packages ...
  # - Benchmark paper: https://pubmed.ncbi.nlm.nih.gov/33338399/ 
  #   "Overall, the DoubletFinder method has the best detection accuracy, 
  #    and the cxds method has the highest computational efficiency"
  # - DoubletFinder: https://github.com/chris-mcginnis-ucsf/DoubletFinder
  #   You need to do the whole processing pipeline and then return to look for doublets
  # - Other R option: https://github.com/plger/scDblFinder - can be run at the beginning of the analysis
  
  # Now using DoubletFinder
  seurat <- FindVariableFeatures(seurat, selection.method = "vst", nfeatures = 2000)
  print(paste0("@RNA: ",length(seurat@assays$RNA@var.features)))
  # [1] 0
  print(paste0("@SCT: ",length(seurat@assays$SCT@var.features)))
  # [1] 2000
        
  # Scale data, regressing out nCount_RNA and percent.mito
  # But, do we need scaled date somewhere in our plots?? :-)
  # seurat <- ScaleData(seurat, vars.to.regress = c("nCount_RNA", "percent.mito"))
  # Already did this using SCTransform !!!  
  # -> You need to specify the assay you working on, see line 34..
  # This happened because, I added the ADT and set the DefaultAssay to 'ADT' in the previous step
  
  # Run PCA
  # The default is to store 50 pc's, increase this to 100
  seurat <- RunPCA(seurat, features = VariableFeatures(seurat), npcs = 100, nfeatures.print = 10)
  
  # NOTE:
  #    JackStraw cannot be run on SCTransform-normalized data.
  # SEE:
  #    https://satijalab.org/seurat/articles/sctransform_vignette.html
  
  # But the ElbowPlot can be made
  pdf(paste(outDir,paste0(newName,"_ElbowPlot_dims100.pdf"), sep="/"))
     p <- ElbowPlot(seurat, ndims=100)
     print(p)
  dev.off()

  # And save
  saveRDS(seurat, paste0(outDir,"/", newName, ".pca.rds"))
  
  # And remove
  rm(seurat)
}

# Decided to use a 10% cut-off for the mitochondrial content (based on previous analysis - 20220623)
# The choice for the number of PCs based on the Elbow plot generated earlier (see 3_doublet_detection_INITIAL) was 30
# 
# NOTE:
#   - Although this might now have changes as we filtered on the number of housekeeping genes ... 
#   - Inspection of the ElbowPlots generated in the previous analysis (using more or less the same data) suggests 35

pcs <- list(Mito10=35)

for ( sc in sampleList){
  for ( i in 1:length(pcs)){
    mitoDir <- paste0(outDir,"/", paste0(sc,"_pct",names(pcs)[i]))
    dir.create(mitoDir)

    # Read the data
    seurat <- readRDS(paste0(outDir,"/",sc,".filt",names(pcs)[i],".SCTransformed.pca.rds"))
    
    # Remove the previous DoubletFinder columns if present
    idx <- grep(c("^pAnn|^DF.classification|^SCT|seurat_clusters"),colnames(seurat@meta.data))
    if ( length(idx) > 0){
      print("Removing previous results of DoubletFinder")
      seurat@meta.data <- seurat@meta.data[,-idx]
    }
    
    # Run TSNE and UMAP
    seurat <- RunTSNE(seurat, reduction = "pca", dims = 1:pcs[i][[1]], dim.embed = 3)
    seurat <- RunUMAP(seurat, reduction = "pca", dims = 1:pcs[i][[1]], verbose = FALSE)

    # Get number of PCs from list
    seurat <- FindNeighbors(seurat, reduction = "pca", dims = 1:pcs[i][[1]])

    # Find clusters using various resolution settings, 0.3 and 1.8, arbitrarily
    # NOTE:
    # - Do we have to change this based on the later clusterings and decisions?
    # - Should we still do an initial doublet detection as we did previously??
    seurat <- FindClusters(seurat, resolution = 0.3)
    # Number of communities: 7
    seurat <- FindClusters(seurat, resolution = 1.8)
    # resolution: 1.8 -> 13 communities
    
    # Save
    saveRDS(seurat, paste0(outDir,"/",sc,".",names(pcs)[i],".clustered.rds"))
    
    # Now, run DoubletFinder 
    # - Perhaps first take out the odd cluster?
    ## As in the tutorial:
    #pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
    sweep.res.list_kidney <- paramSweep_v3(seurat, PCs = 1:pcs[i][[1]], sct = TRUE)  # Put 'sct = TRUE' to work with SCTransform normalized data
    sweep.stats_kidney <- summarizeSweep(sweep.res.list_kidney, GT = FALSE)
    bcmvn_kidney <- find.pK(sweep.stats_kidney)
    
    ## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
    annotations <- seurat@meta.data$SCT_snn_res.0.3  # Or take resolution 1.8??
    homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
    nExp_poi <- round(0.075*nrow(seurat@meta.data))         ## Assuming 7.5% doublet formation rate - tailor for your dataset - We expect 8%?
    nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
    
    ## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
    seurat <- doubletFinder_v3(seurat, PCs = 1:pcs[i][[1]], pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)
    # [1] "Creating 590 artificial doublets...
    idx <- which(grepl("pANN",colnames(seurat@meta.data)))
    seurat <- doubletFinder_v3(seurat, PCs = 1:pcs[i][[1]], pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, 
                               reuse.pANN = colnames(seurat@meta.data)[idx], sct = TRUE)
    
    idx <- which(grepl("DF.classifications",colnames(seurat@meta.data)))
    if (length(idx) >= 2 ){
      idx <- idx[2]
    }
    
    print("Plotting the singlets and doublets")

    pdf(paste(mitoDir,"dimPlot_DoubletFinder.pdf", sep="/"))
      Idents(seurat) <- colnames(seurat@meta.data)[idx]
      p <- DimPlot(seurat, reduction = 'umap', label = FALSE) + ggtitle(paste0(sc," - pct",names(pcs)[i]))
      print(p)
      p <- DimPlot(seurat, reduction = 'tsne', label = FALSE) + ggtitle(paste0(sc," - pct",names(pcs)[i]))
      print(p)
    dev.off()
    
    # Save the object with the singlet/doublet annotation
    saveRDS(seurat, paste0(outDir,"/",sc,".filt",names(pcs)[i],".SCTransformed.pca.rds"))
    
    # Remove the doublets
    # NOTE: Here, I use the 'unadjusted' (for homotypic doublets) classification ... otherwise use idx[2]!!
    Idents(seurat) <- colnames(seurat@meta.data)[idx[1]]
    seurat <- subset(x = seurat, idents = "Singlet")
    dim(seurat)
    # 
    
    # Save intermediate, but not the scaled data - use DietSeurat
    seurat <- DietSeurat(seurat, counts = TRUE, data = TRUE, scale.data = FALSE, features = NULL,
                         assays = Assays(seurat), dimreducs = Reductions(seurat), graphs = Graphs(seurat) )
    saveRDS(seurat, paste0(outDir,"/", sc, ".", names(pcs)[i],".doubletsRemoved.rds"))
    
    # And check whether the ADT assay is still there...
    # print(seurat)  # - And it is!
  }
}

#### Print some results ####
for ( sc in sampleList){
  for ( i in 1:length(pcs)){
    # Read the data
    seurat <- readRDS(paste0(outDir,"/",sc,".filt",names(pcs)[i],".SCTransformed.pca.rds"))
    
    idx <- which(grepl("DF.classifications",colnames(seurat@meta.data)))
    if (length(idx) >= 2 ){
      idx <- idx[2]
    }
    
    print(sc)
    print(table(seurat@meta.data[,idx]))
    
    # NOTE: This is not taking the same column as the one I used to extract Singlets and Doublets !!!
  }
}

# [1] "MOMA17"
# Doublet Singlet 
# 81    1225 

# [1] "MOMA52"
# Doublet Singlet 
# 401    5709 

# [1] "MOMA57"
# Doublet Singlet 
# 609    8565

# [1] "MOMA67"
# Doublet Singlet 
# 310    4659 

# [1] "MOMA68"
# Doublet Singlet 
# 172    2568 

# [1] "MOMA72"
# Doublet Singlet 
# 75    1239 

# [1] "MOMA302"
# Doublet Singlet 
# 176    2785 

# Seems to be rather OK?
# No changes w.r.t. previous analysis! (also no to be expected :-)

## Are cells with both T and B cells (how many remained after filtering step?) now detected as doublets?
#

## MOMA52
seurat <- readRDS(paste0(outDir,"/", "MOMA52.Mito10.doubletsRemoved.rds"))
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
table(seurat@meta.data$TB, seurat@meta.data$DF.classifications_0.25_0.09_458)
# 
#         Doublet Singlet
# B            12     395
# T            17    3749
# T_and_B       2      36
# unknown      26    1472

# @1 Huh, Doublets should have been removed !!?!?!?!
#  - Ah, we have two DF.classifications_0.25_0.09xxx classifications slots and I only look in one of them to
#    make the subset of Singlets
table(seurat@meta.data$TB, seurat@meta.data$DF.classifications_0.25_0.09_458)
#         Singlet
# B           407
# T          3766
# T_and_B      38
# unknown    1498

# @2 14 TB cells have been filtered out as we had 52 to start with 

## MOMA57
seurat <- readRDS(paste0(outDir,"/", "MOMA57.Mito10.doubletsRemoved.rds"))
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
table(seurat@meta.data$TB, seurat@meta.data$DF.classifications_0.25_0.09_688)
#         Doublet Singlet
# B             4     398
# T            24    4815
# T_and_B       0      48
# unknown      51    3225

# Only 4 cells have been filtered out ...

#### First remove TB-cells ####
## Remove the cells with T and B cell receptor sequences before Doublet Finder and run Doublet Finder again ...

# Sample list
sampleList <- c("MOMA52", "MOMA57")

file_names <- dir("./2_filtering_normalization_and_scaling", pattern = "*.rds")
file_names <- file_names[grepl("SCTransformed",file_names)]
file_names <- file_names[grepl(paste(sampleList, collapse="|"),file_names)]

for ( rds in file_names ){
  # Read data
  seurat <- readRDS(paste0("2_filtering_normalization_and_scaling/",rds))

  # Set the DefaultAssay to SCT, i.e. the normalized data
  DefaultAssay(seurat) <- "SCT"
  
  # File name
  newName <- paste0(gsub(".rds","", rds),".noTB_cells")
  
  # Remove the TB cells
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
  
  seurat <- subset(seurat, subset = TB != "T_and_B" )
  
  origDim <- dim(seurat)
  print(paste0(rds, " - number of cells after removing TB cells: ", origDim[2]))
  
  # Detect doublets
  # - https://github.com/EDePasquale/DoubletDecon
  #   Problems installing loads of packages ...
  # - Benchmark paper: https://pubmed.ncbi.nlm.nih.gov/33338399/ 
  #   "Overall, the DoubletFinder method has the best detection accuracy, 
  #    and the cxds method has the highest computational efficiency"
  # - DoubletFinder: https://github.com/chris-mcginnis-ucsf/DoubletFinder
  #   You need to do the whole processing pipeline and then return to look for doublets
  # - Other R option: https://github.com/plger/scDblFinder - can be run at the beginning of the analysis
  
  # Now using DoubletFinder
  seurat <- FindVariableFeatures(seurat, selection.method = "vst", nfeatures = 2000)
  print(paste0("@RNA: ",length(seurat@assays$RNA@var.features)))
  # [1] 0
  print(paste0("@SCT: ",length(seurat@assays$SCT@var.features)))
  # [1] 2000
  
  # Scale data, regressing out nCount_RNA and percent.mito
  # But, do we need scaled date somewhere in our plots?? :-)
  # seurat <- ScaleData(seurat, vars.to.regress = c("nCount_RNA", "percent.mito"))
  # Already did this using SCTransform !!!  
  # -> You need to specify the assay you working on, see line 34..
  # This happened because, I added the ADT and set the DefaultAssay to 'ADT' in the previous step
  
  # Run PCA
  # The default is to store 50 pc's, increase this to 100
  seurat <- RunPCA(seurat, features = VariableFeatures(seurat), npcs = 100, nfeatures.print = 10)
  
  # NOTE:
  #    JackStraw cannot be run on SCTransform-normalized data.
  # SEE:
  #    https://satijalab.org/seurat/articles/sctransform_vignette.html
  
  # But the ElbowPlot can be made
  pdf(paste(outDir,paste0(newName,"_ElbowPlot_dims100.pdf"), sep="/"))
  p <- ElbowPlot(seurat, ndims=100)
  print(p)
  dev.off()
  
  # And save
  saveRDS(seurat, paste0(outDir,"/", newName, ".pca.rds"))
  
  # And remove
  rm(seurat)
}

# Decided to use a 10% cut-off for the mitochondrial content (based on previous analysis - 20220623)
# The choice for the number of PCs based on the Elbow plot generated earlier (see 3_doublet_detection_INITIAL) was 30
# 
# NOTE:
#   - Although this might now have changes as we filtered on the number of housekeeping genes ... 
#   - Inspection of the ElbowPlots generated in the previous analysis (using more or less the same data) suggests 35

pcs <- list(Mito10=35)

for ( sc in sampleList){
  for ( i in 1:length(pcs)){
    mitoDir <- paste0(outDir,"/", paste0(sc,"_noTB_cells_pct",names(pcs)[i]))
    dir.create(mitoDir)
    
    # Read the data
    seurat <- readRDS(paste0(outDir,"/",sc,".filt",names(pcs)[i],".SCTransformed.noTB_cells.pca.rds"))
    
    # Remove the previous DoubletFinder columns if present
    idx <- grep(c("^pAnn|^DF.classification|^SCT|seurat_clusters"),colnames(seurat@meta.data))
    if ( length(idx) > 0){
      print("Removing previous results of DoubletFinder")
      seurat@meta.data <- seurat@meta.data[,-idx]
    }
    
    # Run TSNE and UMAP
    seurat <- RunTSNE(seurat, reduction = "pca", dims = 1:pcs[i][[1]], dim.embed = 3)
    seurat <- RunUMAP(seurat, reduction = "pca", dims = 1:pcs[i][[1]], verbose = FALSE)
    
    # Get number of PCs from list
    seurat <- FindNeighbors(seurat, reduction = "pca", dims = 1:pcs[i][[1]])
    
    # Find clusters using various resolution settings, 0.3 and 1.8, arbitrarily
    # NOTE:
    # - Do we have to change this based on the later clusterings and decisions?
    # - Should we still do an initial doublet detection as we did previously??
    seurat <- FindClusters(seurat, resolution = 0.3)
    # Number of communities: 7
    seurat <- FindClusters(seurat, resolution = 1.8)
    # resolution: 1.8 -> 13 communities
    
    # Save
    saveRDS(seurat, paste0(outDir,"/",sc,".",names(pcs)[i],".noTB_cells.clustered.rds"))
    
    # Now, run DoubletFinder 
    # - Perhaps first take out the odd cluster?
    ## As in the tutorial:
    #pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
    sweep.res.list_kidney <- paramSweep_v3(seurat, PCs = 1:pcs[i][[1]], sct = TRUE)  # Put 'sct = TRUE' to work with SCTransform normalized data
    sweep.stats_kidney <- summarizeSweep(sweep.res.list_kidney, GT = FALSE)
    bcmvn_kidney <- find.pK(sweep.stats_kidney)
    
    ## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
    annotations <- seurat@meta.data$SCT_snn_res.0.3  # Or take resolution 1.8??
    homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
    nExp_poi <- round(0.075*nrow(seurat@meta.data))         ## Assuming 7.5% doublet formation rate - tailor for your dataset - We expect 8%?
    nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
    
    ## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
    seurat <- doubletFinder_v3(seurat, PCs = 1:pcs[i][[1]], pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)
    # [1] "Creating 590 artificial doublets...
    idx <- which(grepl("pANN",colnames(seurat@meta.data)))
    seurat <- doubletFinder_v3(seurat, PCs = 1:pcs[i][[1]], pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, 
                               reuse.pANN = colnames(seurat@meta.data)[idx], sct = TRUE)
    
    idx <- which(grepl("DF.classifications",colnames(seurat@meta.data)))
    if (length(idx) >= 2 ){
      idx <- idx[2]
    }
    
    print("Plotting the singlets and doublets")
    
    pdf(paste(mitoDir,"dimPlot_DoubletFinder.pdf", sep="/"))
    Idents(seurat) <- colnames(seurat@meta.data)[idx]
    p <- DimPlot(seurat, reduction = 'umap', label = FALSE) + ggtitle(paste0(sc," - pct",names(pcs)[i]))
    print(p)
    p <- DimPlot(seurat, reduction = 'tsne', label = FALSE) + ggtitle(paste0(sc," - pct",names(pcs)[i]))
    print(p)
    dev.off()
    
    # Save the object with the singlet/doublet annotation
    saveRDS(seurat, paste0(outDir,"/",sc,".filt",names(pcs)[i],".SCTransformed.noTB_cells.pca.rds"))
    
    # Remove the doublets
    Idents(seurat) <- colnames(seurat@meta.data)[idx[1]]
    seurat <- subset(x = seurat, idents = "Singlet")
    dim(seurat)
    # 
    
    # Save intermediate, but not the scaled data - use DietSeurat
    seurat <- DietSeurat(seurat, counts = TRUE, data = TRUE, scale.data = FALSE, features = NULL,
                         assays = Assays(seurat), dimreducs = Reductions(seurat), graphs = Graphs(seurat) )
    saveRDS(seurat, paste0(outDir,"/", sc, ".", names(pcs)[i],".noTB_cells.doubletsRemoved.rds"))
    
    # And check whether the ADT assay is still there...
    # print(seurat)  # - And it is!
  }
}

#### Print some results - 2 ####
for ( sc in sampleList){
  for ( i in 1:length(pcs)){
    # Read the data
    seurat <- readRDS(paste0(outDir,"/",sc,".filt",names(pcs)[i],".SCTransformed.noTB_cells.pca.rds"))
    
    idx <- which(grepl("DF.classifications",colnames(seurat@meta.data)))
    if (length(idx) >= 2 ){
      idx <- idx[2]
    }
    
    print(sc)
    print(table(seurat@meta.data[,idx]))
  }
}
# [1] "MOMA52"
# Doublet Singlet 
#     402    5656 
#
# [1] "MOMA57"
# Doublet Singlet 
#     605    8517 

# To compare with (see line 203 etc.), i.e. before removal of "TB" cells prior to Doublet detection
# [1] "MOMA52"
# Doublet Singlet 
#     401    5709 

# [1] "MOMA57"
# Doublet Singlet 
#     609    8565

## OBSERVATION:
# - Almost no difference ...


# #### 3D TSNE plot ####
# # For fun ...  ;-)
# source("tSNE_3D_plots.r")
# 
# TSNE_3D_Plot(object=seurat, plotDir=outDir)
# TSNE_3D_GenExprPlot(object=seurat, gene="CD1C", plotDir=outDir)
# 
# # Can we color differently for HTO-1/HTO-2/HTO-3 etc. ?? 

# And clean
rm(seurat)

