## VISION - How To: Working with Seurat: https://yoseflab.github.io/VISION/articles/web_only/Seurat.html

library(ggplot2)
library(Seurat)
library(VISION)

outputDir <- "12b_reclustering_vision"
dir.create(outputDir)

# Hallmark and C8 sets
signatures <- list(hallmark = c("../../../../../MSigDB/v7.4/h.all.v7.4.symbols.gmt"),
                   c8 = c("../../../../../MSigDB/v7.4/c8.all.v7.4.symbols.gmt"))

# Adding Signatures for celltypes
# see mail Yosta 20230206 and Perry 20230206 (and https://yoseflab.github.io/VISION/articles/Signatures.html)
# and see `12a_vision.r`

Classical = c("HLA-DRB5" = 1, "CD68" = 1, "CD14" = 1, "CD36" = 1,"S100A9" = 1,"LYZ" = 1,"VCAN" = 1)
NonClassical = c("HLA-DRB5" = 1, "CD68" = 1, "FCGR3A" = 1, "CX3CR1" = 1, "CCR2" = -1)
TissueRes = c("HLA-DRB5"= 1, "CD68" = 1, "CD74" = 1, "CD81" = 1, "C1QA" = 1)
T_cell = c("CD3E" = 1, "CD4" = 1,"CD8A" = 1, "GZMB" = 1, "GZMK" = 1, "CD27" = 1)
cDC = c("HLA-DRB5" = 1,"CD1C" = 1,"FLT3" = 1,"CLEC10A" = 1, "CD14" = -1)
pDC = c("HLA-DRB5" = 1,"IL3RA" = 1,"CLEC4C" = 1,"GZMB" = 1)
NK = c("NCAM1" = 1, "FCGR3A" = 1, "HLA-DRB5" = 1, "HLA-DRB5" = -1, "CD3E" = -1)

cellTypeList <- list(Classical_monocyte_derived_Macrophages=Classical,
                     Non_classical_monocyte_derived_Macrophages=NonClassical,
                     Tissue_Resident_Macrophages=TissueRes,
                     T_cells=T_cell, cDCs=cDC, pDCs=pDC, NKs=NK)

sig <- createGeneSignature(name = "Classical_monocyte_derived_Macrophages_", sigData = Classical)
sig1 <- createGeneSignature(name = "Non_classical_monocyte_derived_Macrophages_", sigData = NonClassical)
sig2 <- createGeneSignature(name = "Tissue_Resident_Macrophages_", sigData = TissueRes)
sig3 <- createGeneSignature(name = "T_Cells_", sigData = T_cell)
sig4 <- createGeneSignature(name = "cDCs_", sigData = cDC)
sig5 <- createGeneSignature(name = "pDCs_", sigData = pDC)

NK = c("NCAM1" = 1, "FCGR3A" = 1, "HLA-DRB5" = 1, "CD3E" = -1)
sig6 <- createGeneSignature(name = "NKs", sigData = NK)

# Construct a second signature for the NK cells with negative HLA-DRB5
NK_pos = c("NCAM1" = 1, "FCGR3A" = 1, "HLA-DRB5" = 1, "CD3E" = -1)
sig6 <- createGeneSignature(name = "NK_pos_", sigData = NK_pos)
NK_neg = c("NCAM1" = 1, "FCGR3A" = 1, "HLA-DRB5" = -1, "CD3E" = -1)
sig7 <- createGeneSignature(name = "NK_neg_", sigData = NK_neg)

# Finally!
signatures <- c(sig, sig1, sig2, sig3, sig4, sig5, sig6, sig7, "../../../../../MSigDB/v7.4/h.all.v7.4.symbols.gmt", "../../../../../MSigDB/v7.4/c8.all.v7.4.symbols.gmt")

#### Read the data ####
seurat <- readRDS("10_Reclustering_MP_and_MC/selectedClusters.subjectIntegrated.reclustered.kidney_immune_reference.rds")

#### Set the DefaultAssay to 'SCT' ####
DefaultAssay(seurat) <- 'SCT'

#### Slim down and reorder the meta data ####
seurat@meta.data$disease[seurat@meta.data$orig.ident %in% c("MOMA17","MOMA52","MOMA67","MOMA72","MOMA302")] <- "ANCA"
seurat@meta.data$disease[seurat@meta.data$orig.ident %in% c("MOMA68")] <- "SLE"
seurat@meta.data$disease[seurat@meta.data$orig.ident %in% c("MOMA57")] <- "Control"

# See also mail with Yosta 20220919
seurat@meta.data$anca.type[seurat@meta.data$orig.ident %in% c("MOMA17","MOMA52", "MOMA302")] <- "PR3"
seurat@meta.data$anca.type[seurat@meta.data$orig.ident %in% c("MOMA67", "MOMA72")] <- "MPO"

seurat@meta.data$celltype.disease <- factor(with(seurat@meta.data, paste0(predicted.labels.celltypist.after,"_",disease)))
seurat@meta.data$celltype.anca.type <- factor(with(seurat@meta.data, paste0(predicted.labels.celltypist.after,"_",anca.type)))

# 20220825 - AJ: convert character columns to factors
seurat@meta.data[sapply(seurat@meta.data, is.character)] <- lapply(seurat@meta.data[sapply(seurat@meta.data, is.character)], as.factor)

# Remove most of the CellTypist annotation before integration in order to slim down the object
# Otherwise you will run out of memory when analyzing the VISION object
seurat@meta.data <- seurat@meta.data[,!(grepl("celltypist$|^[A-Z].*celltypist\\.after$", colnames(seurat@meta.data), ignore.case = FALSE))]
seurat@meta.data <- seurat@meta.data[,!(grepl("^DF|^pANN", colnames(seurat@meta.data), ignore.case = FALSE))]

# Reorder the metadata  
metadata <- seurat@meta.data
metanames <- colnames(metadata)
metanames.phenotype  <- c("orig.ident", "disease", "anca.type", "celltype.disease", "celltype.anca.type")
metanames.cluster <- c(metanames[grep("integrated",metanames)], metanames[grep("SCT_snn",metanames)], "seurat_clusters")
metanames.celltype <- c(metanames[grep("mature",metanames)], metanames[grep("celltypist",metanames)],
                        metanames[grep("hemapoietic",metanames)], metanames[grep("hpca",metanames)], 
                        metanames[grep("immCell",metanames)], metanames[grep("monaco",metanames)]) 
metanames.vdj <- metanames[grep("TB|vdj",metanames)]
metanames.remainder <- setdiff(metanames, union(metanames.phenotype, union(metanames.cluster, union(metanames.celltype, metanames.vdj))))
metanames <- c(metanames.phenotype, metanames.remainder, metanames.cluster, metanames.celltype, metanames.vdj)
metanames <- metanames[!(metanames %in% c("old.ident", "Cluster1"))]
seurat@meta.data <- metadata[,metanames]

# And clean
rm(metanames, metanames.celltype, metanames.cluster, metanames.phenotype, metanames.remainder,
   metanames.vdj, metadata)

#### Add the gene signatures also as module scores ####
# - But should probably not add the downregulated genes?
Classical = c("HLA-DRB5" = 1, "CD68" = 1, "CD14" = 1, "CD36" = 1,"S100A9" = 1,"LYZ" = 1,"VCAN" = 1)
NonClassical = c("HLA-DRB5" = 1, "CD68" = 1, "FCGR3A" = 1, "CX3CR1" = 1)
TissueRes = c("HLA-DRB5"= 1, "CD68" = 1, "CD74" = 1, "CD81" = 1, "C1QA" = 1)
T_cell = c("CD3E" = 1, "CD4" = 1,"CD8A" = 1, "GZMB" = 1, "GZMK" = 1, "CD27" = 1)
cDC = c("HLA-DRB5" = 1,"CD1C" = 1,"FLT3" = 1,"CLEC10A" = 1)
pDC = c("HLA-DRB5" = 1,"IL3RA" = 1,"CLEC4C" = 1,"GZMB" = 1)
NK = c("NCAM1" = 1, "FCGR3A" = 1, "HLA-DRB5" = 1)

cellTypeList <- list(Classical_monocyte_derived_Macrophages=Classical,
                     Non_classical_monocyte_derived_Macrophages=NonClassical,
                     Tissue_Resident_Macrophages=TissueRes,
                     T_cells=T_cell, cDCs=cDC, pDCs=pDC, NKs=NK)

# see https://github.com/satijalab/seurat/issues/5841 and 1965
# But set the DefaultAssay to "SCT" before, so this is probably redundant here ...
DefaultAssay(seurat) <- "SCT"
for (cellType in names(cellTypeList)){
  seurat <- AddModuleScore(
    object = seurat,
    features = list(c(names(cellTypeList[[cellType]]))),
    ctrl = 100,   # default
    name = cellType
  )
  # And visualise
  pdf(paste(outputDir,paste0("moduleScores_",cellType,".pdf"), sep ="/"), width=14)
  for (clust in colnames(seurat@meta.data)[grep("integrated_snn", colnames(seurat@meta.data))]){
    Idents(seurat) <- clust
    res <- gsub("integrated_snn_r","R", clust)
    p1 <- DimPlot(seurat, reduction = "umap", label = TRUE,  repel = TRUE) +
      guides(colour=guide_legend(override.aes = list(size=2), ncol=1))
    p1[[1]]$layers[[1]]$aes_params$alpha = 0.6
    p1[[1]]$layers[[1]]$mapping$alpha = 0.4
    p2 <- FeaturePlot(object = seurat, features = paste0(cellType,"1"), pt.size = 0.1, reduction = "umap" ) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
    p3 <- VlnPlot(object = seurat, features = paste0(cellType,"1"), pt.size = 0.1, y.max = NULL, same.y.lims = TRUE) + NoLegend() +
      theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
      ggtitle(res)
    p <- p1 | p2
    print(p)
    p <- p1 | p3
    print(p)
  }
  dev.off()
}

#### Prepare protein data - ADT ####
ADT <- t(seurat@assays$ADT@data)
ADTnames <- colnames(ADT)
ADTnames <- sub("anti-human ","", ADTnames)
colnames(ADT) <- ADTnames
# Next line doesn't have any effect on how the ADTs are shown in VISION
# ADT <- ADT[, order(colnames(ADT))]


#### Run Vision ####
# I do not need the `sig_gene_threshold = 0` anymore ... see the discussion above on the signatures
# Still need `min_signature_genes = 3` as some of the celltype signatires are small ...
# - dimRedComponents: Selected 30 components, since this was the number determined using the elbowPlot
# - To instruct VISION to not run any additional visualization projections, set projection_methods = NULL 

vision.obj <- Vision(seurat, signatures = signatures, proteinData = ADT, 
                     dimRedComponents = 30, projection_methods = NULL, 
                     meta=seurat@meta.data,
                     min_signature_genes = 3)
# Using 17254/22769 genes detected in 0.10% of cells for signature analysis.

# Get rid of the large seurat object to save space 
rm(seurat, ADT, sig, sig1, sig2, sig3, sig4, sig5, sig6, sig7, signatures, ADTnames)

#### Analyze ####
vision.obj <- analyze(vision.obj)
# NOTE:
# - If you keep a lot of the meta data columns the object becomes big and this calculation will
#   run out of memory
# - This is why we made a selection (see `Slim down and reorder meta data`)

#### Save ####
saveRDS(vision.obj, paste0(outputDir,"/vision.rds"))

### Visualize ####
vision.obj <- readRDS(paste0(outputDir,"/vision.rds"))
viewResults(vision.obj)

#### Plots ####
head(getSignatureAutocorrelation(vision.obj), n=10)

# Plot signature scores for a signature of interest
umap <- getProjections(vision.obj)[["Seurat_umap"]]
sigScores <- getSignatureScores(vision.obj)[, "LAKE_ADULT_KIDNEY_C5_PROXIMAL_TUBULE_EPITHELIAL_CELLS_STRESS_INFLAM"]

pdf(paste0(outputDir,"/LAKE_ADULT_KIDNEY_C5_PROXIMAL_TUBULE_EPITHELIAL_CELLS_STRESS_INFLAM_PTPRC.pdf"), width=14)
  ggplot() + aes(x=umap[, 1], y=umap[, 2], color=sigScores) + geom_point()
  ggplot() + aes(x=umap[, 1], y=umap[, 2], color=log2(vision.obj@exprData["PTPRC",])) + geom_point() +
    ggtitle("PTPRC")
dev.off()

# Plot signature scores for a signature of interest
umap <- getProjections(vision.obj)[["Seurat_umap"]]
sigScores <- getSignatureScores(vision.obj)[, "HALLMARK_INTERFERON_ALPHA_RESPONSE"]

pdf(paste0(outputDir,"/HALLMARK_INTERFERON_ALPHA_RESPONSE.pdf"), width=14)
ggplot() + aes(x=umap[, 1], y=umap[, 2], color=sigScores) + geom_point()
ggplot() + aes(x=umap[, 1], y=umap[, 2], color=log2(vision.obj@exprData["CD74",])) + geom_point() +
  ggtitle("CD74")
dev.off()
