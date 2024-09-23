# AJ - 20230112

# Some extra stats, figures etc.

#### Set default for a good start ####
source("initAnalysis.r")
init()

# Create a directory to hold results and plots
outDir <- "11_extras"
dir.create(outDir)

#### Library ####
library(Seurat)
library(cowplot)
library(patchwork)
library(ggplot2)
library(dplyr)

#### Read the data ####
seurat <- readRDS("9_kidney_reference_mapping/Mito10.integrated.clustered.annot.kidney_immune_reference.rds")
DefaultAssay(seurat) <- "RNA"

#### Celltype annotation ####
# How are cells labelled 'MNP-d/Tissue macrophage' according to the annotation of Stewart et al labelled using CellTypist annotation?
tissueResMacrophage.seurat <- subset(seurat, subset = predicted.celltype == "MNP-d/Tissue macrophage")
dim(tissueResMacrophage.seurat)
# [1] 22769   195

# Using CellTYpist probability matching followed by majority voting
table(tissueResMacrophage.seurat@meta.data$probMatchMajorityVoting.labels.celltypist.after)
# Classical monocytes                         DC2     Non-classical monocytes Tem/Temra cytotoxic T cells                  Unassigned 
#                   4                           2                           2                           1                         186 

# And similar for the 'best match' approach
table(tissueResMacrophage.seurat@meta.data$bestMatchNoMajorityVoting.labels.celltypist.after)
# Age-associated B cells          Alveolar macrophages           Classical monocytes                            DC                           DC2 
#                      1                            32                             4                            15                            26 
# Double-positive thymocytes Erythrophagocytic macrophages      Intermediate macrophages        Intestinal macrophages                   Macrophages 
#                          1                            39                             1                            19                            26 
# Mono-mac                      NK cells       Non-classical monocytes                  Plasma cells                  Plasmablasts 
#        2                             1                             4                             1                             1 
# Regulatory T cells      Tcm/Naive helper T cells   Tem/Effector helper T cells   Tem/Temra cytotoxic T cells     Tem/Trm cytotoxic T cells 
#                  3                             1                             2                             9                             7

#### Attempt to diagnose integration and batch effects ####
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7994321/
# https://bioconductor.org/books/release/OSCA.multisample/integrating-datasets.html#mnn-correction
library(batchelor)

#### MNN ####
## Read the data ##
MOMA17 <- readRDS("5_removal_tubuluscells/MOMA17.Mito10.doubletsRemoved.annot.rds")
DefaultAssay(MOMA17) <- "RNA"
MOMA52 <- readRDS("5_removal_tubuluscells/MOMA52.Mito10.noTB_cells.doubletsRemoved.annot.rds")
DefaultAssay(MOMA52) <- "RNA"

library(scran)
all.sce <- list(MOMA17 = MOMA17@assays$RNA@data, MOMA52 = MOMA52@assays$RNA@data)
all.dec <- lapply(all.sce, modelGeneVar)

MOMA17.dec <- all.dec$MOMA17
MOMA52.dec <- all.dec$MOMA52

universe <- intersect(rownames(MOMA17.dec), rownames(MOMA52.dec))
length(universe)
#[1] 15572

MOMA17.dec <- MOMA17.dec[universe,]
MOMA52.dec <- MOMA52.dec[universe,]

combined.dec <- combineVar(MOMA17.dec, MOMA52.dec)
chosen.hvgs <- combined.dec$bio > 0
sum(chosen.hvgs)
#[1] 1565  - Relatively low number compared with example

set.seed(1000101001)
MOMA17 <- subset(MOMA17, features=universe)
MOMA52 <- subset(MOMA52, features=universe)

mnn.out <- fastMNN(MOMA17@assays$RNA@data, MOMA52@assays$RNA@data, d=50, k=20, subset.row=chosen.hvgs,
                   BSPARAM=BiocSingular::RandomParam(deferred=TRUE))
mnn.out

head(mnn.out$batch)
dim(reducedDim(mnn.out, "corrected"))
# [1] 6019   50

snn.gr <- buildSNNGraph(mnn.out, use.dimred="corrected")
clusters.mnn <- igraph::cluster_walktrap(snn.gr)$membership
tab.mnn <- table(Cluster=clusters.mnn, Batch=mnn.out$batch)
tab.mnn

library(scater)
set.seed(0010101010)
mnn.out <- runTSNE(mnn.out, dimred="corrected")

mnn.out$batch <- factor(mnn.out$batch)
plotTSNE(mnn.out, colour_by="batch")

## Now run with all samples and let MNN decide what is the wisest way to merge, i.e. auto.merge = TRUE
MOMA17 <- readRDS("5_removal_tubuluscells/MOMA17.Mito10.doubletsRemoved.annot.rds")
DefaultAssay(MOMA17) <- "RNA"
MOMA52 <- readRDS("5_removal_tubuluscells/MOMA52.Mito10.noTB_cells.doubletsRemoved.annot.rds")
DefaultAssay(MOMA52) <- "RNA"
MOMA57 <- readRDS("5_removal_tubuluscells/MOMA57.Mito10.noTB_cells.doubletsRemoved.annot.rds")
DefaultAssay(MOMA57) <- "RNA"
MOMA67 <- readRDS("5_removal_tubuluscells/MOMA67.Mito10.doubletsRemoved.annot.rds")
DefaultAssay(MOMA67) <- "RNA"
MOMA68 <- readRDS("5_removal_tubuluscells/MOMA68.Mito10.doubletsRemoved.annot.rds")
DefaultAssay(MOMA68) <- "RNA"
MOMA72 <- readRDS("5_removal_tubuluscells/MOMA72.Mito10.doubletsRemoved.annot.rds")
DefaultAssay(MOMA72) <- "RNA"
MOMA302 <- readRDS("5_removal_tubuluscells/MOMA302.Mito10.doubletsRemoved.annot.rds")
DefaultAssay(MOMA302) <- "RNA"

library(scran)
# Get the raw counts and let 'scater' normalize etc .. ?
all.sce <- list(MOMA17 = as.SingleCellExperiment(MOMA17), MOMA52 = as.SingleCellExperiment(MOMA17),
                MOMA57 = as.SingleCellExperiment(MOMA57), MOMA67 = as.SingleCellExperiment(MOMA67),
                MOMA68 = as.SingleCellExperiment(MOMA68), MOMA72 = as.SingleCellExperiment(MOMA72),
                MOMA302 = as.SingleCellExperiment(MOMA302))

all.sce <- lapply(all.sce, logNormCounts)
lapply(all.sce, function(x) summary(sizeFactors(x)))
# $MOMA17
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.07676 0.65859 0.84480 1.00000 1.19360 5.33451 
# 
# $MOMA52
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.07676 0.65859 0.84480 1.00000 1.19360 5.33451 
# 
# $MOMA57
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.2243  0.6554  0.8128  1.0000  1.0637  8.1444 
# 
# $MOMA67
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.2074  0.6671  0.8028  1.0000  1.0154  7.9030 
# 
# $MOMA68
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.09956 0.38205 0.70065 1.00000 1.25784 7.03940 
# 
# $MOMA72
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.06211 0.55142 0.70005 1.00000 1.19038 4.31477 
# 
# $MOMA302
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.1138  0.6714  0.8352  1.0000  1.0627  8.5711 

all.dec <- lapply(all.sce, modelGeneVar)

quick.corrected2 <- quickCorrect(`MOMA17`=all.sce$MOMA17, `MOMA52`=all.sce$MOMA52, `MOMA57`=all.sce$MOMA57,
                                 `MOMA67`=all.sce$MOMA67, `MOMA68`=all.sce$MOMA68, `MOMA72`=all.sce$MOMA72, `MOMA302`=all.sce$MOMA302,
                                 precomputed=list(all.dec$MOMA17, all.dec$MOMA52, all.dec$MOMA67, all.dec$MOMA68,
                                                  all.dec$MOMA68, all.dec$MOMA72, all.dec$MOMA302),
                                 PARAM=FastMnnParam(BSPARAM=BiocSingular::RandomParam(), auto.merge=TRUE))

quick.sce2 <- quick.corrected2$corrected

quick.corrected2$corrected@metadata$merge.info
# DataFrame with 6 rows and 6 columns
#                        left  right                               pairs batch.size   skipped                            lost.var
#                      <List> <List>                     <DataFrameList>  <numeric> <logical>                            <matrix>
# 1                   MOMA302 MOMA67 17652:10952,17652:13775,17652:12651   0.439007     FALSE 0.00000000:0.0000000:0.00000000:...
# 2            MOMA302,MOMA67 MOMA57    17660:6317,17660:8103,17660:2689   0.769610     FALSE 0.00000000:0.0000000:0.02176111:...
# 3     MOMA302,MOMA67,MOMA57 MOMA72 17661:16747,17666:17246,17670:16842   0.623978     FALSE 0.00000000:0.0000000:0.00818093:...
# 4 MOMA302,MOMA67,MOMA57,... MOMA68 17653:16120,17658:16089,17660:14736   0.597808     FALSE 0.00000000:0.0000000:0.01577933:...
# 5 MOMA302,MOMA67,MOMA57,... MOMA17       17653:437,17675:449,17675:324   0.495319     FALSE 0.03465704:0.0000000:0.00391166:...
# 6 MOMA302,MOMA67,MOMA57,... MOMA52       17653:903,17659:734,17671:771   0.582561     FALSE 0.00106613:0.0383126:0.00120855:...

set.seed(00101010)
quick.sce2 <- runTSNE(quick.sce2, dimred="corrected")
quick.sce2 <- runUMAP(quick.sce2, dimred="corrected")

pdf(paste(outDir,paste0("plots_MNN_quickCorrectBatch.pdf"), sep="/"), width = 7, height = 7) 
  plotUMAP(quick.sce2, colour_by="batch")
  plotTSNE(quick.sce2, colour_by="batch")
dev.off()

#### LIGER ####
# Just following https://broadinstitute.github.io/2019_scWorkshop/correcting-batch-effects.html
install.packages('rliger')
library(rliger)

# Read the data 
MOMA17 <- readRDS("5_removal_tubuluscells/MOMA17.Mito10.doubletsRemoved.annot.rds")
DefaultAssay(MOMA17) <- "RNA"
MOMA52 <- readRDS("5_removal_tubuluscells/MOMA52.Mito10.noTB_cells.doubletsRemoved.annot.rds")
DefaultAssay(MOMA52) <- "RNA"
MOMA57 <- readRDS("5_removal_tubuluscells/MOMA57.Mito10.noTB_cells.doubletsRemoved.annot.rds")
DefaultAssay(MOMA57) <- "RNA"
MOMA67 <- readRDS("5_removal_tubuluscells/MOMA67.Mito10.doubletsRemoved.annot.rds")
DefaultAssay(MOMA67) <- "RNA"
MOMA68 <- readRDS("5_removal_tubuluscells/MOMA68.Mito10.doubletsRemoved.annot.rds")
DefaultAssay(MOMA68) <- "RNA"
MOMA72 <- readRDS("5_removal_tubuluscells/MOMA72.Mito10.doubletsRemoved.annot.rds")
DefaultAssay(MOMA72) <- "RNA"
MOMA302 <- readRDS("5_removal_tubuluscells/MOMA302.Mito10.doubletsRemoved.annot.rds")
DefaultAssay(MOMA302) <- "RNA"

# Next we create a LIGER object with raw counts data from each batch.
# You need to ,ake the columnames unique ... or 'createLiger' will stumble
MOMA17.counts <- MOMA17@assays$RNA@counts
colnames(MOMA17.counts) <- gsub("-1","_MOMA17", colnames(MOMA17.counts))
MOMA52.counts <- MOMA52@assays$RNA@counts
colnames(MOMA52.counts) <- gsub("-1","_MOMA52", colnames(MOMA52.counts))
MOMA57.counts <- MOMA57@assays$RNA@counts
colnames(MOMA57.counts) <- gsub("-1","_MOMA57", colnames(MOMA57.counts))
MOMA67.counts <- MOMA67@assays$RNA@counts
colnames(MOMA67.counts) <- gsub("-1","_MOMA67", colnames(MOMA67.counts))
MOMA68.counts <- MOMA68@assays$RNA@counts
colnames(MOMA68.counts) <- gsub("-1","_MOMA68", colnames(MOMA68.counts))
MOMA72.counts <- MOMA72@assays$RNA@counts
colnames(MOMA72.counts) <- gsub("-1","_MOMA72", colnames(MOMA72.counts))
MOMA302.counts <- MOMA302@assays$RNA@counts
colnames(MOMA302.counts) <- gsub("-1","_MOMA302", colnames(MOMA302.counts))

ob.list <- list("MOMA17" = MOMA17.counts, "MOMA52" = MOMA52.counts, "MOMA57" = MOMA57.counts,
                "MOMA67" = MOMA67.counts, "MOMA68" = MOMA68.counts, "MOMA72" = MOMA72.counts,
                "MOMA302" = MOMA302.counts)
data.liger <- createLiger(ob.list)

# Normalize gene expression for each batch.
data.liger <- rliger::normalize(data.liger)

# For variable gene selection, we can either use those we identified in the earlier CCA 
# batch correction analysis (genes.use) or we can use the LIGER function selectGenes().
data.liger <- selectGenes(data.liger, var.thresh = 0.1)
#data.liger@var.genes <- genes.use

# Print out the number of variable genes for LIGER analysis.
print(length(data.liger@var.genes))
#[1] 5494

# Scale the gene expression across the datasets. 
# Why does LIGER not center the data? Hint, think about the use of 
# non-negative matrix factorization and the constraints that this imposes.
data.liger <- scaleNotCenter(data.liger)

# These two steps take 10-20 min. Only run them if you finish with the rest of the code.
# Use the `suggestK` function to determine the appropriate number of factors to use.
# Use the `suggestLambda` function to find the smallest lambda for which the alignment metric stabilizes.
k.suggest <- suggestK(data.liger, return.data = TRUE, num.cores = 4)  # returning data frame for ggplot
# How to read the plot, it suggest k-20 to be reasonable?

# Save data for future use?
saveRDS(k.suggest, paste(outDir,"k.suggest.rds"))

lambda.suggest <- suggestLambda(data.liger, k = 20)
# This operation may take several minutes depending on number of values being tested
# Preprocessing for rep 1: optimizing initial factorization with smallest test lambda=0.25
# |==================================================================================================| 100%
# Finished in 16.65647 mins, 24 iterations.
# Max iterations set: 100.
# Final objective delta: 9.758011e-05.
# Best results with seed 1.
# Testing different choices of lambda values
# Error in { : task 15 failed - "missing value where TRUE/FALSE needed"
saveRDS(lambda.suggest, paste(outDir,"lambda.suggest.rds"))

# Use alternating least squares (ALS) to factorize the matrix.
k.suggest <- 20  # with this line, we do not use the suggested k by suggestK()
data.liger <- optimizeALS(data.liger, k = k.suggest, rand.seed = 1) 

# What do matrices H, V, and W represent, and what are their dimensions?
dim(data.liger@H$MOMA17)  #[1] 466  20
dim(data.liger@V$MOMA17)  #[1] 20 5494
dim(data.liger@W)         #[1] 20 5494

# Let's see what the integrated data looks like mapped onto a tSNE visualization.
data.liger <- runTSNE(data.liger, use.raw = T)
p <- plotByDatasetAndCluster(data.liger, return.plots = T)
pdf(paste(outDir,"plotLiger_DatasetAndCluster.pdf", sep="/"), width = 7, height = 7)
  print(p[[1]]) 
dev.off()

# Next, do clustering of cells in shared nearest factor space, and then quantile alignment.
data.liger <- quantileAlignSNF(data.liger, resolution = 1.0) # SNF clustering and quantile alignment

# What are the dimensions of H.norm. What does this represent? 
dim(data.liger@H.norm)  # [1] 25485    20

# Visualize liger batch correction results.
data.liger <- runTSNE(data.liger)
p <- plotByDatasetAndCluster(data.liger, return.plots = T) 
pdf(paste(outDir,"plotLiger_DatasetAndCluster_quantileAligned.pdf", sep="/"), width = 7, height = 7)
  print(p[[1]])  # plot by dataset
  print(p[[2]])  # plot by dataset
  
  plot_grid(p[[1]], p[[2]])
dev.off()

# Let's look to see how the adjusted rand index changed compared to using no batch correction.
tech <- unlist(lapply(1:length(data.liger@H), function(x) { 
  rep(names(data.liger@H)[x], nrow(data.liger@H[[x]]))}))
clusters <- data.liger@clusters
ari <- data.frame("tech" = tech, "clusters" = clusters)
ari$tech <- plyr::mapvalues(ari$tech, from = c("MOMA17", "MOMA52", "MOMA57", "MOMA67", "MOMA68", "MOMA72", "MOMA302"), to = c(0, 1, 2, 3, 4, 5, 6))
adj.rand.index(as.numeric(ari$tech), as.numeric(ari$clusters)) 
# No such function ....

# Use calcARI ??
calcARI(data.liger, ari$clusters)

# Look at proportion of each batch in each cluster, and look at factor loadings across clusters
pdf(paste(outDir,"plotLiger_ClusterProportions.pdf", sep="/"), width = 7, height = 7)
  p <- plotClusterProportions(data.liger)
  print(p)
dev.off()

p <- plotClusterFactors(data.liger, use.aligned = T)
pdf(paste(outDir,"plotLiger_ClusterFactors.pdf", sep="/"), width = 7, height = 7)
  plotClusterFactors(data.liger, use.aligned = T)
dev.off()

# Look at genes that are specific to a dataset and shared across datasets.
# Use the plotWordClouds function and choose 2 datasets.
pdf(paste0(outDir, "/word_clouds_MOMA17_vs_52.pdf"))
   plotWordClouds(data.liger, dataset1 = "MOMA17", dataset2 = "MOMA52")
dev.off()

# Look at factor loadings for each cell using plotFactors. 
pdf(paste(outDir, "plot_factors.pdf", sep="/"), width = 7, height= 7)
  plotFactors(data.liger)
dev.off()

# Identify shared and batch-specific marker genes from liger factorization.
# Use the getFactorMarkers function and choose 2 datasets.
# Then plot some genes of interest using plotGene and plotGeneViolin.
markers <- getFactorMarkers(data.liger, dataset1 = "MOMA17", dataset2 = "MOMA52", num.genes = 10)
pdf(paste(outDir, "plotLiger_Genes_MOMA17_vs_52.pdf", sep="/"), width=7, height=7)
  plotGene(data.liger, gene = "HLA-DPA1")
  plotGeneViolin(data.liger, gene = "HLA-DPA1")
dev.off()

# Save current progress.
saveRDS(data.liger, paste(outDir, "liger.rds", sep="/"))

#### kBET ####
library(devtools)
install_github('theislab/kBET')

library(kBET)

# Read the data 
MOMA17 <- readRDS("5_removal_tubuluscells/MOMA17.Mito10.doubletsRemoved.annot.rds")
DefaultAssay(MOMA17) <- "RNA"
MOMA52 <- readRDS("5_removal_tubuluscells/MOMA52.Mito10.noTB_cells.doubletsRemoved.annot.rds")
DefaultAssay(MOMA52) <- "RNA"
MOMA57 <- readRDS("5_removal_tubuluscells/MOMA57.Mito10.noTB_cells.doubletsRemoved.annot.rds")
DefaultAssay(MOMA57) <- "RNA"
MOMA67 <- readRDS("5_removal_tubuluscells/MOMA67.Mito10.doubletsRemoved.annot.rds")
DefaultAssay(MOMA67) <- "RNA"
MOMA68 <- readRDS("5_removal_tubuluscells/MOMA68.Mito10.doubletsRemoved.annot.rds")
DefaultAssay(MOMA68) <- "RNA"
MOMA72 <- readRDS("5_removal_tubuluscells/MOMA72.Mito10.doubletsRemoved.annot.rds")
DefaultAssay(MOMA72) <- "RNA"
MOMA302 <- readRDS("5_removal_tubuluscells/MOMA302.Mito10.doubletsRemoved.annot.rds")
DefaultAssay(MOMA302) <- "RNA"

allCounts <- merge(MOMA17@assays$RNA@data, MOMA52@assays$RNA@data, by="row.names", all=TRUE)
rownames(allCounts) <- allCounts$Row.names
allCounts$Row.names <- NULL

allCounts <- merge(allCounts, MOMA57@assays$RNA@data, by="row.names", all=TRUE)
rownames(allCounts) <- allCounts$Row.names
allCounts$Row.names <- NULL

allCounts <- merge(allCounts, MOMA67@assays$RNA@data, by="row.names", all=TRUE)
rownames(allCounts) <- allCounts$Row.names
allCounts$Row.names <- NULL

allCounts <- merge(allCounts, MOMA68@assays$RNA@data, by="row.names", all=TRUE)
rownames(allCounts) <- allCounts$Row.names
allCounts$Row.names <- NULL

allCounts <- merge(allCounts, MOMA72@assays$RNA@data, by="row.names", all=TRUE)
rownames(allCounts) <- allCounts$Row.names
allCounts$Row.names <- NULL

allCounts <- merge(allCounts, MOMA302@assays$RNA@data, by="row.names", all=TRUE)
rownames(allCounts) <- allCounts$Row.names
allCounts$Row.names <- NULL

dim(allCounts)
# [1] 22769 25485

# There are rows with only NA ....
allCounts <- allCounts[rowSums(is.na(allCounts)) != ncol(allCounts), ]
dim(allCounts)
#[1] 14134 25485

# Remove rows that only contain zeroes
allCounts <- allCounts[rowSums(allCounts) > 0, ]
dim(allCounts)
#[1] 14134 25485

# Make the batch vector, i.e. a concatenation of the number of cells from each sample with the right sample name (in the right order!!)
# I could also add the sample name to the cells earlier and then extract? Less error-prone?
batch <- c(rep("MOMA17", ncol(MOMA17@assays$RNA@data)), rep("MOMA52", ncol(MOMA52@assays$RNA@data)), rep("MOMA57", ncol(MOMA57@assays$RNA@data)),
           rep("MOMA67", ncol(MOMA67@assays$RNA@data)), rep("MOMA68", ncol(MOMA68@assays$RNA@data)), rep("MOMA72", ncol(MOMA72@assays$RNA@data)),
           rep("MOMA302", ncol(MOMA302@assays$RNA@data)))
length(batch) == ncol(allCounts)
# [1] TRUE

table(batch)
# batch
# MOMA17 MOMA302  MOMA52  MOMA57  MOMA67  MOMA68  MOMA72 
#    466    2747    5553    8359    4583    2546    1231 

## NOTE: Perhaps stratify (as suggested in the kBET GitHub pages)?

# And test
# NOTE: the samples should be in the rows, features in the columns!!!
batch.estimate <- kBET(t(allCounts), batch, plot=FALSE)
# Error: cannot allocate vector of size 4.5 Gb

## Subsample
# put all count data together into one object, features in columns, samples in rows
# But this quickly becomes too big .... so subsample (also suggested by authors)
# "In case of differently sized batches, one should consider stratified sampling in order to keep more samples from smaller batches."

subset_size <- 0.1 #subsample to 10% of the data

subset_id_MOMA17 <- sample.int(n = ncol(MOMA17@assays$RNA@data), size = floor(subset_size * ncol(MOMA17@assays$RNA@data)), replace=FALSE)
subset_id_MOMA52 <- sample.int(n = ncol(MOMA52@assays$RNA@data), size = floor(subset_size * ncol(MOMA52@assays$RNA@data)), replace=FALSE)
subset_id_MOMA57 <- sample.int(n = ncol(MOMA57@assays$RNA@data), size = floor(subset_size * ncol(MOMA57@assays$RNA@data)), replace=FALSE)
subset_id_MOMA67 <- sample.int(n = ncol(MOMA67@assays$RNA@data), size = floor(subset_size * ncol(MOMA67@assays$RNA@data)), replace=FALSE)
subset_id_MOMA68 <- sample.int(n = ncol(MOMA68@assays$RNA@data), size = floor(subset_size * ncol(MOMA68@assays$RNA@data)), replace=FALSE)
subset_id_MOMA72 <- sample.int(n = ncol(MOMA72@assays$RNA@data), size = floor(subset_size * ncol(MOMA72@assays$RNA@data)), replace=FALSE)
subset_id_MOMA302 <- sample.int(n = ncol(MOMA302@assays$RNA@data), size = floor(subset_size * ncol(MOMA302@assays$RNA@data)), replace=FALSE)

allCounts <- merge(MOMA17@assays$RNA@data[,subset_id_MOMA17], MOMA52@assays$RNA@data[,subset_id_MOMA52], by="row.names", all=TRUE)
rownames(allCounts) <- allCounts$Row.names
allCounts$Row.names <- NULL

allCounts <- merge(allCounts, MOMA57@assays$RNA@data[,subset_id_MOMA57], by="row.names", all=TRUE)
rownames(allCounts) <- allCounts$Row.names
allCounts$Row.names <- NULL

allCounts <- merge(allCounts, MOMA67@assays$RNA@data[,subset_id_MOMA67], by="row.names", all=TRUE)
rownames(allCounts) <- allCounts$Row.names
allCounts$Row.names <- NULL

allCounts <- merge(allCounts, MOMA68@assays$RNA@data[,subset_id_MOMA68], by="row.names", all=TRUE)
rownames(allCounts) <- allCounts$Row.names
allCounts$Row.names <- NULL

allCounts <- merge(allCounts, MOMA72@assays$RNA@data[,subset_id_MOMA72], by="row.names", all=TRUE)
rownames(allCounts) <- allCounts$Row.names
allCounts$Row.names <- NULL

allCounts <- merge(allCounts, MOMA302@assays$RNA@data[,subset_id_MOMA302], by="row.names", all=TRUE)
rownames(allCounts) <- allCounts$Row.names
allCounts$Row.names <- NULL

dim(allCounts)
# [1] 22769  2545

# Remove rows that only contain zeroes
allCounts <- allCounts[rowSums(allCounts) > 0, ]

# Make the batch vector, i.e. a concatenation of the number of cells from each sample with the right sample name (in the right order!!)
# I could also add the sample name to the cells earlier and then extract? Less error-prone?
batch <- c(rep("MOMA17", length(subset_id_MOMA17)), rep("MOMA52", length(subset_id_MOMA52)), rep("MOMA57", length(subset_id_MOMA57)),
           rep("MOMA67", length(subset_id_MOMA67)), rep("MOMA68", length(subset_id_MOMA68)), rep("MOMA72", length(subset_id_MOMA72)),
           rep("MOMA302", length(subset_id_MOMA302)))
length(batch) == ncol(allCounts)
# [1] TRUE

table(batch)
# batch
# MOMA17 MOMA302  MOMA52  MOMA57  MOMA67  MOMA68  MOMA72 
#     46     274     555     835     458     254     123 

## NOTE: Perhaps stratify (as suggested in the kBET GitHub pages)?

# And test
# NOTE: the samples should be in the rows, features in the columns!!!
batch.estimate <- kBET(t(allCounts), batch, plot=FALSE)
# Error in svd(x = dataset, nu = dim.comp, nv = 0) : 
#   infinite or missing values in 'x'

# There are rows with only NA ....
m <- allCounts[rowSums(is.na(allCounts)) != ncol(allCounts), ]
dim(m)
#[1] 14131  2545

# And test
# NOTE: the samples should be in the rows, features in the columns!!!
batch.estimate <- kBET(t(m), batch, plot=FALSE)

# And make a plot
plot.data <- data.frame(class=rep(c('observed', 'expected'), 
                                  each=length(batch.estimate$stats$kBET.observed)), 
                        data =  c(batch.estimate$stats$kBET.observed,
                                  batch.estimate$stats$kBET.expected))
g <- ggplot(plot.data, aes(class, data)) + geom_boxplot() + 
  labs(x='Test', y='Rejection rate',title='kBET test results') +
  theme_bw() +  
  scale_y_continuous(limits=c(0,1))

# And shows ~0 for expected and 1 for observed ...
# so, compute the average silhouette width and PCA-based batch-effect measures to explore the degree of the batch effect
pca.data <- prcomp(t(m), center=TRUE) #compute PCA representation of the data
batch.silhouette <- batch_sil(pca.data, batch)
# [1] -0.1205051
batch.pca <- pcRegression(pca.data, batch)

                     
# Or do it on cluster level (after integration??)
seurat <- readRDS("9_kidney_reference_mapping/Mito10.integrated.clustered.annot.kidney_immune_reference.rds")
DefaultAssay(seurat) <- "RNA"

# sample name is in 'orig.ident', we can take cluster annotation at different levels... fx. 'integrated_snn_res.0.6'

subset_id_MOMA17 <- colnames(subset(seurat, orig.ident == "MOMA17"))[sample.int(n = ncol(subset(seurat, orig.ident == "MOMA17")), size = floor(0.5 * ncol(subset(seurat, orig.ident == "MOMA17"))), replace=FALSE)]
subset_id_MOMA52 <- colnames(subset(seurat, orig.ident == "MOMA52"))[sample.int(n = ncol(subset(seurat, orig.ident == "MOMA52")), size = floor(0.05 * ncol(subset(seurat, orig.ident == "MOMA52"))), replace=FALSE)]
subset_id_MOMA57 <- colnames(subset(seurat, orig.ident == "MOMA57"))[sample.int(n = ncol(subset(seurat, orig.ident == "MOMA57")), size = floor(0.025 * ncol(subset(seurat, orig.ident == "MOMA57"))), replace=FALSE)]
subset_id_MOMA67 <- colnames(subset(seurat, orig.ident == "MOMA67"))[sample.int(n = ncol(subset(seurat, orig.ident == "MOMA67")), size = floor(0.05 * ncol(subset(seurat, orig.ident == "MOMA67"))), replace=FALSE)]
subset_id_MOMA68 <- colnames(subset(seurat, orig.ident == "MOMA68"))[sample.int(n = ncol(subset(seurat, orig.ident == "MOMA68")), size = floor(0.1 * ncol(subset(seurat, orig.ident == "MOMA68"))), replace=FALSE)]
subset_id_MOMA72 <- colnames(subset(seurat, orig.ident == "MOMA72"))[sample.int(n = ncol(subset(seurat, orig.ident == "MOMA72")), size = floor(0.2 * ncol(subset(seurat, orig.ident == "MOMA72"))), replace=FALSE)]
subset_id_MOMA302 <- colnames(subset(seurat, orig.ident == "MOMA302"))[sample.int(n = ncol(subset(seurat, orig.ident == "MOMA302")), size = floor(0.1 * ncol(subset(seurat, orig.ident == "MOMA302"))), replace=FALSE)]

allCounts <- seurat@assays$RNA@data[, c(subset_id_MOMA17, subset_id_MOMA52, subset_id_MOMA57, subset_id_MOMA67, subset_id_MOMA68, subset_id_MOMA72, subset_id_MOMA302)]
# Remove rows that only contain zeroes
allCounts <- allCounts[rowSums(allCounts) > 0, ]

# There are rows with only NA ....
m <- as.matrix(t(allCounts[rowSums(is.na(allCounts)) != ncol(allCounts), ]))
dim(m)
# [1] 1721  19663

batch <- c(rep("MOMA17", length(subset_id_MOMA17)), rep("MOMA52", length(subset_id_MOMA52)), rep("MOMA57", length(subset_id_MOMA57)),
           rep("MOMA67", length(subset_id_MOMA67)), rep("MOMA68", length(subset_id_MOMA68)), rep("MOMA72", length(subset_id_MOMA72)),
           rep("MOMA302", length(subset_id_MOMA302)))

length(batch) == nrow(m)
#[1] TRUE

table(batch)
# batch
# MOMA17 MOMA302  MOMA52  MOMA57  MOMA67  MOMA68  MOMA72 
#    233     274     277     208     229     254     246 

# For this more stratified example
batch.estimate <- kBET(m, batch, plot=FALSE)

# And make a plot
plot.data <- data.frame(class=rep(c('observed', 'expected'), 
                                  each=length(batch.estimate$stats$kBET.observed)), 
                        data =  c(batch.estimate$stats$kBET.observed,
                                  batch.estimate$stats$kBET.expected))
g <- ggplot(plot.data, aes(class, data)) + geom_boxplot() + 
  labs(x='Test', y='Rejection rate',title='kBET test results') +
  theme_bw() +  
  scale_y_continuous(limits=c(0,1))
# Still, 0 vs 1 ...

## Go for the cluster approach
# - Do we go back to the original data, i.e. all cells per cluster??
# - Take the clustering at 0.3
#   0.6 gave too many small clusters for which kBET did not work (i.e. gave NA)
clusters <- seurat@meta.data[c(subset_id_MOMA17, subset_id_MOMA52, subset_id_MOMA57, subset_id_MOMA67, subset_id_MOMA68, subset_id_MOMA72, subset_id_MOMA302), "integrated_snn_res.0.3"]
length(batch) == length(clusters)
#[1] TRUE

kBET_result_list <- list()
sum_kBET <- 0
for (cluster_level in unique(clusters)){
  batch_tmp <- batch[clusters == as.numeric(cluster_level)]
  print(paste0("Cluster: ",cluster_level, " has ", length(batch_tmp), " cells"))
  data_tmp <- m[clusters == cluster_level,]
  kBET_tmp <- kBET(df=data_tmp, batch=batch_tmp, plot=FALSE)
  # If the optimal neighbourhood size (k0) is smaller than 10, NA is returned.... damn!!
  if (!is.na(kBET_tmp)){
    print("the optimal neighbourhood size (k0) is larger than 10 ...")
    kBET_result_list[[cluster_level]] <- kBET_tmp
    sum_kBET <- sum_kBET + kBET_tmp$summary$kBET.observed[1]
  }
}
# [1] "Cluster: 6 has 81 cells"
# [1] "Cluster: 0 has 423 cells"
# [1] "the optimal neighbourhood size (k0) is larger than 10 ..."
# [1] "Cluster: 3 has 196 cells"
# [1] "the optimal neighbourhood size (k0) is larger than 10 ..."
# [1] "Cluster: 1 has 353 cells"
# [1] "the optimal neighbourhood size (k0) is larger than 10 ..."
# [1] "Cluster: 8 has 113 cells"
# [1] "the optimal neighbourhood size (k0) is larger than 10 ..."
# [1] "Cluster: 4 has 174 cells"
# [1] "the optimal neighbourhood size (k0) is larger than 10 ..."
# [1] "Cluster: 7 has 68 cells"
# [1] "Cluster: 5 has 71 cells"
# [1] "Cluster: 2 has 184 cells"
# [1] "the optimal neighbourhood size (k0) is larger than 10 ..."
# [1] "Cluster: 9 has 50 cells"
# [1] "Cluster: 11 has 3 cells"
# [1] "Cluster: 10 has 5 cells"

#averaging
mean_kBET = sum_kBET/length(unique(clusters))
mean_kBET
#[1] 0.1783342

## And try with all data ...
allCounts <- seurat@assays$RNA@data

# There are rows with only NA ....
allCounts <- as.matrix(allCounts[rowSums(is.na(allCounts)) != ncol(allCounts), ])
dim(allCounts)
#[1] 22769 25485

# Remove rows that only contain zeroes
allCounts <- allCounts[rowSums(allCounts) > 0, ]
dim(allCounts)
# [1] 22265 25485
allCounts <- t(allCounts)

batch <- seurat@meta.data$orig.ident

length(batch) == nrow(allCounts)
#[1] TRUE

table(batch)
# batch
# MOMA17 MOMA302  MOMA52  MOMA57  MOMA67  MOMA68  MOMA72 
#    466    2747    5553    8359    4583    2546    1231 

clusters <- seurat@meta.data$integrated_snn_res.0.3
length(batch) == length(clusters)
#[1] TRUE

kBET_result_list <- list()
sum_kBET <- 0
for (cluster_level in unique(clusters)){
  batch_tmp <- batch[clusters == as.numeric(cluster_level)]
  print(paste0("Cluster: ",cluster_level, " has ", length(batch_tmp), " cells"))
  data_tmp <- allCounts[clusters == cluster_level,]
  kBET_tmp <- kBET(df=data_tmp, batch=batch_tmp, plot=FALSE)
  # If the optimal neighbourhood size (k0) is smaller than 10, NA is returned.... damn!!
  if (!is.na(kBET_tmp)){
    print("the optimal neighbourhood size (k0) is larger than 10 ...")
    kBET_result_list[[cluster_level]] <- kBET_tmp
    sum_kBET <- sum_kBET + kBET_tmp$summary$kBET.observed[1]
  }
}
# [1] "Cluster: 7 has 1109 cells"
# [1] "the optimal neighbourhood size (k0) is larger than 10 ..."
# [1] "Cluster: 3 has 3083 cells"
# [1] "the optimal neighbourhood size (k0) is larger than 10 ..."
# [1] "Cluster: 2 has 3885 cells"
# [1] "the optimal neighbourhood size (k0) is larger than 10 ..."
# [1] "Cluster: 1 has 5057 cells"
# [1] "the optimal neighbourhood size (k0) is larger than 10 ..."
# [1] "Cluster: 10 has 84 cells"
# [1] "Cluster: 0 has 5756 cells"
# [1] "the optimal neighbourhood size (k0) is larger than 10 ..."
# [1] "Cluster: 5 has 1279 cells"
# [1] "the optimal neighbourhood size (k0) is larger than 10 ..."
# [1] "Cluster: 4 has 2228 cells"
# [1] "the optimal neighbourhood size (k0) is larger than 10 ..."
# [1] "Cluster: 9 has 740 cells"
# [1] "the optimal neighbourhood size (k0) is larger than 10 ..."
# [1] "Cluster: 8 has 1039 cells"
# [1] "the optimal neighbourhood size (k0) is larger than 10 ..."
# [1] "Cluster: 6 has 1162 cells"
# [1] "the optimal neighbourhood size (k0) is larger than 10 ..."
# [1] "Cluster: 11 has 63 cells"

#averaging
mean_kBET = sum_kBET/length(unique(clusters))
mean_kBET
# [1] 0.8065303
# And this is now high (compared to the 0.17 obtained earlier)

# And probably 0.6 also works as we have more cells per cluster ..
clusters <- seurat@meta.data$integrated_snn_res.0.6
length(batch) == length(clusters)
#[1] TRUE

kBET_result_list <- list()
sum_kBET <- 0
for (cluster_level in unique(clusters)){
  batch_tmp <- batch[clusters == as.numeric(cluster_level)]
  print(paste0("Cluster: ",cluster_level, " has ", length(batch_tmp), " cells"))
  data_tmp <- allCounts[clusters == cluster_level,]
  kBET_tmp <- kBET(df=data_tmp, batch=batch_tmp, plot=FALSE)
  # If the optimal neighbourhood size (k0) is smaller than 10, NA is returned.... damn!!
  if (!is.na(kBET_tmp)){
    print("the optimal neighbourhood size (k0) is larger than 10 ...")
    kBET_result_list[[cluster_level]] <- kBET_tmp
    sum_kBET <- sum_kBET + kBET_tmp$summary$kBET.observed[1]
  }
}
# [1] "Cluster: 8 has 1115 cells"
# [1] "the optimal neighbourhood size (k0) is larger than 10 ..."
# [1] "Cluster: 2 has 3078 cells"
# [1] "the optimal neighbourhood size (k0) is larger than 10 ..."
# [1] "Cluster: 0 has 3887 cells"
# [1] "the optimal neighbourhood size (k0) is larger than 10 ..."
# [1] "Cluster: 3 has 3037 cells"
# [1] "the optimal neighbourhood size (k0) is larger than 10 ..."
# [1] "Cluster: 16 has 114 cells"
# [1] "the optimal neighbourhood size (k0) is larger than 10 ..."
# [1] "Cluster: 4 has 2376 cells"
# [1] "the optimal neighbourhood size (k0) is larger than 10 ..."
# [1] "Cluster: 1 has 3166 cells"
# [1] "the optimal neighbourhood size (k0) is larger than 10 ..."
# [1] "Cluster: 6 has 1170 cells"
# [1] "the optimal neighbourhood size (k0) is larger than 10 ..."
# [1] "Cluster: 13 has 648 cells"
# [1] "the optimal neighbourhood size (k0) is larger than 10 ..."
# [1] "Cluster: 10 has 971 cells"
# [1] "the optimal neighbourhood size (k0) is larger than 10 ..."
# [1] "Cluster: 5 has 1281 cells"
# [1] "the optimal neighbourhood size (k0) is larger than 10 ..."
# [1] "Cluster: 17 has 109 cells"
# [1] "the optimal neighbourhood size (k0) is larger than 10 ..."
# [1] "Cluster: 11 has 808 cells"
# [1] "the optimal neighbourhood size (k0) is larger than 10 ..."
# [1] "Cluster: 7 has 1121 cells"
# [1] "the optimal neighbourhood size (k0) is larger than 10 ..."
# [1] "Cluster: 14 has 643 cells"
# [1] "the optimal neighbourhood size (k0) is larger than 10 ..."
# [1] "Cluster: 19 has 61 cells"
# [1] "the optimal neighbourhood size (k0) is larger than 10 ..."
# [1] "Cluster: 18 has 63 cells"
# [1] "Cluster: 12 has 654 cells"
# [1] "the optimal neighbourhood size (k0) is larger than 10 ..."
# [1] "Cluster: 15 has 199 cells"
# [1] "the optimal neighbourhood size (k0) is larger than 10 ..."
# [1] "Cluster: 9 has 984 cells"
# [1] "the optimal neighbourhood size (k0) is larger than 10 ..."
# There were 50 or more warnings (use warnings() to see the first 50)

#averaging
mean_kBET = sum_kBET/length(unique(clusters))
mean_kBET
# [1] 0.6588933
# A little smaller ...

#### Picasso ####
# Just for fun ... https://github.com/pachterlab/picasso
# Yosta is fond of 'tekkels' (dachshund), so perhaps we could put here data in the form of a famous Dutch
# tekkel, i.e. 'Takkie' (by Fien Westendorp) :-)
# Picasso needs the count data and the associated meta data
seurat <- readRDS("7_cellTypist_annotation_after_integration/Mito10.integrated.rds")

# The raw counts are already in "7_cellTypist_annotation_after_integration"
# Now, I also write the associated meta data
DefaultAssay(seurat) <- "integrated"

colsToKeep <- c("orig.ident", "Cluster1", "SCT_snn_res.0.3", "SCT_snn_res.1.8",
                "predicted.labels.celltypist.after")
write.table(as.data.frame(seurat@meta.data[, colsToKeep]), "7_cellTypist_annotation_after_integration/Mito10.integrated.metaData.csv",
            sep = ',', row.names = T, col.names = T, quote = F)

