#### NOTES ####
# - Which final results to use? 14_libra_deg_20230426 or other (see fx. initial aPear trial "Between clusters" )
#   HALLMARK only or selectedSets (HALLMARK, REACTOME, KEGG, BIOCARTA)?
# - Which cut-offs to use?
# - Default aPEAR (i.e. clusterProfiler) or aPEAR with MSigDB results?
# - ....

#### NOTES 20230711 ####
# Code reinstalled monocle3 and SeuratWrappers, but ran into some problems (I probably messed up with the !require() cluase ....)
# Needed a higher version of Seurat (how does this influence everything I did before???)
#     Error in loadNamespace(j <- i[[1L]], c(lib.loc, .libPaths()), versionCheck = vI[[j]]) : 
#     namespace 'Seurat' 4.1.1 is being loaded, but >= 4.2.0 is required
#     Calls: <Anonymous> ... namespaceImportFrom -> asNamespace -> loadNamespace
#     Execution halted
#     ERROR: lazy loading failed for package 'SeuratWrappers'
#     * removing 'D:/R/myRLibs/R-4.1.0/library/SeuratWrappers'
#
# And when I tried to get 'SeuratWrappers' again I got a warning:
# Error: Failed to install 'unknown package' from GitHub:
#   HTTP error 403.
# API rate limit exceeded for 145.117.145.139. (But here's the good news: Authenticated requests get a higher rate limit. Check out the documentation for more details.)
# 
#   Rate limit remaining: 0/60
#   Rate limit reset at: 2023-07-11 07:23:43 UTC
# 
#   To increase your GitHub API rate limit
#   - Use `usethis::create_github_token()` to create a Personal Access Token.
#   - Use `usethis::edit_r_environ()` and add the token as `GITHUB_PAT`.

# In D:\Dropbox\Support\Yosta_Vegting\scRNASeq\Processing\20221220_AllSamples_v4\Code\Sharing\AnalyseDEG_20230317\14_libra_deg
# I found a 'sessionInfo.txt' from an earlier analysis!!! (Saved that one (as a shortcut), 'sessionInfo.15032023.txt')


#### Set default for a good start ####
source("initAnalysis.r")
init()
# The legacy packages maptools, rgdal, and rgeos, underpinning the sp package,
# which was just loaded, will retire in October 2023.
# Please refer to R-spatial evolution reports for details, especially
# https://r-spatial.org/r/2023/05/15/evolution4.html.
# It may be desirable to make the sf package available;
# package maintainers should consider adding sf to Suggests:.
# The sp package is now running under evolution status 2
# (status 2 uses the sf package in place of rgdal)
# Attaching SeuratObject

# Create a directory to hold results and plots
outDir <- "20_trajectory_inference"
dir.create(outDir)


#### Libraries ####
if (!require("SeuratWrappers")){
  if (!require("monocle3")){
    install.packages("devtools")
  }
  devtools::install_github('satijalab/seurat-wrappers')
}

if (!require("monocle3")){
  # https://cole-trapnell-lab.github.io/monocle3/docs/installation/
  BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                         'limma', 'lme4', 'S4Vectors', 'SingleCellExperiment',
                         'SummarizedExperiment', 'batchelor', 'HDF5Array',
                         'terra', 'ggrastr'))
  if (!require("devtools")){
    install.packages("devtools")
  }
  devtools::install_github('cole-trapnell-lab/monocle3')
}

library(SeuratWrappers)
library(monocle3)
library(ggplot2)
library(patchwork)
library(magrittr)

## Follow https://rdrr.io/github/satijalab/seurat-wrappers/f/docs/monocle3.Rmd
## And also see https://satijalab.org/signac/articles/monocle.htmlfor inspiration

#### Original data ####
## Read the data 
seurat <- readRDS("10_Reclustering_MP_and_MC/selectedClusters.subjectIntegrated.reclustered.kidney_immune_reference.rds")
DefaultAssay(seurat) <- "RNA"

pdf(paste0(outDir,"/dimPlot_initial.pdf"))
   DimPlot(object = seurat, group.by = "integrated_snn_res.1")  
dev.off()



#### Trial 1 ####
# Used code from project Tom Groot-Kormelink ()
# so basically as the tutorial, more or less

# Read the data
seurat <- readRDS("10_Reclustering_MP_and_MC/selectedClusters.subjectIntegrated.reclustered.kidney_immune_reference.rds")
DefaultAssay(seurat) <- "RNA"

cds <- SeuratWrappers::as.cell_data_set(seurat)
cds <- cluster_cells(cds, random_seed=42, reduction_method = "UMAP")
# Error in leidenbase::leiden_find_partition(graph_result[["g"]], partition_type = partition_type,  : 
#                                              REAL() can only be applied to a 'numeric', not a 'NULL'
# see https://github.com/cole-trapnell-lab/monocle3/issues/666
cds <- cluster_cells(cds, random_seed=42, reduction_method = "UMAP", cluster_method = 'louvain')
p1 <- plot_cells(cds, show_trajectory_graph = FALSE)
p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE)

pdf(paste0(outDir,"/trial_1.partitioning.pdf"))
  wrap_plots(p1, p2)
dev.off()

rm(seurat)
gc(verbose = FALSE)

seurat.sub <- subset(as.Seurat(cds), monocle3_partitions == 1) # This selects the large cluster/partition
# Error in h(simpleError(msg, call)) : 
#   error in evaluating the argument 'x' in selecting a method for function 'subset': One or more of the assays 
#   you are trying to convert is not in the SingleCellExperiment object
# see https://github.com/satijalab/seurat-wrappers/issues/114
seurat.sub <- subset(as.Seurat(cds, assay = NULL), monocle3_partitions == 1)
cds <- SeuratWrappers::as.cell_data_set(seurat.sub)
cds <- learn_graph(cds)

pdf(paste0(outDir,"/trial_1.learnedGraph.pdf"))
  plot_cells(
    cds,
    label_groups_by_cluster = FALSE,
    label_leaves = FALSE,
    label_branch_points = FALSE,
    label_principal_points = TRUE 
  )
dev.off()

gc(verbose = FALSE)

# Determine the root cell to start from
max.sub <- which.max(unlist(FetchData(seurat.sub, "S100A8")))
max.sub <- colnames(seurat.sub)[max.sub]
# 
cds <- order_cells(cds, root_cells = max.sub)

pdf(paste0(outDir,"/trial_1.orderedCells_maxS100A8.pdf"))
plot_cells(
  cds,
  color_cells_by = "pseudotime",
  label_cell_groups = FALSE,
  label_leaves = FALSE,
  label_branch_points = FALSE
)
dev.off()

#### Trial 2 - more finetuned ####
seurat <- readRDS("10_Reclustering_MP_and_MC/selectedClusters.subjectIntegrated.reclustered.kidney_immune_reference.rds")
DefaultAssay(seurat) <- "RNA"

cds <- as.cell_data_set(seurat)

# A little more finegrained
cds <- cluster_cells(cds, verbose = TRUE, k=5, partition_qval = 0.005, 
                     random_seed = 42, cluster_method = 'louvain')
p1 <- plot_cells(cds, show_trajectory_graph = FALSE)
p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE)
pdf(paste0(outDir,"/trial_2.partitioning.pdf"))
  wrap_plots(p1, p2)
dev.off()

rm(seurat)
gc(verbose = FALSE)

# Hmmm, which clusters to select?? 
# As we already have a subset of clusters ///
#selected.sub <- subset(as.Seurat(cds, assay = NULL), monocle3_clusters %in% c(3,4,)) 
selected.sub <- subset(as.Seurat(cds, assay = NULL), integrated_snn_res.1 %in% c(0,1,4,8, 9, 11)) 
selected.sub
# An object of class Seurat 
# 22769 features across 2133 samples within 1 assay 
# Active assay: RNA (22769 features, 0 variable features)
# 2 dimensional reductions calculated: UMAP, REF.UMAP
#
# It is essential to set the DefaultAssay to 'RNA', otherwise you might only get the 'integrated' features etc.!!

cds.sub <- as.cell_data_set(selected.sub)
cds.sub <- learn_graph(cds.sub, use_partition = FALSE, close_loop = TRUE,
                       learn_graph_control = NULL)
# if 'use_partitions' has its default, i.e. 'FALSE' you will get an error ...
#     Error in if (kmean_res$ifault != 0) { : argument is of length zero

pdf(paste0(outDir,"/trial_2.learnedGraph.pdf"))

plot_cells(
  cds.sub,
  label_groups_by_cluster = FALSE,
  label_leaves = FALSE,
  label_branch_points = FALSE,
  label_principal_points = TRUE
)
dev.off()

# Hmm, what do I actually see/do here??? Very hard to interpret!!
# - Fx, is 'Y_138' a starting node?
# - I only see a connection between the upper left part and the lower left part at 'Y_1173' 
#   (so to the left of these clusters interafces), I do not see a connection betweeb 'Y_992'
#   and 'Y_1082'
# - The path seems to go from 'Y_969' at the right blob of cells, via 'Y_1251' upwards through 
#   the left clusters of cells and then down to the lower left cluster...

gc(verbose = FALSE)

# Determine the root cell to start from
# - does it start at 'Y_138' or somewhere in the mmiddle of the right cluster?
cds.sub <- order_cells(cds.sub, root_pr_nodes = "Y_138", reduction_method = "UMAP", verbose = TRUE)

pdf(paste0(outDir,"/trial_2.orderedCells_node.pdf"))
plot_cells(
  cds.sub,
  color_cells_by = "pseudotime",
  label_cell_groups = FALSE,
  label_leaves = FALSE,
  label_branch_points = FALSE
)
dev.off()

# # Set the assay back as "integrated"
# integrated.sub <- as.Seurat(cds, assay = "integrated")
# FeaturePlot(integrated.sub, "monocle3_pseudotime")

# NOTE: Something wrong with graphing as Seurat simply uses nCount as size factor
# https://github.com/satijalab/seurat/issues/1658


#### Trial 3 - Clusters 0, 1, 4, 8, 9, 11 ####
library(igraph)  # neede for 'get.edgelist'

seurat <- readRDS("10_Reclustering_MP_and_MC/selectedClusters.subjectIntegrated.reclustered.kidney_immune_reference.rds")
seurat
# An object of class Seurat 
# 41653 features across 4194 samples within 5 assays 
# Active assay: integrated (3000 features, 3000 variable features)
# 4 other assays present: RNA, ADT, SCT, prediction.score.celltype
# 4 dimensional reductions calculated: pca, umap, ref.pca, ref.umap

DefaultAssay(seurat) <- "RNA"
# An object of class Seurat 
# 41653 features across 4194 samples within 5 assays 
# Active assay: RNA (22769 features, 0 variable features)
# 4 other assays present: ADT, SCT, integrated, prediction.score.celltype
# 4 dimensional reductions calculated: pca, umap, ref.pca, ref.umap

DimPlot(object = seurat, group.by = "integrated_snn_res.1", label = TRUE)  

selClusters_seurat <- seurat[, seurat$integrated_snn_res.1 %in% c(0,1,4,8,9,11)]
selClusters_seurat
# An object of class Seurat 
# 41653 features across 2356 samples within 5 assays 
# Active assay: RNA (22769 features, 0 variable features)
# 4 other assays present: ADT, SCT, integrated, prediction.score.celltype
# 4 dimensional reductions calculated: pca, umap, ref.pca, ref.umap

selClusters.cds <- as.cell_data_set(selClusters_seurat, assay="RNA")
selClusters.cds
# class: cell_data_set 
# dim: 22769 2356 
# metadata(0):
#   assays(2): counts logcounts
# rownames(22769): AL627309.5 LINC01409 ... GGT5 AC244090.3
# rowData names(0):
#   colnames(2356): AACACGTCAAAGGCGT-1_1 AACCATGAGCCGATTT-1_1 ... TTTGTTGGTAGAGACC-1_7 TTTGTTGTCTGCTTTA-1_7
# colData names(344): orig.ident nCount_RNA ... ident Size_Factor
# reducedDimNames(2): UMAP REF.UMAP
# mainExpName: RNA
# altExpNames(0):
selClusters.cds <- cluster_cells(cds = selClusters.cds, reduction_method = "UMAP", 
                                 random_seed = 42, cluster_method = "louvain")
selClusters.cds <- learn_graph(selClusters.cds, use_partition = TRUE, close_loop = FALSE)

# And plot
pdf(paste0(outDir,"/trial_3.learnedGraph_noClosedLoop.pdf"))
  plot_cells(
    selClusters.cds,
    label_groups_by_cluster = FALSE,
    label_leaves = FALSE,
    label_branch_points = FALSE,
    label_principal_points = TRUE
  )
dev.off()

# order cells
selClusters.cds <- order_cells(selClusters.cds, root_pr_nodes = "Y_384", reduction_method = "UMAP")

# plot trajectories colored by pseudotime
pdf(paste0(outDir,"/trial_3.orderedCells_noClosedLoop.pdf"))
  plot_cells(
    cds = selClusters.cds,
    color_cells_by = "pseudotime",
    show_trajectory_graph = TRUE
  )
dev.off()


# Extract the pseudotime and add to the seurat object
seurat <- AddMetaData(
  object = seurat,
  metadata = selClusters.cds@principal_graph_aux@listData$UMAP$pseudotime,
  col.name = "Selected_Monocle3"
)

# And visualise
pdf(paste0(outDir,"/trial_3.featurePlot.pdf"))
  FeaturePlot(seurat, c("Selected_Monocle3"), pt.size = 0.4) & scale_color_viridis_c()
dev.off()

# How to also plot the trajectory on the FeaturePlot ??
# Trial code ....

# all coordinates are perhaps in selClusters.cds@principal_graph_aux@listData$UMAP$dp_mst
Ycoords.df <- data.frame(x=selClusters.cds@principal_graph_aux@listData$UMAP$dp_mst["UMAP_1",],
                         y=selClusters.cds@principal_graph_aux@listData$UMAP$dp_mst["UMAP_2",])
rownames(Ycoords.df) <- colnames(selClusters.cds@principal_graph_aux@listData$UMAP$dp_mst)

# And now translate the edges list in selClusters.cds@principal_graph$UMAP
Yedges.df <- as.data.frame(get.edgelist(selClusters.cds@principal_graph$UMAP))
colnames(Yedges.df) <- c("from","to")

from.df <- as.data.frame(Ycoords.df[Yedges.df$from,])
colnames(from.df) <- c("x.from","y.from")

to.df <- as.data.frame(Ycoords.df[Yedges.df$to,])
colnames(to.df) <- c("x.to","y.to")

from_to.df <- cbind(from.df,to.df)

pdf(paste0(outDir,"/trial_3.featurePlot.v2.pdf"))
  p <- FeaturePlot(seurat, c("Selected_Monocle3"), pt.size = 0.4) & scale_color_viridis_c()
  p <- p + geom_segment(data=from_to.df,aes(x=x.from, xend = x.to, y=y.from, yend = y.to),
                        lwd=1,colour="blue") +
    geom_point(data=Ycoords.df, aes(x=x, y=y), colour = "red")

  p <- p + geom_label_repel(data = Ycoords.df, mapping=aes(x, y, label=rownames(Ycoords.df)),
                            max.overlaps = 10, label.size = NA, 
                            size = 3,
                            alpha = 0.6, 
                            force = 1,
                            label.padding=.1,
                            min.segment.length = 3,
                            fontface = 'bold', color = "black",
                            #box.padding = 0.80, point.padding = 0.5,
                            na.rm=TRUE,
                            seed = 42) +
    geom_label_repel(data = Ycoords.df, mapping=aes(x, y, label=rownames(Ycoords.df)),
                     max.overlaps = 10, label.size = NA, 
                     size = 3,
                     alpha = 1, 
                     force = 1,
                     label.padding=.1, 
                     min.segment.length = 3,
                     fontface = 'bold', color = "black",
                     na.rm=TRUE, fill = NA,
                     seed = 42)
  print(p)
dev.off()


#### Trial 4 - no clustering by Monocle, use Seurat clusters ####
library(igraph)  # needed for 'get.edgelist'
library(ggrepel)

# Read the data
seurat <- readRDS("10_Reclustering_MP_and_MC/selectedClusters.subjectIntegrated.reclustered.kidney_immune_reference.rds")
DefaultAssay(seurat) <- "RNA"

cds <- SeuratWrappers::as.cell_data_set(seurat)
#Error in validityMethod(as(object, superClass)) : 
#  object 'CsparseMatrix_validate' not found

# I did not have the problems earlier, but I seem to have been messing around with package versions
# and actually have broken my installation ...
# see https://github.com/satijalab/seurat/issues/6715
# install.packages("Matrix", repos="http://R-Forge.R-project.org")

# The sizeFactors are wrong ... as observed by Utkarsh & Perry
# Workaround observed from debugging session.
# A github issue has been already raised along with solution
# https://github.com/satijalab/seurat-wrappers/issues/97
n.count <- colSums(seurat@assays$RNA@counts)
cds@colData@listData$Size_Factor <- n.count / exp(mean(log(n.count)))

# Code from Utkarsh
# Create the placeholder for holding the cluster and partition information
cds <- cluster_cells(cds)
# plot_cells(cds, show_trajectory_graph = F, group_label_size = 4)

# migrate the cluster info from the original seurat object to the monocle3 object
list_cluster <- seurat@meta.data[[sprintf("integrated_snn_res.1")]]
names(list_cluster) <- seurat@assays[["RNA"]]@data@Dimnames[[2]]
cds@clusters@listData[["UMAP"]][["clusters"]] <- list_cluster

## IMPORTANT ##
# Refer https://github.com/cole-trapnell-lab/monocle3/issues/419 for below line
# since the cds doesn't have gene_short_name column it's unable to plot and
# throws error every time you use plot_cells and want genenames ....
rowData(cds)$gene_short_name <- row.names(rowData(cds))

p1 <- plot_cells(cds, show_trajectory_graph = FALSE)
p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE)
pdf(paste0(outDir,"/trial_4.partitioning_originalClusters.pdf"))
  wrap_plots(p1, p2)
dev.off()

## Clear the memory
gc(verbose = FALSE)

# Learn the graph
cds <- learn_graph(cds)

pdf(paste0(outDir,"/trial_4.learnedGraph.pdf"))
plot_cells(
  cds,
  label_groups_by_cluster = FALSE,
  label_leaves = FALSE,
  label_branch_points = FALSE,
  label_principal_points = TRUE,
  label_roots = TRUE
)
dev.off()

## OBSERVATION
# - Much 'cleaner' path then before, less branches..

# Determine the root cell to start from by max expression of a cell
max.sub <- which.max(unlist(FetchData(seurat, "S100A8")))
max.sub <- colnames(seurat)[max.sub]
# 
cds.maxExpr <- order_cells(cds, root_cells = max.sub)

pdf(paste0(outDir,"/trial_4.orderedCells_maxS100A8.pdf"))
plot_cells(
  cds.maxExpr,
  color_cells_by = "pseudotime",
  label_cell_groups = FALSE,
  label_leaves = FALSE,
  label_branch_points = FALSE,
  label_roots = TRUE
)
dev.off()

# Clear memory
gc(verbose = FALSE)

# Determine the root cell to start from
# - does it start at 'Y_183' 
# - or somewhere in the middle/top of the right cluster? i.e. Y_1?
cds.Node <- order_cells(cds, root_pr_nodes = "Y_183", reduction_method = "UMAP", verbose = TRUE)

pdf(paste0(outDir,"/trial_4.orderedCells_nodeY_183.pdf"))
plot_cells(
  cds.Node,
  color_cells_by = "pseudotime",
  label_cell_groups = FALSE,
  label_leaves = FALSE,
  label_branch_points = FALSE,
  label_roots = TRUE
)
dev.off()

# Before subsetting, put in the necessary information in order to be able to find the gene modules later on
# see below 'find_gene_modules'
# Transfer projected feature loadings, I use umap ... but will transfer both 'pca' and 'umap' information
seurat <- ProjectDim(object = seurat, reduction = "pca")
seurat <- ProjectDim(object = seurat, reduction = "umap")
cds@reduce_dim_aux[['PCA']][['model']][['svd_v']] <- seurat@reductions[["pca"]]@feature.loadings.projected
cds@reduce_dim_aux[['UMAP']][['model']][['svd_v']] <- seurat@reductions[["umap"]]@feature.loadings.projected
# However, this slot is empty, '<0 x 0 matrix>' if you have not used ProjectDim(object = seurat, reduction = "pca") !!
# Both for 'umap' and 'pca' ...

# Transfer PCs stdevs
cds@reduce_dim_aux[['PCA']][['model']][['svd_sdev']] <- seurat@reductions$pca@stdev

# Transfer PCs cell embeddings (not necessary if UMAP transferred from Seurat used)
reducedDim(cds, type = "PCA") <-   seurat@reductions$pca@cell.embeddings

# Subset just clusters 0,1,4,8,9,11
clusters_cds <- cds[, colData(cds)$integrated_snn_res.1 %in% c("0","1","4","8","9","11")]
plot_cells(clusters_cds, color_cells_by="integrated_snn_res.1")

pr_graph_test_res <- graph_test(clusters_cds, neighbor_graph="knn", cores=8)
# Warning message:
#   In pbmcapply::pbmclapply(row.names(exprs_mat), FUN = function(x,  :
#     mc.cores > 1 is not supported on Windows due to limitation of mc*apply() functions.
#     mc.core is set to 1.
range(pr_graph_test_res$morans_I, na.rm = T)
# [1] -0.01284827  0.67941084

pr_graph_test_res_qVal <- subset(pr_graph_test_res[order(pr_graph_test_res$morans_I, decreasing = TRUE),], q_value < 0.05)
pr_deg_ids_top25 <- rownames(pr_graph_test_res_qVal)[1:25]
pr_deg_ids_top25
# [1] "C1QC"     "RGS1"     "C1QB"     "S100A8"   "C3"       "S100A9"   "VCAN"     "C1QA"     "FCN1"     "LST1"     "FCGR3A"   "S100A4"  
# [13] "S100A12"  "SMIM25"   "A2M"      "RPS19"    "S100A6"   "CD9"      "ALOX5AP"  "FABP5"    "MS4A6A"   "TREM2"    "LYZ"      "C15orf48"
# [25] "APOE"

## OBSERVATION:
# - Interestingly enough, C1QC, C1QB, C1QA, TREM2 etc pop up!!
#   Confirmation of what we already knew from our Suerat analysis

# Plot some of these genes
pdf(paste0(outDir,"/trial_4.selectedGenes.pdf"))
plot_cells(clusters_cds, 
           genes=c("C1QC", "S100A8","TREM2", "FCGR3A"),
           show_trajectory_graph=TRUE,
           trajectory_graph_color = "red",
           trajectory_graph_segment_size = 0.1,
           label_cell_groups=FALSE,
           label_branch_points = FALSE,
           label_leaves=FALSE,
           label_roots = FALSE,
           cell_size = 0.5,
           scale_to_range = TRUE)   #i.e. the expressions are scaled to percent of maximum expression
                                    # there seem only tobe a cuple of cells that have a high expression and they determine the color scale ?
dev.off()

# # Error in plot_cells(clusters_cds, genes = c("C1QC", "S100A8", "TREM2",  : 
# #         None of the provided genes were found in the cds
# # Huh ?!
# c("C1QC", "S100A8","TREM2", "FCGR3A") %in% rownames(clusters_cds)
# [1] TRUE TRUE TRUE TRUE
#
## Solved by including the code (from Utkarsh) at line 398 !!!

## Find modules of co-regulated genes
#  This essentially runs UMAP on the genes (as opposed to the cells) and then groups 
#  them into modules using Louvain community analysis
gene_module_df <- find_gene_modules(clusters_cds[rownames(pr_graph_test_res_qVal),], resolution=1e-2)

# It is essential to run 'ProjectDim' in order to have the reduced_dim_aux, see line 469-472
# If you havenot done this properly, you will get the errors below:
#
#   Error in cds@reduce_dim_aux[[preprocess_method]][["model"]]$svd_v %*%  : 
#     requires numeric/complex matrix/vector arguments
#
# Probably, bacause I went from a Seurat object to Monocle3 and then did a subset ...
# there now seems to be no clusters_cds@reduce_dim_aux[["UMAP"]][["model"]]$svd_v  slot
#
# see: https://github.com/cole-trapnell-lab/monocle3/issues/623  -> NOT A SOLVED ISSUE YET
# - do this before you subset to the clusters ....
# - However, it seems as if I haven't taken this data along in the Seurat object to begin with 
#   and probably can not recalculate it, as it would probably change all 
# - I do seem to have the loadings for the PCA ...
# - As explained in: https://github.com/satijalab/seurat/issues/6979
#   UMAP is non-linear multidimensional reduction therefore it works on distances calculated from the PCs if 
#   you need feature embeddings you will have to use lower dimensional equivalents like the PCA.
#
# Also, this will not work:
#   gene_module_df <- find_gene_modules(clusters_cds[rownames(pr_graph_test_res_qVal),], resolution=1e-2, preprocess_method = "PCA")
# 
# Error in cds@reduce_dim_aux[[preprocess_method]][["model"]]$svd_v %*%  : 
#   non-conformable arguments

cell_cluster_df <- tibble::tibble(cell=row.names(colData(clusters_cds)), 
                                cluster_seurat=colData(cds)[colnames(clusters_cds), "integrated_snn_res.1"])
agg_mat <- aggregate_gene_expression(clusters_cds, gene_module_df, cell_cluster_df)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
colnames(agg_mat) <- stringr::str_c("Cluster ", colnames(agg_mat))

pdf(paste0(outDir,"/trial_4.geneModules_selectedClusters.pdf"))
  pheatmap::pheatmap(agg_mat, cluster_rows=TRUE, cluster_cols=TRUE,
                     scale="column", clustering_method="ward.D2",
                     fontsize=6)
dev.off()

## Plot some genes over time
interesting_genes <- c("C1QC", "S100A8","TREM2", "FCGR3A")
genes_cds <- clusters_cds[rowData(clusters_cds)$gene_short_name %in% interesting_genes,]
genes_cds <- order_cells(genes_cds, root_pr_nodes = "Y_183", reduction_method = "UMAP", verbose = TRUE)

p <- plot_genes_in_pseudotime(genes_cds,
                         color_cells_by="integrated_snn_res.1",
                         min_expr=0.5)
p <- p + guides(colour = guide_legend(override.aes = list(size = 5)))

pdf(paste0(outDir,"/trial_4.selectedGenes_over_Time.pdf"))
  print(p)
dev.off()

# Save the modules of co-regulated genes
write.table(as.data.frame(gene_module_df[order(gene_module_df$module),]), 
            file=paste0(outDir,"/trial4_gene_modules.txt"), sep="\t", quote=FALSE,
            row.names = FALSE, col.names = TRUE)


#### Trajectory without T cells - version 1 ####
# 'Remove' by putting these clusters in a different partition
library(igraph)  # needed for 'get.edgelist'
library(ggrepel)

# Read the data
seurat <- readRDS("10_Reclustering_MP_and_MC/selectedClusters.subjectIntegrated.reclustered.kidney_immune_reference.rds")
DefaultAssay(seurat) <- "RNA"

# Convert to Monocle ready object
cds <- SeuratWrappers::as.cell_data_set(seurat)

# The sizeFactors are wrong ... as observed by Utkarsh & Perry
# Workaround observed from debugging session.
# A github issue has been already raised along with solution
# https://github.com/satijalab/seurat-wrappers/issues/97
n.count <- colSums(seurat@assays$RNA@counts)
cds@colData@listData$Size_Factor <- n.count / exp(mean(log(n.count)))

# Code from Utkarsh
# Create the placeholder for holding the cluster and partition information
cds <- cluster_cells(cds)

# migrate the cluster info from the original seurat object to the monocle3 object
list_cluster <- seurat@meta.data[[sprintf("integrated_snn_res.1")]]
names(list_cluster) <- seurat@assays[["RNA"]]@data@Dimnames[[2]]
cds@clusters@listData[["UMAP"]][["clusters"]] <- list_cluster
# Alter the partition, i.e. all T cell clusters in partition = 2
partitionVec <- ifelse(cds@clusters@listData[["UMAP"]][["clusters"]] %in% c(0,1,4,8,9,11), 1, 2)
names(partitionVec) <- names(list_cluster)
cds@clusters@listData[["UMAP"]][["partition"]] <- factor(partitionVec)
cds@clusters@listData[["UMAP"]][["partitions"]] <- factor(partitionVec)

# And check
table(cds@clusters@listData[["UMAP"]][["partition"]], cds@clusters@listData[["UMAP"]][["clusters"]])
#       0   1   2   3   4   5   6   7   8   9  10  11  12
#   1 799 783   0   0 364   0   0   0 187 162   0  61   0
#   2   0   0 494 381   0 304 273 246   0   0  82   0  58

## IMPORTANT ##
# Refer https://github.com/cole-trapnell-lab/monocle3/issues/419 for below line
# since the cds doesn't have gene_short_name column it's unable to plot and
# throws error every time you use plot_cells and want genenames ....
rowData(cds)$gene_short_name <- row.names(rowData(cds))

p1 <- plot_cells(cds, show_trajectory_graph = FALSE)
p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE, reduction_method = "UMAP")
pdf(paste0(outDir,"/without_Tcells.partitioning_originalClusters.pdf"))
wrap_plots(p1, p2)
dev.off()

## Clear the memory
gc(verbose = FALSE)

# Learn the graph
# It will learn graphs for each partition
cds <- learn_graph(cds, use_partition = TRUE)

pdf(paste0(outDir,"/without_Tcells.learnedGraph.pdf"))
plot_cells(
  cds,
  label_groups_by_cluster = FALSE,
  label_leaves = FALSE,
  label_branch_points = FALSE,
  label_principal_points = TRUE,
  label_roots = TRUE
)
dev.off()

## OBSERVATION
# - Very interesting!! There still are some cells belonging to cluster 12 that are in the left cluster
#   and Monocle tries to get these into the graph learned for partition 2 ..

## Conclusion:
# This does not work ...

#### Trajectory without T cells - version 2 ####
# 'Remove' T cell clusters by subsetting after conversion to Monocle object
library(igraph)  # needed for 'get.edgelist'
library(ggrepel)

# Read the data
seurat <- readRDS("10_Reclustering_MP_and_MC/selectedClusters.subjectIntegrated.reclustered.kidney_immune_reference.rds")
DefaultAssay(seurat) <- "RNA"

# Convert to Monocle ready object
cds <- SeuratWrappers::as.cell_data_set(seurat)

# The sizeFactors are wrong ... as observed by Utkarsh & Perry
# Workaround observed from debugging session.
# A github issue has been already raised along with solution
# https://github.com/satijalab/seurat-wrappers/issues/97
n.count <- colSums(seurat@assays$RNA@counts)
cds@colData@listData$Size_Factor <- n.count / exp(mean(log(n.count)))

# Code from Utkarsh
# Create the placeholder for holding the cluster and partition information
cds <- cluster_cells(cds)

# migrate the cluster info from the original seurat object to the monocle3 object
list_cluster <- seurat@meta.data[[sprintf("integrated_snn_res.1")]]
names(list_cluster) <- seurat@assays[["RNA"]]@data@Dimnames[[2]]
cds@clusters@listData[["UMAP"]][["clusters"]] <- list_cluster

# And subset to get rid of the T cell clusters (provided by Yosta, get rid of
#   "De voornaamste T cell clusters zijn: 2, 5,6, 10"
# We keep the DC's: 
#   "Aangezien DCs mogelijk ook naar deze LAMs differentieren, laten we die er maar in, dat zijn cluster 3,7,12."
dim(cds)
# [1] 22769  4194
table(colData(cds)$integrated_snn_res.1)
#   0   1   2   3   4   5   6   7   8   9  10  11  12 
# 799 783 494 381 364 304 273 246 187 162  82  61  58

cds <- cds[, colData(cds)$integrated_snn_res.1 %in% c(0,1,3,4,7,8,9,11,12)]
dim(cds)
# [1] 22769  3041

# However, the cluster slot is not subsetted
length(cds@clusters@listData[["UMAP"]][["clusters"]])
# [1] 4194

# So, have to redo the cluster annotation
# Generate the structure
cds <- cluster_cells(cds)

# migrate the cluster info from the original seurat object to the monocle3 object
cds@clusters@listData[["UMAP"]][["clusters"]] <- colData(cds)$integrated_snn_res.1
cds@clusters@listData[["UMAP"]][["clusters"]] <- factor(cds@clusters@listData[["UMAP"]][["clusters"]])

## IMPORTANT ##
# Refer https://github.com/cole-trapnell-lab/monocle3/issues/419 for below line
# since the cds doesn't have gene_short_name column it's unable to plot and
# throws error every time you use plot_cells and want genenames ....
rowData(cds)$gene_short_name <- row.names(rowData(cds))

p1 <- plot_cells(cds, show_trajectory_graph = FALSE)
p2 <- plot_cells(cds, color_cells_by = "cluster", show_trajectory_graph = FALSE, reduction_method = "UMAP")
pdf(paste0(outDir,"/without_Tcells.v2.partitioning_originalClusters.pdf"))
wrap_plots(p1, p2)
dev.off()
# Error: Must request at least one colour from a hue palette.

#### Trajectory without T cells - version 3 ####
# 'Remove' T cell clusters before conversion to Monocle object
library(igraph)  # needed for 'get.edgelist'
library(ggrepel)

seurat <- readRDS("10_Reclustering_MP_and_MC/selectedClusters.subjectIntegrated.reclustered.kidney_immune_reference.rds")
DefaultAssay(seurat) <- "RNA"

seurat <- subset(seurat, subset = integrated_snn_res.1 %in% c(0,1,3,4,7,8,9,11,12))

# Convert to Monocle ready object
cds <- SeuratWrappers::as.cell_data_set(seurat)

# The sizeFactors are wrong ... as observed by Utkarsh & Perry
# Workaround observed from debugging session.
# A github issue has been already raised along with solution
# https://github.com/satijalab/seurat-wrappers/issues/97
n.count <- colSums(seurat@assays$RNA@counts)
cds@colData@listData$Size_Factor <- n.count / exp(mean(log(n.count)))

# Code from Utkarsh
# Create the placeholder for holding the cluster and partition information
cds <- cluster_cells(cds)

# migrate the cluster info from the original seurat object to the monocle3 object
list_cluster <- seurat@meta.data[[sprintf("integrated_snn_res.1")]]
names(list_cluster) <- seurat@assays[["RNA"]]@data@Dimnames[[2]]
cds@clusters@listData[["UMAP"]][["clusters"]] <- list_cluster

## IMPORTANT ##
# Refer https://github.com/cole-trapnell-lab/monocle3/issues/419 for below line
# since the cds doesn't have gene_short_name column it's unable to plot and
# throws error every time you use plot_cells and want genenames ....
rowData(cds)$gene_short_name <- row.names(rowData(cds))

p1 <- plot_cells(cds, show_trajectory_graph = FALSE)
p2 <- plot_cells(cds, color_cells_by = "cluster", show_trajectory_graph = FALSE, reduction_method = "UMAP")
pdf(paste0(outDir,"/without_Tcells.v3.partitioning_originalClusters.pdf"))
wrap_plots(p1, p2)
dev.off()

# And this works ... !!

## NOTE
# In the code above I used subsetting as well (subset(as.Seurat)), did this now work??? CHECK?

# See some spurious cells, probably going to mess up the trajectory?
# Should we remove these beforehand?

## Clear the memory
gc(verbose = FALSE)

# Learn the graph
# It will learn graphs for each partition
cds <- learn_graph(cds)

pdf(paste0(outDir,"/without_Tcells.v3.learnedGraph.pdf"))
plot_cells(
  cds,
  label_groups_by_cluster = FALSE,
  label_leaves = FALSE,
  label_branch_points = FALSE,
  label_principal_points = TRUE,
  label_roots = TRUE,
  cell_size = 0.6           # increased to better see the different colors of the clusters
)
dev.off()

## OBSERVATION
# This leads to a complete different trajectory, it seems...

### Start node/cell
# also see http://cole-trapnell-lab.github.io/monocle-release/monocle3/
# They have a function for detecting the correct root node: get_correct_root_state

# Determine the root cell to start from by max expression of a cell
# Mail Yosta 20230718:
# "Als we een startpunt moeten kiezen zou je de cel met de meeste CD14 expressie kunnen nemen (dat zal in cluster 1 zijn inderdaad) , 
#  want dat is dan de meest 'classical' monocyte like cell en in principe obv literatuur gaan classicals over naar non-classicals
max.sub <- which.max(unlist(FetchData(seurat, "CD14")))  # This uses the counts
max.sub <- colnames(seurat)[max.sub]
max.sub
# [1] "TACTGCCAGCATGTTC-1_3"
cds.maxExpr <- order_cells(cds, root_cells = max.sub)

pdf(paste0(outDir,"/without_Tcells.v3.orderedCells_maxCD14.pdf"))
plot_cells(
  cds.maxExpr,
  color_cells_by = "pseudotime",
  trajectory_graph_color = "red", trajectory_graph_segment_size = 0.9,
  label_cell_groups = FALSE,
  label_leaves = FALSE,
  label_branch_points = FALSE,
  label_roots = TRUE,
  cell_size = 0.6
)
dev.off()

# NOTE
# The cell with the maximum value of CD14 does not lie in cluster 1 ... :-)

cds_cd14 <- cds[rowData(cds)$gene_short_name %in% c("CD14"), cds@clusters$UMAP$clusters %in% c(1)]
max.sub <- colnames(cds_cd14)[which.max(as.vector(counts(cds_cd14)))]
cds.maxExpr <- order_cells(cds, root_cells = max.sub)

# # Just to check, whatt does this code find if you leave out the constraint of membership of cluster1
# cds_cd14 <- cds[rowData(cds)$gene_short_name %in% c("CD14"), ]
# max.sub <- colnames(cds_cd14)[which.max(as.vector(counts(cds_cd14)))]
# max.sub
# [1] "TACTGCCAGCATGTTC-1_3"  - same as the 'seurat FetchData' code, line 812!!)

pdf(paste0(outDir,"/without_Tcells.v3.orderedCells_maxCD14_cluster1.pdf"))
plot_cells(
  cds.maxExpr,
  color_cells_by = "pseudotime",
  trajectory_graph_color = "red", trajectory_graph_segment_size = 0.9,
  label_cell_groups = FALSE,
  label_leaves = FALSE,
  label_branch_points = FALSE,
  label_roots = TRUE,
  cell_size = 0.6
)
dev.off()


# Plot some of the interesting genes
pdf(paste0(outDir,"/without_Tcells.v3.selectedGenes.pdf"))
plot_cells(cds.maxExpr, 
           genes=c("CD14","FCGR3A", "CD36", "HLA-DRB5"),
           show_trajectory_graph=TRUE,
           trajectory_graph_color = "red",
           trajectory_graph_segment_size = 0.1,
           label_cell_groups=FALSE,
           label_branch_points = FALSE,
           label_leaves=FALSE,
           label_roots = FALSE,
           cell_size = 0.5,
           scale_to_range = TRUE)   #i.e. the expressions are scaled to percent of maximum expression
                                    # there seem only to be a couple of cells that have a high expression and they determine the color scale ?
plot_cells(cds.maxExpr, 
           genes=c("C1QC", "S100A8","TREM2", "CD1C" ),
           show_trajectory_graph=TRUE,
           trajectory_graph_color = "red",
           trajectory_graph_segment_size = 0.1,
           label_cell_groups=FALSE,
           label_branch_points = FALSE,
           label_leaves=FALSE,
           label_roots = FALSE,
           cell_size = 0.5,
           scale_to_range = TRUE)
dev.off()

## Plot some genes over time
interesting_genes <- c("CD14","FCGR3A", "CD36", "HLA-DRB5", "C1QC", "S100A8","TREM2", "CD1C")
genes_cds <- cds.maxExpr[rowData(cds.maxExpr)$gene_short_name %in% interesting_genes,]
genes_cds <- order_cells(genes_cds, root_cells = max.sub, reduction_method = "UMAP", verbose = TRUE)
genes_cds <- genes_cds[interesting_genes,]

p <- plot_genes_in_pseudotime(genes_cds[1:4,],
                              color_cells_by="integrated_snn_res.1",
                              min_expr=0.5)
p <- p + guides(colour = guide_legend(override.aes = list(size = 5)))

p2 <- plot_genes_in_pseudotime(genes_cds[5:8,],
                              color_cells_by="integrated_snn_res.1",
                              min_expr=0.5)
p2 <- p2 + guides(colour = guide_legend(override.aes = list(size = 5)))


pdf(paste0(outDir,"/without_Tcells.v3.selectedGenes_over_Time.pdf"))
  print(p)
  print(p2)
dev.off()

#### Trajectory without T cells - version 3 :Stability test trajectory ####
# iterating over different sub-samples of cells to understand the stability of
# trajectory inferred by monocle3 algorithm

# Gleaned code from Utkarsh Mahamune
# D:\surfdrive2\Shared\UtkarshMahamune\20220815_GRN\Processing\2_Trajectory_inference\Code\Monocle3_traj_inf_stability.R

# 'Remove' T cell clusters before conversion to Monocle object
library(igraph)  # needed for 'get.edgelist'
library(ggrepel)
library(cowplot)

seurat <- readRDS("10_Reclustering_MP_and_MC/selectedClusters.subjectIntegrated.reclustered.kidney_immune_reference.rds")
DefaultAssay(seurat) <- "RNA"

seurat <- subset(seurat, subset = integrated_snn_res.1 %in% c(0,1,3,4,7,8,9,11,12))
seurat@meta.data$integrated_snn_res.1 <- factor(seurat@meta.data$integrated_snn_res.1)

# Convert to Monocle ready object
cds <- SeuratWrappers::as.cell_data_set(seurat)
rowData(cds)$gene_short_name <- row.names(rowData(cds))
# migrate the cluster info from the original seurat object to the monocle3 object
list_cluster <- seurat@meta.data[[sprintf("integrated_snn_res.1")]]
names(list_cluster) <- seurat@assays[["RNA"]]@data@Dimnames[[2]]
cds@clusters@listData[["UMAP"]][["clusters"]] <- list_cluster

cds_cd14 <- cds[rowData(cds)$gene_short_name %in% c("CD14"), cds@clusters$UMAP$clusters %in% c(1)]
max.sub <- colnames(cds_cd14)[which.max(as.vector(counts(cds_cd14)))]
print(max.sub)
# "ACCAGTACAGCGTCCA-1_2"

# define a step to increase the percent
stp <- 5

pdf(file = paste0(outDir, "/without_Tcells.v3.Monocle3_traj_stability_", stp, "percent.pdf"),
    width = 19.2, height = 10.8)

for (pl in 1:10) {
  
  plts <- lapply(1: (100/stp), function(i) {
    message(i)
    tmp.cells <- c()
    for (j in 1:length(unique(seurat$integrated_snn_res.1))) {
      Idents(seurat) <- "integrated_snn_res.1"
      t <- WhichCells(seurat, idents = unique(seurat$integrated_snn_res.1)[j])
      set.seed(stp*pl)
      random.cells <- sample(t, floor(length(t)/100 * stp*i))
      # print(head(random.cells, 3))
      
      tmp.cells <- append(tmp.cells, random.cells)
    }
    
    # adding fixed starting/root cell, i.e cell with maximum expression of CD14 in cluster 1 to cells list
    tmp.cells <- union(tmp.cells, max.sub)
    #print(tail(tmp.cells))
    
    tmp.sc <- seurat[, tmp.cells]
    # print(table(tmp.sc$integrated_snn_res.1))
    
    tmp.cds <- SeuratWrappers::as.cell_data_set(tmp.sc)
    
    tmp.cds <- cluster_cells(tmp.cds)
    
    # migrate the cluster info from the original seurat object to the monocle3 object
    list_cluster <- tmp.sc@meta.data[[sprintf("integrated_snn_res.1")]]
    names(list_cluster) <- tmp.sc@assays[["RNA"]]@data@Dimnames[[2]]
    tmp.cds@clusters@listData[["UMAP"]][["clusters"]] <- list_cluster
    #print(tmp.cds@clusters@listData[["UMAP"]][["clusters"]][max.sub])
    
    # Monocle3 specifically looks for gene_short_name column
    rowData(tmp.cds)$gene_short_name <- row.names(rowData(tmp.cds))
    
    tmp.cds <- learn_graph(tmp.cds, use_partition = F)
    
    #tmp.cds <- order_cells(tmp.cds, root_cells = "CCTCTCCAGATGGTCG-1_7")
    tmp.cds <- order_cells(tmp.cds, root_cells = "ACCAGTACAGCGTCCA-1_2")
    plot_cells(tmp.cds,
               color_cells_by = "integrated_snn_res.1",
               label_groups_by_cluster = F,
               trajectory_graph_color = "black",
               label_leaves = F,
               label_branch_points = F,
               alpha = 0.7,
               group_label_size = 3) + 
      ggtitle(paste0(stp*i, " percent cells per cluster - ", pl))
  })
  
  
  # png(file = paste0(Results, "Monocle3_traj_stability_", stp, "percent.png"),
  #     res = 450, width = 19.2, height = 10.8, units = "in")
  p <- plot_grid(plts[[1]], plts[[2]], plts[[3]], plts[[4]], plts[[5]],
                 plts[[6]], plts[[7]], plts[[8]], plts[[9]], plts[[10]],
                 plts[[11]], plts[[12]], plts[[13]], plts[[14]], plts[[15]],
                 plts[[16]], plts[[17]], plts[[18]], plts[[19]], plts[[20]],
                 nrow = 4)
  print(p)
  # print(plts[[1]])
  # dev.off()
  rm(plts, p)
  
}

dev.off()

## Fixed percentages
# and 'closed_loop = FALSE``

# define a step to increase the percent
stp <- 5

pdf(file = paste0(outDir, "/without_Tcells.v3.Monocle3_traj_stability_noClosedLoop", stp, "percent.pdf"),
    width = 19.2, height = 10.8)

for (pl in 1:10) {
  
  plts <- lapply(1: (75/stp), function(i) {
    message(i)
    tmp.cells <- c()
    for (j in 1:length(unique(seurat$integrated_snn_res.1))) {
      Idents(seurat) <- "integrated_snn_res.1"
      t <- WhichCells(seurat, idents = unique(seurat$integrated_snn_res.1)[j])
      set.seed(stp*pl)
      random.cells <- sample(t, floor(length(t)/100 * (25 + stp*i)))
      # print(head(random.cells, 3))
      
      tmp.cells <- append(tmp.cells, random.cells)
    }
    
    # adding fixed starting/root cell, i.e cell with maximum expression of CD14 in cluster 1 to cells list
    tmp.cells <- union(tmp.cells, max.sub)
    #print(tail(tmp.cells))
    
    tmp.sc <- seurat[, tmp.cells]
    # print(table(tmp.sc$integrated_snn_res.1))
    
    tmp.cds <- SeuratWrappers::as.cell_data_set(tmp.sc)
    
    tmp.cds <- cluster_cells(tmp.cds)
    
    # migrate the cluster info from the original seurat object to the monocle3 object
    list_cluster <- tmp.sc@meta.data[[sprintf("integrated_snn_res.1")]]
    names(list_cluster) <- tmp.sc@assays[["RNA"]]@data@Dimnames[[2]]
    tmp.cds@clusters@listData[["UMAP"]][["clusters"]] <- list_cluster
    #print(tmp.cds@clusters@listData[["UMAP"]][["clusters"]][max.sub])
    
    # Monocle3 specifically looks for gene_short_name column
    rowData(tmp.cds)$gene_short_name <- row.names(rowData(tmp.cds))
    
    tmp.cds <- learn_graph(tmp.cds, use_partition = F, close_loop = FALSE)
    
    tmp.cds <- order_cells(tmp.cds, root_cells = "ACCAGTACAGCGTCCA-1_2")
    plot_cells(tmp.cds,
               color_cells_by = "integrated_snn_res.1",
               label_groups_by_cluster = F,
               trajectory_graph_color = "black",
               label_leaves = F,
               label_branch_points = F,
               alpha = 0.7,
               group_label_size = 3) + 
      ggtitle(paste0((25 + stp*i), " percent cells per cluster - ", pl))
  })
  
  
  # png(file = paste0(Results, "Monocle3_traj_stability_", stp, "percent.png"),
  #     res = 450, width = 19.2, height = 10.8, units = "in")
  p <- plot_grid(plts[[1]], plts[[2]], plts[[3]], plts[[4]], plts[[5]],
                 plts[[6]], plts[[7]], plts[[8]], plts[[9]], plts[[10]],
                 plts[[11]], plts[[12]], plts[[13]], plts[[14]], plts[[15]],
                 nrow = 3)
  print(p)
  # print(plts[[1]])
  # dev.off()
  rm(plts, p)
  
}

dev.off()

# And do the same, but now coloring using pseudotime
pdf(file = paste0(outDir, "/without_Tcells.v3.Monocle3_traj_stability_noClosedLoop", stp, "percent_pseudoTimeColored.pdf"),
    width = 19.2, height = 10.8)

for (pl in 1:10) {
  plts <- list()
  plts_pt <- list()
  
  for (i in 1: (75/stp)) {
    message(i)
    tmp.cells <- c()
    for (j in 1:length(unique(seurat$integrated_snn_res.1))) {
      Idents(seurat) <- "integrated_snn_res.1"
      t <- WhichCells(seurat, idents = unique(seurat$integrated_snn_res.1)[j])
      set.seed(stp*pl)
      random.cells <- sample(t, floor(length(t)/100 * (25 + stp*i)))
      # print(head(random.cells, 3))
      
      tmp.cells <- append(tmp.cells, random.cells)
    }
    
    # adding fixed starting/root cell, i.e cell with maximum expression of CD14 in cluster 1 to cells list
    tmp.cells <- union(tmp.cells, max.sub)
    #print(tail(tmp.cells))
    
    tmp.sc <- seurat[, tmp.cells]
    # print(table(tmp.sc$integrated_snn_res.1))
    
    tmp.cds <- SeuratWrappers::as.cell_data_set(tmp.sc)
    
    tmp.cds <- cluster_cells(tmp.cds)
    
    # migrate the cluster info from the original seurat object to the monocle3 object
    list_cluster <- tmp.sc@meta.data[[sprintf("integrated_snn_res.1")]]
    names(list_cluster) <- tmp.sc@assays[["RNA"]]@data@Dimnames[[2]]
    tmp.cds@clusters@listData[["UMAP"]][["clusters"]] <- list_cluster
    #print(tmp.cds@clusters@listData[["UMAP"]][["clusters"]][max.sub])
    
    # Monocle3 specifically looks for gene_short_name column
    rowData(tmp.cds)$gene_short_name <- row.names(rowData(tmp.cds))
    
    tmp.cds <- learn_graph(tmp.cds, use_partition = F, close_loop = FALSE)
    
    tmp.cds <- order_cells(tmp.cds, root_cells = "ACCAGTACAGCGTCCA-1_2")
    
    plts[[i]] <- plot_cells(tmp.cds,
               color_cells_by = "integrated_snn_res.1",
               label_groups_by_cluster = F,
               trajectory_graph_color = "black",
               label_leaves = F,
               label_branch_points = F,
               alpha = 0.7,
               group_label_size = 3) + 
      ggtitle(paste0((25 + stp*i), " percent cells per cluster - ", pl))
    
    plts_pt[[i]] <- plot_cells(tmp.cds,
               color_cells_by = "pseudotime",
               label_groups_by_cluster = F,
               trajectory_graph_color = "black",
               label_leaves = F,
               label_branch_points = F,
               alpha = 0.7,
               group_label_size = 3) 
    
    
  }
  
  # png(file = paste0(Results, "Monocle3_traj_stability_", stp, "percent.png"),
  #     res = 450, width = 19.2, height = 10.8, units = "in")
  p1 <- plot_grid(plts[[1]]/plts_pt[[1]], plts[[2]]/plts_pt[[2]], plts[[3]]/plts_pt[[3]], plts[[4]]/plts_pt[[4]], plts[[5]]/plts_pt[[5]], nrow = 1, ncol=5)
  p2 <- plot_grid(plts[[6]]/plts_pt[[6]], plts[[7]]/plts_pt[[7]], plts[[8]]/plts_pt[[8]], plts[[9]]/plts_pt[[9]], plts[[10]]/plts_pt[[10]], nrow = 1, ncol=5)
  p3 <- plot_grid(plts[[11]]/plts_pt[[11]], plts[[12]]/plts_pt[[12]], plts[[13]]/plts_pt[[13]], plts[[14]]/plts_pt[[14]], plts[[15]]/plts_pt[[15]],nrow = 1, ncol=5)
  print(p1)
  print(p2)
  print(p3)

  # And clena
  rm(plts, plts_pt, p1, p2, p3)
  
}

dev.off()


####
# Extract the pseudotime and add to the seurat object
seurat <- AddMetaData(
  object = seurat,
  metadata = cds.Node@principal_graph_aux@listData$UMAP$pseudotime,
  col.name = "pseudoTime_Monocle3"
)

# And visualise - make it a function
plotPseudotime_Seurat <- function(seuratObject = seurat, monocleObject = cds, pseudotimeCol = "pseudoTime_Monocle3", 
                                  outDir = "resultsDir", fileName = "plotPseudotime_Seurat"){
  if (!dir.exists(outDir)){
    dir.create(outDir)
  }
  
  pdf(paste(outDir,paste0(fileName,".pdf"), sep="/"))
    p <- FeaturePlot(seuratObject, pseudotimeCol, pt.size = 0.4) & scale_color_viridis_c()
    print(p)
  dev.off()
  
  # How to also plot the trajectory on the FeaturePlot ??
  # Trial code ....
  
  # all coordinates are perhaps in monocleObject@principal_graph_aux@listData$UMAP$dp_mst
  Ycoords.df <- data.frame(x=monocleObject@principal_graph_aux@listData$UMAP$dp_mst["UMAP_1",],
                           y=monocleObject@principal_graph_aux@listData$UMAP$dp_mst["UMAP_2",])
  rownames(Ycoords.df) <- colnames(monocleObject@principal_graph_aux@listData$UMAP$dp_mst)
  
  # And now translate the edges list in monocleObject@principal_graph$UMAP
  Yedges.df <- as.data.frame(get.edgelist(monocleObject@principal_graph$UMAP))
  colnames(Yedges.df) <- c("from","to")
  
  from.df <- as.data.frame(Ycoords.df[Yedges.df$from,])
  colnames(from.df) <- c("x.from","y.from")
  
  to.df <- as.data.frame(Ycoords.df[Yedges.df$to,])
  colnames(to.df) <- c("x.to","y.to")
  
  from_to.df <- cbind(from.df,to.df)
  
  # # Select the endpoints and branchpoints (basically, this is also offered by monocle3 ...)
  # # So, they either occur ones or occur 3 or more times
  # interestingPoints <- table(sub("\\..*", "", rownames(from_to.df)))
  # interestingPoints <- names(interestingPoints[table(sub("\\..*", "", rownames(from_to.df))) >= 3 | table(sub("\\..*", "", rownames(from_to.df))) == 1])
  # Ycoords.df <- Ycoords.df[interestingPoints,]
  
  pdf(paste(outDir,paste0(fileName,".v2.pdf"), sep="/"))
  p <- FeaturePlot(seuratObject, pseudotimeCol, pt.size = 0.4) & scale_color_viridis_c()
  p <- p + geom_segment(data=from_to.df,aes(x=x.from, xend = x.to, y=y.from, yend = y.to),
                        lwd=1,colour="blue") +
    geom_point(data=Ycoords.df, aes(x=x, y=y), colour = "red")
  
  p <- p + geom_label_repel(data = Ycoords.df, mapping=aes(x, y, label=rownames(Ycoords.df)),
                            max.overlaps = 10, label.size = NA, 
                            size = 3,
                            alpha = 0.6, 
                            force = 1,
                            label.padding=.1,
                            min.segment.length = 3,
                            fontface = 'bold', color = "black",
                            #box.padding = 0.80, point.padding = 0.5,
                            na.rm=TRUE,
                            seed = 42) +
    geom_label_repel(data = Ycoords.df, mapping=aes(x, y, label=rownames(Ycoords.df)),
                     max.overlaps = 10, label.size = NA, 
                     size = 3,
                     alpha = 1, 
                     force = 1,
                     label.padding=.1, 
                     min.segment.length = 3,
                     fontface = 'bold', color = "black",
                     na.rm=TRUE, fill = NA,
                     seed = 42)
  print(p)
  
  # And a plot with just the root node (if present)
  # Get the root node if specified
  if (!is.null(monocleObject@principal_graph_aux@listData$UMAP$root_pr_nodes)){
    rootNode.df <- Ycoords.df[monocleObject@principal_graph_aux@listData$UMAP$root_pr_nodes,]
    p2 <- FeaturePlot(seuratObject, pseudotimeCol, pt.size = 0.4) & scale_color_viridis_c()
    p2 <- p2 + geom_segment(data=from_to.df,aes(x=x.from, xend = x.to, y=y.from, yend = y.to),
                        lwd=1,colour="blue") +
          geom_point(data=rootNode.df, aes(x=x, y=y), colour = "red")
  
    p2 <- p2 + geom_label_repel(data = rootNode.df, mapping=aes(x, y, label=rownames(rootNode.df)),
                            max.overlaps = 10, label.size = NA, 
                            size = 3,
                            alpha = 0.6, 
                            force = 1,
                            label.padding=.1,
                            min.segment.length = 3,
                            fontface = 'bold', color = "red",
                            #box.padding = 0.80, point.padding = 0.5,
                            na.rm=TRUE,
                            seed = 42) +
      geom_label_repel(data = rootNode.df, mapping=aes(x, y, label=rownames(rootNode.df)),
                     max.overlaps = 10, label.size = NA, 
                     size = 3,
                     alpha = 1, 
                     force = 1,
                     label.padding=.1, 
                     min.segment.length = 3,
                     fontface = 'bold', color = "red",
                     na.rm=TRUE, fill = NA,
                     seed = 42)
  
    print(p2)
  
  }
  dev.off()
}

#
plotPseudotime_Seurat(seuratObject = seurat, monocleObject = cds.Node, pseudotimeCol = "pseudoTime_Monocle3", outDir = outDir, fileName = "plotPseudotime_Seurat")


#### RNA Velocity ####
# see https://github.com/basilkhuder/Seurat-to-RNA-Velocity#integrating-loom-file-and-meta-data
seurat <- readRDS("10_Reclustering_MP_and_MC/selectedClusters.subjectIntegrated.reclustered.kidney_immune_reference.rds")

# Write out the cell ID's
write.csv(Cells(seurat), file = paste0(outDir,"/cellID_integrated.csv"), row.names = FALSE)

# Write the embeddings
write.csv(Embeddings(seurat, reduction = "umap"), file = paste0(outDir,"/cell_embeddings.csv"))

# And extract the clusters
write.csv(data.frame('Cell ID'=rownames(seurat@meta.data), cluster=seurat@meta.data$integrated_snn_res.0.3), file = paste0(outDir,"/clusters_integrated_snn_res.0.3.csv"), row.names = FALSE)
write.csv(data.frame('Cell ID'=rownames(seurat@meta.data), cluster=seurat@meta.data$integrated_snn_res.0.6), file = paste0(outDir,"/clusters_integrated_snn_res.0.6.csv"), row.names = FALSE)
write.csv(data.frame('Cell ID'=rownames(seurat@meta.data), cluster=seurat@meta.data$integrated_snn_res.1), file = paste0(outDir,"/clusters_integrated_snn_res.1.csv"), row.names = FALSE)
table(seurat@meta.data$orig.ident)
# 
# MOMA17 MOMA302  MOMA52  MOMA57  MOMA67  MOMA68  MOMA72 
# 107     201     276    1517     292    1538     263 

sum(grepl("-1_1", rownames(seurat@meta.data)))
# [1] 276

sum(grepl("-1_2", rownames(seurat@meta.data)))
# [1] 1517

sum(grepl("-1_3", rownames(seurat@meta.data)))
# [1] 107

sum(grepl("-1_4", rownames(seurat@meta.data)))
# [1] 201

sum(grepl("-1_5", rownames(seurat@meta.data)))
# [1] 292

sum(grepl("-1_6", rownames(seurat@meta.data)))
# [1] 1538

sum(grepl("-1_7", rownames(seurat@meta.data)))
# [1] 263


#### Monocle - 2 starting points ####
# Yosta/ Koen Prange 20230828: Choose two starting points, one in the classical monocytes, one in the resident macrophages
# see also:
#  https://github.com/cole-trapnell-lab/monocle3/issues/328 for automagically choosiing startiing cells

library(igraph)  # needed for 'get.edgelist'
library(ggrepel)
library(cowplot)

seurat <- readRDS("10_Reclustering_MP_and_MC/selectedClusters.subjectIntegrated.reclustered.kidney_immune_reference.rds")
DefaultAssay(seurat) <- "RNA"

seurat <- subset(seurat, subset = integrated_snn_res.1 %in% c(0,1,3,4,7,8,9,11,12))
seurat@meta.data$integrated_snn_res.1 <- factor(seurat@meta.data$integrated_snn_res.1)

# Convert to Monocle ready object
cds <- SeuratWrappers::as.cell_data_set(seurat)
rowData(cds)$gene_short_name <- row.names(rowData(cds))
# migrate the cluster info from the original seurat object to the monocle3 object
list_cluster <- seurat@meta.data[[sprintf("integrated_snn_res.1")]]
names(list_cluster) <- seurat@assays[["RNA"]]@data@Dimnames[[2]]
cds@clusters@listData[["UMAP"]][["clusters"]] <- list_cluster

cds_cd14 <- cds[rowData(cds)$gene_short_name %in% c("CD14"), cds@clusters$UMAP$clusters %in% c(1)]
max.sub.CMDM <- colnames(cds_cd14)[which.max(as.vector(counts(cds_cd14)))]
print(max.sub.CMDM)
# Huh, now it is [1] "CCTCTCCAGATGGTCG-1_7" (as it was previously ....), last time it came with "ACCAGTACAGCGTCCA-1_2" 

cds_c1qc <- cds[rowData(cds)$gene_short_name %in% c("C1QC"), cds@clusters$UMAP$clusters %in% c(4)]
max.sub.ResM <- colnames(cds_c1qc)[which.max(as.vector(counts(cds_c1qc)))]
print(max.sub.ResM)
# "GTGGGTCCATTAGCCA-1_1"

# # Or use the helper function supplied by the tutorial 
# # http://cole-trapnell-lab.github.io/monocle-release/monocle3/#tutorial-1-learning-trajectories-with-monocle-3
# # a helper function to identify the root principal points:
# get_correct_root_state <- function(cds, cell_phenotype, root_type){
#   cell_ids <- which(pData(cds)[, cell_phenotype] == root_type)
#   
#   closest_vertex <-
#     cds@principal_graph_aux$UMAP$pr_graph_cell_proj_closest_vertex
#   closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
#   root_pr_nodes <-
#     V(cds@principal_graph_aux$UMAP$pr_graph_cell_proj_tree)$name[as.numeric(names
#                                            (which.max(table(closest_vertex[cell_ids,]))))]
#   
#   root_pr_nodes
# }

# Determine middle point of cells belonging to cluster
centerPoint <- function(cds, cell_phenotype, root_type){
  cell_ids <- which(pData(cds)[, cell_phenotype] == root_type)
  m <- cds@int_colData$reducedDims$UMAP[cell_ids,]
  cnt = c(mean(m[,1]),mean(m[,2]))
  dist_to_cnt <- apply(m,1,function(x,cnt) {(sqrt((x[1] - cnt[1])^2+(x[2]-cnt[2])^2))},cnt)
  df <- data.frame(CellID = rownames(m), dist=dist_to_cnt)
  
  idx <- which.min(df$dist)
  centerCell <- df[idx,'CellID']

  centerCell
}
# Just run 'learnGraph' once on the complete object to deduce the cells we want as root nodes
tmp.cds <- SeuratWrappers::as.cell_data_set(seurat)

tmp.cds <- cluster_cells(tmp.cds)

# migrate the cluster info from the original seurat object to the monocle3 object
list_cluster <- seurat@meta.data[[sprintf("integrated_snn_res.1")]]
names(list_cluster) <- seurat@assays[["RNA"]]@data@Dimnames[[2]]
tmp.cds@clusters@listData[["UMAP"]][["clusters"]] <- list_cluster
#print(tmp.cds@clusters@listData[["UMAP"]][["clusters"]][max.sub])

# Monocle3 specifically looks for gene_short_name column
rowData(tmp.cds)$gene_short_name <- row.names(rowData(tmp.cds))

tmp.cds <- learn_graph(tmp.cds, use_partition = T, close_loop = FALSE)

#cluster1_ids = get_correct_root_state(tmp.cds, cell_phenotype = 'integrated_snn_res.1', "1")
#cluster4_ids = get_correct_root_state(tmp.cds, cell_phenotype = 'integrated_snn_res.1', "4")
cluster1_ids = centerPoint(tmp.cds, cell_phenotype = 'integrated_snn_res.1', "1")
cluster4_ids = centerPoint(tmp.cds, cell_phenotype = 'integrated_snn_res.1', "4")

## Fixed percentages
# and 'closed_loop = FALSE``

# define a step to increase the percent
stp <- 5

pdf(file = paste0(outDir, "/without_Tcells.v4.Monocle3_traj_stability_noClosedLoop", stp, "percent.pdf"),
    width = 19.2, height = 10.8)

for (pl in 1:10) {
  
  plts <- lapply(1: (75/stp), function(i) {
    message(i)
    tmp.cells <- c()
    for (j in 1:length(unique(seurat$integrated_snn_res.1))) {
      Idents(seurat) <- "integrated_snn_res.1"
      t <- WhichCells(seurat, idents = unique(seurat$integrated_snn_res.1)[j])
      set.seed(stp*pl)
      random.cells <- sample(t, floor(length(t)/100 * (25 + stp*i)))
      # print(head(random.cells, 3))
      
      tmp.cells <- append(tmp.cells, random.cells)
    }
    
    # adding fixed starting/root cell, i.e cells closest to the principal nodes iin cluster 1 and 4
    tmp.cells <- union(tmp.cells, c(cluster1_ids, cluster4_ids))
    #print(tail(tmp.cells))
    
    tmp.sc <- seurat[, tmp.cells]
    # print(table(tmp.sc$integrated_snn_res.1))
    
    tmp.cds <- SeuratWrappers::as.cell_data_set(tmp.sc)
    
    tmp.cds <- cluster_cells(tmp.cds)
    
    # migrate the cluster info from the original seurat object to the monocle3 object
    list_cluster <- tmp.sc@meta.data[[sprintf("integrated_snn_res.1")]]
    names(list_cluster) <- tmp.sc@assays[["RNA"]]@data@Dimnames[[2]]
    tmp.cds@clusters@listData[["UMAP"]][["clusters"]] <- list_cluster
    #print(tmp.cds@clusters@listData[["UMAP"]][["clusters"]][max.sub])
    
    # Monocle3 specifically looks for gene_short_name column
    rowData(tmp.cds)$gene_short_name <- row.names(rowData(tmp.cds))
    
    tmp.cds <- learn_graph(tmp.cds, use_partition = F, close_loop = FALSE)
    
    tmp.cds <- order_cells(tmp.cds, root_cells = c(cluster1_ids, cluster4_ids))
    plot_cells(tmp.cds,
               color_cells_by = "integrated_snn_res.1",
               label_groups_by_cluster = F,
               trajectory_graph_color = "black",
               label_leaves = F,
               label_branch_points = F,
               alpha = 0.7,
               group_label_size = 3) + 
      ggtitle(paste0((25 + stp*i), " percent cells per cluster - ", pl))
  })
  
  
  # png(file = paste0(Results, "Monocle3_traj_stability_", stp, "percent.png"),
  #     res = 450, width = 19.2, height = 10.8, units = "in")
  p <- plot_grid(plts[[1]], plts[[2]], plts[[3]], plts[[4]], plts[[5]],
                 plts[[6]], plts[[7]], plts[[8]], plts[[9]], plts[[10]],
                 plts[[11]], plts[[12]], plts[[13]], plts[[14]], plts[[15]],
                 nrow = 3)
  print(p)
  # print(plts[[1]])
  # dev.off()
  rm(plts, p)
  
}

dev.off()

# And do the same, but now coloring using pseudotime
pdf(file = paste0(outDir, "/without_Tcells.v4.Monocle3_traj_stability_noClosedLoop", stp, "percent_pseudoTimeColored.pdf"),
    width = 19.2, height = 10.8)

for (pl in 1:10) {
  plts <- list()
  plts_pt <- list()
  
  for (i in 1: (75/stp)) {
    message(i)
    tmp.cells <- c()
    for (j in 1:length(unique(seurat$integrated_snn_res.1))) {
      Idents(seurat) <- "integrated_snn_res.1"
      t <- WhichCells(seurat, idents = unique(seurat$integrated_snn_res.1)[j])
      set.seed(stp*pl)
      random.cells <- sample(t, floor(length(t)/100 * (25 + stp*i)))
      # print(head(random.cells, 3))
      
      tmp.cells <- append(tmp.cells, random.cells)
    }
    
    # adding fixed starting/root cell, i.e cells closest to the principal nodes iin cluster 1 and 4
    tmp.cells <- union(tmp.cells, c(cluster1_ids, cluster4_ids))
    #print(tail(tmp.cells))
    
    tmp.sc <- seurat[, tmp.cells]
    # print(table(tmp.sc$integrated_snn_res.1))
    
    tmp.cds <- SeuratWrappers::as.cell_data_set(tmp.sc)
    
    tmp.cds <- cluster_cells(tmp.cds)
    
    # migrate the cluster info from the original seurat object to the monocle3 object
    list_cluster <- tmp.sc@meta.data[[sprintf("integrated_snn_res.1")]]
    names(list_cluster) <- tmp.sc@assays[["RNA"]]@data@Dimnames[[2]]
    tmp.cds@clusters@listData[["UMAP"]][["clusters"]] <- list_cluster
    #print(tmp.cds@clusters@listData[["UMAP"]][["clusters"]][max.sub])
    
    # Monocle3 specifically looks for gene_short_name column
    rowData(tmp.cds)$gene_short_name <- row.names(rowData(tmp.cds))
    
    tmp.cds <- learn_graph(tmp.cds, use_partition = F, close_loop = FALSE)
    
    tmp.cds <- order_cells(tmp.cds, root_cells = c(cluster1_ids, cluster4_ids))
    
    plts[[i]] <- plot_cells(tmp.cds,
                            color_cells_by = "integrated_snn_res.1",
                            label_groups_by_cluster = F,
                            trajectory_graph_color = "black",
                            label_leaves = F,
                            label_branch_points = F,
                            alpha = 0.7,
                            group_label_size = 3) + 
      ggtitle(paste0((25 + stp*i), " percent cells per cluster - ", pl))
    
    plts_pt[[i]] <- plot_cells(tmp.cds,
                               color_cells_by = "pseudotime",
                               label_groups_by_cluster = F,
                               trajectory_graph_color = "black",
                               label_leaves = F,
                               label_branch_points = F,
                               alpha = 0.7,
                               group_label_size = 3) 
    
    
  }
  
  # png(file = paste0(Results, "Monocle3_traj_stability_", stp, "percent.png"),
  #     res = 450, width = 19.2, height = 10.8, units = "in")
  p1 <- plot_grid(plts[[1]]/plts_pt[[1]], plts[[2]]/plts_pt[[2]], plts[[3]]/plts_pt[[3]], plts[[4]]/plts_pt[[4]], plts[[5]]/plts_pt[[5]], nrow = 1, ncol=5)
  p2 <- plot_grid(plts[[6]]/plts_pt[[6]], plts[[7]]/plts_pt[[7]], plts[[8]]/plts_pt[[8]], plts[[9]]/plts_pt[[9]], plts[[10]]/plts_pt[[10]], nrow = 1, ncol=5)
  p3 <- plot_grid(plts[[11]]/plts_pt[[11]], plts[[12]]/plts_pt[[12]], plts[[13]]/plts_pt[[13]], plts[[14]]/plts_pt[[14]], plts[[15]]/plts_pt[[15]],nrow = 1, ncol=5)
  print(p1)
  print(p2)
  print(p3)
  
  # And clean
  rm(plts, plts_pt, p1, p2, p3)
  
}

dev.off()


#### Monocle - no DC's, 2 starting points ####
# Yosta/ Koen Prange 20230828: Choose two starting points, one in the classical monocytes, one in the resident macrophages
#- And leave out the cDC clusters (i.e. cluster, 3, 7 and 12??)
# see also:
#  https://github.com/cole-trapnell-lab/monocle3/issues/328 for automagically choosing starting cells
# 
# For visualization, plot all trajectories at once as a kind of 'trumpet', and perhaps use the density as some kind of measure?
#

library(igraph)  # needed for 'get.edgelist'
library(ggrepel)
library(cowplot)

seurat <- readRDS("10_Reclustering_MP_and_MC/selectedClusters.subjectIntegrated.reclustered.kidney_immune_reference.rds")
DefaultAssay(seurat) <- "RNA"

seurat <- subset(seurat, subset = integrated_snn_res.1 %in% c(0,1,4,8,9,11))
seurat@meta.data$integrated_snn_res.1 <- factor(seurat@meta.data$integrated_snn_res.1)

# Convert to Monocle ready object
cds <- SeuratWrappers::as.cell_data_set(seurat)
rowData(cds)$gene_short_name <- row.names(rowData(cds))
# migrate the cluster info from the original seurat object to the monocle3 object
list_cluster <- seurat@meta.data[[sprintf("integrated_snn_res.1")]]
names(list_cluster) <- seurat@assays[["RNA"]]@data@Dimnames[[2]]
cds@clusters@listData[["UMAP"]][["clusters"]] <- list_cluster

cds_cd14 <- cds[rowData(cds)$gene_short_name %in% c("CD14"), cds@clusters$UMAP$clusters %in% c(1)]
max.sub.CMDM <- colnames(cds_cd14)[which.max(as.vector(counts(cds_cd14)))]
print(max.sub.CMDM)
# Huh, now it is [1] "CCTCTCCAGATGGTCG-1_7" (as it was previously ....), last time it came with "ACCAGTACAGCGTCCA-1_2" 

cds_c1qc <- cds[rowData(cds)$gene_short_name %in% c("C1QC"), cds@clusters$UMAP$clusters %in% c(4)]
max.sub.ResM <- colnames(cds_c1qc)[which.max(as.vector(counts(cds_c1qc)))]
print(max.sub.ResM)
# "GTGGGTCCATTAGCCA-1_1"

# # Or use the helper function supplied by the tutorial 
# # http://cole-trapnell-lab.github.io/monocle-release/monocle3/#tutorial-1-learning-trajectories-with-monocle-3
# # a helper function to identify the root principal points:
# get_correct_root_state <- function(cds, cell_phenotype, root_type){
#   cell_ids <- which(pData(cds)[, cell_phenotype] == root_type)
#   
#   closest_vertex <-
#     cds@principal_graph_aux$UMAP$pr_graph_cell_proj_closest_vertex
#   closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
#   root_pr_nodes <-
#     V(cds@principal_graph_aux$UMAP$pr_graph_cell_proj_tree)$name[as.numeric(names
#                                            (which.max(table(closest_vertex[cell_ids,]))))]
#   
#   root_pr_nodes
# }

# Determine middle point of cells belonging to cluster
centerPoint <- function(cds, cell_phenotype, root_type){
  cell_ids <- which(pData(cds)[, cell_phenotype] == root_type)
  m <- cds@int_colData$reducedDims$UMAP[cell_ids,]
  cnt = c(mean(m[,1]),mean(m[,2]))
  dist_to_cnt <- apply(m,1,function(x,cnt) {(sqrt((x[1] - cnt[1])^2+(x[2]-cnt[2])^2))},cnt)
  df <- data.frame(CellID = rownames(m), dist=dist_to_cnt)
  
  idx <- which.min(df$dist)
  centerCell <- df[idx,'CellID']
  
  centerCell
}
# Just run 'learnGraph' once on the complete object to deduce the cells we want as root nodes
tmp.cds <- SeuratWrappers::as.cell_data_set(seurat)

tmp.cds <- cluster_cells(tmp.cds)

# migrate the cluster info from the original seurat object to the monocle3 object
list_cluster <- seurat@meta.data[[sprintf("integrated_snn_res.1")]]
names(list_cluster) <- seurat@assays[["RNA"]]@data@Dimnames[[2]]
tmp.cds@clusters@listData[["UMAP"]][["clusters"]] <- list_cluster
#print(tmp.cds@clusters@listData[["UMAP"]][["clusters"]][max.sub])

# Monocle3 specifically looks for gene_short_name column
rowData(tmp.cds)$gene_short_name <- row.names(rowData(tmp.cds))

tmp.cds <- learn_graph(tmp.cds, use_partition = T, close_loop = FALSE)

#cluster1_ids = get_correct_root_state(tmp.cds, cell_phenotype = 'integrated_snn_res.1', "1")
#cluster4_ids = get_correct_root_state(tmp.cds, cell_phenotype = 'integrated_snn_res.1', "4")
cluster1_ids = centerPoint(tmp.cds, cell_phenotype = 'integrated_snn_res.1', "1")
cluster4_ids = centerPoint(tmp.cds, cell_phenotype = 'integrated_snn_res.1', "4")

## Fixed percentages
# and 'closed_loop = FALSE``

# define a step to increase the percent
stp <- 5

pdf(file = paste0(outDir, "/without_cDC_and_Tcells.v1.Monocle3_traj_stability_noClosedLoop", stp, "percent.pdf"),
    width = 19.2, height = 10.8)

for (pl in 1:10) {
  
  plts <- lapply(1: (75/stp), function(i) {
    message(i)
    tmp.cells <- c()
    for (j in 1:length(unique(seurat$integrated_snn_res.1))) {
      Idents(seurat) <- "integrated_snn_res.1"
      t <- WhichCells(seurat, idents = unique(seurat$integrated_snn_res.1)[j])
      set.seed(stp*pl)
      random.cells <- sample(t, floor(length(t)/100 * (25 + stp*i)))
      # print(head(random.cells, 3))
      
      tmp.cells <- append(tmp.cells, random.cells)
    }
    
    # adding fixed starting/root cell, i.e cells closest to the principal nodes iin cluster 1 and 4
    tmp.cells <- union(tmp.cells, c(cluster1_ids, cluster4_ids))
    #print(tail(tmp.cells))
    
    tmp.sc <- seurat[, tmp.cells]
    # print(table(tmp.sc$integrated_snn_res.1))
    
    tmp.cds <- SeuratWrappers::as.cell_data_set(tmp.sc)
    
    tmp.cds <- cluster_cells(tmp.cds)
    
    # migrate the cluster info from the original seurat object to the monocle3 object
    list_cluster <- tmp.sc@meta.data[[sprintf("integrated_snn_res.1")]]
    names(list_cluster) <- tmp.sc@assays[["RNA"]]@data@Dimnames[[2]]
    tmp.cds@clusters@listData[["UMAP"]][["clusters"]] <- list_cluster
    #print(tmp.cds@clusters@listData[["UMAP"]][["clusters"]][max.sub])
    
    # Monocle3 specifically looks for gene_short_name column
    rowData(tmp.cds)$gene_short_name <- row.names(rowData(tmp.cds))
    
    tmp.cds <- learn_graph(tmp.cds, use_partition = F, close_loop = FALSE)
    
    tmp.cds <- order_cells(tmp.cds, root_cells = c(cluster1_ids, cluster4_ids))
    plot_cells(tmp.cds,
               color_cells_by = "integrated_snn_res.1",
               label_groups_by_cluster = F,
               trajectory_graph_color = "black",
               label_leaves = F,
               label_branch_points = F,
               alpha = 0.7,
               group_label_size = 3) + 
      ggtitle(paste0((25 + stp*i), " percent cells per cluster - ", pl))
  })
  
  
  # png(file = paste0(Results, "Monocle3_traj_stability_", stp, "percent.png"),
  #     res = 450, width = 19.2, height = 10.8, units = "in")
  p <- plot_grid(plts[[1]], plts[[2]], plts[[3]], plts[[4]], plts[[5]],
                 plts[[6]], plts[[7]], plts[[8]], plts[[9]], plts[[10]],
                 plts[[11]], plts[[12]], plts[[13]], plts[[14]], plts[[15]],
                 nrow = 3)
  print(p)
  # print(plts[[1]])
  # dev.off()
  rm(plts, p)
  
}

dev.off()

# And do the same, but now coloring using pseudotime
pdf(file = paste0(outDir, "/without_cDC_and_Tcells.v1.Monocle3_traj_stability_noClosedLoop", stp, "percent_pseudoTimeColored.pdf"),
    width = 19.2, height = 10.8)

for (pl in 1:10) {
  plts <- list()
  plts_pt <- list()
  
  for (i in 1: (75/stp)) {
    message(i)
    tmp.cells <- c()
    for (j in 1:length(unique(seurat$integrated_snn_res.1))) {
      Idents(seurat) <- "integrated_snn_res.1"
      t <- WhichCells(seurat, idents = unique(seurat$integrated_snn_res.1)[j])
      set.seed(stp*pl)
      random.cells <- sample(t, floor(length(t)/100 * (25 + stp*i)))
      # print(head(random.cells, 3))
      
      tmp.cells <- append(tmp.cells, random.cells)
    }
    
    # adding fixed starting/root cell, i.e cells closest to the principal nodes iin cluster 1 and 4
    tmp.cells <- union(tmp.cells, c(cluster1_ids, cluster4_ids))
    #print(tail(tmp.cells))
    
    tmp.sc <- seurat[, tmp.cells]
    # print(table(tmp.sc$integrated_snn_res.1))
    
    tmp.cds <- SeuratWrappers::as.cell_data_set(tmp.sc)
    
    tmp.cds <- cluster_cells(tmp.cds)
    
    # migrate the cluster info from the original seurat object to the monocle3 object
    list_cluster <- tmp.sc@meta.data[[sprintf("integrated_snn_res.1")]]
    names(list_cluster) <- tmp.sc@assays[["RNA"]]@data@Dimnames[[2]]
    tmp.cds@clusters@listData[["UMAP"]][["clusters"]] <- list_cluster
    #print(tmp.cds@clusters@listData[["UMAP"]][["clusters"]][max.sub])
    
    # Monocle3 specifically looks for gene_short_name column
    rowData(tmp.cds)$gene_short_name <- row.names(rowData(tmp.cds))
    
    tmp.cds <- learn_graph(tmp.cds, use_partition = F, close_loop = FALSE)
    
    tmp.cds <- order_cells(tmp.cds, root_cells = c(cluster1_ids, cluster4_ids))
    
    plts[[i]] <- plot_cells(tmp.cds,
                            color_cells_by = "integrated_snn_res.1",
                            label_groups_by_cluster = F,
                            trajectory_graph_color = "black",
                            label_leaves = F,
                            label_branch_points = F,
                            alpha = 0.7,
                            group_label_size = 3) + 
      ggtitle(paste0((25 + stp*i), " percent cells per cluster - ", pl))
    
    plts_pt[[i]] <- plot_cells(tmp.cds,
                               color_cells_by = "pseudotime",
                               label_groups_by_cluster = F,
                               trajectory_graph_color = "black",
                               label_leaves = F,
                               label_branch_points = F,
                               alpha = 0.7,
                               group_label_size = 3) 
    
    
  }
  
  # png(file = paste0(Results, "Monocle3_traj_stability_", stp, "percent.png"),
  #     res = 450, width = 19.2, height = 10.8, units = "in")
  p1 <- plot_grid(plts[[1]]/plts_pt[[1]], plts[[2]]/plts_pt[[2]], plts[[3]]/plts_pt[[3]], plts[[4]]/plts_pt[[4]], plts[[5]]/plts_pt[[5]], nrow = 1, ncol=5)
  p2 <- plot_grid(plts[[6]]/plts_pt[[6]], plts[[7]]/plts_pt[[7]], plts[[8]]/plts_pt[[8]], plts[[9]]/plts_pt[[9]], plts[[10]]/plts_pt[[10]], nrow = 1, ncol=5)
  p3 <- plot_grid(plts[[11]]/plts_pt[[11]], plts[[12]]/plts_pt[[12]], plts[[13]]/plts_pt[[13]], plts[[14]]/plts_pt[[14]], plts[[15]]/plts_pt[[15]],nrow = 1, ncol=5)
  print(p1)
  print(p2)
  print(p3)
  
  # And clean
  rm(plts, plts_pt, p1, p2, p3)
  
}

dev.off()

## Run 1000 times taking 90% of cells per cluster with different seed
## And see what the found path is?

# Also get the correct colors
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

Idents(seurat) <- "integrated_snn_res.1"

clustColors <- gg_color_hue(13)
clustColors <- clustColors[c(0,1,4,8,9,11)+1]

# tmp.sc <- seurat
# # print(table(tmp.sc$integrated_snn_res.1))
# 
# tmp.cds <- SeuratWrappers::as.cell_data_set(tmp.sc)
# 
# tmp.cds <- cluster_cells(tmp.cds)
# 
# # migrate the cluster info from the original seurat object to the monocle3 object
# list_cluster <- tmp.sc@meta.data[[sprintf("integrated_snn_res.1")]]
# names(list_cluster) <- tmp.sc@assays[["RNA"]]@data@Dimnames[[2]]
# tmp.cds@clusters@listData[["UMAP"]][["clusters"]] <- list_cluster
# #print(tmp.cds@clusters@listData[["UMAP"]][["clusters"]][max.sub])
# 
# # Monocle3 specifically looks for gene_short_name column
# rowData(tmp.cds)$gene_short_name <- row.names(rowData(tmp.cds))

# Ran with 75, 90 and 100
pdf(file = paste0(outDir, "/without_cDC_and_Tcells.v2.Monocle3_traj_stability_noClosedLoop_75percent.pdf"),
    width = 19.2, height = 10.8)

p <- DimPlot(seurat, reduction="umap", pt.size = 0.4, cols = clustColors)

for (i in 1:1000) {
  
  message(i)
  set.seed(i)
  
  tmp.cells <- c()
  for (j in 1:length(unique(seurat$integrated_snn_res.1))) {
    Idents(seurat) <- "integrated_snn_res.1"
    t <- WhichCells(seurat, idents = unique(seurat$integrated_snn_res.1)[j])
    random.cells <- sample(t, floor((length(t)/100) * 75))   # changed to 75, 90 , 100
    # print(head(random.cells, 3))
    
    tmp.cells <- append(tmp.cells, random.cells)
  }
  
  # adding fixed starting/root cell, i.e cells closest to the principal nodes in cluster 1 and 4
  tmp.cells <- union(tmp.cells, c(cluster1_ids, cluster4_ids))
  #print(tail(tmp.cells))
  
  tmp.sc <- seurat[, tmp.cells]
  # print(table(tmp.sc$integrated_snn_res.1))
  
  tmp.cds <- SeuratWrappers::as.cell_data_set(tmp.sc)
  
  tmp.cds <- cluster_cells(tmp.cds)
  
  # migrate the cluster info from the original seurat object to the monocle3 object
  list_cluster <- tmp.sc@meta.data[[sprintf("integrated_snn_res.1")]]
  names(list_cluster) <- tmp.sc@assays[["RNA"]]@data@Dimnames[[2]]
  tmp.cds@clusters@listData[["UMAP"]][["clusters"]] <- list_cluster
  #print(tmp.cds@clusters@listData[["UMAP"]][["clusters"]][max.sub])
  
  # Monocle3 specifically looks for gene_short_name column
  rowData(tmp.cds)$gene_short_name <- row.names(rowData(tmp.cds))
  
  monocleObject <- learn_graph(tmp.cds, use_partition = F, close_loop = FALSE)
  monocleObject <- order_cells(monocleObject, root_cells = c(cluster1_ids, cluster4_ids))
  
  # all coordinates are perhaps in monocleObject@principal_graph_aux@listData$UMAP$dp_mst
  Ycoords.df <- data.frame(x=monocleObject@principal_graph_aux@listData$UMAP$dp_mst["UMAP_1",],
                           y=monocleObject@principal_graph_aux@listData$UMAP$dp_mst["UMAP_2",])
  rownames(Ycoords.df) <- colnames(monocleObject@principal_graph_aux@listData$UMAP$dp_mst)
  
  # And now translate the edges list in monocleObject@principal_graph$UMAP
  Yedges.df <- as.data.frame(get.edgelist(monocleObject@principal_graph$UMAP))
  colnames(Yedges.df) <- c("from","to")
  
  from.df <- as.data.frame(Ycoords.df[Yedges.df$from,])
  colnames(from.df) <- c("x.from","y.from")
  
  to.df <- as.data.frame(Ycoords.df[Yedges.df$to,])
  colnames(to.df) <- c("x.to","y.to")
  
  from_to.df <- cbind(from.df,to.df)
  
  if (i == 1){
    rootNode.df <- Ycoords.df[monocleObject@principal_graph_aux@listData$UMAP$root_pr_nodes,]
    p <- p + geom_segment(data=from_to.df,aes(x=x.from, xend = x.to, y=y.from, yend = y.to),
                          lwd=1,colour="lightblue") +
      geom_point(data=rootNode.df, aes(x=x, y=y), colour = "red")
    
    p <- p + geom_label_repel(data = rootNode.df, mapping=aes(x, y, label=rownames(rootNode.df)),
                              max.overlaps = 10, label.size = NA, 
                              size = 3,
                              alpha = 0.6, 
                              force = 1,
                              label.padding=.1,
                              min.segment.length = 3,
                              fontface = 'bold', color = "red",
                              #box.padding = 0.80, point.padding = 0.5,
                              na.rm=TRUE,
                              seed = 42) +
      geom_label_repel(data = rootNode.df, mapping=aes(x, y, label=rownames(rootNode.df)),
                       max.overlaps = 10, label.size = NA, 
                       size = 3,
                       alpha = 1, 
                       force = 1,
                       label.padding=.1, 
                       min.segment.length = 3,
                       fontface = 'bold', color = "red",
                       na.rm=TRUE, fill = NA,
                       seed = 42)
  } else {
    p <- p + geom_segment(data=from_to.df,aes(x=x.from, xend = x.to, y=y.from, yend = y.to),
                          lwd=1,colour="lightblue")
      
  }
}
print(p)


dev.off()

# There is no randomness involved if you aways have the same cluster info from Seurat.
# Both 'learn_graph' and 'order_cells' do not have a seed setting?
