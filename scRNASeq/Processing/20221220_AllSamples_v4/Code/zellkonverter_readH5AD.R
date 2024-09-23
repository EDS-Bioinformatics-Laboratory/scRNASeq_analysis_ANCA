## AJ - 20231124
# This script was used on an old MacBook Pro (2010) using R4.1 and the following versions
# (as obtained from /Library/Frameworks/R.framework/Versions/4.1/Resources/library)
# zellkonverter: v1.2.1
# singleCellTK : v2.2.0
# The script was dated July 10th 2022 (and used for 20220623_AllSamples_InitialAnalysis)

library(zellkonverter)

url <- "https://cellgeni.cog.sanger.ac.uk/kidneycellatlas/Mature_Immune_v2.1.h5ad"
curl::curl_download(url, basename(url))

sce <- readH5AD("Mature_Immune_v2.1.h5ad", reader = "python")
# Error in adata$uns$keys() : attempt to apply non-function

sce <- readH5AD("Mature_Immune_v2.1.h5ad", reader = "R", )
# Warning messages:
#   1: In value[[3L]](cond) :
#   setting 'colData' failed for 'Mature_Immune_v2.1.h5ad': HDF5. Links. Can't get
#   value.
# 2: In value[[3L]](cond) :
#   setting 'rowData' failed for 'Mature_Immune_v2.1.h5ad': HDF5. Links. Can't get
# value.
# 3: In value[[3L]](cond) :
#   setting 'reducedDims' failed for 'Mature_Immune_v2.1.h5ad': error in evaluating
# the argument 'x' in selecting a method for function 't': HDF5. Links. Can't get
#   value.

sce
# class: SingleCellExperiment 
# dim: 33694 7803 
# metadata(4): Experiment_categories Project_categories celltype_categories
# compartment_categories
# assays(1): X
# rownames: NULL
# rowData names(0):
#   colnames: NULL
# colData names(0):
#   reducedDimNames(0):
#   mainExpName: NULL
# altExpNames(0):

#### singleCellTK ####
# https://camplab.net/sctk/v2.5.1/articles/import_data.html
# Need to install python etc.
library('reticulate')
library('anndata')
py_install('scrublet', pip = TRUE)
py_install('scanorama', pip = TRUE)
py_install('bbknn', pip = TRUE)

library(singleCellTK)
sce <- importAnnData(sampleDirs = "./", sampleNames = "Mature_Immune_v2.1")
# Warning: unable to add 'X_umap' from .obsm AnnData slot to SCE metadata. Skipping. 
# sys:1: DeprecationWarning: Use dict(obj) instead of as_dict, as_dict will be removed in the future.

saveRDS(sce, "Mature_Immune_v2.1.noUMAP.sce")
