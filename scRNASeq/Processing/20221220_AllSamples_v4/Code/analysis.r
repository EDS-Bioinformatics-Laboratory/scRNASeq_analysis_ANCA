# AJ - 20221117
#

# Extended analysis of the scRNASeq data of Yosta Vegting / Marc Hilhorst
# 
#  For an initial analysis of MOMA17, 52,57, 67 and 68, see 20220623_AllSamples_InitialAnalysis
#   MOMA17  - ANCA/PR3 - 3' sequencing  (lots of tubulus cells)
#   MOMA52  - ANCA/PR3 - 5' with B and T cell hashtag
#   MOMA57  - Tumor controle - 5' with B and T cell hashtag
#   MOMA67  - ANCA/MPO - 3'
#   MOMA68  - SLE - 3'

# For the analysis including two more sample (MOMA72 and MOMA302, see 20220916_AllSamples_v2)

# In this analysis we try to learn from the previous analyses and finalize the analysis
# - There is a cluster, that contains a low number of expressed housekeepinig genes 
#   (see cluster 4 in 20220916_AllSamples_v2/Code/6_clustering/Mito10_clusters/violinPlot_qcStats_res.0.3.pdf)
#   as I seemingly missed these while filtering
# - We make decisions on the annotations to be used
#   - No SignacX
#   - SingleR before integration 
#   - CellTypist before annotation: bestMatch withou majorityVoting and probMatch+majorityVoting
#   - CellTypist after integration: probMatch+majorityVoting (as the integration can have influence on the majority voting)

# Integration of datasets, we decide for SCTransform
# see https://satijalab.org/seurat/articles/integration_introduction.html#performing-integration-on-datasets-normalized-with-sctransform-1

#### Set default for a good start ####
source("initAnalysis.r")
init()

# Rename the ../Data/NameOfDataset_1 folder to scRNASeq (via copy and remove)
if (dir.exists("../Data/NameOfDataset_1")){
  R.utils::copyDirectory("../Data/NameOfDataset_1", "../Data/scRNASeq")
  unlink("../Data/NameOfDataset_1", recursive=TRUE)
}

#### Analysis workflow ####
# Preprocessing - read in data separately, also read in ADT and VDJ info
source("1_preprocessing_and_QC.r")

# Filtering and normalization - per dataset
#
# Schroeder et al:
# "we expanded the final threshold to keep cells with <80% mitochondrial content and >200 genes per cell in an effort to 
#  maximize our cellular population. Additionally, cells with >6,000 genes per cell were excluded to eliminate the likelihood
#  of bias introduced by the presence of doublets."
#
# Filter on 'subset = nFeature_RNA >= 200 & nFeature_RNA < 6000 & percent.mito < pctMito'
#
# AJ - 20221117 - Also filter on the housekeeping genes!!
#  - subset = nFeature_RNA >= 200 & nFeature_RNA < 6000 & percent.mito < pctMito & n.exp.hkgenes > 55 
#  - We no longer explore cut-offs for mitochondrial content -> we use 10% (as that is what we finally used in 20220916_AllSamples_v2)

# Normalization use SCTransform
source("2_filtering_normalization_and_scaling.r")

## OBSERVATION:
# Mainly MOMA17 is now more heavily filtered, n.exp.hkgenes >55 might be to harsh, although the violin plots from step 1
# seem to suggest that it would be OK ...

# Doublet detection and removal
# Only for case where we filtered everything above 10% mitochondrial content
# Used xxxx PC when clustering ...
# Ran doublet detection with and without cells containing both B and T cell receptor sequences.
source("3_doublet_detection.r")

# Need to do annotation before integration
source("4_cell_annotation_before_integration.r")
# You can actually not run the script at once, as it also involves a
# step to annotate using CellTypist (and that runs in Python, a Jupyter Notebook)

# Remove tubulus cells
source("5_removal_tubuluscells.r")

# Integration
# - Integrate using data without cells containing both B and T cell receptor sequences for MOMA52 and 57.
source("6_integration.r")

# Annotation
source("7_cellTypist_annotation_after_integration.r")
# You can actually not run the script at once, as it also involves a
# step to annotate using CellTypist (and that runs in Python, a Jupyter Notebook)

# Clustering
source("8_clustering.r")

# Reference mapping
source("9_kidney_reference_mapping.r")

# Recluster macrophageand moncyte clusters
source("10_reclustering_MP_and_MC.r")
# You can actually not run the script at once, as it also involves a
# step to annotate using CellTypist (and that runs in Python, a Jupyter Notebook)


