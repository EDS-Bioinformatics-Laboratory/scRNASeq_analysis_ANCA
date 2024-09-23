# Just initiate VISION with the object
# see ' 12b_reclustering_vision.r' for the generation of the VISION object

# Install VISION
if (!requireNamespace("BiocManager") && tools:::.BioC_version_associated_with_R_version() >= 3.8 ){
  install.packages("BiocManager")
}

if(require("VISION") == FALSE) {
  BiocManager::install("devtools")
  library(devtools)
  install_github("YosefLab/VISION")
}
  
# Library
library(VISION)

outputDir <- "12b_reclustering_vision"
vision.obj <- readRDS(paste0(outputDir,"/vision.rds"))
viewResults(vision.obj)
