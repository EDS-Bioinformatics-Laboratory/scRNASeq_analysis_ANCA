# Start analysis 

#### Set locations, dropbox, helper directories and scripts ####
init <- function() {
  # Create directory to hold scripts
  scriptDir <<- "applied_scripts"
  if ( !file.exists(scriptDir)){
    dir.create(scriptDir)
  }
  
  # Create directory to hold results
  tablesDir <<- "resultTables"
  if ( !file.exists(tablesDir)){
    dir.create(tablesDir)
  }
  
  # Copy helper script to directory and source from there
  if ( !file.exists(paste0(scriptDir,"/generalFunctions.r"))){
    # Read the Dropbox location for additional scripts
    dropbox = ifelse(Sys.info()['sysname'][[1]] == "Darwin", "/Users/Aldo/Desktop/AMC/Dropbox", ifelse(Sys.info()['sysname'][[1]] == "Linux", "/home/ajongejan/Dropbox", "D:/Dropbox"))
    
    file.copy(paste0(dropbox,"/Support/R/generalFunctions.r"),scriptDir)
    source(paste(scriptDir,"generalFunctions.r", sep="/"))
  } else {
    source(paste(scriptDir,"generalFunctions.r", sep="/"))
  }
  
  # Set the dropbox location
  setDropboxLocation()
  
  #### Seurat ####
  # see: https://satijalab.org/seurat/v3.0/hashing_vignette.html
  #      https://satijalab.org/seurat/v3.1/multimodal_vignette.html
  library(Seurat)
}