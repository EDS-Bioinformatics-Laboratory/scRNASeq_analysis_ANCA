# Set dropbox location
setDropboxLocation <- function(){
  if ( Sys.info()['user'] == 'ajongejan' ){
    # Read the Dropbox location for additional scripts
    dropbox <<- ifelse(Sys.info()['sysname'][[1]] == "Darwin", "/Users/Aldo/Dropbox", ifelse(Sys.info()['sysname'][[1]] == "Linux", "/home/ajongejan/Dropbox", "D:/Dropbox"))
  } else if ( Sys.info()['user'] == 'pdmoerland' ){
    dropbox <<- ifelse(Sys.info()['sysname'][[1]] == "Darwin", "/Users/Aldo/Dropbox", ifelse(Sys.info()['sysname'][[1]] == "Linux", "/home/ajongejan/Dropbox", "D:/Dropbox"))
  }
}

# Function to catch the version of R and either use 'source(http://bioconductor.org/biocLite.R)' or 'BiocManager'
biocInstall <- function(x){
  if ( version$major >= 4 || version$minor > 5.1 ){
    BiocManager::install(x)
    } else {
      source("http://bioconductor.org/biocLite.R")
      biocLite(x, suppressUpdates=TRUE)
    }
}

# Function to write out Bioconductor version and sessionInfo (include R version)
writeVersions <- function(sessionDir="../Documentation"){
  write(paste0("Bioconductor version ", capture.output(tools:::.BioC_version_associated_with_R_version()),"\n"), paste0(sessionDir,"/sessionInfo.txt"))
  write(capture.output(sessionInfo()), paste0(sessionDir,"/sessionInfo.txt"), append=TRUE)
}

# Function checking whether to make a link to the dropbox for additional scripts
loadScript <- function(script, scriptDir.= "applied_scripts"){
  if ( Sys.info()['user'] == 'ajongejan' ){
    # Read the Dropbox location for additional scripts
    dropbox = ifelse(Sys.info()['sysname'][[1]] == "Darwin", "/Users/Aldo/Dropbox", ifelse(Sys.info()['sysname'][[1]] == "Linux", "/home/ajongejan/Dropbox", "D:/Dropbox"))
    file.copy(paste(dropbox,paste0("Support/R/",x),sep="/"),scriptDir.)
    source(paste(scriptDir.,x, sep="/"))
  } else if ( Sys.info()['user'] == 'pdmoerland' ){
    dropbox = ifelse(Sys.info()['sysname'][[1]] == "Darwin", "/Users/Aldo/Dropbox", ifelse(Sys.info()['sysname'][[1]] == "Linux", "/home/ajongejan/Dropbox", "D:/Dropbox"))
    file.copy(paste(dropbox,paste0("Support/R/",x),sep="/"),scriptDir.)
    source(paste(scriptDir.,x, sep="/"))
  } else {
    source(paste(scriptDir.,x, sep="/"))
  }
}

# File permissions
# https://stat.ethz.ch/R-manual/R-devel/library/base/html/files2.html
# - Set the permissions afterwards for a directory and all files?
# - put this into the "final"/"cleaning up" function? (fx. extend "writeVersions"??, make it into a "finishUpAnalysis" function?)
#


# Test files
# - test whether the results files really matches the corresponding object they originate from
# - test numbers Venn diagram and csv files
