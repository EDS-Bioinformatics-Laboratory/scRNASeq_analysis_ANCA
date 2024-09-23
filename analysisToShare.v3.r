# AJ - 30082022
#
# Script to extract all files from an analysis performed in our FSS directory structure
# and put it into the Sharing directory to share it with coworkers
# - make it self-contained
# - scripts should work
# - organize files
# - ...
#

# IMPORTANT QUESTIONS: 
# - Do we want out collaborators to run their own analysis, overwriting ours?
# - How could we prevent that?
# - We should allow calculation of Venn-diagrams (genes and genesets) with different cut-offs
#   Close everything except the venn and venn_gs folder?
# - Include the count files, even if they are not made by us (fx. see project Andy Hsiao) in the sharing directory
# - Include a RNASeq Workflow figure with all steps (see GitHub repository), perhaps add numbers in front of directories,
#   draft a text/markdown/pdf 'report'.
# - Extract all kind of information and put in an 'report' (see efforts Bas_Boukens/Human_fetal/Analysis_test)

# FSS 2.0
# <ProjectName> - Processing - DataAnalysisX - Code                        - Meta      - targets.txt, BroadSets.rds etc.
#                                            - Data ------------ Dataset_1 - Processed - .rds objects/ count tables
#                                            - Documentation - ExpInfo.txt/biomaRt versions etc
#                                            - Notebooks
#                                            - Results --- Figures/Tables/SummaryTables
#                                            - Settings
# 
#                                               - Meta
#                            - Data - Dataset_1 - Processed
#                                               - Raw
#
#               - Sharing
#

# 1. Within 'Sharing' make directory with name of Analysis directory in 'Processing'
# 2. make directories to hold scripts, results, r objects
# 3. within 'results' directory make directory for figures, tables etc.

analysisToShare <- function(projectDir, dirToShare, dataSetToShare = "Dataset_1", analysisDataToShare = "Dataset_1"){
  # projectDir - the main directory, i.e. <ProjectName>
  # dirToShare - the analysis you want to share, i.e. <ProjectName>/Processing/DataAnalysisX or just DataAnalysisX
  #              the first way allows you to just follow the path... less error-prone
  
  # You should be in the directory 'above' the <YYYYMMDD_ProjectName> directory
  baseDir <- getwd()  
  if ( !(gsub("\\/","",projectDir) %in% dir(baseDir))){
    print(paste0("Make sure you are in a directory one level above the Project directory!!"))
    break()
  }
  
  # Make sure to only get the name of the directory you want to share
  dirToShare <- basename(dirToShare)
  
  # Create the directory that is going to hold all information
  dir.create(paste0(baseDir,"/",projectDir,"/Sharing/",dirToShare), showWarnings = FALSE)
  shareDir <- paste0(baseDir,"/",projectDir,"/Sharing/",dirToShare)
  analysisDir <- paste0(baseDir,"/", projectDir,"/Processing/",dirToShare)
  
  # Create the directories holding the scripts, results, Robjects
  for (i in c("Scripts","CountTables","Data","R_Objects",
              "ResultTables","VersionInfo","Figures", "QualityChecks")){
    dir.create(paste0(projectDir,"/Sharing/",dirToShare,"/",i), showWarnings = FALSE)
  }
  
  # Go to the analysis directory 
  setwd(analysisDir)
  
  #### And copy
  # scripts
  toBeCopied <- list.files("Code", pattern="*.r|*.R", full.names = TRUE)
  file.copy(toBeCopied,to=paste0(shareDir,"/Scripts"))
  cat("... scripts copied\n")
  # CountTables
  toBeCopied <- list.files(paste0(analysisDir,"/Data/", analysisDataToShare,"/Processed"), pattern="*.tab", full.names = TRUE)
  file.copy(toBeCopied, to=paste0(shareDir,"/CountTables"))
  cat("... count table copied\n")
  # Data - FastQC,dupRadar
  file.copy(paste0("../../Data/", dataSetToShare,"/Processed/"), to=paste0(shareDir,"/Data"), recursive = TRUE)
  cat("... fastqc, dupradar data copied\n")
  # R_objects
  toBeCopied <- list.files(paste0(analysisDir,"/Data/", analysisDataToShare,"/Processed"), pattern="*.rds", full.names = TRUE)
  # If you have multiple datasets and/or with different names:
  # list.files(paste0(list.dirs(paste0(analysisDir,"/Data/"), recursive = FALSE),"/Processed"), pattern=".rds", full.names = TRUE)
  file.copy(toBeCopied, to=paste0(shareDir,"/R_objects"))
  toBeCopied <- list.files(paste0(analysisDir,"/Data/", analysisDataToShare,"/Meta"), pattern="*.rds", full.names = TRUE)
  file.copy(toBeCopied, to=paste0(shareDir,"/R_objects"))
  cat("... R objects copied\n")
  # GSEA_Reports
  if (dir.exists(paste0(analysisDir,"/Results/Tables/reports"))){
    # Filenames were too long, so I first move to the directory and then back again ...
    pwd <- getwd()
    setwd(paste0(analysisDir,"/Results/Tables/reports"))
    #R.utils::copyDirectory(paste0(analysisDir,"/Results/SummaryTables/reports"), paste0(shareDir,"/GSEA_Reports"))
    R.utils::copyDirectory(".", paste0(shareDir,"/GSEA_Reports"))
    
    # Go back to the main working directory
    setwd(pwd)
  }
  cat("... GSEA results copied\n")
  #file.copy(paste0(analysisDir,"/Results/Tables/reports"), to=shareDir, recursive = TRUE)
  #Sys.sleep(30)
  #file.rename(from=paste0(shareDir,"/reports"), to=paste0(shareDir,"/GSEA_Reports"))
  # Venn_Genes
  file.copy(paste0(analysisDir,"/Results/Figures/venn"), to=shareDir, recursive = TRUE)
  file.rename(from=paste0(shareDir,"/venn"), to=paste0(shareDir,"/Venn_Genes"))
  cat("... Venn Genes copied\n")
  # Venn_GeneSets
  file.copy(paste0(analysisDir,"/Results/Figures/venn_gs"), to=shareDir, recursive = TRUE)
  file.rename(from=paste0(shareDir,"/venn_gs"), to=paste0(shareDir,"/Venn_GeneSets"))
  cat("... Venn Genesets copied\n")
  # ResultTables
  toBeCopied <- list.files(paste0(analysisDir,"/Results/Tables"), pattern="*.tab|*.txt|*.csv|*.pptx", full.names=TRUE)
  file.copy(toBeCopied, to=paste0(shareDir,"/ResultTables"))
  # If we performed ROMER, just leave it in the ResultsTables directory
  if (dir.exists(paste0(analysisDir,"/Results/SummaryTables/romerReports"))){
    file.copy(paste0(analysisDir,"/Results/SummaryTables/romerReports"), to=paste0(shareDir,"/ResultTables/"), recursive = TRUE)
  }
  toBeCopied <- list.files(paste0(analysisDir,"/Results/Tables"), pattern="*.rds", full.names=TRUE)
  file.copy(toBeCopied, to=paste0(shareDir,"/R_Objects"))
  cat("... Result tables copied\n")
  # VersionInfo
  toBeCopied <- list.files(paste0(analysisDir,"/Documentation"), pattern="*.versions", full.names=TRUE)
  file.copy(toBeCopied, to=paste0(shareDir,"/VersionInfo"))
  file.copy(paste0(analysisDir,"/Documentation/sessionInfo.txt"), to=paste0(shareDir,"/VersionInfo"))
  file.copy(paste0(analysisDir,"/Documentation/ExpInfo.txt"), to=shareDir)
  cat("... VersionInfo copied\n")
  # Voom_Reports
  file.copy(paste0(analysisDir,"/Results/Figures/voom"), to=shareDir, recursive = TRUE)
  file.rename(from=paste0(shareDir,"/voom"), to=paste0(shareDir,"/Voom_Reports"))
  cat("... Voom reports copied\n")
  # Figures
  toBeCopied <- list.files(paste0(analysisDir,"/Results/Figures"), pattern="*.png|*.tiff|*.html", full.names=TRUE)
  file.copy(toBeCopied, to=paste0(shareDir,"/Figures"))
  dirToCopy <- list.dirs(paste0(analysisDir,"/Results/Figures"), recursive = FALSE)
  dirToCopy <- dirToCopy[!(grepl("venn$|venn_gs$|voom$", dirToCopy)) ]
  file.copy(dirToCopy, to=paste0(shareDir,"/Figures"), recursive=TRUE)
  cat("... Figures copied\n")
  
  # QC_FastQC, QC_dupRadar ?? Separate QC folder? But do I need to put all dupRadar stats/htmls/pdfs here as well then...?
  # Or put everything into figures ...
  toBeCopied <- list.files(paste0(analysisDir,"/Results/SummaryTables"), pattern="QC*|*fastqc*|*.html", full.names=TRUE)
  file.copy(toBeCopied, to=paste0(shareDir,"/QualityChecks"), recursive=TRUE)
  cat("... QC directories copied\n")
  
  # Copy the 'targets' file to the right directory
  file.copy(paste0(analysisDir,"/Data/", analysisDataToShare,"/Meta/targets.txt"), shareDir)
  cat("... Meta data copied\n")
  
  #### Alter scripts with respect to where they should store/get initial data
  cat("\n... Now adapting the paths in the scripts\n")
  
  file.copy(from = paste0(shareDir,"/Scripts/analysis.r"), to=shareDir)
  file.remove(paste0(shareDir,"/Scripts/analysis.r"))
  # Perhaps make a 'codeDir'/'scriptDir' variable in the beginning of the analysis.r, so I can 
  # change that easily ... ?
  
  # Replace the paths in the 'analysis.r' script
  f <- list.files(path=shareDir, pattern = "analysis.r", full.names=TRUE)
  print(f)
  x <- readLines(f)
  x <- gsub("^dropbox","# dropbox", x)
  x <- gsub("Code", dirToShare, x)
  idx <- which(grepl("^scriptDir",x))
  x[idx] <- "scriptDir =\"./Scripts\""
  x <- gsub("source\\(\"", "source\\(\"./Scripts/", x)
  x <- gsub("../Results/Figures", "./Figures", x)
  x <- gsub("\\.\\./Results/SummaryTables", "./ResultTables", x)
  x <- gsub("\\./Results/SummaryTables", "./ResultTables", x)
  x <- gsub(paste0("../Data/", analysisDataToShare,"\\/Processed"), "./R_Objects", x)
  x <- gsub("\\.\\./Results/Tables", "./ResultTables", x)
  idx <- which(grepl("^docDir",x))
  x[idx] <- gsub("=\"\\.\\./Documentation", " = \"./VersionInfo\"\nqualDir = \"./QualityChecks\"\ncountDir = \"./CountTables", x[idx])
  idx <- which(grepl("sessionDir",x))
  x[idx] <- gsub("\\.\\./Documentation", "./VersionInfo", x[idx])
  x <- gsub("docDir)","docDir)\n  dir.create(qualDir)\n  dir.create(countDir)", x )
  cat(x, file=f, sep="\n")

  ## Replace directory paths in all scripts we have used
  scriptDir <- "./Scripts"
  
  filenames <- list.files(path=paste0(shareDir,"/Scripts"), pattern = "^0", full.names=TRUE)
  for( f in filenames ){
    x <- readLines(f)
    
    # General: set 'scriptDir' to the right directory, i.e. "." or "./Scripts ?
    x <- gsub("paste\\(scriptDir", paste0("paste(\"" ,scriptDir,"\""),x)
    
    # Substitute
    x <- gsub("Code",dirToShare, x)  # Check in 01_preprocessing to assure we're in the base directory
    # Set the dropbox variable to the location of the scripts
    idx <- which(grepl("^dropbox",x))
    x[idx] <- "dropbox = './Scripts'"
    # Outcomment all references to dropbox in combination with file.copy.
    idx <- which(grepl("dropbox", x) & grepl("file.copy",x))
    x[idx] <- paste0("# ", x[idx])
    # Check whether there is a directory in Figure with fastqc or dupRadar in the name...
    idx <- which((grepl("[Ff]ast[Qq][Cc]", x) | (grepl("dupRadar",x))) & grepl(paste0("../../../Data/", dataSetToShare,"/Processed"),x))
    # And then replace
    x[idx] <- gsub(paste0("../../../Data/", dataSetToShare,"/Processed"), "./Data/Processed/", x[idx])
    #
    idx <- which((grepl("[Ff]ast[Qq][Cc]", x) | (grepl("dupRadar",x))) & grepl("../Results/SummaryTables",x))
    x[idx] <- gsub("../Results/SummaryTables", "./QualityChecks", x[idx])
    # Check whether there is a directory with the counts
    idx <- which(grepl("*count", x) & grepl(paste0("../../../Data/", dataSetToShare,"/Processed"),x))
    x[idx] <- gsub(paste0("../../../Data/", dataSetToShare,"/Processed"),"./Data/Processed", x[idx])
    idx <- which(grepl("*count|allCounts", x) & grepl(paste0("../Data/", analysisDataToShare,"/Processed"),x))
    x[idx] <- gsub(paste0("../Data/", analysisDataToShare,"/Processed"),"./CountTables/", x[idx])
    idx <- which(grepl(".rds", x))
    x[idx] <- gsub(paste0("../Data/", analysisDataToShare,"/Processed"), "./R_Objects/", x[idx]) 
    idx <- which(grepl("D_chrom.rds", x))
    x[(idx+1)] <- gsub(paste0("../Data/", analysisDataToShare,"/Processed"), "./R_Objects/", x[(idx+1)]) 
    # 
    idx <- which(grepl(".rds", x) & grepl("tablesDir", x))
    x[idx] <- gsub("tablesDir", "\"./R_Objects\"", x[idx]) 
    idx <- which(grepl("write.table", x) & grepl("tablesDir", x))
    x[idx] <- gsub("tablesDir", "\"./ResultTables\"", x[idx])
    idx <- which(grepl(".rds", x) & grepl("../Results/Tables/",x) & !(grepl("tabCameraFile",x)))
    x[idx] <- gsub("../Results/Tables/", "./R_Objects/", x[idx]) 

    # Outcomment the construction of the BroadSets
    if ( grepl("06_", f) ){
      idx1 <- which(grepl("hBroadSets", x))
      idx2 <- which(grepl("saveRDS\\(BroadSets", x))
      for ( i in seq(idx1[1],idx2) ){
        x[i] <- paste0("# ", x[i])
      }
      idx <- which(grepl("readRDS\\(\"BroadSet", x))
      x[idx] <- gsub("readRDS\\(\"BroadSet","readRDS\\(\"\\./R_Objects/BroadSet", x[idx])
    }
    
    # Outcomment dupRadar analysis in 01_preprocessing.r if the directory is not there, but the code is in (should have been removed :-) )
    # Just happens once in a while ..
    if ( grepl("01_", f) & !grepl("dupRadar", list.dirs(paste0(shareDir,"/QualityChecks"))) ){
      idx1 <- which(grepl("dupRadar", x))
      idx2 <- which(grepl("The results from the dupRadar QC analysis", x))
      for ( i in seq(idx1[1],(idx2+1)) ){
        x[i] <- paste0("# ", x[i])
      }
    }
    
    if ( grepl("05_|07_", f) ){
      idx <- which(grepl(".rds", x))
      x[(idx+1)] <- gsub(paste0("../Data/", analysisDataToShare,"/Processed"), "./R_Objects", x[(idx+1)]) 
      x <- gsub("scriptDir(?<=\\S)\\s{0,4}(?=\\S)\\=(?<=\\S)\\s{0,4}(?=\\S)\".\"", "scriptDir\\=\"\\./Scripts\"",x, perl=TRUE)
    }
    
    # And the more general substitutions
    if ( grepl("08_",f) ){
      print(paste0("Adapting ",f))
      # x <- gsub("paste\\(scriptDir", paste0("paste(\"" ,scriptDir,"\""),x)
      x <- gsub("VersionInfo/ExpInfo\\.txt", "./ExpInfo.txt",x)
      x <- gsub("fastqc=\"./Figures", "fastqc=\"./QualityChecks",x)
      x <- gsub("tabCameraFile=\"[A-Za-z0-9=_/\\.]*reports", "tabCameraFile=\"./GSEA_Reports",x)
      x <- gsub(paste0("../Data/", analysisDataToShare,"/Meta"), "./R_Objects",x)
      x <- gsub(paste0("../Data/", analysisDataToShare,"/Processed"), "./R_Objects",x)
      x <- gsub("../Results/Figures", "./Figures",x)
      x <- gsub("../Results/Tables", "./R_Objects",x)
      x <- gsub("../Results/Figures/venn_gs[\\/]?", "./Venn_GeneSets",x)
      x <- gsub("../Results/Figures/venn[\\/]?", "./Venn_Genes",x)
      x <- gsub("scriptDir[A-Za-z0-9 _=/.\"]*",paste0("scriptDir \\= \"",scriptDir,"\""), x)
    } else {
      x <- gsub("../Results/Tables", "./ResultTables",x)
      x <- gsub("../Results/Figures", "./Figures",x)
      
    }
    
    if ( grepl("03_|04_|06_", f) ){
      x <- gsub(paste0("../Data/", analysisDataToShare,"/Meta/"),"./R_Objects/", x)
      x <- gsub(paste0("../Data/", analysisDataToShare,"/Meta"),"./R_Objects/", x)
    } else {
      x <- gsub(paste0("../Data/", analysisDataToShare,"/Meta/"),"", x)
      x <- gsub(paste0("../Data/", analysisDataToShare,"/Meta"),"", x)
    }
    
    x <- gsub("../Documentation", "./VersionInfo", x)
    x <- gsub("VersionInfo/ExpInfo\\.txt", "./ExpInfo.txt",x)
    x <- gsub("ResultTables/reports", "./GSEA_Reports", x )
    x <- gsub("./Figures/venn_gs", "./Venn_GeneSets",x)
    x <- gsub("./Figures/venn", "./Venn_Genes",x)
    
    cat(x, file=f, sep="\n")
  }
  
  f <- list.files(path=shareDir, pattern = "ExpInfo.txt", full.names=TRUE)
  x <- readLines(f)
  x <- gsub("../Results/Figures", "Figures", x)
  x <- gsub("../Results/SummaryTables", "ResultTables", x)
  x <- gsub(paste0("../Data/", analysisDataToShare,"/Processed"), "CountTables", x)
  x <- gsub("../Results/Tables", "ResultTables", x)
  x <- gsub("../Documentation/", "VersionInfo/", x)
  cat(x, file=f, sep="\n")
  
  
  # And execute the writeShiny scripts to get the new 'server.r' and 'ui.r'
  setwd(shareDir)
  source("./Scripts/08_shiny_app.r")
  
  # In the process of rerunning the scripts (again), the ExpInfo.txt file will be
  # relocated to ./ResultTables using 01_preprocessing.r
  # So, change 08_shiny_app.r to put it back ... :-)
  f <- "./Scripts/08_shiny_app.r"
  x <- readLines(f)
  x <- gsub("writeUI_RNASeq", "# First put ExpInfo.txt back to the main directory if it is in \\./ResultTables:
if (!file.exists(\"ExpInfo.txt\")){
  file.copy(\"\\./ResultTables/ExpInfo.txt\", \"\\.\")
  file.remove(\"\\./ResultTables/ExpInfo.txt\")
}\n\n
writeUI_RNASeq",x)
  # Make sure to get some location right
  x <- gsub(paste0("../Data/", analysisDataToShare,"/Processed"), "./R_Objects/", x)
  x <- gsub(paste0("../Data/", analysisDataToShare,"/Meta"), "./R_Objects/", x)
  cat(x, file=f, sep="\n")
  
  # And remove the original 'server.r', 'ui.r' and 'runShinyApp.r' ??
  file.remove(paste0(shareDir,"/Scripts/server.r"))
  file.remove(paste0(shareDir,"/Scripts/ui.r"))
  file.remove(paste0(shareDir,"/Scripts/runShinyApp.r"))
  
  # And remove some empty directories?
  
  # Fix/Lock the directory with the original analysis
  #Sys.chmod(analysisDir, mode = "0555", use_umask = TRUE)
  #system(paste0('chmod -R 0555 ', analysisDir))
  
  # Change some locations in all the helper scripts
  filenames <- list.files(path=paste0(shareDir,"/Scripts"), pattern = "^[a-zA-Z]", full.names=TRUE)
  for( f in filenames ){
    x <- readLines(f)
    if(f == paste0(shareDir,"/Scripts/","annotateD.r")){
      x <- gsub("annotDir=NULL", "annotDir=\"./R_Objects\"", x)
      # And prevent overwriting the biomart/ensembl version file !!
      x <- gsub("write.table", "# write.table ", x)
    }
    x <- gsub("dropbox,\"/Support/R", "\"\\.", x)
    x <- gsub("dropbox,\"Support/R", "\"\\.", x)
    
    # Outcomment dropbox location ..
    x <- gsub("dropbox =", "# dropbox =", x)
    
    idx <- grepl("multiFastQC/FastQC", x)
    x[idx] <- gsub("file.copy", "# file.copy", x[idx])
    x <- gsub("source\\(\"FastQC\\-class\\.R\")", "source(\"\\./Scripts\\/FastQC-class.R\")", x)
    x <- gsub("source\\(\"FastQC_v0\\.10\\-class\\.R\")", "source(\"\\./Scripts\\/FastQC_v0.10-class.R\")", x)
    
    # cameraReports.r has code to copy the MSigDB.Categories object, but it is already there ..
    if(f == paste0(shareDir,"/Scripts/","cameraReports.r")){
      idx <- which(grepl("/Support/MSigDB/v", x))
      for ( i in seq((idx[1]-3), (idx+4)) ){
        x[i] <- paste0("# ", x[i])
      }
      x <- gsub("\\./barcodeplot.R", "\\./Scripts/barcodeplot.R", x)
    }
    
    cat(x, file=f, sep="\n")
  }
  
  setwd(baseDir)
  
  cat("\n... All done!!\n")
  
}