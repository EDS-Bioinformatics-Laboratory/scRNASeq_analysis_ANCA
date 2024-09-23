startNewRepo <- function(dirName, analysisDir="DataAnalysis"){
  # where are we?
  baseDir <- getwd()
  
  # Create the new directory
  dir.create(dirName)
  
  # Cloning a complete repository will work
  system('git clone https://github.com/EDS-Bioinformatics-Laboratory/Reproducibility Reproducibility')
  
  # And we could now restructure it, i.e. remove the top layer by moving the 'ddmmyy_ProjectName' directory to our new ProjectName
  
  # https://www.datanovia.com/en/blog/how-to-easily-manipulate-files-and-directories-in-r/
  #library(fs)
  
  file.copy("./Reproducibility/ExampleFolderStructure/yyyymmdd_ProjectName/.",dirName, recursive = TRUE)
  
  # And get to the source directory of the project
  print("Now moving into the new directory ...")
  setwd(dirName)
  
  # Now move the Data to the right directory in the FSS, i.e. toplevel 'Data'
  if ( dir.exists("../Data") ){
    dataToMove <- list.dirs('../Data', full.names = FALSE, recursive = FALSE)
    setwd("../Data")
    
    for (f in dataToMove){
      folder_old_path = "../Data"
      path_new = paste0("../", dirName,"/Data/Dataset_1/Processed/")
      file.copy(from = folder_old_path, to = path_new, 
                overwrite = TRUE, recursive = TRUE, copy.mode = TRUE)
    }
    setwd(paste0("../", dirName,"/"))
    
    # And remove the original Data folder
    system('rm -rf ../Data')
  }
  
  # Do your analysis
  setwd("Processing/")
  analysisDir <- paste(format(Sys.time(), "%Y%m%d"),analysisDir,sep="_")
  file.rename("yyyymmdd_NameOfDataAnalysis/", analysisDir)
  setwd(paste0(analysisDir,"/Code"))
  
  # Remove the original 'Reproducbility' directory
  succes = 1
  # Remove the original 'Reproducibility' directory
  succes = 1
  while (succes > 0){
    unlink("Reproducibility", recursive=TRUE, force=TRUE)
    if (!("Reproducibility" %in% dir())) {
      print("Removed the initial Reproducibility folder ...")
      succes = 0
    }
  }
  
  
}

startNewAnalysis <- function(analysisDir="DataAnalysis"){
  # Make sure you're in the 'Processing' directory within the main projectDir
  if ( basename(getwd()) != "Processing"){
    print("To create a new analysis directory, you should be in the \'Processing\' directory ...")
    # We assume that the directory structure follows our FSS .. :-)
    # Should we not just create all of this??
    break()
  }
  
  # Clone a complete FSS structure into this directory
  # AJ - 20230117 - Had to add this system('chmod ugo+rwx .') as I was not permitted to write Reproducibility/.git/config
  #    Cloning into 'Reproducibility'...
  #    error: could not write config file D:/Dropbox/Support/Yosta_Vegting/RNASeq/Processing/Reproducibility/.git/config: Permission denied
  #    fatal: could not set 'core.repositoryformatversion' to '0'
  system('chmod ugo+rwx .')
  system('git clone https://github.com/EDS-Bioinformatics-Laboratory/Reproducibility Reproducibility')
  
  # Add the time in the right format
  analysisDir <- paste(format(Sys.time(), "%Y%m%d"),analysisDir,sep="_")
  file.rename("./Reproducibility/ExampleFolderStructure/yyyymmdd_ProjectName/Processing/yyyymmdd_NameOfDataAnalysis",analysisDir)
  
  # And remove the repo
  # Remove the original 'Reproducibility' directory
  succes = 1
  while (succes > 0){
    unlink("Reproducibility", recursive=TRUE, force=TRUE)
    if (!("Reproducibility" %in% dir())) {
      print("Removed the initial Reproducibility folder ...")
      succes = 0
    }
  }

  # And get to the source directory of the project
  setwd(paste0(analysisDir, "/Code"))
  
  # And give some information
  if (succes == 0 & basename(getwd()) == "Code"){
    print("The directory has been succesfully created and the working directory has been set to the \'Code\' directory!!")
  }
  
  
}


pushToRepo <- function(branchName="supportProject", branchComment="A new support project has been finished!"){
  # And push everything to GitHub
  # Have to think about the organization:
  # - I could create a 'Support_RNASeq' repository and then push everything into that one 
  #   (just make the repository once, all support projects will be a seperate directory)
  #   DECIDED with BioLab on 20210218!!
  #   Put every projct in a separate branch!!  
  # - Or create a a separate repository for every project...
  #
  # My preference is for the first one..
  
  # Where do I put my .git?
  # How can I name the directory on GitHub under the 'Support_RNASeq' level?
  # - i.e. if I put the .git in the directory above 'Code', then every support project is going to be called 'Code' ...
  # Only put the R code into GitHub, so get a .gitignore that filters everything else
  
  # First clone the Support_RNASeq repo
  system('git clone https://github.com/EDS-Bioinformatics-Laboratory/Support_RNASeq D://Data/Dropbox/Support_RNASeq_GitHub')
  
  fileConn <- file("./Code/.gitignore")
  myRules <- paste0("# Blacklist files/folders in same directory as the .gitignore file\n/*\n",
                    "# Whitelist some files\n!.gitignore\n!*.r\n!*.R\n")
  writeLines(myRules, fileConn)
  close(fileConn)
  
  system('git init')
  system('git add Code')
  
  system(paste0('git commit -m ', paste0(branchComment, collapse=" ")))
  
  system(paste0('git checkout -b ',branchName))
  # Create a new branch within the repo on GitHub
  # - but you don't want to add the Personal Access Token hardcoded!!
  # So, see: https://www.r-bloggers.com/2020/07/a-better-way-to-manage-your-github-personal-access-tokens/
  library(credentials)
  if ( credentials::git_credential_ask("https://github.com")$username == "PersonalAccessToken" ){
    myPAT <- credentials::git_credential_ask("https://github.com")$password
  }
  system(paste0("curl -H \'Authorization: token ", myPAT,"\' https://api.github.com/user/repos -d \'{\"name\":\"Support_RNASeq.git\"}\'"))
  #
  system('git remote add origin https://github.com/aldojongejan/Support_RNASeq.git')
  system(paste0('git push -u origin ', branchName))
  system(paste0('git push origin ' ,branchName))
}