vennReportsGeneSets <- function(tab="tt.genesets.rds", filePrefix="./tab.camera.v6.1_", 
                                contrasts=my.contrasts, indices=1, myTitle="vennDiagram_genesets_1",
                                reportsDir="venn_gs", p.value=NULL, fdr=NULL, 
                                col=c("red","blue"), scriptDir="./applied_scripts",... ){
  
  # AJ - 02092018
  
  # Load the neccessary packages
  require(limma)
  require(edgeR)
  
  # set the Dropbox location for additional scripts
  dropbox = ifelse(Sys.info()['sysname'][[1]] == "Darwin", "/Users/Aldo/Desktop/AMC/Dropbox", ifelse(Sys.info()['sysname'][[1]] == "Linux", "/home/ajongejan/Dropbox", "D:/Dropbox"))
  
  # Copy if the file isn't there and source code for makeClickableVenn
  if(!file.exists(paste(scriptDir,"makeClickableVenn.R",sep="/"))){
    file.copy(paste(dropbox,"Support/R/makeClickableVenn.R",sep="/"),scriptDir)
  }
  source(paste(scriptDir,"makeClickableVenn.R",sep="/"))
  
  # Test whether a cut-off has been specified
  if(is.null(p.value) & is.null(fdr)){
    stop("You need to specify either a p-value or a FDR cut-off! \n")
  }

  # If both are specified, the FDR will be taken..
  if(!is.null(p.value) & !is.null(fdr)){
    warning("You specified both a p-value and a FDR cut-off, the FDR cut-off will be used with the specified p-value! \n")
#    p.value=NULL
    fdr=p.value
  }
  
  # Create the caption for the picture
  pValString = ifelse( (!is.null(p.value) | is.null(p.value)) & !is.null(fdr),paste0("FDR cut-off ", p.value, " (BH)" ), 
                      ifelse(is.null(fdr), paste0("p-val. cut-off ",p.value), ""))
  cat(pValString,"\n")
  
  # In case, no name for the picture was provided
  if (myTitle == "vennDiagram_genesets_1"){ 
    myTitle = paste0("vennDiagram_genesets_",paste(attr(contrasts,"dimnames")$Contrasts[indices], collapse="_"))
  } 
  
  # Make the complete table holding the data for the genesets for all contrasts
  if (!exists(tab)){
    for (i in 1:length(attr(contrasts,"dimnames")$Contrasts)){
      contrast = attr(contrasts,"dimnames")$Contrasts[i]
      
      tab.i = read.delim(paste0(filePrefix,contrast,".txt"))
      names(tab.i) = paste0(names(tab.i)," (",attr(contrasts,"dimnames")$Contrasts[i],")")
      names(tab.i)[1:2] <- c("GeneSet", "NGenes")
      if (i == 1){
        tab.all <- tab.i
      } else {
        tab.i$NGenes <- NULL
        tab.all <- merge(tab.all, tab.i, by="GeneSet")
      }
    }
    
    rownames(tab.all) <- tab.all$GeneSet
    
    # And save for future use
    saveRDS(tab.all, "tt.genesets.rds")
  } else {
    # We already have an object holding all information on the Gene sets
    tab.all <- readRDS(tab)
  }  
  
  # Now generate the Venn diagram for the specified combination of contrasts
  for (i in indices){
    contrast = attr(contrasts,"dimnames")$Contrasts[i]
    
    if (!is.null(p.value) & is.null(fdr)){
      cat(c(paste0("PValue (", contrast,")"),paste0("Direction (", contrast,")")),"\n")
      
      tab.i = tab.all[,c(paste0("PValue (", contrast,")"),paste0("Direction (", contrast,")"))]
    
      deTags = apply(tab.i, 1, function(x) ifelse(as.numeric(x[paste0("PValue (", contrast,")")]) <= p.value, 
                                                ifelse(x[paste0("Direction (", contrast,")")] == "Up", 1, -1 ), 0))
    } else if (!is.null(fdr) ){
      cat(c(paste0("FDR (", contrast,")"),paste0("Direction (", contrast,")")),"\n")
      
      tab.i = tab.all[,c(paste0("FDR (", contrast,")"),paste0("Direction (", contrast,")"))]
      
      deTags = apply(tab.i, 1, function(x) ifelse(as.numeric(x[paste0("FDR (", contrast,")")]) <= p.value, 
                                                  ifelse(x[paste0("Direction (", contrast,")")] == "Up", 1, -1 ), 0))
    }
    
    deTags = data.frame(deTags)
    if ( nrow(deTags) > 0) {
      rownames(deTags) <- rownames(tab.i)
    }
    colnames(deTags) <- c(contrast)
    #    
    if(i == indices[1]){ 
      testMatrix <- deTags 
    } else { 
      cat(all(match(rownames(testMatrix), rownames(deTags))),"\n")
      testMatrix <- cbind(testMatrix,deTags) }
  }
  
  mySuffix <- ifelse(!is.null(fdr), paste0(gsub("\\.","_",p.value),"_BH"), paste0(gsub("\\.","_",p.value), "_none"))
  png(paste0(paste0(myTitle,"_", mySuffix,".png")), width=960, height=960, pointsize=20)
  vennDiagram( testMatrix[,1:length(indices)], include=c("up","down"),counts.col=col, cex = 1.4, names=colnames(testMatrix) )
  if(!is.null(p.value) & !is.null(fdr) ){
    mtext(paste0("p.value = ",p.value, ", fdr = ",fdr),side=1,line=3,cex=1.2, adj=0)
  } else if (is.null(fdr)){
    mtext(paste0("p.value = ",p.value),side=1,line=3,cex=1.2, adj=0)
  } else if (is.null(p.value) & !is.null(fdr) ){
    mtext(paste0("p.value = ",p.value, ", fdr = ",fdr),side=1,line=3,cex=1.2, adj=0)
  }
  dev.off()
  
  colnames(testMatrix) <- paste(indices,colnames(testMatrix), sep="_")
  makeClickableVenn(as.matrix(testMatrix), tab.all, reportDir=reportsDir, pValString = pValString, myCol=col)
  
}