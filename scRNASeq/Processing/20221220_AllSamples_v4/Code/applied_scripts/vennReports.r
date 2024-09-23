vennReports <- function(fit=fit, contrasts=my.contrasts, indices=1, myTitle="vennDiagram_constrast_1",
                        reportsDir="venn", p.value=0.0001, adjust.method="none", lfc=0,
                        col=c("red","blue"), trend=FALSE, robust=FALSE, scriptDir="./applied_scripts",... ){
  
  # AJ - 01092017
  
  # Load the neccessary packages
  if (!require("pacman")){ install.packages("pacman") }
  pacman::p_load(checkmate)
  require(limma)
  require(edgeR)

  # set the Dropbox location for additional scripts
  dropbox = ifelse(Sys.info()['sysname'][[1]] == "Darwin", "/Users/Aldo/Desktop/AMC/Dropbox", ifelse(Sys.info()['sysname'][[1]] == "Linux", "/home/ajongejan/Dropbox", "D:/Dropbox"))
  
  # Copy if the file isn't there and source code for makeClickableVenn
  if(!file.exists(paste(scriptDir,"makeClickableVenn.R",sep="/"))){
    file.copy(paste(dropbox,"Support/R/makeClickableVenn.R",sep="/"),scriptDir)
  }
  source(paste(scriptDir,"makeClickableVenn.R",sep="/"))
  
  adjustT <- ifelse(adjust.method=="BH","Benjamini-Hochberg","none")
  pValString = paste0("p-val. cut-off ",p.value,"(",adjustT,"), logFC cut-off ",lfc)
  
  # In case, no name for the picture was provided
  if (myTitle == "vennDiagram_constrast_1"){ 
    myTitle = paste0("vennDiagram_contrasts_",paste(attr(contrasts,"dimnames")$Contrasts[indices], collapse="_"))
  } 

  # Test fit object 
  if (grepl("DGEGLM|DGELRT",class(fit)) && attr(class(fit),"package")=="edgeR"){
    # edgeR
    for (i in indices){
      contrast = attr(contrasts,"dimnames")$Contrasts[i]
      if (!is.null(fit$df.residual.zeros)) {
        lrt = glmQLFTest(fit,contrast=contrasts[,i])
      } else {
        lrt = glmLRT(fit,contrast=contrasts[,i])
      }
      tt.i = topTags(lrt,n=Inf,sort="none")$table
      indx_logFC = which( colnames(tt.i)=="logFC" )
      indx_logCPM = which( colnames(tt.i)=="logCPM" )
      tt.names = names(tt.i)[c(1:(indx_logFC-1),indx_logCPM)]
      names(tt.i) = paste0(names(tt.i)," (",attr(contrasts,"dimnames")$Contrasts[i],")")
      names(tt.i)[c(1:(indx_logFC - 1),indx_logCPM)] = tt.names
      deTags <- decideTestsDGE(lrt,adjust.method=adjust.method, p.value = p.value, lfc = lfc)
      colnames(deTags) <- c(contrast)
      if(i == indices[1]){ testMatrix <- deTags} else {testMatrix <- cbind(testMatrix,deTags) }
      if(i == indices[1]){ tt <- tt.i[,c(1:(indx_logFC-1),indx_logCPM,indx_logFC,(indx_logCPM+1):dim(tt.i)[2])]} 
      else { tt <- cbind(tt,tt.i[,c(indx_logFC,(indx_logCPM+1):dim(tt.i)[2])]) }
    }
  } else if (class(fit)[1]=="MArrayLM" && attr(class(fit),"package")=="limma"){
    # Voom
    for (i in indices){
      contrast = attr(contrasts,"dimnames")$Contrasts[i]
      fit2 = contrasts.fit(fit,contrast=contrasts[,i])
      fit.eb = eBayes(fit2, trend=trend, robust=robust)
      tt.i = topTable(fit.eb,n=Inf,sort="none")
      indx_logFC = which( colnames(tt.i)=="logFC" )
      indx_AveExpr = which( colnames(tt.i)=="AveExpr" )
      tt.names = names(tt.i)[c(1:(indx_logFC-1),indx_AveExpr)]
      names(tt.i) = paste0(names(tt.i)," (",attr(contrasts,"dimnames")$Contrasts[i],")")
      names(tt.i)[c(1:(indx_logFC - 1),indx_AveExpr)] = tt.names
      deTags <- decideTests(fit.eb, adjust.method=adjust.method, p.value = p.value, lfc = lfc)
      colnames(deTags) <- c(contrast)
      if(i == indices[1]){ testMatrix <- deTags} else {testMatrix <- cbind(testMatrix,deTags) }
      if(i == indices[1]){ tt <- tt.i[,c(1:(indx_logFC-1),indx_AveExpr,indx_logFC,(indx_AveExpr+1):dim(tt.i)[2])]} 
      else { tt <- cbind(tt,tt.i[,c(indx_logFC,(indx_AveExpr+1):dim(tt.i)[2])]) }
    }
  } else {
    # No fit object
    stop("You need a fit object for this function to work ... \n")
  }

  FDR <- ifelse(adjustT == "none", "none", "BH")
  mySuffix <- paste0(gsub("\\.","_",p.value),"_",gsub("\\.","_",lfc),"_",FDR)
  
  png(paste0(paste0(myTitle,"_", mySuffix,".png")), width=960, height=960, pointsize=20)
  vennDiagram( testMatrix[,1:length(indices)], include=c("up","down"),counts.col=col, cex = 1.4, names=colnames(testMatrix) )
  mtext(paste0("p.value = ",p.value,"(",adjustT,"), logFC = ",lfc),side=1,line=3,cex=1.2, adj=0)
  dev.off()

  colnames(testMatrix) <- paste(indices,colnames(testMatrix), sep="_")
  makeClickableVenn(testMatrix, tt, reportDir=reportsDir, pValString = pValString, myCol=col)

}