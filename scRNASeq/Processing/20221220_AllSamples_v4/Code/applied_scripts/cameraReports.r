cameraReport <- function(v.=v, fit.=fit, sets.=sets, design.=design, contrasts=my.contrasts, wantedContrast="all",
                         species="human", 
                         reportsDir="reports", interGeneCor=0.01, p.value=0.0001, max.pathways=250, msigdb.version="5.1",
                         trend=trend, robust=robust,... ){

  # AJ - 13092016
  
  # Load the neccessary packages
  require(limma)
  require(hwriter)
  require(ReportingTools)
  
  # set the Dropbox location for additional scripts
  dropbox = ifelse(Sys.info()['sysname'][[1]] == "Darwin", "/Users/Aldo/Desktop/AMC/Dropbox", ifelse(Sys.info()['sysname'][[1]] == "Linux", "/home/ajongejan/Dropbox", "D:/Dropbox"))

  # Copy the R object holding the MSigDB categories to the directory
  if (!file.exists(paste0("MSigDB_Categories.v", msigdb.version,".rds"))){
    # Try to get it from the DropBox location (only works for Aldo and Perry....)
    if (exists("dropbox")){
      file.copy(paste0(dropbox,"/Support/MSigDB/v",msigdb.version, "/MSigDB_Categories.v", msigdb.version,".rds"),".")
    } else {
      cat("Do not have a file with the MSigDB_Categories available\n") 
    }
  }

  # Source our hacked code for the generation of the barcodeplot
  source(paste(dropbox,"/Support/R/barcodeplot.R", sep=""))
  
  # Set the species
  if (species=="human"){
    species="Homo_sapiens"
  } else if (species == "mouse"){
    species="Mus_musculus"
  } else if (species == "c_elegans"){
    species="Caenorhabditis_elegans"
  } else {
    warning("\tEither choose \"human\", \"mouse\" or \"c_elegans\", other species are not (yet?) implemented...")
  }
  
  # Set the directory where to store the results
  dir.create(reportsDir)
  
  # CAMERA using inter.gene.correlation set to 0.01 as CAMERA is conservative when using small genesets
  # (https://support.bioconductor.org/p/70005/)
  if ( wantedContrast == "all") { wantedContrast = 1:ncol(contrasts) }
  
  for (i in wantedContrast){ 
    tab.camera = camera(v., sets., design., contrast=contrasts[,i],inter.gene.cor=interGeneCor, use.ranks=TRUE)
    contrast = attr(contrasts,"dimnames")$Contrasts[i]
    write.table(data.frame(geneset=rownames(tab.camera),tab.camera),file=paste0(reportsDir,paste0("/tab.camera.v",msigdb.version,"_"),colnames(contrasts)[i],".txt"),row.names=FALSE,quote=FALSE,sep="\t")
    
    if ((class(v.)[1]=="DGEList" || class(v.)[1]=="DGEGLM") && attr(class(v.),"package")=="edgeR"){
      if (!is.null(fit$df.residual.zeros)) {
        lrt = glmQLFTest(fit.,contrast=my.contrasts[,i])
      } else {
        lrt = glmLRT(fit.,contrast=my.contrasts[,i])
      }
      tt  = topTags(lrt,n=Inf,sort="none")$table
      y = -log10(tt$PValue)
    } 
    if (class(v.)[1] %in% c("EList","ExpressionSet") && attr(class(v.),"package") %in% c("limma","Biobase")){
      fit2 = contrasts.fit(fit.,contrasts[,i])
      fit.eb = eBayes(fit2, trend=trend, robust=robust)
      tt = topTable(fit.eb,n=Inf,sort.by="none")
      y = -log10(tt$P.Value)
    }
    
    x = tt$logFC
    x_min = round(min(x))
    x_max = round(max(x))
    y_max = round(max(y)) + 1
    
    tab.all = read.delim(paste0(reportsDir,paste0("/tab.camera.v",msigdb.version,"_"),contrast,".txt"))
    tab = tab.all[which(tab.all$PValue < p.value & !(is.na(tab.all$PValue))),]  
    ## AJ - added clause to remove "NA"s as these give problems later on!
    
    cat("\n",contrast,"\n----------------------------------------------------------------------\n")
    
    if (nrow(tab) > 0){
      addEGIDLink <- function(object, ...){        
        object$entrezgene <- ifelse(!is.na(as.character(object$entrezgene)),
                                    hwrite(as.character(object$entrezgene),link = paste0("http://www.ncbi.nlm.nih.gov/gene/",as.character(object$entrezgene)), table = FALSE),
                                    hwrite("NA",table=FALSE)
        )                           
        object$description <- ifelse(!is.na(as.character(object$description)),
                                     hwrite(as.character(object$description), table = FALSE),
                                     hwrite("NA",table=FALSE)
        )                           
        object$ensembl_gene_id <- ifelse(!is.na(as.character(object$ensembl_gene_id)),
                                         hwrite(as.character(object$ensembl_gene_id), link = paste0("http://www.ensembl.org/",species,"/Gene/Summary?db=core;g=",as.character(object$ensembl_gene_id)), table = FALSE),
                                         hwrite("NA",table=FALSE)
        )                           
        return(object)
      }
      addNGenesLink = function(object, ...){
        filename <- c()
        for (i in 1:nrow(object)){
          filename[i] = paste0(as.character(object$geneset[i]),".html")
          cat(i,filename[i],"\n")
          indx = sets[[as.character(object$geneset[i])]]  
          
          basefilename = as.character(object$geneset[i])
          
          #     Some pathway names are too long to be acceptable by R:
          if (nchar(basefilename) > 80){
            basefilename = gsub("\\.", "", basefilename)
            basefilename = strsplit(basefilename, '')[[1]]
            indxSplit = c(1,which(basefilename == "_") + 1)
            basefilename = paste(basefilename[indxSplit[1:length(indxSplit)-1]],collapse = "_")
            basefilename = gsub(" ", "", basefilename)
            filename[i] = paste0(basefilename,".html")
            cat("Pathway abbreviated to: ", basefilename," (due to length of name!!)\n")
          }
          plotfilename=paste0(basefilename,".png")
          png(filename=paste0(reportsDir,"/",contrast,"/",plotfilename))
          plot(x,y,xlab="log2 fold change",ylab="-log10 p-value", pch = 16, cex = 0.3,xlim=c(x_min-1,x_max+1),ylim=c(0,y_max+1),cex.lab=1.2)
          points(x[indx], y[indx], pch = 16, cex = 1.25, col = "red")
          abline(h=-log10(0.05),lty=2)
          text(-2.3,1.2, labels = "p = 0.05", cex = 0.8, col="black")
          abline(v=0,lty=2)
          legend(x = x_min-1, y = y_max+1, legend = as.character(object$geneset[i]), pch=16, col = c("red"), cex = 0.9)
          dev.off()
          
          plotfilename2 = paste0(basefilename,"_barcode.png")
          png(filename=paste0(reportsDir,"/",contrast,"/",plotfilename2))
          barcodeplot(-tt$logFC,indx,labels=c("DOWN","UP"))
          legend("top",legend=object$geneset[i],lty=1)
          dev.off()
          
          textfilename = paste0(basefilename,".txt")
          write.table(tt[indx,],file=paste0(reportsDir,"/",contrast,"/",textfilename),quote=FALSE,sep="\t",row.names=FALSE)
          
          if (nrow(tt[indx,]) > 0){
            genesetReport = HTMLReport(shortName=basefilename,title=as.character(object$geneset[i]),reportDirectory=paste0("./",reportsDir,"/",contrast))
            himg  = hwriteImage(plotfilename,link=plotfilename)
            himg2 = hwriteImage(plotfilename2,link=plotfilename2)
            publish(hwrite("Volcano plot"), genesetReport)
            publish(hwrite(himg, br=TRUE),genesetReport,name="Volcano plot") 
            publish(hwrite("Barcode plot"), genesetReport)
            publish(hwrite(himg2, br=TRUE),genesetReport,name="Barcode plot") 
            publish(hwrite("Top table",link=textfilename),genesetReport)
            publish(tt[indx,],genesetReport,.modifyDF=list(addEGIDLink),name="Table")
            finish(genesetReport) 
          }
        }  
        object$NGenes = hwrite(as.character(object$NGenes),link=filename,table = FALSE)
        return(object)
      }
      addMSigDBLink = function(object, ...){           
        object$geneset = 
          hwrite(as.character(object$geneset),link=paste0("http://www.broadinstitute.org/gsea/msigdb/geneset_page.jsp?geneSetName=",as.character(gsub("^C[1-8]_|^H_","",as.character(object$geneset)))) , table = 
                   FALSE)
        return(object)
      }
      cameraReport = HTMLReport(shortName="CAMERA",title=paste0("CAMERA (contrast: ",contrast,", geneset cut-off: p  &lt;",p.value,")"),reportDirectory=paste0("./",reportsDir,"/",contrast))
      if (dim(tab)[1] > max.pathways){
        publish(tab[1:max.pathways,],cameraReport,.modifyDF=list(addNGenesLink,addMSigDBLink))
      } else {
        publish(tab,cameraReport,.modifyDF=list(addNGenesLink,addMSigDBLink))
      }
      finish(cameraReport)
    }
  }
  
  rm(cameraReport)
  
}