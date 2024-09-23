voomReports <- function(fit, contrasts=my.contrasts, my.dir="./voom/", p.value=0.01, fdr=0.005, eBayesTrend=FALSE, eBayesRobust=FALSE, ...){

  # AJ - 13092016
  
  # Load necessarry packages
  require(limma)
  require(hwriter)
  require(ReportingTools)
  
  # You will need an object with 'entrezgene', 'ensembl_gene_id' and 'description'..
  # Maybe it can be made more generic??
  
  fit2 = contrasts.fit(fit, contrasts)
  fit.eb = eBayes(fit2, trend=eBayesTrend, robust=eBayesRobust)
  
  for (i in 1:ncol(contrasts)){   
    contrast = attr(contrasts,"dimnames")$Contrasts[i]
    tab.all = topTable(fit.eb,coef=contrast,n=Inf,sort.by="none")
    tab = tab.all[tab.all$P.Value < p.value,]
    
    voomReport = HTMLReport(shortName="voom",title=paste0("voom (contrast: ",contrast,")"),reportDirectory=paste0(my.dir,contrast))
    
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
                                       hwrite(as.character(object$ensembl_gene_id), table = FALSE),
                                       hwrite("NA",table=FALSE)
      )                           
      return(object)
    }
    
    #detags = rownames(D)[which(p.adjust(tab.all$P.Value,"BH")<0.05)]
    detags = ifelse(p.adjust(tab.all$P.Value,"BH")< fdr,"de","non")
    voomDir = paste0(my.dir,contrast)
    plotfilename = paste0(contrast,".png")
    png(paste0(voomDir,"/",plotfilename))
    #limma::plotMA(fit.eb)
    limma::plotMA(fit.eb, coef=i, values=c("non","de"), status=detags, col=c("black","red"), cex=c(0.1,1), main=paste0("MA plot - ", contrast))
    mtext(paste0("p.adj (BH) = ",fdr),side=3, line=0.2, cex=0.8, col="black")
    dev.off()
    rm(detags)
    
    himg  = hwriteImage(plotfilename,link=plotfilename)
    #  publish(hwrite("Smear plot", style='text-align: center; font-weight: bold; font-style: italic;font-size:150%; color: blue'), voomReport)
    publish(hwrite(himg, br=TRUE),voomReport,name="MA plot")
    if (nrow(tab) > 0){
      publish(tab,voomReport,.modifyDF=list(addEGIDLink))
    }
    finish(voomReport)
  }
  
  rm(voomReport)

}