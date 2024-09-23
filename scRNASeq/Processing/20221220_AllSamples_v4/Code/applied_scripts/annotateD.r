annotateD <- function(D, geneSet, species="human", geneFilter="hgnc_symbol", host="www.ensembl.org", 
                      annotDir=NULL, ensemblVersion=NULL, extraFields=NULL, ...){

  # AJ - 12092016
  # AJ - 16062022 - Added 'ExtraFields' option to be able to fx. add the 'gene_type' information form Ensembl
  
  # Load necessarry packages
  require(biomaRt)
  require(plyr)

  #### Annotation ####
  # This is the point to add the gene info (Entrez ID, Ensembl ID, description)
  #  geneSet <- rownames(D$counts)
  # or:
  #  geneSet <- sapply(strsplit(rownames(D$counts),split=":"),'[',1)
  #  transcripts <- sapply(strsplit(rownames(D$counts),split=":"),'[',3)
  #  rownames(D$counts) <- paste(sapply(rownames(total),function(x) strsplit(x,":")[[1]][1]),sapply(rownames(total),function(x) strsplit(x,":")[[1]][2]),sep=":")
  #  geneSet <- sapply(strsplit(geneSet,split="_dup"),'[',1)
  # The latter has been used in combination with Lifescope (i.e. SOLiD) pipeline. The 'transcript' column is dropped here...
  
  # https://support.bioconductor.org/p/74322/
  # Biomart doesn't provide the service any longer...
  # "www.ensembl.org" reflects the latest version
  # NOTE that this means the latest GRCh, GRCm, hg or mm as well and that this might be of importance when
  # using coordinates etc. Here we only use identifiers, so we should be save..
  #
  if (species == "human"){
    if (!is.null(ensemblVersion)){
      ensembl = useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host=host, version=ensemblVersion)
    } else {
      ensembl = useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host=host)
    }
    # Write out the versions used...
    write.table(listMarts(host=host),paste0("biomaRt.",species,".versions"), sep="\t", quote=FALSE)
    symbol = "hgnc_symbol"
  } else if (species == "mouse"){
    if (!is.null(ensemblVersion)){
      ensembl = useMart("ENSEMBL_MART_ENSEMBL",dataset="mmusculus_gene_ensembl", host=host, version=ensemblVersion)
    } else {
       ensembl = useMart("ENSEMBL_MART_ENSEMBL",dataset="mmusculus_gene_ensembl", host=host)
    }
    # Write out the versions used...
    write.table(listMarts(host=host),paste0("biomaRt.",species,".versions"), sep="\t", quote=FALSE)
    symbol = "mgi_symbol"
  } else {
    warning("\tEither choose \"human\" or \"mouse\", other species are not (yet?) implemented...")
  }
  
  # 05072019 - Ensembl has renamed their "entrezgene" identifier to "entrezgene_id"...
  entrezGene = "entrezgene"
  if ( "entrezgene_id" %in% listAttributes(ensembl)$name ){
    entrezGene = "entrezgene_id"
  }
  
  if ( is.null(annotDir) ){
    # Get annotation
    annot.biomart = getBM(attributes = c(symbol, entrezGene, "ensembl_gene_id","description", extraFields), filters = geneFilter, values = geneSet, mart = ensembl, uniqueRows=FALSE)
    
    # getBM returns "" in case of missing symbol (not NA)
    annot.biomart[annot.biomart[,symbol]=="",symbol] <- NA
    
    saveRDS(annot.biomart,file="annot.biomart.rds")
    rm(ensembl)
    annot.biomart = readRDS("annot.biomart.rds")
  } else {
    annot.biomart = readRDS(paste0(annotDir,"/annot.biomart.rds"))
  }
  
  ## remove duplicated combinations of GeneSymbol, Entrez Gene ID and Ensembl Gene ID
  annot.biomart = subset(annot.biomart,!duplicated(annot.biomart[,c(symbol, entrezGene, "ensembl_gene_id")]))
  
  ## concatenate probes mapping to multiple Entrez Gene IDs
#  df = aggregate(annot.biomart,list(annot.biomart[,geneFilter]),function(x) paste(unlist(x),collapse="//"))
  df = aggregate(annot.biomart,list(annot.biomart[,geneFilter]),function(x) paste(unique(x),collapse="//"))
  
  # If entrezgene etc. is a concatenation of two (or more) the same ids, uniguefy and save, otherwise just leave it ..
#  df$entrezgene = unlist(lapply(df$entrezgene,function(x){a = unique(unlist(strsplit(x,"//")));ifelse(length(a) == 1,a,x)}))
#  df$ensembl_gene_id = unlist(lapply(df$ensembl_gene_id,function(x){a = unique(unlist(strsplit(x,"//")));ifelse(length(a) == 1,a,x)}))
#  df$description = unlist(lapply(df$description,function(x){a = unique(unlist(strsplit(x,"//")));ifelse(length(a) == 1,a,x)}))
#  df[,symbol] = unlist(lapply(df[,symbol],function(x){a = unique(unlist(strsplit(x,"//")));ifelse(length(a) == 1,a,paste(a,collapse="//"))}))
  
  # Rename "entrezgene_id" to "entrezgene"
  names(df)[names(df) == entrezGene] <- 'entrezgene'
  
  # Replace by NA
  df[df$entrezgene=="NA","entrezgene"] <- NA
  
  # Remove "entrezgene_id"
  df$entrezgene_id <- NULL
  
  idx <- which(colnames(df)==geneFilter)
  df <- df[,-c(idx)]
  colnames(df)[1] <- geneFilter
  
  if ( "counts" %in% colnames(D)){
    genes <- as.data.frame(rownames(D$counts))
  } else {
    genes <- as.data.frame(rownames(D))
  }
  genes <- cbind(genes,geneSet)
  colnames(genes) <- c("original",geneFilter)
  
  # Make sure to use 'join' from plyr package
  genesListNew <- plyr::join(genes,df,by=geneFilter,type="left")
  
  # Replace the old genenames with the new annotation information
  D$genes <- genesListNew
  
  # Replace <NA> by "NA" 
  D$genes[is.na(D$genes)] <- "NA"

  # Give some stats back
  nrSymbol <- sum(D$genes[,symbol] != "NA")
  nrEnsembl <- sum(D$genes$ensembl_gene_id != "NA")
  nrEntrez <- sum(D$genes$entrezgene != "NA")
  cat("\nAnnotation:\n+++++++++++\n")
  cat("  #genes:\t",dim(D$genes)[1],"\n\nRetrieved:\n")
  cat("  #symbols    :\t",nrSymbol,"(",(nrSymbol/dim(D$genes)[1])*100,"%)\n")
  cat("  #ensembl ids:\t",nrEnsembl,"(",(nrEnsembl/dim(D$genes)[1])*100,"%)\n")
  cat("  #entrez  ids:\t",nrEntrez,"(",(nrEntrez/dim(D$genes)[1])*100,"%)\n\n")
  
  # Clean
  rm(annot.biomart,genes,geneSet,genesListNew,df, nrSymbol, nrEntrez, nrEnsembl)
  
  # Return object
  return(D)
}
