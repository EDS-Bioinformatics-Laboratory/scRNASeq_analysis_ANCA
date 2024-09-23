# Cell annotation with SingleR and CellTypist

#### Set default for a good start ####
source("initAnalysis.r")
init()

# Create a directory to hold results and plots
outDir <- "4_cell_annotation_before_integration"
dir.create(outDir)

#### Library ####
library(Seurat)
library(ggplot2)
library(cowplot)

dataDir <- "3_doublet_detection"

#### SingleR - individual databases ####
library(SingleR)
library(SingleCellExperiment)
library(scater)
library(pheatmap)
library(ggplot2)
library(IntEREst)
library(grid)
library(corrplot)

# Load the Immunecell, MonacoImmune and NovershternHematopoietic databases
immCellExpr <- DatabaseImmuneCellExpressionData()
monacoImmData <- MonacoImmuneData()
hemapoieticData <- NovershternHematopoieticData()
hpca.se <- HumanPrimaryCellAtlasData()
hpca.se

# Get the specific annotation for the mature kidney from https://www.kidneycellatlas.org/
# Stewart BJ. et al. Science 2019: https://science.sciencemag.org/content/365/6460/1461
# see '5a_cell_annotation_emptyDroples_retain.r' in '20210630_Analysis_MitoContent/Code'
# I copied the object to this directory 4_cell_annotation_before_integration

# NOTE:
# - Also do this for the mature_immune database?

matureKidney.sce <- readRDS(paste0(outDir,"/matureKidney.sce.rds"))

# List of databases
dbList <- list(immCellExpr=immCellExpr, monaco=monacoImmData, hemapoietic=hemapoieticData, hpca = hpca.se, matureKidney = matureKidney.sce)

# List of samples
file_names <- dir("./3_doublet_detection/", pattern = "*doubletsRemoved.rds")
file_names <- gsub(".doubletsRemoved.rds","",file_names)

for ( i in 1:length(file_names) ){
  print(file_names[i])
  annotDir <- paste0(outDir,"/", paste0(gsub("\\.","_",file_names[i]), "_SingleR_annotation"))
  dir.create(annotDir)
  
    # Read data
    if ( file.exists(paste0(outDir,"/",file_names[i],".doubletsRemoved.annot.rds"))) {
      seurat <- readRDS(paste0(outDir,"/",file_names[i],".doubletsRemoved.annot.rds"))
    } else {
      seurat <- readRDS(paste0(dataDir,"/",file_names[i],".doubletsRemoved.rds"))
    }
    
    # Convert the Seurat object to a 'SingleCellExperiment' format
    seurat.sce <- as.SingleCellExperiment(seurat, assay = "RNA")
    
    for (grain in c("main","fine","ont")){
      print(paste0("... ",grain))
      for (j in 1:length(dbList)){
        db <- dbList[j]
        # And extract the features we are going to use for cell identification
        common <- intersect(rownames(seurat.sce), rownames(db[[1]]))
        seurat.db <- seurat.sce[common,]
        
        # Test and reference sets should always be log normalized. The included reference sets are already normalized. 
        seurat.db <- scater::logNormCounts(seurat.db)
        
        # use the more specific fine-grained cell type labels
        if (grain == "main"){
          pred.db <- SingleR(test = seurat.db, ref = db[[1]][common,], 
                             labels = db[[1]][common,]$label.main, assay.type.ref = "logcounts")
        } else if (grain == "fine"){
          pred.db <- SingleR(test = seurat.db, ref = db[[1]][common,], 
                             labels = db[[1]][common,]$label.fine, assay.type.ref = "logcounts")
        } else {
          if (names(dbList[j]) == "matureKidney"){
            pred.db <- SingleR(test = seurat.db, ref = db[[1]][common,], 
                               labels = db[[1]][common,]$label.compartment, assay.type.ref = "logcounts")
          } else {
            pred.db <- SingleR(test = seurat.db, ref = db[[1]][common,], 
                               labels = db[[1]][common,]$label.ont, assay.type.ref = "logcounts")
          }
        }
        df <- as.data.frame(table(pred.db$labels))
        colnames(df) <- c("CellType", "Freq")
        write.table(df, file=paste(annotDir,paste("cellTypesPredicted", names(db), grain, "txt", sep="."), sep="/"), sep="\t", row.names = FALSE, quote = FALSE)
        
        # Annotate the original seurat object with the pruned labels
        if (names(dbList[j]) == "matureKidney" & grain == "ont"){ 
          seurat@meta.data[,paste(names(db), "compartment","before", sep=".")] <- pred.db$pruned.labels
        } else {
          seurat@meta.data[,paste(names(db), grain, "before",sep=".")] <- pred.db$pruned.labels
        }
      }
      
      # And now we have to put some classes/types together for clarity...
      pdf(paste(annotDir,paste("cellTyping", grain,file_names[i],"beforeIntegration.pdf",sep="."), sep="/"))
      for (db in c("immCellExpr","monaco","hemapoietic","hpca", "matureKidney")){
        # Need this to get the legend
        if (grain %in% c("fine","ont")){ my_ncol = 3 } else { my_ncol = 2 }
        if (db == "matureKidney" & grain == "ont"){grain <- "compartment"}
        
        p <- DimPlot(seurat, reduction = "tsne", group.by = paste(db, grain, "before", sep="."), label = TRUE, repel = TRUE) +
          guides(color = guide_legend(override.aes = list(size = 2), ncol = my_ncol)) + 
          ggtitle(paste0(paste(db, grain, sep="."), " before integration" )) +
          theme(legend.text = element_text(size = 9), legend.spacing.y = unit(0, 'cm'), legend.spacing.x = unit(0, 'cm'),
                legend.margin=margin(0,0,0,0)) 
        
        # Using the cowplot package
        legend <- cowplot::get_legend(p)
        p <- DimPlot(seurat, reduction = "tsne", group.by = paste(db, grain, "before", sep="."), label = TRUE, repel = TRUE) +
          NoLegend() + ggtitle(paste0(paste(db, grain, sep="."), " before integration" ))
        
        print(p)
        grid.newpage()
        grid.draw(legend)
        #
        p <- DimPlot(seurat, reduction = "umap", group.by = paste(db, grain, "before", sep="."), label = TRUE, repel = TRUE) +
          NoLegend() + ggtitle(paste0(paste(db, grain, sep="."), " before integration" ))
        
        print(p)
        grid.newpage()
        grid.draw(legend)
        
      }
      dev.off()
      
    }
    
    # Save the object with the new cell type annotation
    saveRDS(seurat, paste0(outDir,"/",file_names[i],".doubletsRemoved.annot.rds"))
}

# And clean
rm(p, df, legend, pred.db, dbList, seurat, seurat.sce, seurat.db)
rm(immCellExpr, monacoImmData, hemapoieticData, hpca.se, matureKidney.sce)

#### CellTypist - before integration ####
dataDir <- "3_doublet_detection"

## Read the data ##
file_names <- dir(dataDir, pattern = "*doubletsRemoved.rds")
file_names <- gsub(".doubletsRemoved.rds","",file_names)

for ( i in 1:length(file_names) ){
  print(file_names[i])
  
  # Read data
  if ( file.exists(paste0(outDir,"/",file_names[i],".doubletsRemoved.annot.rds"))) {
    seurat <- readRDS(paste0(outDir,"/",file_names[i],".doubletsRemoved.annot.rds"))
  } else {
    seurat <- readRDS(paste0(dataDir,"/",file_names[i],".doubletsRemoved.rds"))
  }
  
  #### Get the raw counts ####
  write.table(as.data.frame(seurat@assays$RNA@counts), paste0(outDir,"/",gsub("\\.","_",file_names[i]), "_countsBeforeIntegration.csv"),
              sep = ',', row.names = T, col.names = T, quote = F)
  
}

#### CellTypist - Run ####
# see ../NoteBooks/runCellTypist_beforeIntegration.ipynb
#
# I tried different scenario's
# 1. Just starting from the raw counts and without majority voting, but with 'best match'
#    - This will be the same before and after integration
# 2. Apply majority voting
#    - This has a dependency on neighbouring cells, so will be different before and after integration

#### CellTypist - Read predicted labels ####
# see ../NoteBooks/runCellTypist_beforeIntegration.ipynb
# 
# see https://www.celltypist.org/tutorials/onlineguide
#
# I ran both scenario's on every subject
# - saved in their own subdirectory with the sample directory
#   i.e. <SAMPLE>/majorityvoting
#        <SAMPLE>/bestmatch

file_names <- dir(outDir, pattern = "*doubletsRemoved.annot.rds")
file_names <- gsub(".doubletsRemoved.annot.rds", "",file_names)

for ( i in 1:length(file_names) ){
  print(file_names[i])
  
  sampleDir <- paste(outDir, paste0(gsub("\\.", "_",file_names[i]),"_CellTypist_annotation"), sep="/" )
  
  # Read data
  if ( file.exists(paste0(outDir,"/",file_names[i],".doubletsRemoved.annot.rds"))) {
    seurat <- readRDS(paste0(outDir,"/",file_names[i],".doubletsRemoved.annot.rds"))
  } else {
    seurat <- readRDS(paste0(dataDir,"/",file_names[i],".doubletsRemoved.rds"))
  }
  
  if (length(grep("celltypist", colnames(seurat@meta.data))) > 0){
    # Make a fresh start, remove all columns containing 'celltypist'in their name ...
    print("Removing all columns holding celltypist information ...")
    seurat@meta.data <- seurat@meta.data[,-c(grep("celltypist", colnames(seurat@meta.data)))]
  }
  # The 'probability_matrix.csv' contains the probabilities per cell of to which celltype it would fit.
  # NOTE: the 'probability_matrix' is of course the same for both the 'best match' and 'majority voting' scenario
  annot.bestMatch.df <- read.csv(paste0(sampleDir,"/bestMatch_noMajorityVoting/predicted_labels.csv"))
  annot.majorityVoting.df <- read.csv(paste0(sampleDir,"/probMatch_majorityVoting/predicted_labels.csv"))
  probMatrix.df <- read.csv(paste0(sampleDir,"/bestMatch_noMajorityVoting/probability_matrix.csv"), check.names = FALSE)
  
  # check whether cells match
  all(colnames(seurat) == annot.bestMatch.df$X)
  all(colnames(seurat) == annot.majorityVoting.df$X)
  
  # # Make sure to keep the rownames !!
  seurat@meta.data$barcodes.orig <- rownames(seurat@meta.data)
  
  # # now, merge by rowname and celltypistName
  colnames(annot.bestMatch.df) <- c("cellName","bestMatchNoMajorityVoting.labels.celltypist.before")
  if (ncol(annot.majorityVoting.df) > 2){
    colnames(annot.majorityVoting.df) <- c("cellName","predicted.labels.celltypist", "over_clustering.celltypist", "probMatchMajorityVoting.labels.celltypist.before")
  } else {
    colnames(annot.majorityVoting.df) <- c("cellName","probMatchMajorityVoting.labels.celltypist.before")
  } 
  seurat@meta.data <- merge(seurat@meta.data, annot.bestMatch.df, by.x="row.names",by.y="cellName", all.x = TRUE, sort = FALSE )
  seurat@meta.data <- merge(seurat@meta.data, annot.majorityVoting.df, by.x="Row.names",by.y="cellName", all.x = TRUE, sort = FALSE )
  # 
  
  # And add the probability scores ... becoming a very large object :-)
  colnames(probMatrix.df) <- paste0(colnames(probMatrix.df),".celltypist")
  celltypist.max_probability <- 0
  celltypist.max_probability <- apply(probMatrix.df[,-1], 1, function(x) max(x))
  
  seurat@meta.data <- merge(seurat@meta.data, probMatrix.df, by.x="Row.names",by.y=".celltypist", all.x = TRUE, sort = FALSE )
  
  # Put the rownames back
  rownames(seurat@meta.data) <- seurat@meta.data$barcodes.orig
  seurat@meta.data$Row.names <- NULL
  seurat@meta.data$barcodes.orig <- NULL
  # And add the column with the maximum probability
  seurat@meta.data$celltypist.max_probability <- celltypist.max_probability
  
  # And save
  saveRDS(seurat, paste0(outDir,"/",file_names[i],".doubletsRemoved.annot.rds"))
} 

#### CellTypist - Plots and tables ####
file_names <- dir(outDir, pattern = "*doubletsRemoved.annot.rds")
file_names <- gsub(".doubletsRemoved.annot.rds", "",file_names)

for ( i in 1:length(file_names) ){
  print(file_names[i])
  
  sampleDir <- paste(outDir, paste0(gsub("\\.", "_",file_names[i]),"_CellTypist_annotation"), sep="/" )
  
  # Read data
  if ( file.exists(paste0(outDir,"/",file_names[i],".doubletsRemoved.annot.rds"))) {
    seurat <- readRDS(paste0(outDir,"/",file_names[i],".doubletsRemoved.annot.rds"))
  } else {
    seurat <- readRDS(paste0(dataDir,"/",file_names[i],".doubletsRemoved.rds"))
  }
  
  # Make plots for bestMatchNoMajorityVoting and probMatchMajorityVoting
  for (scenario in c("bestMatchNoMajorityVoting", "probMatchMajorityVoting")){
    pdf(paste(sampleDir, paste0("umap_celltypist_", scenario,"_annot.", file_names[i],".pdf"), sep="/"), width = 14, height = 7)
    p1 <- DimPlot(object = seurat, reduction = "umap", group.by = paste0(scenario,".labels.celltypist.before"), label = TRUE) + NoLegend()
    p2 <- DimPlot(object = seurat, reduction = "umap", group.by = "SCT_snn_res.0.3", label = TRUE) + NoLegend()
    print((plot_grid(p1, p2)))
    # Get the legends on a separate page using the cowplot package
    p1 <- DimPlot(object = seurat, reduction = "umap", group.by = paste0(scenario,".labels.celltypist.before"), label = TRUE) +
      guides(color = guide_legend(override.aes = list(size=3), nrow = 30) ) +
      theme(legend.text = element_text(size=9))
    legend_p1 <- cowplot::get_legend(p1)
    p2 <- DimPlot(object = seurat, reduction = "umap", group.by = "SCT_snn_res.0.3", label = TRUE)
    legend_p2 <- cowplot::get_legend(p2)
    print(plot_grid(NULL, legend_p1, NULL, NULL, legend_p2,ncol = 5))
    dev.off()
    
    pdf(paste(sampleDir, paste0("umap_celltypist_maxProbability.", file_names[i],".pdf"), sep="/"), width = 14, height = 7)
    p1 <- FeaturePlot(object = seurat, features = "celltypist.max_probability", reduction = "umap")
    print(plot_grid(p1))
    dev.off()
    
    pdf(paste(sampleDir, paste0("violin_celltypist_probabilities.", file_names[i],".pdf"), sep="/"), width = 14, height = 14)
    p1 <- VlnPlot(object = seurat, features = "celltypist.max_probability", group.by = paste0(scenario,".labels.celltypist.before")) + NoLegend()
    print(plot_grid(p1))
    dev.off()
    
    pdf(paste(sampleDir, paste0(scenario,"_celltypist_annotation.", file_names[i],".pdf"), sep="/"), width = 14, height =14)
    tab <- table(seurat$SCT_snn_res.0.3, seurat$bestMatchNoMajorityVoting.labels.celltypist.before)
    corrplot(tab/rowSums(tab), is.corr = FALSE, col = 
               colorRampPalette(c("#FFFFFF", "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC", "#053061"))(200))
    dev.off()
    
    # And get some tables 
    # - of the distribution of the different celltypes over the different clusters
    
    tab <- table(seurat$SCT_snn_res.0.3, seurat@meta.data[,paste0(scenario,".labels.celltypist.before")])
    write.table(tab, paste(sampleDir,paste0(scenario,"_celltypes_per_cluster_res03.", file_names[i],".txt"), sep="/"), quote=FALSE, sep="\t", row.names = TRUE, col.names = TRUE)
    
    df <- as.data.frame(tab)
    colnames(df) <- c("cluster", "celltype", "frequency")
    
    df$celltype <- factor(df$celltype)
    nrColors <- length(levels(factor(df$celltype)))
    rainbowColors = rainbow(nrColors)
    df$colors <- rainbowColors[factor(df$celltype)]
    
    pdf(paste(sampleDir, paste0("distribution_per_cluster_stack.", scenario,".", file_names[i],".pdf"), sep="/"), width = 14, height =14)
    p <- ggplot(df, aes(x=cluster, y=frequency, fill=celltype)) +
      geom_bar(stat='identity', position="stack")+
      scale_fill_manual(values = rainbowColors) +
      guides(fill=guide_legend(ncol =1))
    print(p)
    legend_p <- cowplot::get_legend(p)
    p <- ggplot(df, aes(x=cluster, y=frequency, fill=celltype)) +
      geom_bar(stat='identity', position="stack")+
      scale_fill_manual(values = rainbowColors) + 
      NoLegend()
    print(p)
    print(plot_grid(legend_p))
    dev.off()
  }
}

#### Annotation TB cells? ####
file_names <- dir(outDir, pattern = "*doubletsRemoved.annot.rds")
file_names <- gsub(".doubletsRemoved.annot.rds", "",file_names)
file_names <- file_names[grepl("MOMA52|MOMA57", file_names)]

for ( i in 1:length(file_names) ){
  print(file_names[i])
  
  sampleDir <- paste(outDir, paste0(gsub("\\.", "_",file_names[i]),"_CellTypist_annotation"), sep="/" )
  
  # Read data
  if ( file.exists(paste0(outDir,"/",file_names[i],".doubletsRemoved.annot.rds"))) {
    seurat <- readRDS(paste0(outDir,"/",file_names[i],".doubletsRemoved.annot.rds"))
  } else {
    seurat <- readRDS(paste0(dataDir,"/",file_names[i],".doubletsRemoved.rds"))
  }
  
  # Identify the TB cells
  idx <- which(!is.na(seurat@meta.data$umis_vdj_t))
  selCells_T <- rownames(seurat@meta.data)[idx]
  idx <- which(!is.na(seurat@meta.data$umis_vdj_b))
  selCells_B <- rownames(seurat@meta.data)[idx]
  T_and_B <- intersect(selCells_B, selCells_T)
  # Identify the cells present in both T and B cell sets from the original ones
  selCells_B <- setdiff(selCells_B, selCells_T)
  selCells_B <- selCells_B[!(selCells_B %in% T_and_B)]
  selCells_T <- setdiff(selCells_T, selCells_B)
  selCells_T <- selCells_T[!(selCells_T %in% T_and_B)]
  
  seurat@meta.data$TB <- ifelse(rownames(seurat@meta.data) %in% selCells_B, "B", 
                                ifelse(rownames(seurat@meta.data) %in% selCells_T, "T",
                                       ifelse(rownames(seurat@meta.data) %in% T_and_B, "T_and_B", "unknown")))
  
  #seurat.tb <- subset(seurat, subset = TB == "T_and_B" )
  
  # Make tables for bestMatchNoMajorityVoting and probMatchMajorityVoting
  for (scenario in c("bestMatchNoMajorityVoting", "probMatchMajorityVoting")){
    print(scenario)
    annot <- paste0(scenario,".labels.celltypist.before")
    tab <- table(seurat@meta.data$TB, seurat@meta.data[,annot])
    print(tab)
    write.table(tab, paste(sampleDir,paste0(scenario,"_celltypes_VDJ.", file_names[i],".txt"), sep="/"), quote=FALSE, sep="\t", row.names = TRUE, col.names = TRUE)
    
    cat("\n")
  }
}

#### Concatenate tables with cellType frquencies ####
file_names <- dir(paste0(outDir,"/", dir(outDir, pattern = "*Mito10_SingleR_annotation")), pattern = "*matureKidney.main.txt", full.names = TRUE)

for ( i in 1:length(file_names) ){
  print(file_names[i])
  s <- gsub("4_cell_annotation_before_integration/|_Mito10_SingleR_annotation/cellTypesPredicted.matureKidney.main.txt","", file_names[i])
  
  f <- read.table(file_names[i], header=TRUE, sep="\t")
  names(f)[2] <- s
    
  if (i==1){
    df <- f
  } else {
    df <- merge(df,f, by="CellType", all=TRUE)
  }
}

rownames(df) <- df$CellType
df$CellType <- NULL

write.table(df, file=paste(outDir,"predictedCellTypes_matureKidney_main.tsv", sep="/"), sep="\t", 
            quote=FALSE, row.names = TRUE, col.names = NA )
  