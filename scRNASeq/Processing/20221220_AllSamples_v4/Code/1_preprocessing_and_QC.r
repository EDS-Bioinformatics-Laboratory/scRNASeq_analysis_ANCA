# Preprocessing and QC 
#
# Previously decided for MOMA17:
#   20210716 - Decided: Cut off for mitochondrial genes 20, 50, 80% on filtered CellRanger data

#### Set default for a good start ####
source("initAnalysis.r")
init()

# Create a directory to hold results and plots
outDir <- "1_preprocessing_and_QC"
dir.create(outDir)

#### Library ####
library(viridis)
library(ggplot2)
library(patchwork)
library(data.table)

#### Read data ####
# I have data sets with different modalities (Gene expression, B and T cell hashtags, ...)
# First concentrate on the gene expression

# Previously for 20210628_Patient1
# # For output from CellRanger >= 3.0 with multiple data types, see end of vignette!!
# data_dir <- '../../../Data/Dataset_1/Raw/AMC-SC-S79/YV-20210401-1/outs/filtered_feature_bc_matrix/'
# list.files(data_dir) # Should show barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz - and it does!!
# data <- Read10X(data.dir = data_dir)
# 
# And turn it into a Seurat object
sampleList <- c("MOMA17", "MOMA52", "MOMA57", "MOMA67", "MOMA68", "MOMA72", "MOMA302")

# Function to concatenate entries (see https://www.biostars.org/p/384640/)
data_concater <- function(x){
  x <- levels(factor(x))
  paste(x, collapse = "__")
}

for ( sc in sampleList){
  print(paste0("Reading data for ", sc))
  if ( sc %in% c("MOMA52","MOMA57")){
    data_dir <- paste0("../../../Data/Dataset_1/Processed/", sc, "/sample_feature_bc_matrix/")
    data <- Read10X(data.dir = data_dir)
    #
    print(paste0("... contains following data: ", paste(names(data), collapse=" and "), " - taking both!"))
    
    seurat = CreateSeuratObject(counts = data$`Gene Expression`, min.cells = 0, min.features = 0, project = sc)
    seurat[['ADT']] <- CreateAssayObject(counts = data$`Antibody Capture`, min.cells = 0, min.features = 0)
    #
    print(paste0("... ", ncol(seurat), " cells were found"))
    
    # Add VDJ (only for MOMA52 and MOMA 57)
    vdj_b <- read.csv(paste0("../../../Data/Dataset_1/Processed/",sc,"/vdj_b/filtered_contig_annotations.csv"))
    vdj_t <- read.csv(paste0("../../../Data/Dataset_1/Processed/",sc,"/vdj_t/filtered_contig_annotations.csv"))
    
    # Change some column names
    colnames(vdj_b)[2:ncol(vdj_b)] <- paste0(colnames(vdj_b)[2:ncol(vdj_b)],"_vdj_b")  
    colnames(vdj_t)[2:ncol(vdj_t)] <- paste0(colnames(vdj_t)[2:ncol(vdj_t)],"_vdj_t")
    
    # 'is_cell' and 'high_confidence' could in principle be removed (all 'true')
    # Concatenate
    vdj_b.collapsed <- as.data.table(vdj_b)[, lapply(.SD, data_concater) , by=barcode]
    vdj_t.collapsed <- as.data.table(vdj_t)[, lapply(.SD, data_concater) , by=barcode]
    rownames(vdj_b.collapsed) <- vdj_b.collapsed$barcode
    rownames(vdj_t.collapsed) <- vdj_t.collapsed$barcode
    vdj_b.collapsed$barcode <- NULL
    vdj_t.collapsed$barcode <- NULL
    # And add to the meta data
    seurat <- AddMetaData(seurat, metadata = vdj_b.collapsed)
    seurat <- AddMetaData(seurat, metadata = vdj_t.collapsed)
    
    #
    print(paste0("... for ", nrow(vdj_b.collapsed), " cells B cell receptor information was added"))
    print(paste0("... for ", nrow(vdj_t.collapsed), " cells T cell receptor information was added"))
    
  } else {
    data_dir <- paste0("../../../Data/Dataset_1/Processed/", sc, "/filtered_feature_bc_matrix/")
    data <- Read10X(data.dir = data_dir)
    #
    print(paste0("... contains following data: ", paste(names(data), collapse=" and "), " - taking both!"))
    
    seurat = CreateSeuratObject(counts = data$`Gene Expression`, min.cells = 0, min.features = 0, project = sc)
    seurat[['ADT']] <- CreateAssayObject(counts = data$`Antibody Capture`, min.cells = 0, min.features = 0)
    #
    print(paste0("... ", ncol(seurat)," cells were found"))
    
  } 
  
  # And save
  saveRDS(seurat, paste0("../Data/scRNASeq/Processed/", sc, ".raw.rds"))
}

# [1] "Reading data for MOMA17"
# 10X data contains more than one type and is being returned as a list containing matrices of each type.
# [1] "... contains following data: Gene Expression and Antibody Capture - taking both!"
# [1] "... 7242 cells were found"
# [1] "Reading data for MOMA52"
# 10X data contains more than one type and is being returned as a list containing matrices of each type.
# [1] "... contains following data: Gene Expression and Antibody Capture - taking both!"
# [1] "... 6695 cells were found"
# [1] "... for 494 cells B cell receptor information was added"
# [1] "... for 4069 cells T cell receptor information was added"
# [1] "Reading data for MOMA57"
# 10X data contains more than one type and is being returned as a list containing matrices of each type.
# [1] "... contains following data: Gene Expression and Antibody Capture - taking both!"
# [1] "... 9605 cells were found"
# [1] "... for 486 cells B cell receptor information was added"
# [1] "... for 5140 cells T cell receptor information was added"
# [1] "Reading data for MOMA67"
# 10X data contains more than one type and is being returned as a list containing matrices of each type.
# [1] "... contains following data: Gene Expression and Antibody Capture - taking both!"
# [1] "... 5752 cells were found"
# [1] "Reading data for MOMA68"
# 10X data contains more than one type and is being returned as a list containing matrices of each type.
# [1] "... contains following data: Gene Expression and Antibody Capture - taking both!"
# [1] "... 3128 cells were found"
# [1] "Reading data for MOMA72"
# 10X data contains more than one type and is being returned as a list containing matrices of each type.
# [1] "... contains following data: Gene Expression and Antibody Capture - taking both!"
# [1] "... 1663 cells were found"
# [1] "Reading data for MOMA302"
# 10X data contains more than one type and is being returned as a list containing matrices of each type.
# [1] "... contains following data: Gene Expression and Antibody Capture - taking both!"
# [1] "... 3549 cells were found"

# And clean
rm(data, seurat)

#### Add extra info - mitochondrial fraction, housekeeping genes and cell cycle markers ####

## Housekeeping genes ##
# Load the the list of house keeping genes
file.copy(paste0(dropbox, "/Support/Data/Tirosh_2016_Housekeeping_Genes.txt"),"../Data/scRNASeq/Meta/Tirosh_2016_Housekeeping_Genes.txt")
hkgenes <- read.table("../Data/scRNASeq//Meta/Tirosh_2016_Housekeeping_Genes.txt", skip = 1)
hkgenes <- as.vector(hkgenes$V1)

for ( sc in sampleList){
  seurat <- readRDS(paste0("../Data/scRNASeq/Processed/", sc, ".raw.rds"))
  
  ## Mitochondrial genes ##
  # The number of genes and UMIs (nGene and nUMI) are automatically calculated for every object by Seurat.  For non-UMI
  # data, nUMI represents the sum of the non-normalized values within a cell We calculate the percentage of mitochondrial
  # genes here and store it in percent.mito using AddMetaData.  We use object@raw.data since this represents
  # non-transformed and non-log-normalized counts The % of UMI mapping to MT-genes is a common scRNA-seq QC metric.
  mito.genes <- grep(pattern = "^MT-", x = rownames(x = seurat@assays$RNA), value = TRUE)
  percent.mito <- (Matrix::colSums(seurat@assays$RNA[mito.genes, ])/Matrix::colSums(seurat@assays$RNA)) * 100
  
  # AddMetaData adds columns to object@meta.data, and is a great place to stash QC stats.  This also allows us to plot the
  # metadata values using the Seurat's VlnPlot().
  head(seurat@meta.data)  # Before adding
  
  seurat <- AddMetaData(object = seurat, metadata = percent.mito, col.name = "percent.mito")
  head(seurat@meta.data) # After adding
  #                       orig.ident nCount_RNA nFeature_RNA percent.mito
  # AAACCCAAGCACCAGA-1 SeuratProject       4831         1028     56.17884
  # AAACCCAAGGACAAGA-1 SeuratProject       6641         2226      7.72474
  # AAACCCAAGTCACAGG-1 SeuratProject       4705          249     94.32519
  # AAACCCAAGTGCGACA-1 SeuratProject       7444         2088     45.88931
  # AAACCCACAAGCCATT-1 SeuratProject       8819         2428     26.90781
  # AAACCCACATCGTCCT-1 SeuratProject       7837          729     86.25750
  
  # Save as table
  myTable.df <- as.data.frame(table(round(seurat@meta.data$percent.mito,2)))
  colnames(myTable.df) <- c("Percentage", "Freq")
  write.table(myTable.df, file = paste(outDir,paste0(sc,"_pctMitoGenesExpressed.txt"), sep="/"), row.names = FALSE, col.names = TRUE, quote=FALSE)
  
  # remove hkgenes that were not found
  hkgenes.found <- which(toupper(rownames(seurat@assays$RNA)) %in% hkgenes)
  
  # Add to the metadata and write
  n.expressed.hkgenes <- Matrix::colSums(seurat@assays$RNA[hkgenes.found, ] > 0)
  seurat <- AddMetaData(object = seurat, metadata = n.expressed.hkgenes, col.name = "n.exp.hkgenes")
  myTable.df <- as.data.frame(table(seurat@meta.data$n.exp.hkgenes))
  colnames(myTable.df) <- c("NrGenes", "Freq")
  write.table(myTable.df, file = paste(outDir,paste0(sc,"_numberHouskeepingGenes.txt"),sep="/"), row.names = FALSE, col.names = TRUE, quote=FALSE)
  
  ## Cell cycle markers ##
  # Read in a list of cell cycle markers, from Tirosh et al, 2015.
  # We can segregate this list into markers of G2/M phase and markers of S phase.
  #cc.genes <- readLines(paste0(dropbox,"/Support/Data/regev_lab_cell_cycle_genes.txt"))
  # There already is a list of these genes present!! -> cc.genes.updated.2019
  s.genes <- cc.genes.updated.2019$s.genes
  g2m.genes <- cc.genes.updated.2019$g2m.genes
  
  seurat <- CellCycleScoring(object = seurat, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE, nbin=12)
  # https://github.com/satijalab/seurat/issues/1227
  
  # Genes upregulated during dissociation of tissue into single cells. ( https://broadinstitute.github.io/2019_scWorkshop/data-wrangling-scrnaseq.html)
  genes.dissoc <- c("ATF3", "BTG2", "CEBPB", "CEBPD", "CXCL3", "CXCL2", "CXCL1", "DNAJA1", "DNAJB1", "DUSP1", 
                    "EGR1", "FOS", "FOSB", "HSP90AA1", "HSP90AB1", "HSPA1A", "HSPA1B", "HSPA1A", "HSPA1B", "HSPA8", 
                    "HSPB1", "HSPE1", "HSPH1", "ID3", "IER2", "JUN", "JUNB", "JUND", "MT1X", "NFKBIA", "NR4A1", 
                    "PPP1R15A", "SOCS3", "ZFP36")
  seurat <- AddModuleScore(seurat, features = list(genes.dissoc), ctrl.size = 20, enrich.name = "genes_dissoc", nbin=12)
  # This adds a "S.Score" column to meta.data
  
  # And save
  saveRDS(seurat, paste0("../Data/scRNASeq/Processed/", sc, ".raw.annot.rds"))
}

#### Plots - initial ####
# Save plots per sample in their own directory ?

for ( sc in sampleList){
  seurat <- readRDS(paste0("../Data/scRNASeq/Processed/", sc, ".raw.annot.rds"))
  
  # Extract some needed values
  counts_per_cell <- Matrix::colSums(seurat@assays$RNA@counts)
  counts_per_gene <- Matrix::rowSums(seurat@assays$RNA@counts)
  genes_per_cell <- Matrix::colSums(seurat@assays$RNA@counts > 0) # count gene only if it has non-zero reads mapped.
  cells_per_gene <- Matrix::rowSums(seurat@assays$RNA@counts > 0) # count number of cells a gene is present in, ie has at least 1 count
  
  # And plot (see Luecken & Theis, Figure 2)
  pdf(paste(outDir,paste0(sc, "_basicInitialStats_Luecken_Figure2.pdf"), sep="/"))
  # A
  hist(counts_per_cell, main=sc, breaks=100, col='cornflowerblue', xlab = "Count depth")
  hist(counts_per_cell, main=sc, xlim=c(0,5000), breaks=1000, col='cornflowerblue', xlab = "Count depth")
  
  # A, but then using ggplot2 to place one figure into another
  counts_per_cell.df <- as.data.frame(counts_per_cell)

  histLarge.df <- hist(counts_per_cell, main=sc, breaks=100, col='cornflowerblue', xlab = "Count depth", plot = FALSE)
  binWidthLarge <- histLarge.df$breaks[2]-histLarge.df$breaks[1]
  breaksLarge <- histLarge.df$breaks
  maxLarge <- max(histLarge.df$counts)
  
  histDetail.df <- hist(counts_per_cell, main=sc, xlim=c(0,5000), breaks=1000, col='cornflowerblue', xlab = "Count depth", plot=FALSE)
  binWidthDetail <- histDetail.df$breaks[2]-histDetail.df$breaks[1]
  breaksDetail <- histDetail.df$breaks
  maxDetail <- max(histDetail.df$counts)
  
  p.overview <- ggplot(counts_per_cell.df,aes(x=counts_per_cell)) +
    geom_histogram(binwidth = binWidthLarge, color='black', fill='cornflowerblue') +
#    scale_x_continuous(breaks = breaksLarge) +
    labs(x="Count depth", y="Frequency") +
    theme_classic()
  p.detail <- ggplot(counts_per_cell.df,aes(x=counts_per_cell)) +
    geom_histogram(binwidth = binWidthDetail, color='black', fill='cornflowerblue') +
    xlim(0,5000) +
#    scale_x_continuous(breaks = breaksDetail) +
    labs(x="Count depth", y="Frequency") +
    theme_classic()
  p.all <- p.overview + annotation_custom(ggplotGrob(p.detail),
                                          xmin = max(counts_per_cell.df)/2, ymin = maxLarge/2,
                                          xmax = max(counts_per_cell.df)) + ggtitle(sc)
  print(p.all)
  rm(counts_per_cell.df)
  
  # B
  hist(genes_per_cell, main=sc, breaks=100, col='cornflowerblue', xlab="Number of genes", xaxt='n')
  axis(side=1, at=c(seq(0,999,250),seq(1000,6000,500)), labels=c(seq(0,999,250),seq(1000,6000,500)), las=2)
  # Q: Where to filter? at 500, 750, 800
  abline(v=c(500,750,800), col="red", lwd=2)
  
  # C
  plot(sort(counts_per_cell, decreasing = TRUE), log='y', main=sc, xlab="Cell rank", ylab="Count depth",
       yaxt='n')
  axis(2, at=c(seq(0,900,100),seq(1000,9000, 1000),seq(10000,70000, 10000)), labels=FALSE)
  labels <- sapply(c(3,4), function(i) as.expression(bquote(10^ .(i))))
  axis(2, at=c(1000, 10000), labels=labels)
  abline(h=1500, col="red", lwd=2) # as in Luecken & Theis
  
  # D
  #Create a function to generate a continuous color palette
  #rbPal <- colorRampPalette(c('red','blue'))
  
  #This adds a column of color values based on the y values
  range(seurat@meta.data$percent.mito)
  # [1] 
  # hist(seurat@meta.data$percent.mito) -> 
  gradCol <- seurat@meta.data$percent.mito
  myCol <- viridis(10)[as.numeric(cut(gradCol,breaks = seq(0,100,10)))]
  layout(matrix(1:2,ncol=2), width = c(3,1),height = c(1,1))
  par(mar=c(5.1,4.1,4.1,0.1))
  plot(counts_per_cell, genes_per_cell, main=sc, log='xy', col=myCol, xlab="Count depth", ylab="Number of genes", las=2, pch=16, cex=0.6)
  legend_image <- as.raster(matrix(rev(viridis(10)), ncol=1))
  par(mar=c(5.1,1,4.1,2.1))
  plot(c(0,2),c(0,1), type = 'n', axes = F, xaxt='n', xlab = '', ylab='', main = '')
  axis(4, at=seq(0,100,10), labels=FALSE, pos=0.5)
  text(x=0.9, y = seq(0,1,0.1), labels = seq(0,100,10), adj=0, cex = 0.6)
  mtext("Percentage mitochondrial counts", side=4, line=-1)
  rasterImage(legend_image, -1, 0, 0.5, 1)
  # coloring from 90-100 is wrong ... white...
  
  # And reset plotting margins,
  par(mar=c(5.1,4.1,4.1,2.1))
  par(mfrow=c(1,1))
  
  # Extra
  plot(sort(genes_per_cell), log='y', main=paste0(sc, ' - genes per cell (ordered)'), xlab="Cell rank", ylab="Number of genes", las=2)
  
  # Mitochondrial genes: > 20% will be filtered out, <= will be left in
  myTable.df <- as.data.frame(table(seurat@meta.data$percent.mito))
  myTable.df$Var1 <- as.numeric(levels(myTable.df$Var1))[myTable.df$Var1]
  p <- hist(seurat@meta.data$percent.mito, main=paste0(sc, ' - Percentage mitochondrial genes expressed'), xlab="Percentage of expression", col='cornflowerblue')
  p
  abline(v=20, col="red", lwd=2)
  abline(v=50, col="red", lwd=2)
  abline(v=80, col="red", lwd=2)
  filteredOut <- sum(myTable.df[myTable.df$Var1 > 20, "Freq"])
  kept <- sum(myTable.df[myTable.df$Var1 <= 20, "Freq"])
  text(25, max(p$counts)-(max(p$counts)*0.1), paste0("----> ", filteredOut), col = "red")
  text(13.50, max(p$counts)-(max(p$counts)*0.1), paste0(kept, "<--"), col = "blue")
  text(30, max(p$counts)/2, "Cutoff: <= 20% should be kept", col = "black", cex = 1.15, adj = 0)
  
  # Housekeeping genes: 55 or more hk genes should be present
  myTable.df <- as.data.frame(table(seurat@meta.data$n.exp.hkgenes))
  myTable.df$Var1 <- as.numeric(levels(myTable.df$Var1))[myTable.df$Var1]
  p <- hist(seurat@meta.data$n.exp.hkgenes, main=paste0(sc, ' - Number of housekeeping genes expressed'), xlab="Number of housekeeping genes expressed", col='cornflowerblue')
  p
  abline(v=55, col="red", lwd=2)
  filteredOut <- sum(myTable.df[myTable.df$Var1 < 55, "Freq"])
  kept <- sum(myTable.df[myTable.df$Var1 >= 55, "Freq"])
  text(62, max(p$counts), paste0("--> ", kept), col = "blue")
  text(48, max(p$counts), paste0(filteredOut, "<----"), col = "red")
  text(48, max(p$counts)/2, "Cutoff: >= 55 should be kept", col = "black", cex = 1.15, adj = 1)
  
  # Cell cycle markers
  p <- FeatureScatter(object = seurat, feature1 = "S.Score", feature2 = "nFeature_RNA") + ggtitle(sc)
  print(p)
  
  dev.off()
  
  # And plot
  pdf(paste(outDir,paste0(sc, "_violinPlot_initial.pdf"), sep="/"))
  Idents(seurat) <- "hash.ID"
  p1 <- VlnPlot(object = seurat, features = c("nFeature_RNA"), pt.size = 0.1) + ggtitle(paste0(sc, " - nFeature_RNA"))
  p2 <- VlnPlot(object = seurat, features = c("nCount_RNA"), pt.size = 0.1)
  #CombinePlots(plots=list(p1,p2), ncol=2, legend = "bottom")
  p <- p1|p2
  print(p)
  p3 <- VlnPlot(object = seurat, features = c("percent.mito"), pt.size = 0.1) + ggtitle(paste0(sc, " - percent mito"))
  p4 <- VlnPlot(object = seurat, features = c("n.exp.hkgenes"), pt.size = 0.1)
  p <- p3|p4
  print(p)
  
  dev.off()
  
  # and clean
  rm(n.expressed.hkgenes, mito.genes, percent.mito, p1,p2,p3,p4,p,p.overview, p.detail, p.all, hkgenes, hkgenes.found, myTable.df)
  rm(counts_per_cell, counts_per_gene, genes_per_cell, cells_per_gene)
  rm(genes.dissoc, g2m.genes, s.genes)
  
  pdf(paste(outDir,paste0(sc, "_nrGenes_vs_Counts_initial.pdf"), sep="/"))
  Idents(seurat) <- "hash.ID"
  p <- FeatureScatter(object = seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + ggtitle(sc)
  print(p)
  dev.off()
  
}

## OBSERVATION
# - For MOMA72 and 302, most cells are in the G1 phase 
#   - Colors jump around, ordering of legend is not the same -> alter in script?
# - MOMA72, 302, for most cells less than 20% of expression is of mitochondrial origin
# 

#### Plots - VDJ ####
# Only MOMA52 and MOMA57 have VDJ data
sampleList <- c("MOMA52", "MOMA57")

dir.create(paste0(outDir,"/tiff"))

for (i in 1:length(pcs)){
  for (sc in sampleList){
    seurat <- readRDS(paste0("../Data/scRNASeq/Processed/", sc, ".raw.annot.rds"))
    
    # Just do a quick and dirty normalizatin etc. to get a PCA or UMAP
    seurat <- SCTransform(seurat, vars.to.regress = "percent.mito")
    seurat <- RunPCA(seurat)
    seurat <- RunUMAP(seurat, dims=1:10)
    
    # Just show all cells that have data, i.e. do not have <NA> in fx. umis_vdj_t or umis_vdj_b
    idx <- which(!is.na(seurat@meta.data$umis_vdj_t))
    selCells_T <- rownames(seurat@meta.data)[idx]
    pdf(paste(outDir, paste0("VDJ_T_cells_dimPlot_rawData.", sc,".pdf"), sep="/"), width = 14, height = 7)
      p <- DimPlot(seurat, cells.highlight = selCells_T, cols.highlight = "blue", cols = "lightgrey", pt.size = 0.5)
      print(p)
    dev.off()
    
    # And as TIFF
    ggsave(paste(outDir, "tiff", paste0("VDJ_T_cells_featurePlot_rawData.", sc,".tiff"), sep="/"), p, scale=2, width=5, height=5, dpi=360,compression="lzw")
    
    # And same for B cells
    idx <- which(!is.na(seurat@meta.data$umis_vdj_b))
    selCells_B <- rownames(seurat@meta.data)[idx]
    pdf(paste(outDir, paste0("VDJ_B_cells_dimPlot_rawData.", sc,".pdf"), sep="/"), width = 14, height = 7)
      p <- DimPlot(seurat, cells.highlight = selCells_B, cols.highlight = "blue", cols = "lightgrey", pt.size = 0.5)
      print(p)
    dev.off()
    
    # And as TIFF
    ggsave(paste(outDir, "tiff", paste0("VDJ_B_cells_featurePlot_rawData.", sc,".tiff"), sep="/"), p, scale=2, width=5, height=5, dpi=360,compression="lzw")
    
    # And combine?
    T_and_B <- intersect(selCells_B, selCells_T)
    # Remove the cells present in both T and B cell sets from the original ones
    selCells_B <- setdiff(selCells_B, selCells_T)
    selCells_B <- selCells_B[!(selCells_B %in% T_and_B)]
    selCells_T <- setdiff(selCells_T, selCells_B)
    selCells_T <- selCells_T[!(selCells_T %in% T_and_B)]
    
    seurat@meta.data$TB <- ifelse(rownames(seurat@meta.data) %in% selCells_B, "B", 
                                  ifelse(rownames(seurat@meta.data) %in% selCells_T, "T",
                                         ifelse(rownames(seurat@meta.data) %in% T_and_B, "T_and_B", "unknown")))
    
    print(sc)
    print(table(seurat@meta.data$TB))
    write.table(as.data.frame(table(seurat@meta.data$TB)), file=paste(outDir, paste0("distribution_T_and_B_cells.",sc,".txt"), sep="/"), sep="\t", quote=FALSE,
                col.names = NA)
    
    # Print the cell names of the cells that have both T and B cells pesent (for later)
    write.table(as.data.frame(rownames(seurat@meta.data[seurat@meta.data$TB == "T_and_B",])), 
                file=paste(outDir, paste0("T_and_B_cells.",sc,".txt"), sep="/"), sep="\t", quote=FALSE,
                col.names = NA)
    
    # and plot, but see https://github.com/satijalab/seurat/issues/3750
    pdf(paste(outDir, paste0("VDJ_All_cells_dimPlot_rawData.", sc,".pdf"), sep="/"), width = 14, height = 7)
      # p <- DimPlot(seurat, cells.highlight = list("T cells only"=selCells_T,"B cells only"=selCells_B, "T and B cells"=T_and_B), cols.highlight = c("blue", "red", "purple"), cols = "lightgrey", pt.size = 0.5) 
      p <- DimPlot(seurat, group.by = "TB", order=c("unknown","B","T","T_and_B"), pt.size = 0.5) 
      p <- p + ggtitle(s)
      print(p)
    dev.off()
    
    # And as TIFF
    ggsave(paste(outDir, "tiff", paste0("VDJ_All_cells_dimPlot_rawData.", sc,".tiff"), sep="/"), p, scale=2, width=5, height=5, dpi=360,compression="lzw")
  }
}

# [1] "MOMA52"
#   B       T T_and_B unknown 
# 442    4017      52    2184 
#
# [1] "MOMA57"
#   B       T T_and_B unknown 
# 434    5088      52    4031 

