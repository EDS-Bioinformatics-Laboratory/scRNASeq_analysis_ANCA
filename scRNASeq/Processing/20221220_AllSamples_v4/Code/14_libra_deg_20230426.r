#### See 14_libra_deg.r ####
# Here, we also include cluster 0 (see mail 20230414)
# And put all comparisons into one directory etc

#### Set default for a good start ####
source("initAnalysis.r")
init()

# Create a directory to hold results and plots
outDir <- "14_libra_deg_20230426"
dir.create(outDir)

#### Library ####
library(Seurat)
library(ggplot2)
library(grid)
library(cowplot)
library(corrplot)
library(tidyverse)
library(edgeR)
library(limma)
library(RColorBrewer)
library(gplots)

if (!require("Libra")){
  # Install Libra
  # https://github.com/neurorestore/Libra
  devtools::install_github("neurorestore/Libra")
  # Had to reinstall 'edgeR' as the DLL could not be loaded:
  #BiocManager::install("edgeR", INSTALL_opts = c("--no-multiarch"), force = TRUE)
  # Obtained the hint of using 'INSTALL_opts' from https://github.com/rstudio/renv/issues/162
  # I did not update any packages when installing Libra
}

library(Libra)
# Loading required package: Libra
# Warning message:
#   In checkMatrixPackageVersion() : Package version inconsistency detected.
# TMB was built with Matrix version 1.3.3
# Current Matrix version is 1.5.1
# Please re-install 'TMB' from source using install.packages('TMB', type = 'source') or ask CRAN for a binary version of 'TMB' matching CRAN's 'Matrix' package

# 20230309 - AJ:  So, reinstalled TMB and consequently also glmmTMB ... gRRR
# install.packages("TMB", type = "source")
# install.packages("glmmTMB", type = "source")


#### Only use limma/voom ####
# We decide (09032023) to opt for limma/voom (and not edgeR), bacause that is what we generally use
# All comparisons will be made:
# - between HC, ANCA, SLE and
# - between ANCA_MPO and ANCA_PR3 (as performed above)
# - And between clusters (within sample, between samples???)
#
# Can (and should) we do all at once??
# - Get all cluster info in the sample names etc...
# - And we use 'duplicateCorrelation' to inlcude the pairing bewteen samples accross the clusters

# Most handy to do this using a loop?
library(limma)
library(RColorBrewer)
library(gplots)

seurat <- readRDS("10_Reclustering_MP_and_MC/selectedClusters.subjectIntegrated.reclustered.kidney_immune_reference.rds")
DefaultAssay(seurat) <- "RNA"

# Extract the pseudobulk data
matrices = to_pseudobulk(seurat, replicate_col = "orig.ident", cell_type_col = "integrated_snn_res.1", label_col = "Condition",
                         min_reps = 1)

# Get pseudobulk matrices for the following cell types
names(matrices)
#  [1] "0"  "1"  "10" "2"  "3"  "4"  "5"  "6"  "7"  "8"  "9" 

# We are interested in clusters 0, 1, 4 and 8
counts_cluster0 <- matrices$`0`
colnames(counts_cluster0) <- paste0(colnames(counts_cluster0), "_0")
counts_cluster1 <- matrices$`1`
colnames(counts_cluster1) <- paste0(colnames(counts_cluster1), "_1")
counts <- cbind(counts_cluster0, counts_cluster1)
counts_cluster4 <- matrices$`4`
colnames(counts_cluster4) <- paste0(colnames(counts_cluster4), "_4")
counts <- cbind(counts, counts_cluster4)
counts_cluster8 <- matrices$`8`
colnames(counts_cluster8) <- paste0(colnames(counts_cluster8), "_8")
counts <- cbind(counts, counts_cluster8)

# And clean
rm(counts_cluster0, counts_cluster1, counts_cluster4, counts_cluster8, matrices)

# Create the new sample info table from the old seurat@meta.data
colData <- seurat@meta.data[, c("orig.ident","Condition","subType")]
colData <- colData[!duplicated(colData),]
colData <- bind_rows(replicate(4, colData, simplify = FALSE))
colData$cluster <- c(rep(0, 7), rep(1, 7), rep(4, 7), rep(8, 7))
rownames(colData) <- paste(paste(colData$orig.ident,colData$Condition, sep=":"), colData$cluster, sep="_")
colData <- colData[colnames(counts),]

y = DGEList(counts=counts,group=with(colData, paste(orig.ident, Condition, subType, cluster,sep="_")))
saveRDS(y,paste(outDir,"y.raw.rds", sep="/"))

# Put in some annotation
y$targets <- colData

geneSet <- rownames(y$counts)

# Copy helper script to directory and source from there
file.copy(paste(dropbox,"Support/R/annotateD.r",sep="/"),scriptDir)
source(paste(scriptDir,"annotateD.r", sep="/"))

# Get annotation
y <- annotateD(y,geneSet,species="human",geneFilter="hgnc_symbol")

# # Annotation:
# +++++++++++
#   #genes:	 22769 
#   
#   Retrieved:
#   #symbols    :	 22769 ( 100 %)
#   #ensembl ids:	 17440 ( 76.6 %)
#   #entrez  ids:	 16511 ( 72.5 %)

# And save
saveRDS(y,paste(outDir,"y.annot.rds", sep="/"))

#### Design, filtering and normalization ####
# Read the data
y <- readRDS(paste(outDir,"y.annot.rds", sep="/"))
dim(y)
# [1] 22769    28

# Design
condition = paste(factor(y$targets$Condition),factor(y$targets$subType), factor(y$targets$cluster), sep = "_")
subjects = factor(y$targets$orig.ident)

#
# design <- model.matrix(~ 0 + condition + subject)
# This will not work as we only have one sample for HC and SLE
# - > duplicateCorrelation??
design <- model.matrix(~ 0 + condition)

# Clean up the columnnames of the design
colnames(design) = gsub("condition","",colnames(design))

# The 'design' is used to determine the MinSampleSize:
#     h <- hat(design)           # 'hat' is a function of 'stats' and part of base-R
#     MinSampleSize <- 1/max(h)
# MinSampleSize
#[1] 1

keep <- filterByExpr(y, design, min.count = 2)
y <- y[keep,]
dim(y)
# [1] 12990    28

# Normalize
y = calcNormFactors(y)

# And save
saveRDS(y,paste(outDir,"y.annot.filt.norm.rds", sep="/"))

#### LibSizes, Heatmap and MDS plot ####
write.table(y$samples, file = paste(outDir,"normFactors.txt", sep="/"), sep="\t", col.names=NA, row.names=TRUE, quote=FALSE)

png(paste(outDir,"libSizes_merged.png", sep="/"),width=720,height=720,pointsize=24)
par(las=2) # make label text perpendicular to axis
par(mar=c(8,5,4,2)) # increase y-axis margin.
cols = ifelse(y$samples$lib.size > 2e+06, "grey","red")
p <- barplot(y$samples$lib.size, col=cols, main="LibSize Distribution", las=2, names.arg=y$samples$group, cex.names=0.4)
p
abline(h=2e+06, col="red")
abline(h=median(y$samples$lib.size), col="red", lty = 2)
text(p[1]+0.5,median(y$samples$lib.size),expression(italic("median")),srt=0.2,pos=3, cex=0.7)
dev.off()

opar <- par(no.readonly = TRUE)# Save original "par" setting

treatmentFactors = factor(paste0(y$targets$Condition, "_", y$targets$subType), levels = c("ANCA_PR3", "ANCA_MPO", "SLE_NONE", "Control_NONE"))
treatmentColors = c("cyan", "blue", "orange", "red")

colors <- treatmentColors[treatmentFactors]

par(mar=c(13.1, 4.1, 4.1, 8.1), xpd=TRUE)
par(mgp = c(3, 1, 0))
png(paste(outDir,"plotMDS_initial.png", sep="/"), width=720,height=720,pointsize=24)
myPlot <- plotMDS(y,top=500,labels="",col=colors, cex=1.0)
plot(myPlot,xlab="Leading logFC dim 1", ylab="Leading logFC dim 2",
     xlim=c(min(myPlot$x)-0.5,max(myPlot$x)+1),ylim=c(min(myPlot$y)-0.5,max(myPlot$y)+1),pch=16,col=colors, cex=1.0)
# add legend
legend("topright", col=treatmentColors, 
       legend=levels(treatmentFactors), 
       pch = 16, cex = 0.6, title=expression(bold(Condition)))
dev.off()
par(opar)

opar <- par(no.readonly = TRUE)# Save original "par" setting

par(mar=c(13.1, 4.1, 4.1, 8.1), xpd=TRUE)
par(mgp = c(3, 1, 0))
png(paste(outDir,"plotMDS_initial.labelled.png", sep="/"), width=720,height=720,pointsize=24)
myPlot <- plotMDS(y,top=500,labels=y$targets$Replicate,col=colors, cex=1.0)
plot(myPlot,xlab="Leading logFC dim 1", ylab="Leading logFC dim 2",
     xlim=c(min(myPlot$x)-0.5,max(myPlot$x)+1),ylim=c(min(myPlot$y)-0.5,max(myPlot$y)+1),pch=16,col=colors, cex=1.0)
# add labels
for ( i in 1:length(myPlot$x)){
  #  if (myPlot$x[i] < -3){
  text(myPlot$x[i],myPlot$y[i],labels=paste0(y$targets$cluster[i],":",y$targets$orig.ident[i]), pos=4, cex=0.6, col="black")
  #  }
}
# add legend
legend("topright", col=treatmentColors, 
       legend=levels(treatmentFactors), 
       pch = 16, cex = 0.6, title=expression(bold(Condition)))
dev.off()
par(opar)

## Heatmap ##
# Construct a heatmap of individual RNA-seq samples
# - the edgeR user's guide suggests using moderated log-counts-per-million 
cpm.y = cpm(y, prior.count=2, log=TRUE)
hmcol = colorRampPalette(brewer.pal(9, "GnBu"))(100)
dists = dist(t(cpm.y))
mat   = as.matrix(dists)
rownames(mat) = colnames(mat) = with(y$targets,paste(orig.ident, Condition, cluster ,sep=":"))
png(paste(outDir,"heatmap_initial.png", sep="/"), width=960, height=720)
heatmap.2(mat, trace="none", col = rev(hmcol), margin=c(15,15),cexRow = 1.0, cexCol = 1.0 )
dev.off()

# Use NMF package for heatmaps
library(NMF)
library(grid)
Condition <- factor(y$targets$Condition, levels=c("ANCA", "SLE", "Control"))
SubType <- factor(y$targets$subType, levels=c("PR3", "MPO", "NONE"))
Cluster <- factor(y$targets$cluster, levels=c("0", "1", "4", "8"))

annotation <- data.frame(Condition, SubType, Cluster)

conditionCol <- c("blue","red", "black")
subtypeCol <- c("cyan","green", "gray")
clusterCol <- c("forestgreen","magenta", "orange", "yellow")

annColors = list(Condition = conditionCol, SubType = subtypeCol, Cluster = clusterCol )
# Width and height are in inches !!
# Default Powerpoint size are 10 * 7.5 (or 960 * 720 in pixels)
aheatmap(mat,col = rev(hmcol), width=10, height=15.0, annCol = annotation,
         annColors = annColors, gp = gpar(cex = 1.1),
         filename = paste(outDir,"aheatmap_initial.png", sep="/"))

#### Run Voom with duplicateCorrelation ####
v <- voom(y, design, plot=TRUE, save.plot=TRUE)

# plot
png(paste(outDir, "All_voomMeanVarianceTrend_beforeDupCor.png", sep="/"), width=720, height=720)
plot(v$voom.xy$x,v$voom.xy$y, xlab=v$voom.xy$xlab,ylab=v$voom.xy$ylab, cex=0.25, pch=16)
title("voom: Mean-variance trend")
lines(v$voom.line$x,v$voom.line$y, col="red")
dev.off()

# duplicate correlation to make the most of the pairing
corfit <- duplicateCorrelation(v, design, block=subjects)
corfit$consensus.correlation
# [1] 0.116515  

saveRDS(corfit,paste(outDir,"corfit.AllComparisons.rds", sep="/"))

# And rerun Voom
v <- voom(y, design, block=subjects, correlation=corfit$consensus.correlation, plot=TRUE, save.plot=TRUE)

# plot
png(paste(outDir, "All_voomMeanVarianceTrend.png", sep="/"), width=720, height=720)
plot(v$voom.xy$x,v$voom.xy$y, xlab=v$voom.xy$xlab,ylab=v$voom.xy$ylab, cex=0.25, pch=16)
title("voom: Mean-variance trend")
lines(v$voom.line$x,v$voom.line$y, col="red")
dev.off()
# 

# And save for the geneset enrichment analysis 
v$targets <- cbind(v$targets, colData)
v$distanceMatrix <- mat
v$design <- design

# And save
saveRDS(v, paste(outDir,"v.AllComparisons.rds", sep="/"))

y.estDisp <- estimateDisp(y, design, robust=TRUE)
y.estDisp$common.dispersion
# [1] 0.126

png(paste(outDir, "plotBCV.png", sep="/"), width=720, height=720)
plotBCV(y.estDisp)
dev.off()

fit <- glmQLFit(y.estDisp,design,robust=TRUE)

png(paste(outDir, "plotQLDisp.png", sep="/"), width=720, height=720)
plotQLDisp(fit)
# Warning message:
#   In sqrt(s2) : NaNs produced
dev.off()

# Perform the fit and save
fit <- lmFit(v, design,block = subjects, correlation = corfit$consensus.correlation)
saveRDS(fit, paste(outDir,"fit.AllComparisons.rds",sep="/"))

fit.eb <- eBayes(fit)
options(digits=3)

topTable(fit.eb,coef=2,n=50)
# Lot of RPL, MT genes, MALAT1, EEF1A1, FTH1 ... :-( ?

fit.eb$df.prior
# [1] 14.2 -- Large ...

summary(decideTests(fit.eb))
#        ANCA_MPO_0 ANCA_MPO_1 ANCA_MPO_4 ANCA_MPO_8 ANCA_PR3_0 ANCA_PR3_1 ANCA_PR3_4 ANCA_PR3_8 Control_NONE_0 Control_NONE_1
# Down            0          0          0          0          0          0          0          0              0              0
# NotSig       2288       2480       1605       2044         18       1676        995       1060           3175           3236
# Up          10702      10510      11385      10946      12972      11314      11995      11930           9815           9754
#        Control_NONE_4 Control_NONE_8 SLE_NONE_0 SLE_NONE_1 SLE_NONE_4 SLE_NONE_8
# Down                0              0         25          0          0          0
# NotSig           2685           3227       2857       2755       2752       3116
# Up              10305           9763      10108      10235      10238       9874

# Only Up ....

rm(fit.eb)

#### Contrasts and topTable ####
my.contrasts = makeContrasts(
  # Within cluster
  C0_ANCA_vs_CONTROL = (ANCA_MPO_0 + ANCA_PR3_0)/2 - Control_NONE_0,
  C0_SLE_vs_CONTROL = SLE_NONE_0 - Control_NONE_0,
  C0_SLE_vs_ANCA = SLE_NONE_0 - (ANCA_MPO_0 + ANCA_PR3_0)/2,
  C0_ANCA_MPO_vs_ANCA_PR3 = ANCA_MPO_0 - ANCA_PR3_0,
  C1_ANCA_vs_CONTROL = (ANCA_MPO_1 + ANCA_PR3_1)/2 - Control_NONE_1,
  C1_SLE_vs_CONTROL = SLE_NONE_1 - Control_NONE_1,
  C1_SLE_vs_ANCA = SLE_NONE_1 - (ANCA_MPO_1 + ANCA_PR3_1)/2,
  C1_ANCA_MPO_vs_ANCA_PR3 = ANCA_MPO_1 - ANCA_PR3_1,
  C4_ANCA_vs_CONTROL = (ANCA_MPO_4 + ANCA_PR3_4)/2 - Control_NONE_4,
  C4_SLE_vs_CONTROL = SLE_NONE_4 - Control_NONE_4,
  C4_SLE_vs_ANCA = SLE_NONE_4 - (ANCA_MPO_4 + ANCA_PR3_4)/2,
  C4_ANCA_MPO_vs_ANCA_PR3 = ANCA_MPO_4 - ANCA_PR3_4,
  C8_ANCA_vs_CONTROL = (ANCA_MPO_8 + ANCA_PR3_8)/2 - Control_NONE_8,
  C8_SLE_vs_CONTROL = SLE_NONE_8 - Control_NONE_8,
  C8_SLE_vs_ANCA = SLE_NONE_8 - (ANCA_MPO_8 + ANCA_PR3_8)/2,
  C8_ANCA_MPO_vs_ANCA_PR3 = ANCA_MPO_8 - ANCA_PR3_8,
  # Between clusters
  # Correct for their baseline, i.e. CONTROL??
  # Are we really interested in these??
  ANCA_C1_vs_ANCA_C0 = (ANCA_MPO_1 + ANCA_PR3_1)/2 - (ANCA_MPO_0 + ANCA_PR3_0)/2,
  ANCA_C4_vs_ANCA_C0 = (ANCA_MPO_4 + ANCA_PR3_4)/2 - (ANCA_MPO_0 + ANCA_PR3_0)/2,
  ANCA_C8_vs_ANCA_C0 = (ANCA_MPO_8 + ANCA_PR3_8)/2 - (ANCA_MPO_0 + ANCA_PR3_0)/2,
  ANCA_C4_vs_ANCA_C1 = (ANCA_MPO_4 + ANCA_PR3_4)/2 - (ANCA_MPO_1 + ANCA_PR3_1)/2,
  ANCA_C8_vs_ANCA_C1 = (ANCA_MPO_8 + ANCA_PR3_8)/2 - (ANCA_MPO_1 + ANCA_PR3_1)/2,
  ANCA_C8_vs_ANCA_C4 = (ANCA_MPO_8 + ANCA_PR3_8)/2 - (ANCA_MPO_4 + ANCA_PR3_4)/2,
  SLE_C1_vs_SLE_C0 = SLE_NONE_1 - SLE_NONE_0,
  SLE_C4_vs_SLE_C0 = SLE_NONE_4 - SLE_NONE_0,
  SLE_C8_vs_SLE_C0 = SLE_NONE_8 - SLE_NONE_0,
  SLE_C4_vs_SLE_C1 = SLE_NONE_4 - SLE_NONE_1,
  SLE_C8_vs_SLE_C1 = SLE_NONE_8 - SLE_NONE_1,
  SLE_C8_vs_SLE_C4 = SLE_NONE_8 - SLE_NONE_4,
  CONTROL_C1_vs_CONTROL_C0 = Control_NONE_1 - Control_NONE_0,
  CONTROL_C4_vs_CONTROL_C0 = Control_NONE_4 - Control_NONE_0,
  CONTROL_C8_vs_CONTROL_C0 = Control_NONE_8 - Control_NONE_0,
  CONTROL_C4_vs_CONTROL_C1 = Control_NONE_4 - Control_NONE_1,
  CONTROL_C8_vs_CONTROL_C1 = Control_NONE_8 - Control_NONE_1,
  CONTROL_C8_vs_CONTROL_C4 = Control_NONE_8 - Control_NONE_4,
  ANCA_MPO_C1_vs_ANCA_MPO_C0 = ANCA_MPO_1 - ANCA_MPO_0,
  ANCA_MPO_C4_vs_ANCA_MPO_C0 = ANCA_MPO_4 - ANCA_MPO_0,
  ANCA_MPO_C8_vs_ANCA_MPO_C0 = ANCA_MPO_8 - ANCA_MPO_0,
  ANCA_MPO_C4_vs_ANCA_MPO_C1 = ANCA_MPO_4 - ANCA_MPO_1,
  ANCA_MPO_C8_vs_ANCA_MPO_C1 = ANCA_MPO_8 - ANCA_MPO_1,
  ANCA_MPO_C8_vs_ANCA_MPO_C4 = ANCA_MPO_8 - ANCA_MPO_4,
  ANCA_PR3_C1_vs_ANCA_PR3_C0 = ANCA_PR3_1 - ANCA_PR3_0,
  ANCA_PR3_C4_vs_ANCA_PR3_C0 = ANCA_PR3_4 - ANCA_PR3_0,
  ANCA_PR3_C8_vs_ANCA_PR3_C0 = ANCA_PR3_8 - ANCA_PR3_0,
  ANCA_PR3_C4_vs_ANCA_PR3_C1 = ANCA_PR3_4 - ANCA_PR3_1,
  ANCA_PR3_C8_vs_ANCA_PR3_C1 = ANCA_PR3_8 - ANCA_PR3_1,
  ANCA_PR3_C8_vs_ANCA_PR3_C4 = ANCA_PR3_8 - ANCA_PR3_4,
  # And the complete clusters relative to each other
  C1_vs_C0 = (ANCA_MPO_1 + ANCA_PR3_1 + SLE_NONE_1 + Control_NONE_1 )/4 - (ANCA_MPO_0 + ANCA_PR3_0 + SLE_NONE_0 + Control_NONE_0)/4,
  C4_vs_C0 = (ANCA_MPO_4 + ANCA_PR3_4 + SLE_NONE_4 + Control_NONE_4 )/4 - (ANCA_MPO_0 + ANCA_PR3_0 + SLE_NONE_0 + Control_NONE_0)/4,
  C8_vs_C0 = (ANCA_MPO_8 + ANCA_PR3_8 + SLE_NONE_8 + Control_NONE_8 )/4 - (ANCA_MPO_0 + ANCA_PR3_0 + SLE_NONE_0 + Control_NONE_0)/4,
  C4_vs_C1 = (ANCA_MPO_4 + ANCA_PR3_4 + SLE_NONE_4 + Control_NONE_4 )/4 - (ANCA_MPO_1 + ANCA_PR3_1 + SLE_NONE_1 + Control_NONE_1)/4,
  C8_vs_C1 = (ANCA_MPO_8 + ANCA_PR3_8 + SLE_NONE_8 + Control_NONE_8 )/4 - (ANCA_MPO_1 + ANCA_PR3_1 + SLE_NONE_1 + Control_NONE_1)/4,
  C8_vs_C4 = (ANCA_MPO_8 + ANCA_PR3_8 + SLE_NONE_8 + Control_NONE_8 )/4 - (ANCA_MPO_4 + ANCA_PR3_4 + SLE_NONE_4 + Control_NONE_4)/4,
  # And between cluster - 53:56
  C0_vs_C1_C4andC8 = (ANCA_MPO_0 + ANCA_PR3_0 + SLE_NONE_0 + Control_NONE_0)/4 - (ANCA_MPO_1 + ANCA_PR3_1 + SLE_NONE_1 + Control_NONE_1 + 
                                                                                    ANCA_MPO_4 + ANCA_PR3_4 + SLE_NONE_4 + Control_NONE_4 + 
                                                                                    ANCA_MPO_8 + ANCA_PR3_8 + SLE_NONE_8 + Control_NONE_8)/12,
  C1_vs_C0_C4andC8 = (ANCA_MPO_1 + ANCA_PR3_1 + SLE_NONE_1 + Control_NONE_1)/4 - (ANCA_MPO_0 + ANCA_PR3_0 + SLE_NONE_0 + Control_NONE_0 +
                                                                                    ANCA_MPO_4 + ANCA_PR3_4 + SLE_NONE_4 + Control_NONE_4 + 
                                                                                    ANCA_MPO_8 + ANCA_PR3_8 + SLE_NONE_8 + Control_NONE_8)/12,
  C4_vs_C0_C1andC8 = (ANCA_MPO_4 + ANCA_PR3_4 + SLE_NONE_4 + Control_NONE_4)/4 - (ANCA_MPO_0 + ANCA_PR3_0 + SLE_NONE_0 + Control_NONE_0 +
                                                                                    ANCA_MPO_1 + ANCA_PR3_1 + SLE_NONE_1 + Control_NONE_1 + 
                                                                                    ANCA_MPO_8 + ANCA_PR3_8 + SLE_NONE_8 + Control_NONE_8)/12,
  C8_vs_C0_C1andC4 = (ANCA_MPO_8 + ANCA_PR3_8 + SLE_NONE_8 + Control_NONE_8)/4 - (ANCA_MPO_0 + ANCA_PR3_0 + SLE_NONE_0 + Control_NONE_0 +
                                                                                    ANCA_MPO_1 + ANCA_PR3_1 + SLE_NONE_1 + Control_NONE_1 +
                                                                                    ANCA_MPO_4 + ANCA_PR3_4 + SLE_NONE_4 + Control_NONE_4)/12,
  # 57:58
  C0_vs_C4andC8 = (ANCA_MPO_0 + ANCA_PR3_0 + SLE_NONE_0 + Control_NONE_0)/4 - (ANCA_MPO_4 + ANCA_PR3_4 + SLE_NONE_4 + Control_NONE_4 + ANCA_MPO_8 + ANCA_PR3_8 + SLE_NONE_8 + Control_NONE_8)/8,
  C1_vs_C4andC8 = (ANCA_MPO_1 + ANCA_PR3_1 + SLE_NONE_1 + Control_NONE_1)/4 - (ANCA_MPO_4 + ANCA_PR3_4 + SLE_NONE_4 + Control_NONE_4 + ANCA_MPO_8 + ANCA_PR3_8 + SLE_NONE_8 + Control_NONE_8)/8,
  # 59:60
  C4_vs_C0andC1 = (ANCA_MPO_4 + ANCA_PR3_4 + SLE_NONE_4 + Control_NONE_4)/4 - (ANCA_MPO_0 + ANCA_PR3_0 + SLE_NONE_0 + Control_NONE_0 + ANCA_MPO_1 + ANCA_PR3_1 + SLE_NONE_1 + Control_NONE_1)/8,
  C8_vs_C0andC1 = (ANCA_MPO_8 + ANCA_PR3_8 + SLE_NONE_8 + Control_NONE_8)/4 - (ANCA_MPO_0 + ANCA_PR3_0 + SLE_NONE_0 + Control_NONE_0 + ANCA_MPO_1 + ANCA_PR3_1 + SLE_NONE_1 + Control_NONE_1)/8,
   # 61
  C0andC1_vs_C4andC8 = (ANCA_MPO_0 + ANCA_PR3_0 + SLE_NONE_0 + Control_NONE_0 + ANCA_MPO_1 + ANCA_PR3_1 + SLE_NONE_1 + Control_NONE_1)/8 - 
    (ANCA_MPO_4 + ANCA_PR3_4 + SLE_NONE_4 + Control_NONE_4 + ANCA_MPO_8 + ANCA_PR3_8 + SLE_NONE_8 + Control_NONE_8)/8,
  
  levels = design
)

fit2 = contrasts.fit(fit,my.contrasts)
fit.eb = eBayes(fit2)
options(digits=3)

## Construct topTable and write results ##
# The AveExpr is the same for every contrast, so add it here and delete it when adding the data of every contrast to the table:
tt   = topTable(fit.eb,n=Inf,sort.by="none")

for (i in 1:ncol(my.contrasts)){
  tt.i = topTable(fit.eb,coef=i,n=Inf,sort.by="none")[,c(6,8:11)]
  names(tt.i) = paste0(names(tt.i)," (",colnames(my.contrasts)[i],")")
  tt = cbind(tt,tt.i)
}
saveRDS(tt,paste(outDir,"tt.AllComparisons.rds", sep="/"))

# And add the CPM values
cpm.y = cpm(y, log=TRUE)
colnames(cpm.y) = paste0("CPM_",y$samples$group)
tt = data.frame(tt,cpm.y,check.names=FALSE)

saveRDS(tt,paste(outDir,"tt_all.AllComparisons.rds", sep="/"))
write.table(tt, paste(outDir,"topTable.AllComparisons.txt", sep="/"), quote=FALSE, sep="\t", row.names = TRUE, col.names = NA)

# And save the voom object with the contrasts
v$contrasts <- my.contrasts
v$design <- design
saveRDS(v, paste(outDir,"v.AllComparisons.rds", sep="/"))

## Voom reports
file.copy(paste(dropbox,"Support/R/voomReports.r",sep="/"), scriptDir)
source(paste(scriptDir,"voomReports.r",sep="/"))

v <- readRDS(paste(outDir,"v.AllComparisons.rds", sep="/"))
targets <- v$targets
design <- v$design
my.contrasts <- v$contrasts

fit <- readRDS(paste(outDir,"fit.AllComparisons.rds", sep="/"))

voomReports(fit, contrasts = my.contrasts)

# Relocate 'voom' directory
file.rename("voom", to=paste0(outDir,"/voom"))

#### Genes - Venn Diagrams ####
file.copy(paste(dropbox,"Support/R/vennReports.r",sep="/"),scriptDir)
source(paste(scriptDir,"vennReports.r", sep="/"))

# Generate a Venn diagram - 
# Below is example code from other scripts - Does this work? TEST!!!
# Which comparisons do we want to see??
v <- readRDS(paste(outDir,"v.AllComparisons.rds", sep="/"))
tt <- readRDS(paste(outDir,"tt_all.AllComparisons.rds", sep="/"))
targets <- v$targets
design <- v$design
my.contrasts <- v$contrasts
fit <- readRDS(paste(outDir,"fit.AllComparisons.rds", sep="/"))

pVal = 0.05
lfcVal = 0.0
adjustMethode = "BH"

vennReports(fit,contrasts = my.contrasts,reportsDir = "venn", myTitle="vennDiagram_Cluster0",
            p.value = pVal, lfc = lfcVal, adjust.method = adjustMethode,
            indices = c(1:4), col = c("red","blue"), scriptDir = ".")

vennReports(fit,contrasts = my.contrasts,reportsDir = "venn", myTitle="vennDiagram_Cluster1",
            p.value = pVal, lfc = lfcVal, adjust.method = adjustMethode,
            indices = c(5:8), col = c("red","blue"), scriptDir = ".")

vennReports(fit,contrasts = my.contrasts,reportsDir = "venn", myTitle="vennDiagram_Cluster4",
            p.value = pVal, lfc = lfcVal, adjust.method = adjustMethode,
            indices = c(9:12), col = c("red","blue"), scriptDir = ".")

vennReports(fit,contrasts = my.contrasts,reportsDir = "venn", myTitle="vennDiagram_Cluster8",
            p.value = pVal, lfc = lfcVal, adjust.method = adjustMethode,
            indices = c(13:16), col = c("red","blue"), scriptDir = ".")

vennReports(fit,contrasts = my.contrasts,reportsDir = "venn", myTitle="vennDiagram_ANCA_MPO_vs_PR3",
            p.value = pVal, lfc = lfcVal, adjust.method = adjustMethode,
            indices = c(4,8,12,16), col = c("red","blue"), scriptDir = ".")

vennReports(fit,contrasts = my.contrasts,reportsDir = "venn", myTitle="vennDiagram_ANCA_betweenClusters_part1",
            p.value = pVal, lfc = lfcVal, adjust.method = adjustMethode,
            indices = c(17:19), col = c("red","blue"), scriptDir = ".")

vennReports(fit,contrasts = my.contrasts,reportsDir = "venn", myTitle="vennDiagram_ANCA_betweenClusters_part2",
            p.value = pVal, lfc = lfcVal, adjust.method = adjustMethode,
            indices = c(20:22), col = c("red","blue"), scriptDir = ".")

vennReports(fit,contrasts = my.contrasts,reportsDir = "venn", myTitle="vennDiagram_SLE_betweenClusters_part1",
            p.value = pVal, lfc = lfcVal, adjust.method = adjustMethode,
            indices = c(23:25), col = c("red","blue"), scriptDir = ".")

vennReports(fit,contrasts = my.contrasts,reportsDir = "venn", myTitle="vennDiagram_SLE_betweenClusters_part2",
            p.value = pVal, lfc = lfcVal, adjust.method = adjustMethode,
            indices = c(26:28), col = c("red","blue"), scriptDir = ".")

vennReports(fit,contrasts = my.contrasts,reportsDir = "venn", myTitle="vennDiagram_CONTROLE_betweenClusters_part1",
            p.value = pVal, lfc = lfcVal, adjust.method = adjustMethode,
            indices = c(29:31), col = c("red","blue"), scriptDir = ".")

vennReports(fit,contrasts = my.contrasts,reportsDir = "venn", myTitle="vennDiagram_CONTROLE_betweenClusters_part2",
            p.value = pVal, lfc = lfcVal, adjust.method = adjustMethode,
            indices = c(32:34), col = c("red","blue"), scriptDir = ".")

vennReports(fit,contrasts = my.contrasts,reportsDir = "venn", myTitle="vennDiagram_ANCA_MPO_betweenClusters_part1",
            p.value = pVal, lfc = lfcVal, adjust.method = adjustMethode,
            indices = c(35:37), col = c("red","blue"), scriptDir = ".")

vennReports(fit,contrasts = my.contrasts,reportsDir = "venn", myTitle="vennDiagram_ANCA_MPO_betweenClusters_part2",
            p.value = pVal, lfc = lfcVal, adjust.method = adjustMethode,
            indices = c(38:40), col = c("red","blue"), scriptDir = ".")

vennReports(fit,contrasts = my.contrasts,reportsDir = "venn", myTitle="vennDiagram_ANCA_PR3_betweenClusters_part1",
            p.value = pVal, lfc = lfcVal, adjust.method = adjustMethode,
            indices = c(41:43), col = c("red","blue"), scriptDir = ".")

vennReports(fit,contrasts = my.contrasts,reportsDir = "venn", myTitle="vennDiagram_ANCA_PR3_betweenClusters_part2",
            p.value = pVal, lfc = lfcVal, adjust.method = adjustMethode,
            indices = c(44:46), col = c("red","blue"), scriptDir = ".")

vennReports(fit,contrasts = my.contrasts,reportsDir = "venn", myTitle="vennDiagram_Clusters_part1",
            p.value = pVal, lfc = lfcVal, adjust.method = adjustMethode,
            indices = c(47:49), col = c("red","blue"), scriptDir = ".")

vennReports(fit,contrasts = my.contrasts,reportsDir = "venn", myTitle="vennDiagram_Clusters_part1",
            p.value = pVal, lfc = lfcVal, adjust.method = adjustMethode,
            indices = c(50:52), col = c("red","blue"), scriptDir = ".")

vennReports(fit,contrasts = my.contrasts,reportsDir = "venn", myTitle="vennDiagram_Clusters_One_vs_Rest",
            p.value = pVal, lfc = lfcVal, adjust.method = adjustMethode,
            indices = c(53:56), col = c("red","blue"), scriptDir = ".")

vennReports(fit,contrasts = my.contrasts,reportsDir = "venn", myTitle="vennDiagram_Cluster_C0_or_C1_vs_C4andC8",
            p.value = pVal, lfc = lfcVal, adjust.method = adjustMethode,
            indices = c(57:58), col = c("red","blue"), scriptDir = ".")

vennReports(fit,contrasts = my.contrasts,reportsDir = "venn", myTitle="vennDiagram_Cluster_C4_or_C8_vs_C0andC1",
            p.value = pVal, lfc = lfcVal, adjust.method = adjustMethode,
            indices = c(59:60), col = c("red","blue"), scriptDir = ".")

vennReports(fit,contrasts = my.contrasts,reportsDir = "venn", myTitle="vennDiagram_Cluster_C0andC1_vs_C4andC8",
            p.value = pVal, lfc = lfcVal, adjust.method = adjustMethode,
            indices = c(61), col = c("red","blue"), scriptDir = ".")

# Relocate 'venn' directory
file.rename("venn", to=paste0(outDir,"/venn"))

# Relocate figures and stats
toBeMoved <- list.files(".", pattern="*.png|*.html", full.names = TRUE)
file.copy(toBeMoved, to=outDir)
file.remove(toBeMoved)

#### Enrichment analysis #####
library(GSEABase)

v <- readRDS(paste(outDir,"v.AllComparisons.rds", sep="/"))
targets <- v$targets
design <- v$design
my.contrasts <- v$contrasts
fit <- readRDS(paste(outDir,"fit.AllComparisons.rds", sep="/"))

## read the MSigDB database v2023.1.Hs genesets
hBroadSets = getGmt(con=paste(dropbox,"Support/MSigDB/v2023.1.Hs/h.all.v2023.1.Hs.entrez.gmt",sep="/"),geneIdType=EntrezIdentifier(), collectionType=BroadCollection(category="h"))
c1BroadSets = getGmt(con=paste(dropbox,"Support/MSigDB/v2023.1.Hs/c1.all.v2023.1.Hs.entrez.gmt",sep="/"),geneIdType=EntrezIdentifier(), collectionType=BroadCollection(category="c1"))
c2BroadSets = getGmt(con=paste(dropbox,"Support/MSigDB/v2023.1.Hs/c2.all.v2023.1.Hs.entrez.gmt",sep="/"),geneIdType=EntrezIdentifier(), collectionType=BroadCollection(category="c2"))
c3BroadSets = getGmt(con=paste(dropbox,"Support/MSigDB/v2023.1.Hs/c3.all.v2023.1.Hs.entrez.gmt",sep="/"),geneIdType=EntrezIdentifier(), collectionType=BroadCollection(category="c3"))
c5BroadSets = getGmt(con=paste(dropbox,"Support/MSigDB/v2023.1.Hs/c5.all.v2023.1.Hs.entrez.gmt",sep="/"),geneIdType=EntrezIdentifier(), collectionType=BroadCollection(category="c5"))
c6BroadSets = getGmt(con=paste(dropbox,"Support/MSigDB/v2023.1.Hs/c6.all.v2023.1.Hs.entrez.gmt",sep="/"),geneIdType=EntrezIdentifier(), collectionType=BroadCollection(category="c6"))
c7BroadSets = getGmt(con=paste(dropbox,"Support/MSigDB/v2023.1.Hs/c7.all.v2023.1.Hs.entrez.gmt",sep="/"),geneIdType=EntrezIdentifier(), collectionType=BroadCollection(category="c7"))
c8BroadSets = getGmt(con=paste(dropbox,"Support/MSigDB/v2023.1.Hs/c8.all.v2023.1.Hs.entrez.gmt",sep="/"),geneIdType=EntrezIdentifier(), collectionType=BroadCollection(category="c7"))
H = geneIds(hBroadSets)
names(H) = paste0("H_",names(H))
C1 = geneIds(c1BroadSets)
names(C1) = paste0("C1_",names(C1))
C2 = geneIds(c2BroadSets)
names(C2) = paste0("C2_",names(C2))
C3 = geneIds(c3BroadSets)
names(C3) = paste0("C3_",names(C3))
C5 = geneIds(c5BroadSets)
names(C5) = paste0("C5_",names(C5))
C6 = geneIds(c6BroadSets)
names(C6) = paste0("C6_",names(C6))
C7 = geneIds(c7BroadSets)
names(C7) = paste0("C7_",names(C7))
C8 = geneIds(c8BroadSets)
names(C8) = paste0("C8_",names(C8))

BroadSets = c(H,C1,C2,C3,C5,C6,C7,C8)

rm(H,C1,C2,C3,C5,C6,C7,C8)

saveRDS(BroadSets, file="BroadSets_v2023.1.Hs.rds")

BroadSets <- readRDS("BroadSets_v2023.1.Hs.rds")

# This will couple the Entrez IDs to the BroadSets defined earlier)
sets = ids2indices(BroadSets,unlist(v$genes$entrezgene))

file.copy(paste(dropbox,"Support/R/barcodeplot.r",sep="/"),scriptDir)
source(paste(scriptDir,"barcodeplot.r", sep="/"))

file.copy(paste(dropbox,"Support/R/cameraReports.r",sep="/"),scriptDir)
source(paste(scriptDir,"cameraReports.r", sep="/"))

cameraReport(v.=v, fit.=fit, constrasts=my.constrasts, trend=FALSE, robust=FALSE, msigdb="2023.1.Hs", species="human", max.pathways=100)

# Relocate 'reports' directory
file.rename("reports/", to=paste0(outDir,"/reports"))

# And relocate the used BroadSets
toBeMoved <- list.files(".", pattern="*.rds", full.names = TRUE)
file.copy(toBeMoved, to=outDir)
file.remove(toBeMoved)

#### Genesets - Venn Diagrams ####
file.copy(paste(dropbox,"Support/R/vennReportsGeneSets.r",sep="/"),scriptDir)
source(paste(scriptDir,"vennReportsGeneSets.r", sep="/"))

v <- readRDS(paste(outDir,"v.AllComparisons.rds", sep="/"))
targets <- v$targets
design <- v$design
my.contrasts <- v$contrasts
fit <- readRDS(paste(outDir,"fit.AllComparisons.rds", sep="/"))

vennReportsGeneSets(p.value=0.05, fdr="BH", filePrefix = paste0(outDir,"/reports/tab.camera.v2023.1.Hs_"),
                    contrasts = my.contrasts,myTitle="vennGeneSet_Cluster0",
                    indices = c(1:4), scriptDir = ".")

vennReportsGeneSets(p.value=0.05, fdr="BH", filePrefix = paste0(outDir,"/reports/tab.camera.v2023.1.Hs_"),
                    contrasts = my.contrasts,myTitle="vennGeneSet_Cluster1",
                    indices = c(5:8), scriptDir = ".")

vennReportsGeneSets(p.value=0.05, fdr="BH", filePrefix = paste0(outDir,"/reports/tab.camera.v2023.1.Hs_"),
                    contrasts = my.contrasts, myTitle="vennGeneSet_Cluster4",
                    indices = c(9:12), scriptDir = ".")

vennReportsGeneSets(p.value=0.05, fdr="BH", filePrefix = paste0(outDir,"/reports/tab.camera.v2023.1.Hs_"),
                    contrasts = my.contrasts, myTitle="vennGeneSet_Cluster8",
                    indices = c(13:16), scriptDir = ".")

vennReportsGeneSets(p.value=0.05, fdr="BH", filePrefix = paste0(outDir,"/reports/tab.camera.v2023.1.Hs_"),
                    contrasts = my.contrasts, indices = c(4,8,12,16), scriptDir = ".",
                    myTitle="vennGeneSet_ANCA_MPO_vs_PR3")

vennReportsGeneSets(p.value=0.05, fdr="BH", filePrefix = paste0(outDir,"/reports/tab.camera.v2023.1.Hs_"),
                    contrasts = my.contrasts, indices = c(17:19), scriptDir = ".",
                    myTitle="vennGeneSet_ANCA_betweenClusters_part1")

vennReportsGeneSets(p.value=0.05, fdr="BH", filePrefix = paste0(outDir,"/reports/tab.camera.v2023.1.Hs_"),
                    contrasts = my.contrasts, indices = c(20:22), scriptDir = ".",
                    myTitle="vennGeneSet_ANCA_betweenClusters_part2")

vennReportsGeneSets(p.value=0.05, fdr="BH", filePrefix = paste0(outDir,"/reports/tab.camera.v2023.1.Hs_"),
                    contrasts = my.contrasts, indices = c(23:25), scriptDir = ".",
                    myTitle="vennGeneSet_SLE_betweenClusters_part1")

vennReportsGeneSets(p.value=0.05, fdr="BH", filePrefix = paste0(outDir,"/reports/tab.camera.v2023.1.Hs_"),
                    contrasts = my.contrasts, indices = c(26:28), scriptDir = ".",
                    myTitle="vennGeneSet_SLE_betweenClusters_part2")

vennReportsGeneSets(p.value=0.05, fdr="BH", filePrefix = paste0(outDir,"/reports/tab.camera.v2023.1.Hs_"),
                    contrasts = my.contrasts, indices = c(29:31), scriptDir = ".",
                    myTitle="vennGeneSet_CONTROLE_betweenClusters_part1")

vennReportsGeneSets(p.value=0.05, fdr="BH", filePrefix = paste0(outDir,"/reports/tab.camera.v2023.1.Hs_"),
                    contrasts = my.contrasts, indices = c(32:34), scriptDir = ".",
                    myTitle="vennGeneSet_CONTROLE_betweenClusters_part2")

vennReportsGeneSets(p.value=0.05, fdr="BH", filePrefix = paste0(outDir,"/reports/tab.camera.v2023.1.Hs_"),
                    contrasts = my.contrasts, indices = c(35:37), scriptDir = ".",
                    myTitle="vennGeneSet_ANCA_MPO_betweenClusters_part1")

vennReportsGeneSets(p.value=0.05, fdr="BH", filePrefix = paste0(outDir,"/reports/tab.camera.v2023.1.Hs_"),
                    contrasts = my.contrasts, indices = c(38:40), scriptDir = ".",
                    myTitle="vennGeneSet_ANCA_MPO_betweenClusters_part2")

vennReportsGeneSets(p.value=0.05, fdr="BH", filePrefix = paste0(outDir,"/reports/tab.camera.v2023.1.Hs_"),
                    contrasts = my.contrasts, indices = c(41:43), scriptDir = ".",
                    myTitle="vennGeneSet_ANCA_PR3_betweenClusters_part1")

vennReportsGeneSets(p.value=0.05, fdr="BH", filePrefix = paste0(outDir,"/reports/tab.camera.v2023.1.Hs_"),
                    contrasts = my.contrasts, indices = c(44:46), scriptDir = ".",
                    myTitle="vennGeneSet_ANCA_PR3_betweenClusters_part2")

vennReportsGeneSets(p.value=0.05, fdr="BH", filePrefix = paste0(outDir,"/reports/tab.camera.v2023.1.Hs_"),
                    contrasts = my.contrasts, indices = c(47:49), scriptDir = ".",
                    myTitle="vennGeneSet_Clusters_part1")

vennReportsGeneSets(p.value=0.05, fdr="BH", filePrefix = paste0(outDir,"/reports/tab.camera.v2023.1.Hs_"),
                    contrasts = my.contrasts, indices = c(50:52), scriptDir = ".",
                    myTitle="vennGeneSet_Clusters_part2")

vennReportsGeneSets(p.value=0.05, fdr="BH", filePrefix = paste0(outDir,"/reports/tab.camera.v2023.1.Hs_"),
                    contrasts = my.contrasts, indices = c(53:56), scriptDir = ".",
                    myTitle="vennGeneSet_Clusters_One_vs_Rest")

vennReportsGeneSets(p.value=0.05, fdr="BH", filePrefix = paste0(outDir,"/reports/tab.camera.v2023.1.Hs_"),
                    contrasts = my.contrasts, indices = c(57:58), scriptDir = ".",
                    myTitle="vennGeneSet_Cluster_C0_or_C1_vs_C4andC8")

vennReportsGeneSets(p.value=0.05, fdr="BH", filePrefix = paste0(outDir,"/reports/tab.camera.v2023.1.Hs_"),
                    contrasts = my.contrasts, indices = c(59:60), scriptDir = ".",
                    myTitle="vennGeneSet_Cluster_C4_or_C8_vs_C0andC1")

vennReportsGeneSets(p.value=0.05, fdr="BH", filePrefix = paste0(outDir,"/reports/tab.camera.v2023.1.Hs_"),
                    contrasts = my.contrasts, indices = c(61), scriptDir = ".",
                    myTitle="vennGeneSet_Clusters_C0andC1_vs_C4andC8")

# Relocate 'venn_gs' directory
file.rename("venn_gs", to=paste0(outDir,"/venn_gs"))

# Relocate figures and stats
toBeMoved <- list.files(".", pattern="*.png|*.html", full.names = TRUE)
file.copy(toBeMoved, to=outDir)
file.remove(toBeMoved)

# And relocate the tt.genesets.rds
toBeMoved <- list.files(".", pattern="*.rds", full.names = TRUE)
file.copy(toBeMoved, to=outDir)
file.remove(toBeMoved)

#### Extra - Get function of clusters? ####

# Read the used genesets
BroadSets <- readRDS(paste0(outDir,"/BroadSets_v2023.1.Hs.rds"))

v <- readRDS(paste(outDir,"v.AllComparisons.rds", sep="/"))
targets <- v$targets
design <- v$design
my.contrasts <- v$contrasts
fit <- readRDS(paste(outDir,"fit.AllComparisons.rds", sep="/"))

extraComps <- colnames(my.contrasts)[grep("and", colnames(my.contrasts))]

for (comp in c("C1_vs_C0", "C4_vs_C0", "C8_vs_C0", "C4_vs_C1", "C8_vs_C1", "C8_vs_C4", extraComps)){
  # Read in the cluster comparison
  df <- read.table(paste0(outDir, paste0("/reports/tab.camera.v2023.1.Hs_", comp,".txt")), header=TRUE)
  
  # Select different categories or top X (use "topX")
  for (category in c("HALLMARK", "REACTOME", "KEGG", "BIOCARTA", "top50") ){
    if (grepl("top", category)){
      nrSets <- as.numeric(gsub("top", "", category))
      df.selected <- df[1:nrSets, ]
    } else {
      # Select only category, 50 sets can nicely be displayed ...
      if (category == "HALLMARK"){
        df.selected <- df[grep(category, df$geneset), ]
      } else {
        df.selected <- df[grep(category, df$geneset), ][1:50,]
      }
    }
    print(paste0(comp,": ", category," - ",nrow(df.selected)))
    df.selected$Description <- gsub("_", " ",gsub("^C[1-8]+_|^H_", "", df.selected$geneset))
    
    # Get top up- and down regulated
    df.selected.up <- df.selected[df.selected$Direction == "Up", ]
    df.selected.down <- df.selected[df.selected$Direction == "Down", ]
    
    df.selected <- rbind(df.selected.up, df.selected.down)
    
    # ggplot knows ordering if you specify the levels ...
    df.selected$Description <- factor(df.selected$Description, levels=df.selected$Description)
    df.selected$Direction <- factor(df.selected$Direction, levels=c("Up", "Down"))
    
    p <- ggplot(data = df.selected, aes(x = -log10(`PValue`), y = Description, fill = Direction)) + 
      geom_bar(stat="identity") +
      ylab("") + 
      xlab("log10(p value)") + 
      scale_y_discrete(limits = rev(levels(df.selected$Description))) +
      scale_fill_manual(values=c("#F53218", "#1895F5")) +
      scale_x_continuous(name="-log10(p value)", breaks=c(0,1,2,3,4,5,6,7,8,9) , labels=waiver()) +
      theme_bw() + 
      theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1)) +
      theme(axis.text.y = element_text(size = 7.5)) +
      theme(legend.key.height=unit(1.2, "cm")) +
      ggtitle(paste0("Gene set enrichment analysis - ", gsub("andC", " and ", gsub("_vs_C", " vs. ", gsub("^C", "Cluster ", contrast))) )) +
      theme(plot.title = element_text(hjust = 1))
    
    pdf(paste(outDir,paste0("Fig_barPlot_GenesetEnrichment_", category,"_", comp,".pdf"), sep="/"), 
        width = 7, height = 7, pointsize = 1, useDingbats=FALSE)
    print(p)
    dev.off()
    
    #ggsave(paste(outDir,paste0("Fig_barPlot_GenesetEnrichment_", category, "_", comp,".tiff"), sep="/"), p, scale=0.8, width=10, height=10, dpi=360,compression="lzw")
  }
}
# Some geneset names are horribly long, so for the final figures I have to find a solution ...

########################## ##
#### Prepare Shiny input ####
########################## ##
file.copy(paste(dropbox,"Support/R/writeShinyScripts_iHeatMap3.R",sep="/"),scriptDir)
source(paste(scriptDir,"writeShinyScripts_iHeatMap3.R",sep="/"))

writeUI_RNASeq(fileOut = paste0(outDir,"/ui.r"), tt="tt_all.AllComparisons.rds",v="v.AllComparisons.rds", 
               expinfo="ExpInfo.txt", scriptDir="applied_scripts/", venn = "venn", 
               venn_gs = "venn_gs/")
writeServer_RNASeq(fileOut = paste0(outDir,"/server.r"), tt="tt_all.AllComparisons.rds", v="v.AllComparisons.rds", 
                   fit="fit.AllComparisons.rds",
                   fastqc="../Results/SummaryTables/QC_fastqc_final", aheatmap="aheatmap_initial.png",
                   BroadSets="BroadSets_v2023.1.Hs.rds", tabCameraFile="reports/tab.camera.v2023.1.Hs_",
                   speciesSymbol="hgnc_symbol", scriptDir="applied_scripts/", venn = "venn/", venn_gs = "venn_gs/")
write_runShinyApp(fileOut = paste0(outDir,"/runShinyApp.r"))

#### Session info ####
writeVersions <- function(sessionDir=outDir){
  write(paste0("Bioconductor version ", capture.output(tools:::.BioC_version_associated_with_R_version()),"\n"), 
        paste0(sessionDir,"/sessionInfo.txt"))
  write(capture.output(sessionInfo()), paste0(sessionDir,"/sessionInfo.txt"), append=TRUE)
}

writeVersions()

print(paste0("Bioconductor version ", capture.output(tools:::.BioC_version_associated_with_R_version()),"\n"))
sessionInfo()