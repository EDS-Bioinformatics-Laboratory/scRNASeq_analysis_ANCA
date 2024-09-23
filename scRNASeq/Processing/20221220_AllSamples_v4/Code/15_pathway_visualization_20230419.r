#### NOTES ####
# - Which final results to use? 14_libra_deg_20230426 or other (see fx. initial aPear trial "Between clusters" )
#   HALLMARK only or selectedSets (HALLMARK, REACTOME, KEGG, BIOCARTA)?
# - Which cut-offs to use?
# - Default aPEAR (i.e. clusterProfiler) or aPEAR with MSigDB results?
# - ....

#### Set default for a good start ####
source("initAnalysis.r")
init()

# Create a directory to hold results and plots
outDir <- "15_pathway_visualization_20230419"
dir.create(outDir)

#### Library ####
library(devtools)
install_github('ievaKer/aPEAR')

library(aPEAR)
library(plotly)
library(htmlwidgets)

# # Example code aPEAR
# library(clusterProfiler)
# library(org.Hs.eg.db)
# library(DOSE)
# data(geneList)
# 
# enrich_base <- gseGO(geneList, OrgDb = org.Hs.eg.db, ont = 'CC')
# p <- enrichmentNetwork(enrich_base@result, drawEllipses = TRUE, fontSize = 2.5)
# 
# # And interactive using plotly
# ggplotly(p, tooltip=c('ID', 'Cluster', 'Cluster size'))

#### Between clusters ####
# aPEAR seems to need enrichr$results ... https://github.com/ievaKer/aPEAR

# Use enrichDF2enrichResult from https://rdrr.io/github/jmw86069/jamenrich/man/enrichDF2enrichResult.html
# install.packages("remotes")
# remotes::install_github("jmw86069/jamenrich")
library(multienrichjam)
library(stringr)

resultsDir <- "15_pathway_visualization_20230419/aPEAR_MSigDB_14_libra_deg_20230418"
dir.create(resultsDir)

## Read the used genesets
BroadSets <- readRDS("14_libra_deg_20230418/BroadSets_v2023.1.Hs.rds")
my.contrasts <- readRDS("14_libra_deg_20230418/contrasts_betweenClusters.rds")

# Should we only select certain combinations of clusters ??

## Read the data ##
for ( comp in colnames(my.contrasts)){
  df <- read.table(paste0("14_libra_deg_20230418/reports_betweenClusters/tab.camera.v2023.1.Hs_", comp,".txt"), header=TRUE)
  df$Description <- gsub("_", " ",gsub("C[1-8]+_|H_", "", df$geneset))
  
  # TO DO:
  # - Select only subsets?? i.e. REACTOME, HALLMARK, KEGG, BIOCARTA, WP, PID
  
  # aPEAR expects a column called 'pathwayGenes' having all genes in the pathway seperated by a "/"
  df$pathwayGenes <- apply(df, 1, function(x) paste(BroadSets[[x["geneset"]]], collapse="/"))
  
  enrich <- enrichDF2enrichResult(df, keyColname = c('geneset'), pvalueColname = c("PValue"), geneColname = c("pathwayGenes") )
  enrich@result$Count <- enrich@result$NGenes
  enrich@result$setSize <- apply(enrich@result, 1, function(x) length(BroadSets[[x["ID"]]]))
  enrich@result$GeneRatio <- paste(enrich@result$Count,enrich@result$setSize, sep="/") 
  # 'enrichmentNetwork seems to expect a 'core_enrichment' column containing all genes separated by a "/"
  # so, rename the "geneID" column
  #colnames(enrich@result)[which(names(enrich@result) == "geneID")] <- "core_enrichment"
  
  ## Observations
  # - Taking the top 1000 takes some time ... and leads to a pretty messy figure
  # - enrichmentNetwork will automatically put the label of the color legend to either "NES" or "p-value"
  #   depending whether you leave 'colorType' on default or 'pval'
  # - You have to change the label of the colour legend afterwards if it is not the 'NES' or 'pval'
  
  # Find a nice color scale also showing the direction, Up or Down
  enrich@result$log10FDR <- -log10(enrich@result$FDR)
  enrich@result$dirFDR <- ifelse(enrich@result$Direction == "Down", -1 * enrich@result$log10FDR, enrich@result$log10FDR )
  # Hmmm, why is the maximum of the color scale not the maximum of the values??
  # look at aPEAR:::enrichmentNetwork.plot ??
  # 
  myColorFactor = "dirFDR"
  p <- enrichmentNetwork(enrich@result[1:250,], colorBy = myColorFactor, nodeSize = "setSize", drawEllipses = TRUE, fontSize = 2.5,
                         colorType = 'nes', repelLabels = TRUE)
  p$labels$colour <- "-log10(FDR) * Direction\n (Up=1, Down=-1)"
  p <- p + ggtitle(comp)
  p
  
  # And save as figure
  pdf(paste(resultsDir,paste0("GenesetEnrichment_", comp, ".pdf"), sep="/"), width = 7, height = 7, pointsize = 1, useDingbats=FALSE)
  print(p)
  dev.off()
  
  ggsave(paste(resultsDir,paste0("GenesetEnrichment_", comp, ".tiff"), sep="/"), p, scale=0.8, width=10, height=10, dpi=360,compression="lzw")
  
  # And interactive using plotly
  p_plotly <- ggplotly(p, tooltip=c('ID', 'Cluster', 'Cluster size'))
  
  # And save
  htmlwidgets::saveWidget(
    widget = p_plotly, #the plotly object
    file = paste(resultsDir,paste0("GenesetEnrichment_", comp, ".html"), sep="/"), #the path & file name
    selfcontained = TRUE #creates a single html file
  )
  
}

#### ANCA's vs HC ####
# These can be found in the "14_libra_deg_20230418/reports/" directory
v <- readRDS("14_libra_deg_20230418/v.AllComparisons.rds")
my.contrasts <- v$contrasts

selectedContrasts <- colnames(my.contrasts)[grep("ANCA_vs_CONTROL", colnames(my.contrasts))]

for (comp in selectedContrasts){
  df <- read.table(paste0("14_libra_deg_20230418/reports/tab.camera.v2023.1.Hs_", comp,".txt"), header=TRUE)
  df$Description <- gsub("_", " ",gsub("C[1-8]+_|H_", "", df$geneset))
  
  # TO DO:
  # - Select only subsets?? i.e. REACTOME, HALMMARK, KEGG, BIOCARTA, WP, PID
  
  # aPEAR expects a column called 'pathwayGenes' having all genes in the pathway seperated by a "/"
  df$pathwayGenes <- apply(df, 1, function(x) paste(BroadSets[[x["geneset"]]], collapse="/"))
  
  enrich <- enrichDF2enrichResult(df, keyColname = c('geneset'), pvalueColname = c("PValue"), geneColname = c("pathwayGenes") )
  enrich@result$Count <- enrich@result$NGenes
  enrich@result$setSize <- apply(enrich@result, 1, function(x) length(BroadSets[[x["ID"]]]))
  enrich@result$GeneRatio <- paste(enrich@result$Count,enrich@result$setSize, sep="/") 
  # 'enrichmentNetwork seems to expect a 'core_enrichment' column containing all genes separated by a "/"
  # so, rename the "geneID" column
  #colnames(enrich@result)[which(names(enrich@result) == "geneID")] <- "core_enrichment"
  
  ## Observations
  # - Taking the top 1000 takes some time ... and leads to a pretty messy figure
  # - enrichmentNetwork will automatically put the label of the color legend to either "NES" or "p-value"
  #   depending whether you leave 'colorType' on default or 'pval'
  # - You have to change the label of the colour legend afterwards if it is not the 'NES' or 'pval'
  
  # Find a nice color scale also showing the direction, Up or Down
  enrich@result$log10FDR <- -log10(enrich@result$FDR)
  enrich@result$dirFDR <- ifelse(enrich@result$Direction == "Down", -1 * enrich@result$log10FDR, enrich@result$log10FDR )
  # Hmmm, why is the maximum of the color scale not the maximum of the values??
  # look at aPEAR:::enrichmentNetwork.plot ??
  # 
  myColorFactor = "dirFDR"
  p <- enrichmentNetwork(enrich@result[1:250,], colorBy = myColorFactor, nodeSize = "setSize", drawEllipses = TRUE, fontSize = 2.5,
                         colorType = 'nes', repelLabels = TRUE)
  p$labels$colour <- "-log10(FDR) * Direction\n (Up=1, Down=-1)"
  p <- p + ggtitle(comp)
  p
  
  # And save as figure
  pdf(paste(resultsDir,paste0("GenesetEnrichment_", comp, ".pdf"), sep="/"), width = 7, height = 7, pointsize = 1, useDingbats=FALSE)
  print(p)
  dev.off()
  
  ggsave(paste(resultsDir,paste0("GenesetEnrichment_", comp, ".tiff"), sep="/"), p, scale=0.8, width=10, height=10, dpi=360,compression="lzw")
  
  # And interactive using plotly
  p_plotly <- ggplotly(p, tooltip=c('ID', 'Cluster', 'Cluster size'))
  
  # And save
  htmlwidgets::saveWidget(
    widget = p_plotly, #the plotly object
    file = paste(resultsDir,paste0("GenesetEnrichment_", comp, ".html"), sep="/"), #the path & file name
    selfcontained = TRUE #creates a single html file
  )
}

#### Cytoscape: EnrichmentMap ####
# - Download Cytoscape 3.10.0 from http://www.cytoscape.org/
# - Or on Windows 7 go back in time and download an older version that
#   doesn't crash from https://cytoscape.org/download_old_versions.html.
#   For me 3.8.0 worked.
# - Open Cytoscape and install the EnrichmentMap Pipeline Collection app via Apps - App Manager
# - Apps - EnrichmentMap 
# - Common Files - GMT File: MSigDB_v2023.1_Hs.smallSet_Hs_EG_unique.gmt
# - Create Enrichment Map by clicking on the + sign
# - Analysis Type: Generic/...
# - Enrichments: EnrTable_0_vs148_Generic.txt
# ...
# - FDR Q-value Cutoff: 1
# - Build

# gmt file
BroadSets <- readRDS("14_libra_deg_GSEA_selectedSets_20230508/BroadSets_v2023.1.Hs.smallSet.rds")

cat(paste(names(BroadSets),names(BroadSets),sapply(BroadSets,function(x) paste(unique(x),collapse="\t")),sep="\t"),
     sep="\n",file=paste0(outDir,"/MSigDB_v2023.1_Hs.smallSet_Hs_EG_unique.gmt"))

camera.0vs148 <- read.delim("14_libra_deg_GSEA_selectedSets_20230508/reports_selectedSets/tab.camera.v2023.1.Hs_C0_vs_C1_C4andC8.txt")
camera.1vs048 <- read.delim("14_libra_deg_GSEA_selectedSets_20230508/reports_selectedSets/tab.camera.v2023.1.Hs_C1_vs_C0_C4andC8.txt")
camera.4vs018 <- read.delim("14_libra_deg_GSEA_selectedSets_20230508/reports_selectedSets/tab.camera.v2023.1.Hs_C4_vs_C0_C1andC8.txt")
camera.8vs014 <- read.delim("14_libra_deg_GSEA_selectedSets_20230508/reports_selectedSets/tab.camera.v2023.1.Hs_C8_vs_C0_C1andC4.txt")

# Cluster 0 vs. 1,4,8
df  <- data.frame(ID=camera.0vs148$geneset,Description=NA,Pvalue=camera.0vs148$PValue,FDR=camera.0vs148$FDR,Phenotype=ifelse(camera.0vs148$Direction=="Up",1,-1))
write.table(df,file=paste0(outDir,"/EnrTable_0_vs148_Generic.txt"),row.names=FALSE,quote=FALSE,sep="\t")

# Cluster 1 vs. 0,4,8
df  <- data.frame(ID=camera.1vs048$geneset,Description=NA,Pvalue=camera.1vs048$PValue,FDR=camera.1vs048$FDR,Phenotype=ifelse(camera.1vs048$Direction=="Up",1,-1))
write.table(df,file=paste0(outDir,"/EnrTable_1_vs048_Generic.txt"),row.names=FALSE,quote=FALSE,sep="\t")

# Cluster 4 vs. 0,1,8
df  <- data.frame(ID=camera.1vs048$geneset,Description=NA,Pvalue=camera.1vs048$PValue,FDR=camera.1vs048$FDR,Phenotype=ifelse(camera.1vs048$Direction=="Up",1,-1))
write.table(df,file=paste0(outDir,"/EnrTable_4_vs018_Generic.txt"),row.names=FALSE,quote=FALSE,sep="\t")

# Cluster 8 vs. 0,1,4
df  <- data.frame(ID=camera.8vs014$geneset,Description=NA,Pvalue=camera.8vs014$PValue,FDR=camera.8vs014$FDR,Phenotype=ifelse(camera.8vs014$Direction=="Up",1,-1))
write.table(df,file=paste0(outDir,"/EnrTable_8_vs014_Generic.txt"),row.names=FALSE,quote=FALSE,sep="\t")

#### aPEAR default using clusterProfiler ####
# Look at tutorial
library(clusterProfiler)
library(org.Hs.eg.db)

resultsDir <- "15_pathway_visualization_20230419/aPEAR_clusterProfiler"
dir.create(resultsDir)

v <- readRDS("14_libra_deg_20230426/v.AllComparisons.rds")
tt <- readRDS("14_libra_deg_20230426/tt_all.AllComparisons.rds")
my.contrasts <- v$contrasts

selectedContrasts <- c("C0_vs_C1_C4andC8", "C1_vs_C0_C4andC8", "C4_vs_C0_C1andC8", "C8_vs_C0_C1andC4" )

for (comp in selectedContrasts){
  colName <- paste0("P.Value (", comp,")")
  tt.i <- tt[, c(1:5,grep(comp, colnames(tt)))]
  #tt.i <- tt.i %>% dplyr::arrange((glue("P.Value ({comp})")))
  tt.i <- tt.i[order(tt[,colName], decreasing = TRUE),]
  
  set.seed(42)
  geneList <- tt.i[, paste0("P.Value (", comp,")")]
  # Default ENTREZGENE_ID is used ...
  names(geneList) <- tt.i[, "hgnc_symbol"]
  # For "C0_vs_C1_C4andC8":
  # enrich <- gseGO(geneList, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont="CC", pvalueCutoff = 0.05, seed = TRUE)
  # dim(enrich)
  # # [1]  1 11
  # # Only 1 hit .... probably due to the p-valueCutoff and  ....
  # enrich <- gseGO(geneList, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont="CC", pvalueCutoff = 0.1, seed = TRUE)
  # dim(enrich)
  # # [1]  1 11
  # enrich <- gseGO(geneList, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont="CC", pvalueCutoff = 0.1, pAdjustMethod = "none", seed = TRUE)
  # dim(enrich)
  # #[1] 53 11
  
  # Use the Biological Processes from GO, i.e. ont="BP"
  enrich <- gseGO(geneList, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont="BP", pvalueCutoff = 0.05, pAdjustMethod = "none", seed = TRUE,
                  minGSSize = 25)
  cat(paste0(comp, ":\t", dim(enrich)[1], " Biological Processes pass the cut-off\n\n"))
  
  p <- enrichmentNetwork(enrich@result, nodeSize = "setSize", drawEllipses = TRUE, fontSize = 2.5,
                         colorType = 'nes', repelLabels = TRUE)
  #p$labels$colour <- "-log10(FDR) * Direction\n (Up=1, Down=-1)"
  p <- p + ggtitle(comp)
  p
  
  # And save as figure
  pdf(paste(resultsDir,paste0("GenesetEnrichment_", comp, ".pdf"), sep="/"), width = 7, height = 7, pointsize = 1, useDingbats=FALSE)
    print(p)
  dev.off()
  
  ggsave(paste(resultsDir,paste0("GenesetEnrichment_", comp, ".tiff"), sep="/"), p, scale=0.8, width=10, height=10, dpi=360,compression="lzw")
  
  # And interactive using plotly
  p_plotly <- ggplotly(p, tooltip=c('ID', 'Cluster', 'Cluster size'))
  
  # And save
  htmlwidgets::saveWidget(
    widget = p_plotly, #the plotly object
    file = paste(resultsDir,paste0("GenesetEnrichment_", comp, ".html"), sep="/"), #the path & file name
    selfcontained = TRUE #creates a single html file
  )
}

# If you interrupt the function, you might get the following rror message:
# Error in serialize(data, node$con) : error writing to connection
# In addition: Warning messages:
#   1: In fgseaMultilevel(...) :
#   There were 38 pathways for which P-values were not calculated properly due to unbalanced (positive and negative) gene-level statistic values. For such pathways pval, padj, NES, log2err are set to NA. You can try to increase the value of the argument nPermSimple (for example set it nPermSimple = 10000)
# 2: In serialize(data, node$con) :
#   'package:stats' may not be available when loading


# Biological Processes (BP)
# p value = 0.05/0.1, BH adjusted
# - C0_vs_C1_C4andC8:	0 Biological Processes pass the cut-off
# - and stops ....

# p value = 0.01, not adjusted
# - C0_vs_C1_C4andC8:	28 Biological Processes pass the cut-off
# - C1_vs_C0_C4andC8:	35 Biological Processes pass the cut-off
# - C4_vs_C0_C1andC8:	21 Biological Processes pass the cut-off
# - C8_vs_C0_C1andC4:	20 Biological Processes pass the cut-off

# p value = 0.05, not adjusted
# - C0_vs_C1_C4andC8:	111 Biological Processes pass the cut-off
# - C1_vs_C0_C4andC8:	158  Biological Processes pass the cut-off
# - C4_vs_C0_C1andC8:	102 Biological Processes pass the cut-off
# - C8_vs_C0_C1andC4:	144 Biological Processes pass the cut-off

#### aPEAR with MSigDB - 14_libra_deg_GSEA_selectedSets_20230508 ####
# To comply with the Cytoscape Enrichment analysis performed by Perry
library(multienrichjam)
library(stringr)
library(clusterProfiler)

# Without "NoLabels" figures ....
#resultsDir <- "15_pathway_visualization_20230419/aPEAR_MSigDB_14_libra_deg_GSEA_selectedSets_20230508"
resultsDir <- "15_pathway_visualization_20230419/aPEAR_MSigDB_14_libra_deg_GSEA_selectedSets_20230508_v2"
dir.create(resultsDir)

## Read the used genesets
BroadSets <- readRDS("14_libra_deg_GSEA_selectedSets_20230508/BroadSets_v2023.1.Hs.smallSet.rds")

selectedContrasts <- c("C0_vs_C1_C4andC8", "C1_vs_C0_C4andC8", "C4_vs_C0_C1andC8", "C8_vs_C0_C1andC4" )

## Read the data ##
for ( comp in selectedContrasts){
  df <- read.table(paste0("14_libra_deg_GSEA_selectedSets_20230508/reports_selectedSets/tab.camera.v2023.1.Hs_", comp,".txt"), header=TRUE)
  df$Description <- gsub("_", " ",gsub("C[1-8]+_|H_", "", df$geneset))
  
  # aPEAR expects a column called 'pathwayGenes' having all genes in the pathway seperated by a "/"
  df$pathwayGenes <- apply(df, 1, function(x) paste(BroadSets[[x["geneset"]]], collapse="/"))
  
  enrich <- enrichDF2enrichResult(df, keyColname = c('geneset'), pvalueColname = c("PValue"), geneColname = c("pathwayGenes") )
  enrich@result$Count <- enrich@result$NGenes
  enrich@result$setSize <- apply(enrich@result, 1, function(x) length(BroadSets[[x["ID"]]]))
  enrich@result$GeneRatio <- paste(enrich@result$Count,enrich@result$setSize, sep="/") 
  # 'enrichmentNetwork seems to expect a 'core_enrichment' column containing all genes separated by a "/"
  # so, rename the "geneID" column
  #colnames(enrich@result)[which(names(enrich@result) == "geneID")] <- "core_enrichment"
  
  ## Observations
  # - Taking the top 1000 takes some time ... and leads to a pretty messy figure
  # - enrichmentNetwork will automatically put the label of the color legend to either "NES" or "p-value"
  #   depending whether you leave 'colorType' on default or 'pval'
  # - You have to change the label of the colour legend afterwards if it is not the 'NES' or 'pval'
  
  # Find a nice color scale also showing the direction, Up or Down
  enrich@result$log10FDR <- -log10(enrich@result$FDR)
  enrich@result$dirFDR <- ifelse(enrich@result$Direction == "Down", -1 * enrich@result$log10FDR, enrich@result$log10FDR )
  # Hmmm, why is the maximum of the color scale not the maximum of the values??
  # look at aPEAR:::enrichmentNetwork.plot ??
  # 
  myColorFactor = "dirFDR"
  p <- enrichmentNetwork(enrich@result[1:250,], colorBy = myColorFactor, nodeSize = "setSize", drawEllipses = TRUE, fontSize = 2.5,
                         colorType = 'nes', repelLabels = TRUE)
  p$labels$colour <- "-log10(FDR) * Direction\n (Up=1, Down=-1)"
  p <- p + ggtitle(comp)
  #p
  
  # And save as figure
  pdf(paste(resultsDir,paste0("GenesetEnrichment_", comp, ".pdf"), sep="/"), width = 7, height = 7, pointsize = 1, useDingbats=FALSE)
  print(p)
  dev.off()
  
  ggsave(paste(resultsDir,paste0("GenesetEnrichment_", comp, ".tiff"), sep="/"), p, scale=0.8, width=10, height=10, dpi=360,compression="lzw")
  
  # And interactive using plotly
  p_plotly <- ggplotly(p, tooltip=c('ID', 'Cluster', 'Cluster size'))
  
  # And save
  htmlwidgets::saveWidget(
    widget = p_plotly, #the plotly object
    file = paste(resultsDir,paste0("GenesetEnrichment_", comp, ".html"), sep="/"), #the path & file name
    selfcontained = TRUE #creates a single html file
  )
  
  # And without labels - only in "aPEAR_MSigDB_14_libra_deg_GSEA_selectedSets_20230508_v2" 
  # But then it plots how many sets are in there??
  # Remove the 'labels', they are in the last layer (I checked for one figure...so, I presumed this
  # to hld true for all figures..)
  p$layers[[length(p$layers)]] <- NULL
  
  # And save as figure
  pdf(paste(resultsDir,paste0("GenesetEnrichment_", comp, "_NoLabels.pdf"), sep="/"), width = 7, height = 7, pointsize = 1, useDingbats=FALSE)
    print(p)
  dev.off()
  
  ggsave(paste(resultsDir,paste0("GenesetEnrichment_", comp, "_NoLabels.tiff"), sep="/"), p, scale=0.8, width=10, height=10, dpi=360,compression="lzw")
  
}

#### aPEAR with MSigDB - 14_libra_deg_GSEA_selectedSets_20230508 - adapted color scheme ####
# To comply with the Cytoscape Enrichment analysis performed by Perry
library(multienrichjam)
library(stringr)
library(clusterProfiler)
library(dplyr)    # needed for "between()"
library(ggforce)  # needed for "geom_link0()"

# "NoLabels" figures ....
resultsDir <- "15_pathway_visualization_20230419/aPEAR_MSigDB_14_libra_deg_GSEA_selectedSets_20230508_v3"
dir.create(resultsDir)

## Read the used genesets
BroadSets <- readRDS("14_libra_deg_GSEA_selectedSets_20230508/BroadSets_v2023.1.Hs.smallSet.rds")

selectedContrasts <- c("C0_vs_C1_C4andC8", "C1_vs_C0_C4andC8", "C4_vs_C0_C1andC8", "C8_vs_C0_C1andC4" )

# Own function to have influence on the coloring scheme used by aPEAR
# - Original code uses a symmetric scheme
# - See remark on aPEAR:::enrichmentNetwork.plot
# - I will have to change enrichmentNewtork to use my own version of aPEAR:::enrichmentNetwork.plot
my_enrichmentNetworkPlot <- function (dt, sim, clust, innerCutoff = 0.1, outerCutoff = 0.5, 
          colorType = c("nes", "pval"), pCutoff = -10, 
          drawEllipses = FALSE, fontSize = 5, repelLabels = FALSE) 
{
  colorType <- match.arg(colorType)
  graph <- aPEAR:::enrichmentNetwork.connect(sim, clust, innerCutoff = innerCutoff, 
                                     outerCutoff = outerCutoff)
  coordinates <- merge(graph$coordinates, dt, by.x = "ID", 
                       by.y = "ID")
  lines <- graph$edges
  lines <- lines[from %in% coordinates[, ID] & to %in% coordinates[, 
                                                                   ID]]
  if (colorType == "nes") {
    # AJ - 20230703: changed original code:
    #     range <- max(abs(coordinates[, color]))
    #     colors <- scale_color_distiller(limits = c(-range, range), palette = "Spectral")
    rangeMax <- max(coordinates[, color])
    rangeMin <- min(coordinates[, color])
    colors <- scale_color_distiller(limits = c(rangeMin, rangeMax), 
                                    palette = "Spectral")
    colorTitle <- "NES"
  }
  if (colorType == "pval") {
    coordinates[, `:=`(color, log(color))] %>% .[color < 
                                                   pCutoff, `:=`(color, pCutoff)]
    colors <- scale_color_distiller(limits = c(pCutoff, 0), 
                                    direction = -1, palette = "OrRd")
    colorTitle <- "p-value"
  }
  plot <- ggplot()
  if (drawEllipses) {
    plot <- aPEAR:::enrichmentNetwork.addEllipses(plot, coordinates)
  }
  plot <- plot + geom_link0(data = lines, aes(x = xStart, y = yStart, 
                                              xend = xEnd, yend = yEnd), size = 0.25, alpha = 0.2) + 
    geom_point(data = coordinates, aes(x = x, y = y, ID = ID, 
                                       color = color, size = size, Cluster = Cluster, `Cluster size` = `Cluster size`)) + 
    theme(axis.title.x = element_blank(), axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), axis.title.y = element_blank(), 
          axis.text.y = element_blank(), axis.ticks.y = element_blank(), 
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_blank(), legend.position = "right", 
          legend.key = element_rect(fill = "white")) + 
    labs(color = colorTitle, size = "Pathway size") + 
    coord_fixed() + aPEAR:::enrichmentNetwork.clusterLabels(coordinates, 
                                                    fontSize, repelLabels) + colors
  plot
}

# And adapt the main code
my_enrichmentNetwork <- function (enrichment, simMethod = c("jaccard", "cosine","cor"), 
                                  clustMethod = c("markov", "hier",  "spectral"), 
                                  clustNameMethod = c("pagerank", "hits", "none"), 
                                  colorBy = NULL, nodeSize = NULL, innerCutoff = 0.1, outerCutoff = 0.5, 
                                  colorType = c("nes", "pval"), pCutoff = -10, drawEllipses = FALSE, fontSize = 5, 
                                  repelLabels = FALSE, minClusterSize = 2, plotOnly = TRUE, 
                                  verbose = FALSE) 
{
  if (class(enrichment) != "data.frame") {
    stop("An object of class data.frame is expected.")
  }
  simMethod <- match.arg(simMethod)
  clustMethod <- match.arg(clustMethod)
  clustNameMethod <- match.arg(clustNameMethod)
  colorType <- match.arg(colorType)
  if (!dplyr::between(innerCutoff, 0, 1)) {
    stop("innerCutoff must be between 0 and 1.")
  }
  if (!between(outerCutoff, 0, 1)) {
    stop("outerCutoff must be between 0 and 1.")
  }
  if (minClusterSize < 2) {
    stop("Currently supported minClusterSize >= 2")
  }
  params <- aPEAR:::validateEnrichment(enrichment, colorBy = colorBy, 
                               nodeSize = nodeSize, verbose = verbose)
  sim <- aPEAR:::pathwaySimilarity(enrichment, geneCol = params$genesCol, method = simMethod)
  clusters <- aPEAR:::findClusters(sim, method = clustMethod, nameMethod = clustNameMethod, 
                           minClusterSize = minClusterSize, verbose = verbose)
  enrichClust <- aPEAR:::enrichmentNetwork.prepareEnrichmentClusters(enrichment, clusters, params)
  
  # Standard code in aPEAR:::enrichmentNetwork.plot uses an symmetric color legend
  plot <- my_enrichmentNetworkPlot(enrichClust, sim, clusters, 
                                 innerCutoff = innerCutoff, outerCutoff = outerCutoff, 
                                 colorType = colorType, pCutoff = pCutoff, drawEllipses = drawEllipses, 
                                 fontSize = fontSize, repelLabels = repelLabels)
  if (plotOnly) {
    return(plot)
  }
  else {
    return(list(plot = plot, clusters = enrichClust[, list(ID, Cluster)]))
  }
}


## Read the data ##
for ( comp in selectedContrasts){
  df <- read.table(paste0("14_libra_deg_GSEA_selectedSets_20230508/reports_selectedSets/tab.camera.v2023.1.Hs_", comp,".txt"), header=TRUE)
  df$Description <- gsub("_", " ",gsub("C[1-8]+_|H_", "", df$geneset))
  
  # aPEAR expects a column called 'pathwayGenes' having all genes in the pathway seperated by a "/"
  df$pathwayGenes <- apply(df, 1, function(x) paste(BroadSets[[x["geneset"]]], collapse="/"))
  
  enrich <- enrichDF2enrichResult(df, keyColname = c('geneset'), pvalueColname = c("PValue"), geneColname = c("pathwayGenes") )
  enrich@result$Count <- enrich@result$NGenes
  enrich@result$setSize <- apply(enrich@result, 1, function(x) length(BroadSets[[x["ID"]]]))
  enrich@result$GeneRatio <- paste(enrich@result$Count,enrich@result$setSize, sep="/") 
  # 'enrichmentNetwork seems to expect a 'core_enrichment' column containing all genes separated by a "/"
  # so, rename the "geneID" column
  #colnames(enrich@result)[which(names(enrich@result) == "geneID")] <- "core_enrichment"
  
  ## Observations
  # - Taking the top 1000 takes some time ... and leads to a pretty messy figure
  # - enrichmentNetwork will automatically put the label of the color legend to either "NES" or "p-value"
  #   depending whether you leave 'colorType' on default or 'pval'
  # - You have to change the label of the colour legend afterwards if it is not the 'NES' or 'pval'
  
  # Find a nice color scale also showing the direction, Up or Down
  enrich@result$log10FDR <- -log10(enrich@result$FDR)
  enrich@result$dirFDR <- ifelse(enrich@result$Direction == "Down", -1 * enrich@result$log10FDR, enrich@result$log10FDR )
  # 
  # And plto using the custom method that takes the min and max of the range into account instead of the original symmetric range...
  myColorFactor = "dirFDR"
  p <- my_enrichmentNetwork(enrich@result[1:250,], colorBy = myColorFactor, nodeSize = "setSize", drawEllipses = TRUE, fontSize = 2.5,
                         colorType = 'nes', repelLabels = TRUE)
  p$labels$colour <- "-log10(FDR) * Direction\n (Up=1, Down=-1)"
  p <- p + ggtitle(comp)
  #p
  
  # And save as figure
  pdf(paste(resultsDir,paste0("GenesetEnrichment_", comp, ".pdf"), sep="/"), width = 7, height = 7, pointsize = 1, useDingbats=FALSE)
  print(p)
  dev.off()
  
  ggsave(paste(resultsDir,paste0("GenesetEnrichment_", comp, ".tiff"), sep="/"), p, scale=0.8, width=10, height=10, dpi=360,compression="lzw")
  
  # And interactive using plotly
  p_plotly <- ggplotly(p, tooltip=c('ID', 'Cluster', 'Cluster size'))
  
  # And save
  htmlwidgets::saveWidget(
    widget = p_plotly, #the plotly object
    file = paste(resultsDir,paste0("GenesetEnrichment_", comp, ".html"), sep="/"), #the path & file name
    selfcontained = TRUE #creates a single html file
  )
  
  # And without labels - only in "aPEAR_MSigDB_14_libra_deg_GSEA_selectedSets_20230508_v2" 
  # But then it plots how many sets are in there??
  # Remove the 'labels', they are in the last layer (I checked for one figure...so, I presumed this
  # to hld true for all figures..)
  p$layers[[length(p$layers)]] <- NULL
  
  # And save as figure
  pdf(paste(resultsDir,paste0("GenesetEnrichment_", comp, "_NoLabels.pdf"), sep="/"), width = 7, height = 7, pointsize = 1, useDingbats=FALSE)
  print(p)
  dev.off()
  
  ggsave(paste(resultsDir,paste0("GenesetEnrichment_", comp, "_NoLabels.tiff"), sep="/"), p, scale=0.8, width=10, height=10, dpi=360,compression="lzw")
  
}

#### aPEAR with MSigDB - 14_libra_deg_GSEA_selectedSets_20230508 - adapted color scheme + selected sets ####
# To comply with the Cytoscape Enrichment analysis performed by Perry
library(multienrichjam)
library(stringr)
library(clusterProfiler)
library(dplyr)    # needed for "between()"
library(ggforce)  # needed for "geom_link0()"
library(colorspace)
library(plotly)

# "NoLabels" figures ....
resultsDir <- "15_pathway_visualization_20230419/aPEAR_v4_20231027"
dir.create(resultsDir)

## Read the used genesets
BroadSets <- readRDS("14_libra_deg_GSEA_selectedSets_20230508/BroadSets_v2023.1.Hs.smallSet.rds")

selectedContrasts <- c("C0_vs_C1_C4andC8", "C1_vs_C0_C4andC8", "C4_vs_C0_C1andC8", "C8_vs_C0_C1andC4" )

## IMPORTANT:  Read in my own functions (see line 410 above)

# But in order to fix the representation ... adapt the plotting function etc
my_enrichmentNetworkPlot <- function (dt, sim, clust, innerCutoff = 0.1, outerCutoff = 0.5, 
                                      colorType = c("nes", "pval"), pCutoff = -10, 
                                      drawEllipses = FALSE, fontSize = 5, repelLabels = FALSE,
                                      graph = NULL) 
{
  colorType <- match.arg(colorType)
  
  if (is.null(graph)){
    print("Creating a graph ...")
    graph <- aPEAR:::enrichmentNetwork.connect(sim, clust, innerCutoff = innerCutoff, 
                                              outerCutoff = outerCutoff)
  }
  # And else just use a predefined graph ....
  coordinates <- merge(graph$coordinates, dt, by.x = "ID", 
                       by.y = "ID")
  lines <- graph$edges
  lines <- lines[from %in% coordinates[, ID] & to %in% coordinates[, 
                                                                   ID]]
  if (colorType == "nes") {
    # AJ - 20230703: changed original code:
    #     range <- max(abs(coordinates[, color]))
    #     colors <- scale_color_distiller(limits = c(-range, range), palette = "Spectral")
    rangeMax <- max(coordinates[, color])
    rangeMin <- min(coordinates[, color])
    # colors <- scale_color_distiller(limits = c(rangeMin, rangeMax), 
    #                                 palette = "Spectral", mid = 0)
    colors <- scale_color_continuous_divergingx(limits = c(rangeMin, rangeMax),
                                               palette = "RdYlBu", mid = 0, rev = TRUE)
    colorTitle <- "NES"
  }
  if (colorType == "pval") {
    coordinates[, `:=`(color, log(color))] %>% .[color < 
                                                   pCutoff, `:=`(color, pCutoff)]
    #colors <- scale_color_distiller(limits = c(pCutoff, 0), 
    #                                direction = -1, palette = "OrRd", mid = 0)
    colors <- scale_color_continuous_divergingx(limits = c(pCutoff, 0),
                                               rev = TRUE, palette = "OrRd", mid = )
    
    colorTitle <- "p-value"
  }
  plot <- ggplot()
  if (drawEllipses) {
    plot <- aPEAR:::enrichmentNetwork.addEllipses(plot, coordinates)
  }
  plot <- plot + geom_link0(data = lines, aes(x = xStart, y = yStart, 
                                              xend = xEnd, yend = yEnd), size = 0.25, alpha = 0.2) + 
    geom_point(data = coordinates, aes(x = x, y = y, ID = ID, 
                                       color = color, size = size, Cluster = Cluster, `Cluster size` = `Cluster size`)) + 
    theme(axis.title.x = element_blank(), axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), axis.title.y = element_blank(), 
          axis.text.y = element_blank(), axis.ticks.y = element_blank(), 
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_blank(), legend.position = "right", 
          legend.key = element_rect(fill = "white")) + 
    labs(color = colorTitle, size = "Pathway size") + 
    coord_fixed() + aPEAR:::enrichmentNetwork.clusterLabels(coordinates, 
                                                            fontSize, repelLabels) + colors
  return(list(plot = plot, graph = graph))
}

# And adapt the main code
my_enrichmentNetwork <- function (enrichment, simMethod = c("jaccard", "cosine","cor"), 
                                  clustMethod = c("markov", "hier",  "spectral"), 
                                  clustNameMethod = c("pagerank", "hits", "none"), 
                                  colorBy = NULL, nodeSize = NULL, innerCutoff = 0.1, outerCutoff = 0.5, 
                                  colorType = c("nes", "pval"), pCutoff = -10, drawEllipses = FALSE, fontSize = 5, 
                                  repelLabels = FALSE, minClusterSize = 2, plotOnly = TRUE, 
                                  verbose = FALSE, myGraph = NULL) 
{
  if (class(enrichment) != "data.frame") {
    stop("An object of class data.frame is expected.")
  }
  simMethod <- match.arg(simMethod)
  clustMethod <- match.arg(clustMethod)
  clustNameMethod <- match.arg(clustNameMethod)
  colorType <- match.arg(colorType)
  if (!dplyr::between(innerCutoff, 0, 1)) {
    stop("innerCutoff must be between 0 and 1.")
  }
  if (!between(outerCutoff, 0, 1)) {
    stop("outerCutoff must be between 0 and 1.")
  }
  if (minClusterSize < 2) {
    stop("Currently supported minClusterSize >= 2")
  }
  params <- aPEAR:::validateEnrichment(enrichment, colorBy = colorBy, 
                                       nodeSize = nodeSize, verbose = verbose)
  sim <- aPEAR:::pathwaySimilarity(enrichment, geneCol = params$genesCol, method = simMethod)
  clusters <- aPEAR:::findClusters(sim, method = clustMethod, nameMethod = clustNameMethod, 
                                   minClusterSize = minClusterSize, verbose = verbose)
  enrichClust <- aPEAR:::enrichmentNetwork.prepareEnrichmentClusters(enrichment, clusters, params)
  
  # Standard code in aPEAR:::enrichmentNetwork.plot uses an symmetric color legend
  result <- my_enrichmentNetworkPlot(enrichClust, sim, clusters, 
                                   innerCutoff = innerCutoff, outerCutoff = outerCutoff, 
                                   colorType = colorType, pCutoff = pCutoff, drawEllipses = drawEllipses, 
                                   fontSize = fontSize, repelLabels = repelLabels, graph = myGraph)
  if (plotOnly) {
    return(result$plot)
  }
  else {
    return(list(plot = result$plot, clusters = enrichClust[, list(ID, Cluster)], graph = result$graph))
  }
}


# Now, first determine the genesets that are significant per contrast (with a certain cut-off)
# And then make a common list of these that will be used when plotting every contrast

commonGeneSetList <- list()
fdrVal = 0.1

for ( comp in selectedContrasts){
  df <- read.table(paste0("14_libra_deg_GSEA_selectedSets_20230508/reports_selectedSets/tab.camera.v2023.1.Hs_", comp,".txt"), header=TRUE)
  df$Description <- gsub("_", " ",gsub("C[1-8]+_|H_", "", df$geneset))
  
  # aPEAR expects a column called 'pathwayGenes' having all genes in the pathway seperated by a "/"
  df$pathwayGenes <- apply(df, 1, function(x) paste(BroadSets[[x["geneset"]]], collapse="/"))
  
  enrich <- enrichDF2enrichResult(df, keyColname = c('geneset'), pvalueColname = c("PValue"), geneColname = c("pathwayGenes") )
  enrich@result$Count <- enrich@result$NGenes
  enrich@result$setSize <- apply(enrich@result, 1, function(x) length(BroadSets[[x["ID"]]]))
  enrich@result$GeneRatio <- paste(enrich@result$Count,enrich@result$setSize, sep="/") 
  
  selectedSets <- rownames(enrich@result[enrich@result$FDR < fdrVal,  ])
  
  cat(paste0(comp, ": #genesets = ", length(selectedSets)),"\n")
  
  commonGeneSetList <- c(commonGeneSetList, selectedSets)
} 
# C0_vs_C1_C4andC8: #genesets = 42 
# C1_vs_C0_C4andC8: #genesets = 177 
# C4_vs_C0_C1andC8: #genesets = 240 
# C8_vs_C0_C1andC4: #genesets = 45 
  
commonGeneSetList <- unique(unlist(commonGeneSetList))
print(length(commonGeneSetList))
# [1] 391

## Read the data ##
fdrVal <- gsub("\\.","_", fdrVal)  # So, fdrVal can be nicely used as text
i = 0

for ( comp in selectedContrasts){
  df <- read.table(paste0("14_libra_deg_GSEA_selectedSets_20230508/reports_selectedSets/tab.camera.v2023.1.Hs_", comp,".txt"), header=TRUE)
  df$Description <- gsub("_", " ",gsub("C[1-8]+_|H_", "", df$geneset))
  
  # aPEAR expects a column called 'pathwayGenes' having all genes in the pathway seperated by a "/"
  df$pathwayGenes <- apply(df, 1, function(x) paste(BroadSets[[x["geneset"]]], collapse="/"))
  
  enrich <- enrichDF2enrichResult(df, keyColname = c('geneset'), pvalueColname = c("PValue"), geneColname = c("pathwayGenes") )
  enrich@result$Count <- enrich@result$NGenes
  enrich@result$setSize <- apply(enrich@result, 1, function(x) length(BroadSets[[x["ID"]]]))
  enrich@result$GeneRatio <- paste(enrich@result$Count,enrich@result$setSize, sep="/") 
  # 'enrichmentNetwork seems to expect a 'core_enrichment' column containing all genes separated by a "/"
  # so, rename the "geneID" column
  #colnames(enrich@result)[which(names(enrich@result) == "geneID")] <- "core_enrichment"
  
  ## Observations
  # - Taking the top 1000 takes some time ... and leads to a pretty messy figure
  # - enrichmentNetwork will automatically put the label of the color legend to either "NES" or "p-value"
  #   depending whether you leave 'colorType' on default or 'pval'
  # - You have to change the label of the colour legend afterwards if it is not the 'NES' or 'pval'
  
  # Find a nice color scale also showing the direction, Up or Down
  enrich@result$log10FDR <- -log10(enrich@result$FDR)
  enrich@result$dirFDR <- ifelse(enrich@result$Direction == "Down", -1 * enrich@result$log10FDR, enrich@result$log10FDR )
  # 
  # And plot using the custom method that takes the min and max of the range into account instead of the original symmetric range...
  myColorFactor = "dirFDR"
  if (i == 0){
    res <- my_enrichmentNetwork(enrich@result[commonGeneSetList,], colorBy = myColorFactor, nodeSize = "setSize", drawEllipses = TRUE, fontSize = 2.5,
                              colorType = 'nes', repelLabels = TRUE, plotOnly = FALSE)
    i = 1
  } else {
    print('Using a predefined graph ...')
    res <- my_enrichmentNetwork(enrich@result[commonGeneSetList,], colorBy = myColorFactor, nodeSize = "setSize", drawEllipses = TRUE, fontSize = 2.5,
                                colorType = 'nes', repelLabels = TRUE, plotOnly = FALSE, myGraph = g)
    
  }
  p <- res$plot
  g <- res$graph
  p$labels$colour <- "-log10(FDR) * Direction\n (Up=1, Down=-1)"
  p <- p + ggtitle(comp)
  #p
  
  # And save as figure
  pdf(paste(resultsDir,paste0("GenesetEnrichment_", comp, "_FDR_",fdrVal,  ".pdf"), sep="/"), width = 7, height = 7, pointsize = 1, useDingbats=FALSE)
    print(p)
  dev.off()
  
  ggsave(paste(resultsDir,paste0("GenesetEnrichment_", comp,"_FDR_",fdrVal,   ".tiff"), sep="/"), p, scale=0.8, width=10, height=10, dpi=360,compression="lzw")
  
  # And interactive using plotly
  p_plotly <- ggplotly(p, tooltip=c('ID', 'Cluster', 'Cluster size'))
  
  # And save
  htmlwidgets::saveWidget(
    widget = p_plotly, #the plotly object
    file = paste(resultsDir,paste0("GenesetEnrichment_", comp,"_FDR_",fdrVal, ".html"), sep="/"), #the path & file name
    selfcontained = TRUE #creates a single html file
  )
  
  # And without labels - only in "aPEAR_MSigDB_14_libra_deg_GSEA_selectedSets_20230508_v2" 
  # But then it plots how many sets are in there??
  # Remove the 'labels', they are in the last layer (I checked for one figure...so, I presumed this
  # to hld true for all figures..)
  p$layers[[length(p$layers)]] <- NULL
  
  # And save as figure
  pdf(paste(resultsDir,paste0("GenesetEnrichment_", comp, "_FDR_",fdrVal,"_NoLabels.pdf"), sep="/"), width = 7, height = 7, pointsize = 1, useDingbats=FALSE)
  print(p)
  dev.off()
  
  ggsave(paste(resultsDir,paste0("GenesetEnrichment_", comp, "_FDR_",fdrVal,"_NoLabels.tiff"), sep="/"), p, scale=0.8, width=10, height=10, dpi=360,compression="lzw")
  
}


#### Alternative visualizations ####
# enrichplot
# hypeR
# netGSA
# EnrichmentBrowser

# install.packages("remotes")
# remotes::install_github("jmw86069/jamenrich")
library(enrichplot)
library(ggplot2)
library(DOSE)
library(limma)
library(org.Hs.eg.db)
library(multienrichjam)
library(stringr)
library(clusterProfiler)
library(cowplot)
library(ggupset)
library(GOSemSim)
library(ggridges)
library(ggnewscale)
library(qpdf)
library(dplyr)

resultsDir <- "15_pathway_visualization_20230419/enrichplot"
dir.create(resultsDir)

# Read in the data
v <- readRDS("14_libra_deg_20230426/v.AllComparisons.rds")
tt <- readRDS("14_libra_deg_20230426/tt_all.AllComparisons.rds")
my.contrasts <- v$contrasts

## Read the used genesets
BroadSets <- readRDS("14_libra_deg_GSEA_selectedSets_20230508/BroadSets_v2023.1.Hs.smallSet.rds")

selectedContrasts <- c("C0_vs_C1_C4andC8", "C1_vs_C0_C4andC8", "C4_vs_C0_C1andC8", "C8_vs_C0_C1andC4" )

# Enrichment Map
d <- godata('org.Hs.eg.db', ont="BP")

# Generate a PDF for every plot and then combine into one
# using qpdf

# Create temporary directory (locally)
tmpDir <- tempdir()
dir.create(tmpDir)

for (comp in selectedContrasts){
  colName <- paste0("P.Value (", comp,")")
  tt.i <- tt[, c(1:5,grep(comp, colnames(tt)))]
  tt.i <- tt.i[order(tt[,colName], decreasing = TRUE),]
  
  set.seed(42)
  geneList <- tt.i[, paste0("P.Value (", comp,")")]
  # enrichDGN uses ENTREZGENE_IDs ...
  names(geneList) <- tt.i[, "entrezgene"]
  de <- names(geneList)[geneList < 0.01]

  # Over Representation Analysis
  # - enrichDGN has a lot of other parameters for which it is not clear what the default
  #   means:
  #   - universe, i.e. the background genes (but what should I provide)
  #   - pvalueCutoff (of what??)
  #   - readable - do I want that??
  edo <- enrichDGN(de)
  
  # GSEA
  # - Need to 'up' the pvalueCutoff to 0.1 in order to find something ... 
  edo2 <- gseNCG(geneList, pvalueCutoff = 0.1, pAdjustMethod = "none", seed = TRUE)
  
  # But taking our own results
  df <- read.table(paste0("14_libra_deg_GSEA_selectedSets_20230508/reports_selectedSets/tab.camera.v2023.1.Hs_", comp,".txt"), header=TRUE)
  df$Description <- gsub("_", " ",gsub("C[1-8]+_|H_", "", df$geneset))
  
  # aPEAR expects a column called 'pathwayGenes' having all genes in the pathway separated by a "/"
  df$pathwayGenes <- apply(df, 1, function(x) paste(BroadSets[[x["geneset"]]], collapse="/"))
  
  enrich <- enrichDF2enrichResult(df, keyColname = c('geneset'), pvalueColname = c("PValue"), geneColname = c("pathwayGenes") )
  enrich@result$Count <- enrich@result$NGenes
  enrich@result$setSize <- apply(enrich@result, 1, function(x) length(BroadSets[[x["ID"]]]))
  enrich@result$GeneRatio <- paste(enrich@result$Count,enrich@result$setSize, sep="/") 
  
  # Barplot
  pdf(paste(tmpDir,"barplot.pdf", sep="/"), width=14, height = 10)
    barplot(enrich, showCategory=20, cex.names=0.6)
  dev.off()
  
  # Dotplot
  p1 <- dotplot(enrich, showCategory=30) + ggtitle("dotplot for CAMERA-to-enrichResult")
  p2 <- dotplot(edo, showCategory=30) + ggtitle("dotplot for enrichDGN (ORA)")
  pdf(paste(tmpDir,"dotplot.pdf", sep="/"), width=20, height = 14)
    print(plot_grid(p1, p2, ncol=2))
  dev.off()
  
  # Gene-Concept Network
  pdf(paste(tmpDir,"cnetplot.pdf", sep="/"), width=10, height = 10)
    cnetplot(enrich, showCategory=20, node_label="category", 
           cex.params=list(foldChange=NULL, category_node=1, category_label=0.5, gene_node=0.1, gene_label=0.0), 
           color.params=list(foldChange=NULL, edge=TRUE, category="blue", gene="gray"),
           max.overlaps=Inf)
  dev.off()
  
  # On Windows get this error  -> due to difference in R and package version!!
  # Error in UseMethod("rescale") : 
  #     no applicable method for 'rescale' applied to an object of class "AsIs"
  
  pdf(paste(tmpDir,"cnetplot_circular.pdf", sep="/"), width=10, height = 10)
    cnetplot(enrich, showCategory=20, circular = TRUE, colorEdge = TRUE)
  dev.off()
  
  # UpSet plot
  pdf(paste(tmpDir,"upsetplot.pdf", sep="/"), width=15, height = 10)
    upsetplot(enrich, n=25)
  dev.off()
  
  # Heatplot
  pdf(paste(tmpDir,"heatplot.pdf", sep="/"), width=14, height = 14)
    p <- heatplot(enrich, showCategory = 30)
    p <- p + geom_tile(color = 'black') +
      theme(axis.text.x = element_text(angle=90, hjust = 1, vjust=0.5))
    print(p)
    
    # Analogous to https://rdrr.io/bioc/enrichplot/src/R/heatplot.R
    # I try to incorporate the pvalues (in the geneList variable)
    n <- enrichplot:::update_n(enrich, showCategory = 30)
    geneSets <- enrichplot:::extract_geneSets(enrich, n)
    #geneList <- enrichplot:::fc_readable(enrich, geneList)
    geneInfo.df <- enrichplot:::list2df(geneSets)
    
    geneInfo.df$pVal <- geneList[as.character(geneInfo.df[,2])]
    ## palette <- fc_palette(geneInfo.df$pVal)
    p_extra <- ggplot(geneInfo.df, aes(x = Gene, y = categoryID, fill = pVal)) +
      geom_tile() +
      scale_fill_continuous(low="blue", high="red", name = "p value", na.value="white") +
      theme(axis.text.x = element_text(angle=90, hjust = 1, vjust=0.5))
      ##scale_fill_gradientn(name = "p value", colors = palette)
    print(p_extra)
    #heatplot(enrich, pvalue=geneList)
  dev.off()
  
  # Tree plot
  enrich.extra <- pairwise_termsim(enrich, method="JC", semData=d)
  pdf(paste(tmpDir,"treeplot.pdf", sep="/"), width=20, height=14)
    p1 <- treeplot(enrich.extra)
    p2 <- treeplot(enrich.extra, hclust_method = "average")
    print(plot_grid(p1, p2, ncol=2))
  dev.off()
  

  # Enrichment Map
  # d <- godata('org.Hs.eg.db', ont="BP")
  pdf(paste(tmpDir,"emapplot.pdf", sep="/"), width=10, height=10)
    emapplot(enrich.extra,
             cex.params=list(category_node=1, category_label=0.5, line=0.5, pie2axis=1, label_group=0.5))
  dev.off()
  
  # Ridgeline plot - needs a gseaResult object (fgsea or somthing similar)
  # I could also convert my MSigDB results to this class??
  # For now, just generate a gseaResult object
  geneList <- tt.i[, paste0("P.Value (", comp,")")]
  names(geneList) <- tt.i[, "hgnc_symbol"]
  
  enrich.gseaResult <- gseGO(geneList, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont="BP", pvalueCutoff = 0.05, pAdjustMethod = "none", seed = TRUE,
                  minGSSize = 25)
  pdf(paste(tmpDir,"ridgeplot.pdf", sep="/"), width=10, height = 10)
    ridgeplot(enrich.gseaResult)
  dev.off()
  
  # And combine all PDFs into one
  pdfList <- list.files(path = tmpDir, pattern = "*.pdf", full.names = TRUE)
  qpdf::pdf_combine(input = pdfList,
                    output = paste(resultsDir,paste0("enrichplot_", comp, ".pdf"), sep="/"))
  

  # In the for-loop it only incorporates 3 PDFs, dotplot, heatplot and treeplot, although the others are also made...
  # Why??
  # If you run it just from the code it all works...
  # Even if you save to a temporary directory, it will only merge the three PDFs
  
  # And delete temporary files to make sure they will not end up in the wrong figure
  unlink(pdfList)
}

#### Session info ####
writeVersions <- function(sessionDir=outDir){
  write(paste0("Bioconductor version ", capture.output(tools:::.BioC_version_associated_with_R_version()),"\n"), 
        paste0(sessionDir,"/enrichplot_sessionInfo.txt"))
  write(capture.output(sessionInfo()), paste0(sessionDir,"/enrichplot_sessionInfo.txt"), append=TRUE)
}

writeVersions(sessionDir = resultsDir)

print(paste0("Bioconductor version ", capture.output(tools:::.BioC_version_associated_with_R_version()),"\n"))
sessionInfo()

#### Heatmap like plot of all comparisons ####
library(dplyr)
library(pheatmap)

resultsDir <- "15_pathway_visualization_20230419/enrichplot"
dir.create(resultsDir)

# Read in the data
v <- readRDS("14_libra_deg_20230426/v.AllComparisons.rds")
tt <- readRDS("14_libra_deg_20230426/tt_all.AllComparisons.rds")
my.contrasts <- v$contrasts

## Read the used genesets
BroadSets <- readRDS("14_libra_deg_GSEA_selectedSets_20230508/BroadSets_v2023.1.Hs.smallSet.rds")

# 20231024 - As Yosta now wants more and more of them I maybe should work in a loop :-) (see code below)
# see mails Wensenlijst 3.0 etc.
# And then again a new request on 26-10-2023 to also include PR3 vs MPO ... :-)

contrastSet <- list(One_vs_Rest = c("C0_vs_C1_C4andC8", "C1_vs_C0_C4andC8", "C4_vs_C0_C1andC8", "C8_vs_C0_C1andC4" ),
                    One_vs_Rest_TopX = c("C0_vs_C1_C4andC8", "C1_vs_C0_C4andC8", "C4_vs_C0_C1andC8", "C8_vs_C0_C1andC4" ),
                    ANCA_vs_CONTROL = c("C0_ANCA_vs_CONTROL", "C1_ANCA_vs_CONTROL", "C4_ANCA_vs_CONTROL", "C8_ANCA_vs_CONTROL" ),
                    ANCA_vs_CONTROL_TopX = c("C0_ANCA_vs_CONTROL", "C1_ANCA_vs_CONTROL", "C4_ANCA_vs_CONTROL", "C8_ANCA_vs_CONTROL" ),
                    SLE_vs_CONTROL = c("C0_SLE_vs_CONTROL", "C1_SLE_vs_CONTROL", "C4_SLE_vs_CONTROL", "C8_SLE_vs_CONTROL"),
                    SLE_vs_CONTROL_TopX = c("C0_SLE_vs_CONTROL", "C1_SLE_vs_CONTROL", "C4_SLE_vs_CONTROL", "C8_SLE_vs_CONTROL"),
                    SLE_vs_ANCA = c("C0_SLE_vs_ANCA", "C1_SLE_vs_ANCA", "C4_SLE_vs_ANCA", "C8_SLE_vs_ANCA"),
                    SLE_vs_ANCA_TopX = c("C0_SLE_vs_ANCA", "C1_SLE_vs_ANCA", "C4_SLE_vs_ANCA", "C8_SLE_vs_ANCA"),
                    MPO_ANCA_vs_PR3_ANCA = c("C0_ANCA_MPO_vs_ANCA_PR3", "C1_ANCA_MPO_vs_ANCA_PR3", "C4_ANCA_MPO_vs_ANCA_PR3",
                                        "C8_ANCA_MPO_vs_ANCA_PR3"),
                    MPO_ANCA_vs_PR3_ANCA_TopX = c("C0_ANCA_MPO_vs_ANCA_PR3", "C1_ANCA_MPO_vs_ANCA_PR3", "C4_ANCA_MPO_vs_ANCA_PR3",
                                             "C8_ANCA_MPO_vs_ANCA_PR3"))

# The cut-off for ANCA_vs_CONTROL_TopX
nrTop = 15

### Code does not fully work yet!!!
for ( numContrast in 1:length(contrastSet)){
  figTitle = names(contrastSet)[numContrast]
  selectedContrasts = unlist(contrastSet[numContrast])

  for (compNr in 1:length(selectedContrasts)){
    comp <- selectedContrasts[compNr]
    # Getting results and make the intersection
    df.i <- read.table(paste0("14_libra_deg_GSEA_selectedSets_20230508/reports_selectedSets/tab.camera.v2023.1.Hs_", comp,".txt"), header=TRUE)
    df.i$Description <- gsub("_", " ",gsub("C[1-8]+_|H_", "", df.i$geneset))
  
    df.i <- df.i[order(df.i$Description),] # WHY?? Does this mess up things?? Nope...
    # And move it to a nicer position ... cosmetics :-)
    df.i <- df.i %>% dplyr::relocate(Description, .before=NGenes)
    colnames(df.i)[colnames(df.i)=="PValue"] <- paste0("PValue_", comp)
    colnames(df.i)[colnames(df.i)=="FDR"] <- paste0("FDR_", comp)
    # Transform the values and include the sign
    df.i$log10FDR <- -log10(df.i[,paste0("FDR_", comp)])
    df.i$dirFDR <- ifelse(df.i$Direction == "Down", -1 * df.i$log10FDR, df.i$log10FDR )
    df.i[,paste0("transformed_FDR_", comp)] <- df.i$dirFDR
    if (compNr == 1){
      if (grepl("_TopX", figTitle)){
        # Take only top X and bottom X
        df.tmp <- df.i[order(df.i[,paste0("transformed_FDR_", comp)]),] 
        selectedGenesets <- df.tmp[c(1:nrTop, (nrow(df.tmp)-nrTop):nrow(df.tmp)),"Description"]
      } else {
        df <- df.i
      }
    } else {
      if (grepl("_TopX", figTitle)){
        # Take only top X and bottom X
        df.tmp <- df.i[order(df.i[,paste0("transformed_FDR_", comp)]),] 
        df.tmp <- df.tmp[c(1:nrTop, (nrow(df.tmp)-nrTop):nrow(df.tmp)),]
        selectedGenesets <- c(selectedGenesets,df.tmp[,"Description"])
      } else {
        df <- cbind(df, df.i[,c(paste0("PValue_", comp), paste0("FDR_", comp), paste0("transformed_FDR_", comp))])
      }
    }
  }
  
  if (grepl("_TopX", figTitle)){
    for (compNr in 1:length(selectedContrasts)){
      comp <- selectedContrasts[compNr]
      # Getting results and make the intersection
      df.i <- read.table(paste0("14_libra_deg_GSEA_selectedSets_20230508/reports_selectedSets/tab.camera.v2023.1.Hs_", comp,".txt"), header=TRUE)
      df.i$Description <- gsub("_", " ",gsub("C[1-8]+_|H_", "", df.i$geneset))
      df.i <- df.i[df.i$Description %in% unique(selectedGenesets),]
      df.i <- df.i[order(df.i$Description),] # WHY?? Does this mess up things?? Nope...
      # And move it to a nicer position ... cosmetics :-)
      df.i <- df.i %>% dplyr::relocate(Description, .before=NGenes)
      colnames(df.i)[colnames(df.i)=="PValue"] <- paste0("PValue_", comp)
      colnames(df.i)[colnames(df.i)=="FDR"] <- paste0("FDR_", comp)
      # Transform the values and include the sign
      df.i$log10FDR <- -log10(df.i[,paste0("FDR_", comp)])
      df.i$dirFDR <- ifelse(df.i$Direction == "Down", -1 * df.i$log10FDR, df.i$log10FDR )
      df.i[,paste0("transformed_FDR_", comp)] <- df.i$dirFDR
      if (compNr == 1){
        df <- df.i
      } else {
        df <- cbind(df, df.i[,c(paste0("PValue_", comp), paste0("FDR_", comp), paste0("transformed_FDR_", comp))])
      }
    }
  }
  
  for (pVal in c(0.01, 0.001)){
    
    # Get the row numbers of where a certain cut-off is being met
    criterion <- paste(paste0("df$FDR_", selectedContrasts, " < ",pVal), collapse=" | ")
    criterion
    
    # https://stackoverflow.com/questions/38804816/subset-a-dataframe-using-condition-passed-as-string-subset-dataframe-dynamicall
    df.subset <- df[eval(parse(text=criterion)),]
    rownames(df.subset) <- df.subset$Description
    
    mat <- as.matrix(df.subset[, grepl("transformed_FDR", colnames(df.subset))])
    colnames(mat) <- gsub("transformed_FDR_","", colnames(mat))
    #heatmap(mat)
    
    # # Half the rownames by putting a newline in the first space after a certain length
    # lengths <- lapply(rownames(mat), nchar)
    # longest <- unlist(lengths[which.max(lengths)])
    # rownames(mat) <- lapply(rownames(mat), function(x) paste(strwrap(x, width = longest/2), collapse = "\n"))
    
    if (figTitle == "One_vs_Rest"){
      fontsize_row = 10 - nrow(mat) / 15
    } else {
      if (grepl("_TopX", figTitle)){
        fontsize_row = 4
      } else {
        fontsize_row = 4
      }
    }
    
    myBreaks <- c(seq(max(mat),0, length.out = 51),seq(0,min(mat),length.out = 51))
    # '0' is now twice in it, so remove
    myBreaks <- unique(myBreaks)
    # myBreaks is 101 long, one more than the number of colours...
    # breaks - a sequence of numbers that covers the range of values in mat and 
    #          is one element longer than color vector. 
    
    if (grepl("_TopX", figTitle)){
      figTitleNew <- gsub("X",nrTop, figTitle)
      myTitle <- paste0("GSEA Top ", nrTop, " & FDR < ",pVal)
    } else {
      figTitleNew <- figTitle
      myTitle <- paste0("GSEA FDR < ",pVal)
    }
    
    pdf(paste(resultsDir,paste0("pheatmap_",figTitleNew,"_FDR", gsub("\\.","_",pVal),".pdf"), sep="/"), width=7, height=10)
    pheatmap(mat, breaks = rev(myBreaks),
             #legend_breaks = seq(round(min(mat))-1, round(max(mat))+1, 1), legend_labels = seq(round(min(mat))-1, round(max(mat))+1, 1),
             main=myTitle, cluster_cols=F, 
             fontsize_row=fontsize_row, border_color=NA)
    dev.off()
  }
  
}
