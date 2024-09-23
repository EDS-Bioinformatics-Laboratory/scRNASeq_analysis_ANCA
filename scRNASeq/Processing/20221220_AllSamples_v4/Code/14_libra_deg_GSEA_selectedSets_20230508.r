#### See 14_libra_deg.r ####
# Here, we also include cluster 0 (see mail 20230414)
# And put all comparisons into one directory etc

#### Set default for a good start ####
source("initAnalysis.r")
init()

# Read data from the previous analysis
dataDir <- "14_libra_deg_20230426/"

# Create a directory to hold results and plots
#outDir <- "14_libra_deg_GSEA_selectedSets_20230508"
outDir <- "14_libra_deg_GSEA_selectedSets_20230526"
dir.create(outDir)

#### Library ####
library(edgeR)
library(limma)
library(ggplot2)

#### Enrichment analysis #####
library(GSEABase)

v <- readRDS(paste(dataDir,"v.AllComparisons.rds", sep="/"))
targets <- v$targets
design <- v$design
my.contrasts <- v$contrasts
fit <- readRDS(paste(dataDir,"fit.AllComparisons.rds", sep="/"))

## Read the MSigDB database v2023.1.Hs genesets,
# But we only want the HALLMARK, REACTOME, KEGG, BIOCARTA pathways, so H and C2
hBroadSets = getGmt(con=paste(dropbox,"Support/MSigDB/v2023.1.Hs/h.all.v2023.1.Hs.entrez.gmt",sep="/"),geneIdType=EntrezIdentifier(), collectionType=BroadCollection(category="h"))
c2BroadSets = getGmt(con=paste(dropbox,"Support/MSigDB/v2023.1.Hs/c2.all.v2023.1.Hs.entrez.gmt",sep="/"),geneIdType=EntrezIdentifier(), collectionType=BroadCollection(category="c2"))
H = geneIds(hBroadSets)
names(H) = paste0("H_",names(H))
C2 = geneIds(c2BroadSets)
names(C2) = paste0("C2_",names(C2))

BroadSets = c(H,C2)

# Make the selection on REACTOME, KEGG, BIOCARTA
BroadSets <- BroadSets[grepl("HALLMARK|REACTOME|KEGG|BIOCARTA", names(BroadSets))]
length(BroadSets)
#[1] 2182

saveRDS(BroadSets, file=paste0(outDir,"/BroadSets_v2023.1.Hs.smallSet.rds"))

BroadSets <- readRDS(paste0(outDir,"/BroadSets_v2023.1.Hs.smallSet.rds"))

# This will couple the Entrez IDs to the BroadSets defined earlier)
sets = ids2indices(BroadSets,unlist(v$genes$entrezgene))

file.copy(paste(dropbox,"Support/R/barcodeplot.r",sep="/"),scriptDir)
source(paste(scriptDir,"barcodeplot.r", sep="/"))

file.copy(paste(dropbox,"Support/R/cameraReports.r",sep="/"),scriptDir)
source(paste(scriptDir,"cameraReports.r", sep="/"))

cameraReport(v.=v, fit.=fit, constrasts=my.constrasts, trend=FALSE, robust=FALSE, msigdb="2023.1.Hs", species="human", max.pathways=100,
             reportsDir = "reports_selectedSets")

# Relocate 'reports_selectedSets' directory
file.rename("reports_selectedSets/", to=paste0(outDir,"/reports_selectedSets"))

# And just for the HALLMARK set
BroadSets <- BroadSets[grepl("HALLMARK", names(BroadSets))]
length(BroadSets)
#[1] 50

saveRDS(BroadSets, file=paste0(outDir,"/BroadSets_v2023.1.Hs.HALLMARK.rds"))

BroadSets <- readRDS(paste0(outDir,"/BroadSets_v2023.1.Hs.HALLMARK.rds"))

# This will couple the Entrez IDs to the BroadSets defined earlier)
sets = ids2indices(BroadSets,unlist(v$genes$entrezgene))

cameraReport(v.=v, fit.=fit, constrasts=my.constrasts, trend=FALSE, robust=FALSE, msigdb="2023.1.Hs", species="human", max.pathways=100,
             reportsDir = "reports_onlyHALLMARK")

# Relocate 'reports_reports_onlyHALLMARK' directory
file.rename("reports_onlyHALLMARK/", to=paste0(outDir,"/reports_onlyHALLMARK"))


#### Genesets - Venn Diagrams ####
file.copy(paste(dropbox,"Support/R/vennReportsGeneSets.r",sep="/"),scriptDir)
source(paste(scriptDir,"vennReportsGeneSets.r", sep="/"))

v <- readRDS(paste(dataDir,"v.AllComparisons.rds", sep="/"))
targets <- v$targets
design <- v$design
my.contrasts <- v$contrasts
fit <- readRDS(paste(dataDir,"fit.AllComparisons.rds", sep="/"))

# NOTE
# - As I get an error message that the filename is too long (while 'reading' the file!!!), I will put files into
#   a 'tmp' directory and then later try to move them to the 'right' directory?

for (my_dir_use in c("reports_onlyHALLMARK", "reports_selectedSets")){
  toBeMoved <- list.files(paste0(outDir,"/", my_dir_use), pattern="tab.camera.v2023.1.Hs_*.*", full.names = TRUE)
  dir.create(paste0(outDir, "/tmp"))
  file.copy(toBeMoved, to=paste0(outDir, "/tmp"))
  
  my_dir = "tmp"
  
  vennReportsGeneSets(p.value=0.05, fdr="BH", filePrefix = paste0(outDir,"/", my_dir,"/tab.camera.v2023.1.Hs_"),
                      contrasts = my.contrasts,myTitle="vennGeneSet_Cluster0",
                      indices = c(1:4), scriptDir = ".")
  
  vennReportsGeneSets(p.value=0.05, fdr="BH", filePrefix = paste0(outDir,"/", my_dir,"/tab.camera.v2023.1.Hs_"),
                      contrasts = my.contrasts,myTitle="vennGeneSet_Cluster1",
                      indices = c(5:8), scriptDir = ".")
  
  vennReportsGeneSets(p.value=0.05, fdr="BH", filePrefix = paste0(outDir,"/", my_dir,"/tab.camera.v2023.1.Hs_"),
                      contrasts = my.contrasts, myTitle="vennGeneSet_Cluster4",
                      indices = c(9:12), scriptDir = ".")
  
  vennReportsGeneSets(p.value=0.05, fdr="BH", filePrefix = paste0(outDir,"/", my_dir,"/tab.camera.v2023.1.Hs_"),
                      contrasts = my.contrasts, myTitle="vennGeneSet_Cluster8",
                      indices = c(13:16), scriptDir = ".")
  
  vennReportsGeneSets(p.value=0.05, fdr="BH", filePrefix = paste0(outDir,"/", my_dir,"/tab.camera.v2023.1.Hs_"),
                      contrasts = my.contrasts, indices = c(4,8,12,16), scriptDir = ".",
                      myTitle="vennGeneSet_ANCA_MPO_vs_PR3")
  
  vennReportsGeneSets(p.value=0.05, fdr="BH", filePrefix = paste0(outDir,"/", my_dir,"/tab.camera.v2023.1.Hs_"),
                      contrasts = my.contrasts, indices = c(17:19), scriptDir = ".",
                      myTitle="vennGeneSet_ANCA_betweenClusters_part1")
  
  vennReportsGeneSets(p.value=0.05, fdr="BH", filePrefix = paste0(outDir,"/", my_dir,"/tab.camera.v2023.1.Hs_"),
                      contrasts = my.contrasts, indices = c(20:22), scriptDir = ".",
                      myTitle="vennGeneSet_ANCA_betweenClusters_part2")
  
  vennReportsGeneSets(p.value=0.05, fdr="BH", filePrefix = paste0(outDir,"/", my_dir,"/tab.camera.v2023.1.Hs_"),
                      contrasts = my.contrasts, indices = c(23:25), scriptDir = ".",
                      myTitle="vennGeneSet_SLE_betweenClusters_part1")
  
  vennReportsGeneSets(p.value=0.05, fdr="BH", filePrefix = paste0(outDir,"/", my_dir,"/tab.camera.v2023.1.Hs_"),
                      contrasts = my.contrasts, indices = c(26:28), scriptDir = ".",
                      myTitle="vennGeneSet_SLE_betweenClusters_part2")
  
  vennReportsGeneSets(p.value=0.05, fdr="BH", filePrefix = paste0(outDir,"/", my_dir,"/tab.camera.v2023.1.Hs_"),
                      contrasts = my.contrasts, indices = c(29:31), scriptDir = ".",
                      myTitle="vennGeneSet_CONTROLE_betweenClusters_part1")
  
  vennReportsGeneSets(p.value=0.05, fdr="BH", filePrefix = paste0(outDir,"/", my_dir,"/tab.camera.v2023.1.Hs_"),
                      contrasts = my.contrasts, indices = c(32:34), scriptDir = ".",
                      myTitle="vennGeneSet_CONTROLE_betweenClusters_part2")
  
  vennReportsGeneSets(p.value=0.05, fdr="BH", filePrefix = paste0(outDir,"/", my_dir,"/tab.camera.v2023.1.Hs_"),
                      contrasts = my.contrasts, indices = c(35:37), scriptDir = ".",
                      myTitle="vennGeneSet_ANCA_MPO_betweenClusters_part1")
  
  vennReportsGeneSets(p.value=0.05, fdr="BH", filePrefix = paste0(outDir,"/", my_dir,"/tab.camera.v2023.1.Hs_"),
                      contrasts = my.contrasts, indices = c(38:40), scriptDir = ".",
                      myTitle="vennGeneSet_ANCA_MPO_betweenClusters_part2")
  
  vennReportsGeneSets(p.value=0.05, fdr="BH", filePrefix = paste0(outDir,"/", my_dir,"/tab.camera.v2023.1.Hs_"),
                      contrasts = my.contrasts, indices = c(41:43), scriptDir = ".",
                      myTitle="vennGeneSet_ANCA_PR3_betweenClusters_part1")
  
  vennReportsGeneSets(p.value=0.05, fdr="BH", filePrefix = paste0(outDir,"/", my_dir,"/tab.camera.v2023.1.Hs_"),
                      contrasts = my.contrasts, indices = c(44:46), scriptDir = ".",
                      myTitle="vennGeneSet_ANCA_PR3_betweenClusters_part2")
  
  vennReportsGeneSets(p.value=0.05, fdr="BH", filePrefix = paste0(outDir,"/", my_dir,"/tab.camera.v2023.1.Hs_"),
                      contrasts = my.contrasts, indices = c(47:49), scriptDir = ".",
                      myTitle="vennGeneSet_Clusters_part1")
  
  vennReportsGeneSets(p.value=0.05, fdr="BH", filePrefix = paste0(outDir,"/", my_dir,"/tab.camera.v2023.1.Hs_"),
                      contrasts = my.contrasts, indices = c(50:52), scriptDir = ".",
                      myTitle="vennGeneSet_Clusters_part2")
  
  vennReportsGeneSets(p.value=0.05, fdr="BH", filePrefix = paste0(outDir,"/", my_dir,"/tab.camera.v2023.1.Hs_"),
                      contrasts = my.contrasts, indices = c(53:56), scriptDir = ".",
                      myTitle="vennGeneSet_Clusters_One_vs_Rest")
  
  vennReportsGeneSets(p.value=0.05, fdr="BH", filePrefix = paste0(outDir,"/", my_dir,"/tab.camera.v2023.1.Hs_"),
                      contrasts = my.contrasts, indices = c(57:58), scriptDir = ".",
                      myTitle="vennGeneSet_Cluster_C0_or_C1_vs_C4andC8")
  
  vennReportsGeneSets(p.value=0.05, fdr="BH", filePrefix = paste0(outDir,"/", my_dir,"/tab.camera.v2023.1.Hs_"),
                      contrasts = my.contrasts, indices = c(59:60), scriptDir = ".",
                      myTitle="vennGeneSet_Cluster_C4_or_C8_vs_C0andC1")
  
  vennReportsGeneSets(p.value=0.05, fdr="BH", filePrefix = paste0(outDir,"/", my_dir,"/tab.camera.v2023.1.Hs_"),
                      contrasts = my.contrasts, indices = c(61), scriptDir = ".",
                      myTitle="vennGeneSet_Clusters_C0andC1_vs_C4andC8")
  
  # Relocate 'venn_gs' directory
  file.rename("venn_gs", to=paste0(outDir,"/venn_gs", gsub("reports","", my_dir_use)))
  
  # Relocate figures and stats
  toBeMoved <- list.files(".", pattern="*.png|*.html", full.names = TRUE)
  file.copy(toBeMoved, to=paste0(outDir, "/venn_gs", gsub("reports","", my_dir_use)))
  file.remove(toBeMoved)
  
  # And relocate the tt.genesets.rds
  toBeMoved <- list.files(".", pattern="*.genesets.rds", full.names = TRUE)
  file.copy(toBeMoved, to=paste0(outDir, "/venn_gs", gsub("reports","", my_dir_use)))
  file.remove(toBeMoved)
  
  # And remove the 'tmp' directory
  unlink(paste0(outDir, "/tmp"),recursive=TRUE)
}

#### Extra - Get function of clusters? ####
v <- readRDS(paste(dataDir,"v.AllComparisons.rds", sep="/"))
targets <- v$targets
design <- v$design
my.contrasts <- v$contrasts
fit <- readRDS(paste(dataDir,"fit.AllComparisons.rds", sep="/"))

extraComps <- colnames(my.contrasts)[grep("and", colnames(my.contrasts))]

for (my_dir in c("reports_onlyHALLMARK", "reports_selectedSets")){
  for (comp in c("C1_vs_C0", "C4_vs_C0", "C8_vs_C0", "C4_vs_C1", "C8_vs_C1", "C8_vs_C4", extraComps)){
    # Read in the cluster comparison
    df <- read.table(paste0(outDir, paste0("/", my_dir,"/tab.camera.v2023.1.Hs_", comp,".txt")), header=TRUE)
  
    # Select different categories or top X (use "topX")
    if (my_dir == "reports_selectedSets"){
      df.selected <- df[1:50, ]
    } else {
      # Just get the 50 HALLMARK sets
      df.selected <- df
    }
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
      ggtitle(paste0("Gene set enrichment analysis - ", gsub("andC", " and ", gsub("_vs_C", " vs. ", gsub("^C", "Cluster ", comp))) )) +
      theme(plot.title = element_text(hjust = 1))
    
    pdf(paste(outDir,paste0("/", my_dir,"/Fig_barPlot_GenesetEnrichment_", comp,".pdf"), sep="/"), 
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