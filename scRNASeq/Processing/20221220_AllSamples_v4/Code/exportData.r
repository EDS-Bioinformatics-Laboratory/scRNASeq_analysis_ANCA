# Export data from Seurat object
library(data.table)

# 
dataDir <- "9_kidney_reference_mapping/"
seurat <- readRDS(paste(dataDir, "Mito10.integrated.clustered.annot.kidney_immune_reference.rds", sep="/"))

for (res in c("0.3", "0.6", "1")){
  Idents(seurat) <- paste0("integrated_snn_res.",res)
  res <- gsub("\\.","_", res)
  adt_markers <- FindAllMarkers(seurat, assay = "ADT")
  fwrite(x = adt_markers, row.names=TRUE, file = paste0("ADT_findAllMarkers_beforeReclustering_res",res,".csv"), sep = "\t")

  adt_avgExprs <- data.frame(AverageExpression(seurat, assay = "ADT"))
  colnames(adt_avgExprs) <- gsub("ADT.", "Cluster_", colnames(adt_avgExprs))
  fwrite(x = adt_avgExprs, row.names=TRUE, file = paste0("ADT_avgExprs_beforeReclustering_res", res,".csv"), sep = "\t")
}

# And after reclustering
dataDir <- "10_reclustering_MP_and_MC"
seurat <- readRDS(paste(dataDir, "selectedClusters.subjectIntegrated.reclustered.kidney_immune_reference.rds", sep="/"))

for (res in c("0.3", "0.6", "1")){
  Idents(seurat) <- paste0("integrated_snn_res.",res)
  res <- gsub("\\.","_", res)
  adt_markers <- FindAllMarkers(seurat, assay = "ADT")
  fwrite(x = adt_markers, row.names=TRUE, file = paste0("ADT_findAllMarkers_afterReclustering_res",res,".csv"), sep = "\t")
  
  adt_avgExprs <- data.frame(AverageExpression(seurat, assay = "ADT"))
  colnames(adt_avgExprs) <- gsub("ADT.", "Cluster_", colnames(adt_avgExprs))
  fwrite(x = adt_avgExprs, row.names=TRUE, file = paste0("ADT_avgExprs_afterReclustering_res", res,".csv"), sep = "\t")
}

