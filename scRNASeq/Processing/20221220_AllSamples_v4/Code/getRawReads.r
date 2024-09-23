for ( sc in sampleList){
  seurat <- readRDS(paste0("../Data/scRNASeq/Processed/", sc, ".raw.annot.rds"))

  # Extract some needed values
  counts_per_cell <- Matrix::colSums(seurat@assays$RNA@counts)
  counts_per_gene <- Matrix::rowSums(seurat@assays$RNA@counts)
  genes_per_cell <- Matrix::colSums(seurat@assays$RNA@counts > 0) # count gene only if it has non-zero reads mapped.
  cells_per_gene <- Matrix::rowSums(seurat@assays$RNA@counts > 0)
  
  cat(paste0(sc,":\t #cells:\t", length(counts_per_cell),"\n"))
  cat(paste0("\tMean:\t", mean(counts_per_cell), "\n"))
  cat(paste0("\tMedian:\t", median(counts_per_cell), "\n"))
}