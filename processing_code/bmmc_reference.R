library(Seurat)

obj <- readRDS("data/bmmc_reference/fullref.Rds")
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj)
obj <- RunPCA(obj)
obj <- RunUMAP(
  object = obj,
  reduction = "pca",
  dims = 1:40,
  reduction.name = "umap",
  return.model = TRUE
)
obj$dataset <- 'reference'
saveRDS(obj, "objects/bmmc_reference.rds")
