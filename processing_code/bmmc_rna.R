library(Seurat)
library(SummarizedExperiment)
library(SingleCellExperiment)

obj <- readRDS("data/bmmc_rna/scRNA-Healthy-Hematopoiesis-191120.rds")
obj <- as(object = obj, Class = "SingleCellExperiment")
obj.seu <- as.Seurat(x = obj, counts = "counts", data = NULL)

obj.seu <- FindVariableFeatures(obj.seu)
obj.seu <- NormalizeData(obj.seu)
obj.seu <- ScaleData(obj.seu)
obj.seu <- RunPCA(obj.seu)
obj.seu <- RunUMAP(obj.seu, reduction = "pca", dims = 1:30)

saveRDS(obj.seu, "objects/bmmc_rna.rds")
