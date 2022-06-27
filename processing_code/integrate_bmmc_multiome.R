library(Seurat)
library(Signac)
library(future)
plan("multicore", workers = 2)
options(future.globals.maxSize = 100 * 1024 ^ 3)

bmmc_multiome <- readRDS("objects/bmmc_multiome.rds")
obj <- SplitObject(bmmc_multiome, split.by = "batch")
ia <- FindIntegrationAnchors(
  object.list = obj,
  reduction = "rpca",
  assay = rep("SCT", length(obj)),
  dims = 1:30
)

integrated <- IntegrateEmbeddings(
  anchorset = ia,
  new.reduction.name = "integrated_pca",
  reductions = bmmc_multiome[['pca']]
)

saveRDS(integrated[['integrated_pca']], "./bmmc_multiome_integrated_pca.rds")

integrated <- IntegrateEmbeddings(
  anchorset = ia,
  new.reduction.name = "integrated_lsi",
  reductions = bmmc_multiome[['lsi']]
)

saveRDS(integrated[['integrated_lsi']], "./bmmc_multiome_integrated_lsi.rds")
