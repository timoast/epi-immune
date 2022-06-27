library(Seurat)
library(Signac)
library(future)
plan("multicore", workers = 8)
options(future.globals.maxSize = 100 * 1024 ^ 3)

bmmc_atac <- readRDS("objects/bmmc_atac.rds")

atac <- SplitObject(bmmc_atac, split.by = "Group")
atac <- lapply(atac, function(x) {
  x <- FindTopFeatures(x)
  x <- RunTFIDF(x)
  x <- RunSVD(x)
  return(x)
})

atac_anchors <- FindIntegrationAnchors(
  object.list = atac,
  reduction = "rlsi",
  dims = 2:30
)

atac_integrated <- IntegrateEmbeddings(
  anchorset = atac_anchors,
  new.reduction.name = "integrated_lsi",
  reductions = bmmc_atac[["lsi"]],
  dims.to.integrate = 1:30
)

saveRDS(atac_integrated[['integrated_lsi']], "./bmmc_atac_integrated_lsi.rds")