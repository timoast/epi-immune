library(Signac)
library(Seurat)


args <- commandArgs(trailingOnly = TRUE)
samplename <- as.character(args[[1]])
query <- readRDS(paste0("objects/", samplename, ".rds"))

if (grepl('reference', samplename)) {
  # reference object given as input
  pred <- query[[]]
  write.table(
    x = pred,
    file = paste0("annotations/", samplename, ".tsv"),
    sep = "\t",
    quote = FALSE
  )
} else {
  # map query to reference
  
  is_multiome <- grepl("multiome", samplename)
  is_pbmc <- grepl("pbmc", samplename)
  
  if (is_pbmc) {
    reference <- readRDS("objects/pbmc_reference.rds")
  } else {
    reference <- readRDS("objects/bmmc_reference.rds")
  }
  
  if (is_multiome & is_pbmc) {
    # pbmc multiome
    normalization.method <- "SCT"
    recompute.residuals <- TRUE
    query.assay <- 'RNA'
    reference.assay <- 'SCT'
    reduction <- 'pcaproject'
    weight.reduction <- 'pcaproject'
    wt.dims <- 1:40
  } else if (is_multiome) {
    # bmmc multiome
    normalization.method <- "LogNormalize"
    recompute.residuals <- FALSE
    query.assay <- 'RNA'
    reference.assay <- 'RNA'
    reduction <- 'pcaproject'
    weight.reduction <- 'pcaproject'
    wt.dims <- 1:40
  } else {
    # atac
    normalization.method <- 'LogNormalize'
    recompute.residuals <- FALSE
    query.assay <- 'GA'
    reference.assay <- 'RNA'
    reduction <- 'cca'
    weight.reduction <- query[['lsi']]
    wt.dims <- 2:30
  }
  
  anchors <- FindTransferAnchors(
    reference = reference,
    query = query,
    normalization.method = normalization.method,
    recompute.residuals = recompute.residuals,
    reference.assay = reference.assay,
    query.assay = query.assay,
    reduction = reduction,
    dims = 1:40
  )
  
  pred <- TransferData(
    anchorset = anchors,
    refdata = reference$celltype.l2,
    weight.reduction = weight.reduction,
    dims = wt.dims
  )
  
  write.table(
    x = pred,
    file = paste0("annotations/", samplename, ".tsv"),
    sep = "\t",
    quote = FALSE
  )
}
