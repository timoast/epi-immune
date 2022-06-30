library(Seurat)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
samplename <- as.character(args[[1]])

obj <- readRDS(paste0("objects/", samplename, ".rds"))
if (grepl("reference", samplename)) {
  gb <- 'celltype.l2'
} else {
  annot <- read.table(paste0("annotations/", samplename, ".tsv"), sep = "\t")
  obj <- AddMetaData(obj, annot)
  gb <- 'predicted.id'
}

p <- DimPlot(
  object = obj,
  label = TRUE,
  repel = TRUE,
  group.by = gb,
  pt.size = 0.1
) + NoLegend()

ggsave(
  filename = paste0("plots/", samplename, ".png"),
  plot = p,
  height = 6,
  width = 7
)
