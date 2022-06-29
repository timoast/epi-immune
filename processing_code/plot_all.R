library(Seurat)
library(ggplot2)

obj.list <- list.files("objects", full.names = FALSE)

for (i in seq_along(obj.list)) {
  o <- readRDS(paste0("objects/", obj.list[[i]]))
  p <- DimPlot(o, label = FALSE, pt.size = 0.1) + NoLegend()
  pname <- gsub(pattern = ".rds", replacement = ".png", x = obj.list[[i]])
  ggsave(filename = paste0("plots/", pname), plot = p, height = 4, width = 5)
}
