library(Seurat)
library(SeuratDisk)

# load preprocessed reference
obj <- LoadH5Seurat(file = "data/pbmc_reference/pbmc_multimodal.h5seurat")

# add raw RNA counts
counts <- Matrix::readMM("data/pbmc_reference/GSM5008737_RNA_3P-matrix.mtx.gz")
features <- read.table("data/pbmc_reference/GSM5008737_RNA_3P-features.tsv.gz", sep = "\t")
cells <- read.table("data/pbmc_reference/GSM5008737_RNA_3P-barcodes.tsv.gz", sep = "\t")

rownames(counts) <- features$V1
colnames(counts) <- cells$V1
counts <- as(counts, "dgCMatrix")

obj[["RNA"]] <- CreateAssayObject(counts = counts)
DefaultAssay(obj) <- "RNA"
obj <- FindVariableFeatures(obj)
obj <- NormalizeData(obj)

# compute SCT residuals for reference
obj <- SCTransform(
  object = obj,
  reference.SCT.model = slot(object = obj[['SCT']], name = "SCTModel.list")[[1]],
  residual.features = VariableFeatures(obj[["SCT"]]),
  assay = "RNA",
  new.assay.name = "SCT",
  verbose = FALSE
)
obj$dataset <- "reference"

saveRDS(obj, "objects/pbmc_reference.rds")
