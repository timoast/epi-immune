library(Signac)
library(Seurat)
library(Matrix)


rna_counts <- readMM(file = "data/bmmc_multiome/multiome/rna/counts.mtx")
rna_metadata <- read.table(file = "data/bmmc_multiome/multiome/rna/metadata.csv", sep = ",", header = TRUE, row.names = 1)
genes <- read.table(file = "data/bmmc_multiome/multiome/rna/genes.csv", sep = ",", header = TRUE, row.names = 1)

rna_counts <- t(rna_counts)
rna_counts <- as(rna_counts, "dgCMatrix")
rownames(rna_counts) <- rownames(genes)
colnames(rna_counts) <- rownames(rna_metadata)

multiome <- CreateSeuratObject(counts = rna_counts, project = "multiome", meta.data = rna_metadata)

atac_counts <- readMM(file = "data/bmmc_multiome/multiome/atac/counts.mtx")
atac_metadata <- read.table(file = "data/bmmc_multiome/multiome/atac/metadata.csv", sep = ",", header = TRUE, row.names = 1)
peaks <- read.table(file = "data/bmmc_multiome/multiome/atac/peaks.csv", sep = ",", header = TRUE, row.names = 1)

atac_counts <- t(atac_counts)
atac_counts <- as(atac_counts, "dgCMatrix")
rownames(atac_counts) <- rownames(peaks)
colnames(atac_counts) <- rownames(atac_metadata)

common_cells <- intersect(colnames(multiome), colnames(atac_counts))

atac_counts <- atac_counts[, common_cells]
multiome <- multiome[, common_cells]

multiome[["peaks"]] <- CreateChromatinAssay(counts = atac_counts)

# process RNA
DefaultAssay(multiome) <- "RNA"

multiome <- FindVariableFeatures(multiome)
multiome <- NormalizeData(multiome)
multiome <- SCTransform(multiome)
multiome <- RunPCA(multiome)
multiome <- RunUMAP(
  object = multiome,
  reduction = "pca",
  dims = 1:40,
  reduction.name = "umap.sct",
  return.model = TRUE
)

# process ATAC
DefaultAssay(multiome) <- "peaks"

multiome <- FindTopFeatures(multiome)
multiome <- RunTFIDF(multiome)
multiome <- RunSVD(multiome)
multiome <- RunUMAP(
  object = multiome,
  reduction = "lsi",
  dims = 2:40,
  reduction.name = "umap.atac",
  return.model = TRUE
)

saveRDS(object = multiome, file = "objects/bmmc_multiome.rds")
