library(Seurat)
library(Signac)
library(GenomicRanges)

args = commandArgs(trailingOnly = TRUE)
ncore <- as.numeric(args[1])
annot.path <- args[2]
peak.path <- args[3]
sample_name <- args[4]
ref <- args[5]

annot <- readRDS(annot.path)
peaks <- read.table(file = peak.path, col.names = c("chrom", "start", "end"))
peaks <- makeGRangesFromDataFrame(df = peaks)
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")

message("Using ", ncore, " cores")
if (ncore > 1) {
  library(future)
  plan("multicore", workers = ncore)
  options(future.globals.maxSize = 100 * 1024^3)
}

# create object
countfile <- list.files(paste0("data/pbmc_multiome/", sample_name), pattern = "*.h5$", full.names = TRUE)
frags <- list.files(paste0("data/pbmc_multiome/", sample_name), pattern = "*.tsv.gz$", full.names = TRUE)
md_file <- list.files(paste0("data/pbmc_multiome/", sample_name), patter = "*.csv$", full.names = TRUE)

counts <- Read10X_h5(filename = countfile)
metadata <- read.table(file = md_file, sep = ",", header = TRUE, row.names = 1)

obj <- CreateSeuratObject(counts = counts$`Gene Expression`, assay = "RNA", meta.data = metadata)
obj[["ATAC"]] <- CreateChromatinAssay(counts = counts$Peaks, fragments = frags, genome = "hg38", annotation = annot, sep = c(":", "-"))

# filter low/high count cells
obj <- subset(obj, subset = nCount_ATAC > 2000 &
                nCount_ATAC < 70000 &
                nCount_RNA > 500 & 
                nCount_RNA < 15000)

# process ATAC
DefaultAssay(obj) <- "ATAC"
obj <- NucleosomeSignal(obj)
obj <- TSSEnrichment(obj)

# create assay containing shared peaks
common_counts <- FeatureMatrix(
  fragments = Fragments(obj)[[1]],
  features = peaks,
  cells = colnames(obj)
)
obj[["common"]] <- CreateChromatinAssay(counts = common_counts, fragments = frags, annotation = annot, genome = "hg38")

DefaultAssay(obj) <- "common"
obj <- FindTopFeatures(obj, min.cutoff = "q5")
obj <- RunTFIDF(obj)
obj <- RunSVD(obj)
obj <- RunUMAP(obj, reduction = "lsi", dims = 2:30, reduction.name = "umap.atac", store.model = TRUE)

# process RNA
DefaultAssay(obj) <- "RNA"
obj <- FindVariableFeatures(obj)
obj <- NormalizeData(obj)
obj <- ScaleData(obj)
obj <- RunPCA(obj)
obj <- RunUMAP(obj, reduction = "pca", dims = 1:30, reduction.name = "umap.rna", store.model = TRUE)

# map to reference
reference <- readRDS(ref)

anchors <- FindTransferAnchors(
  reference = reference,
  query = obj,
  reference.assay = "SCT",
  query.assay = "RNA",
  dims = 1:30,
  normalization.method = "SCT",
  recompute.residuals = TRUE,
  reduction = 'cca'
)

pred <- TransferData(
    anchorset = anchors,
    refdata = reference$celltype.l2,
    weight.reduction = obj[['pca']],
    dims = 1:30
)
obj <- AddMetaData(object = obj, metadata = pred)

# save object
saveRDS(object = obj, file = paste0("objects/", sample_name, ".rds"))
