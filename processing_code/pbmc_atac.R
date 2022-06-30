library(Seurat)
library(Signac)
library(GenomicRanges)

args = commandArgs(trailingOnly = TRUE)
ncore <- as.numeric(args[1])
annot.path <- args[2]
peak.path <- args[3]
sample_name <- args[4]

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
countfile <- list.files(paste0("data/pbmc_atac/", sample_name), pattern = "*.h5$", full.names = TRUE)
frags <- list.files(paste0("data/pbmc_atac/", sample_name), pattern = "*.tsv.gz$", full.names = TRUE)
md_file <- list.files(paste0("data/pbmc_atac/", sample_name), patter = "*.csv$", full.names = TRUE)

counts <- Read10X_h5(filename = countfile)
metadata <- read.table(file = md_file, sep = ",", header = TRUE, row.names = 1)

assay <- CreateChromatinAssay(counts = counts, fragments = frags, genome = "hg38", annotation = annot, sep = c(":", "-"))
obj <- CreateSeuratObject(counts = assay, meta.data = metadata, assay = "ATAC")

# process ATAC
DefaultAssay(obj) <- "ATAC"
obj <- NucleosomeSignal(obj)
obj <- TSSEnrichment(obj)

# add gene activity
ga <- GeneActivity(object = obj)
obj[["GA"]] <- CreateAssayObject(counts = ga)
obj <- NormalizeData(object = obj, assay = 'GA')

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

# save object
saveRDS(object = obj, file = paste0("objects/", sample_name, ".rds"))