library(Signac)
library(EnsDb.Hsapiens.v86)
library(GenomeInfoDb)

# extract gene annotations from EnsDb
annotations.hg38 <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)

# change to UCSC style since the data was mapped to hg38
seqlevelsStyle(annotations.hg38) <- 'UCSC'
genome(annotations.hg38) <- "hg38"

# save
saveRDS(object = annotations.hg38, file = "data/annotations_hg38.rds")