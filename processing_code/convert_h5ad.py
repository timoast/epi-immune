import anndata
import sys
import os
from scipy import io
import pandas as pd

basepath = sys.argv[1]

bmmc = anndata.read_h5ad(basepath + "/GSE194122_openproblems_neurips2021_multiome_BMMC_processed.h5ad")
df = bmmc.var
atac_features = df.index[df['feature_types'] == 'ATAC'].tolist()
rna_features = df.index[df['feature_types'] == 'GEX'].tolist()
atac = bmmc[:, atac_features]
rna = bmmc[:, rna_features]

# RNA
if not os.path.isdir(basepath + "/rna"):
    os.makedirs(basepath + "/rna")
rnacounts = rna.layers["counts"]
io.mmwrite(target=basepath + "/rna/counts.mtx", a=rnacounts)
genenames = rna.var
rna_md = rna.obs
genenames.to_csv(basepath + "/rna/genes.csv")
rna_md.to_csv(basepath + "/rna/metadata.csv")

# ATAC
if not os.path.isdir(basepath + "/atac"):
    os.makedirs(basepath + "/atac")
ataccounts = atac.layers["counts"]
io.mmwrite(target=basepath + "/atac/counts.mtx", a=ataccounts)
peaknames = atac.var
atac_md = atac.obs
peaknames.to_csv(basepath + "/atac/peaks.csv")
atac_md.to_csv(basepath + "/atac/metadata.csv")
