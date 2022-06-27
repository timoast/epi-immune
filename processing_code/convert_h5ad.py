import anndata
import sys
import os
from scipy import io
import pandas as pd

basepath = sys.argv[1]

if not os.path.isdir(basepath + "/rna"):
    os.makedirs(basepath + "/rna")

rna = anndata.read_h5ad(basepath + "/multiome_gex_processed_training.h5ad")
rnacounts = rna.layers['counts']
io.mmwrite(target=basepath + "/rna/counts.mtx", a=rnacounts)
genenames = rna.var
rna_md = rna.obs
genenames.to_csv(basepath + "/rna/genes.csv")
rna_md.to_csv(basepath + "/rna/metadata.csv")


if not os.path.isdir(basepath + "/atac"):
    os.makedirs(basepath + "/atac")

atac = anndata.read_h5ad(basepath + "/multiome_atac_processed_training.h5ad")
ataccounts = atac.X
io.mmwrite(target=basepath + "/atac/counts.mtx", a=ataccounts)
peaknames = atac.var
atac_md = atac.obs
peaknames.to_csv(basepath + "/atac/peaks.csv")
atac_md.to_csv(basepath + "/atac/metadata.csv")
