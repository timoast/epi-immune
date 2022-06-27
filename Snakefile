pbmc_atac = [
        "pbmc_atac_10k_chromium",
        "pbmc_atac_10k_chromiumX",
        "pbmc_atac_10k",
        "pbmc_atac_5k",
        "pbmc_atac_1k",
        "pbmc_atac_500"
        ]

pbmc_multiome = [
        "pbmc_multiome_3k_sorted",
        "pbmc_multiome_3k_unsorted",
        "pbmc_multiome_10k_sorted",
        "pbmc_multiome_10k_unsorted",
        "pbmc_multiome_10k_chromium",
        "pbmc_multiome_10k_chromiumX"
        ]

bmmc_atac = [
        "bmmc_atac"
        ]

bmmc_multiome = [
        "bmmc_multiome"
        ]

all_samples = pbmc_atac + pbmc_multiome + bmmc_atac + bmmc_multiome

rule all:
  input:
      expand("objects/{sample}.rds", sample = all_samples)

# ---- Download ---- #

rule download_pbmc_multiome:
    input: "datasets/pbmc_multiome/{dset}.txt"
    output: touch("data/pbmc_multiome/{dset}/download.done")
    message: "Download PBMC multiome datasets"
    threads: 1
    shell:
        """
        wget -i {input} -P data/pbmc_multiome/{wildcards.dset}
        """

rule download_pbmc_atac:
    input: "datasets/pbmc_atac/{dset}.txt"
    output: touch("data/pbmc_atac/{dset}/download.done")
    message: "Download PBMC scATAC-seq datasets"
    threads: 1
    shell:
        """
        wget -i {input} -P data/pbmc_atac/{wildcards.dset}
        """

rule download_bmmc_atac:
    input: "datasets/bmmc_atac/bmmc_atac.txt"
    output: touch("data/bmmc_atac/download.done")
    threads: 1
    message: "Download BMMC scATAC-seq hg38 fragment files"
    shell:
        """
        wget -i {input} -P data/bmmc_atac
        aws s3 sync s3://mpal-hg38/public/ ./data/bmmc_atac/ --request-payer
        """

rule download_bmmc_multiome:
    output:
        "data/bmmc_multiome/multiome/multiome_atac_processed_training.h5ad",
        "data/bmmc_multiome/multiome/multiome_gex_processed_training.h5ad"
    message: "Download BMMC multiome"
    threads: 1
    shell:
        """
        aws s3 sync s3://openproblems-bio/public/explore ./data/bmmc_multiome/ --request-payer
        """

rule download_bmmc_reference:
    output: "objects/bmmc_reference.rds"
    message: "Download BMMC reference object"
    threads: 1
    shell:
        """
        aws s3 cp s3://bmmc-reference/public/fullref.Rds {output} --no-sign-request
        """

rule download_pbmc_reference:
    input: "datasets/pbmc_reference/{dset}.txt"
    output: touch("data/{dset}/download.done")
    message: "Download PBMC reference"
    threads: 1
    shell:
        """
        wget -i {input} -P data/{wildcards.dset}
        """

rule get_peaks:
  output: "data/pbmc_multiome_peaks.bed"
  message: "Downloading peaks"
  threads: 1
  shell: "wget https://cf.10xgenomics.com/samples/cell-arc/2.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_atac_peaks.bed -O {output}"

# ---- Process ---- #

rule convert_bmmc_multiome:
    input:
        "data/bmmc_multiome/multiome/multiome_atac_processed_training.h5ad",
        "data/bmmc_multiome/multiome/multiome_gex_processed_training.h5ad"
    output:
        "data/bmmc_multiome/multiome/atac/counts.mtx",
        "data/bmmc_multiome/multiome/atac/metadata.csv",
        "data/bmmc_multiome/multiome/atac/peaks.csv",
        "data/bmmc_multiome/multiome/rna/counts.mtx",
        "data/bmmc_multiome/multiome/rna/metadata.csv",
        "data/bmmc_multiome/multiome/rna/genes.csv"
    message: "Convert BMMC multiome raw data"
    threads: 1
    shell:
      """
      python processing_code/convert_h5ad.py data/bmmc_multiome/multiome
      """

rule process_bmmc_multiome:
    input:
        "data/bmmc_multiome/multiome/atac/counts.mtx",
        "data/bmmc_multiome/multiome/atac/metadata.csv",
        "data/bmmc_multiome/multiome/atac/peaks.csv",
        "data/bmmc_multiome/multiome/rna/counts.mtx",
        "data/bmmc_multiome/multiome/rna/metadata.csv",
        "data/bmmc_multiome/multiome/rna/genes.csv"
    output: "objects/bmmc_multiome.rds"
    threads: 1
    message: "Generate BMMC multiome Seurat object"
    shell:
        """
        Rscript processing_code/bmmc_multiome.R
        """

rule process_bmmc_atac:
    input:
        "objects/bmmc_multiome.rds",
        "data/bmmc_atac/download.done",
        "data/annotations_hg38.rds"
    output: "objects/bmmc_atac.rds"
    threads: 8
    message: "Generate BMMC scATAC Seurat object"
    shell:
        """
        Rscript processing_code/bmmc_atac.R
        """

rule create_annotations:
    output:
        "data/annotations_hg38.rds",
        "data/annotations_hg19.rds"
    message: "Extracting annotations"
    threads: 1
    shell: "Rscript processing_code/get_annotations.R"

rule process_pbmc_reference:
    input: "data/pbmc_reference/download.done"
    output: "objects/pbmc_reference.rds"
    message: "Construct PBMC scRNA-seq reference"
    threads: 1
    shell: "Rscript processing_code/pbmc_reference.R"

rule process_pbmc_multiome:
    input:
        data="data/pbmc_multiome/{sample}/download.done",
        ref="objects/pbmc_reference.rds",
        annot="data/annotations_hg38.rds",
        peaks="data/pbmc_multiome_peaks.bed"
    output: "objects/{sample}.rds"
    message: "Process PBMC scATAC-seq"
    threads: 8
    shell:
        """
        Rscript processing_code/pbmc_multiome.R \
                {threads} \
                {input.annot} \
                {input.peaks} \
                {wildcards.sample} \
                {input.ref}
        """

rule process_pbmc_atac:
    input:
        data="data/pbmc_atac/{sample}/download.done",
        annot="data/annotations_hg38.rds",
        peaks="data/pbmc_multiome_peaks.bed"
    output: "objects/{sample}.rds"
    message: "Process PBMC scATAC-seq"
    threads: 8
    shell:
        """
        Rscript processing_code/pbmc_atac.R \
                {threads} \
                {input.annot} \
                {input.peaks} \
                {wildcards.sample} 
        """
