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

pbmc_reference = [
        "pbmc_reference"
        ]

bmmc_reference = [
        "bmmc_reference"
        ]

all_samples = pbmc_atac + pbmc_multiome + bmmc_atac + bmmc_multiome + pbmc_reference + bmmc_reference

rule all:
  input:
      expand("objects/{sample}.rds", sample = all_samples),
      expand("annotations/{sample}.tsv.gz", sample = all_samples),
      expand("plots/{sample}.png", sample = all_samples)

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
    input: "datasets/bmmc_multiome/bmmc_multiome.txt"
    output: "data/bmmc_multiome/GSE194122_openproblems_neurips2021_multiome_BMMC_processed.h5ad"
    message: "Download BMMC multiome"
    threads: 1
    shell:
        """
        wget -i {input} -P data/bmmc_multiome
        cd data/bmmc_multiome
        gzip -d GSE194122_openproblems_neurips2021_multiome_BMMC_processed.h5ad.gz
        """

rule download_bmmc_multiome_frags:
    message: "Download BMMC multiome fragment files"
    threads: 1
    output: "data/bmmc_multiome/s1d1/atac_fragments.tsv.gz"
    shell:
      """
      aws s3 sync s3://openproblems-bio/public/post_competition/multiome ./data/bmmc_multiome
      """

rule download_bmmc_reference:
    output: "data/bmmc_reference/fullref.Rds"
    message: "Download BMMC reference object"
    threads: 1
    shell:
        """
        aws s3 cp s3://bmmc-reference/public/fullref.Rds {output} --request-payer
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
        "data/bmmc_multiome/GSE194122_openproblems_neurips2021_multiome_BMMC_processed.h5ad"
    output:
        "data/bmmc_multiome/atac/counts.mtx",
        "data/bmmc_multiome/atac/metadata.csv",
        "data/bmmc_multiome/atac/peaks.csv",
        "data/bmmc_multiome/rna/counts.mtx",
        "data/bmmc_multiome/rna/metadata.csv",
        "data/bmmc_multiome/rna/genes.csv"
    message: "Convert BMMC multiome raw data"
    threads: 1
    shell:
      """
      python processing_code/convert_h5ad.py data/bmmc_multiome
      """

rule process_bmmc_multiome:
    input:
        "data/bmmc_multiome/atac/counts.mtx",
        "data/bmmc_multiome/atac/metadata.csv",
        "data/bmmc_multiome/atac/peaks.csv",
        "data/bmmc_multiome/rna/counts.mtx",
        "data/bmmc_multiome/rna/metadata.csv",
        "data/bmmc_multiome/rna/genes.csv",
        "data/annotations_hg38.rds"
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

rule process_bmmc_reference:
    input: "data/bmmc_reference/fullref.Rds"
    output: "objects/bmmc_reference.rds"
    threads: 1
    message: "Process BMMC reference object"
    shell:
        """
        Rscript processing_code/bmmc_reference.R
        """

rule create_annotations:
    output:
        "data/annotations_hg38.rds"
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

# ---- Map ---- #

rule refmap:
    input:
        "objects/{sample}.rds",
        "objects/pbmc_reference.rds",
        "objects/bmmc_reference.rds"
    output: "annotations/{sample}.tsv"
    message: "Reference mapping {wildcards.sample}"
    threads: 1
    shell:
        """
        Rscript processing_code/refmap.R {wildcards.sample}
        gzip annotations/{wildcards.sample}.tsv
        """

# ---- Plot ---- #

rule plot:
    input: "objects/{sample}.rds", "annotations/{sample}.tsv.gz"
    output: "plots/{sample}.png"
    message: "Plotting {wildcards.sample}"
    threads: 1
    shell:
        """
        Rscript processing_code/plot.R {wildcards.sample}
        """
