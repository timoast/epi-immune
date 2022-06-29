# Immune cell public epigenomic data

Processing pipeline for public single-cell epigenomic immune datasets.

## Installing dependencies

All the required dependencies needed to run the workflow can 
be installed automatically by creating a new conda environment.

First ensure that [conda](https://docs.conda.io/en/latest/miniconda.html)
or [mamba](https://github.com/mamba-org/mamba) is installed and available.

To create a new environment with the dependencies installed, run:

```
# using conda
conda env create -f environment.yaml
```

```
# using mamba
mamba env create -f environment.yaml
```

## Running the workflow

This workflow involves downloading data from AWS using the AWS
command line tools. To enable the download, you will first need
to create an AWS account and set up the AWS command line tools by
running `aws configure`. *Note that some of the data downloaded
may incur charges from AWS*.

To run the Snakemake workflow, first activate the conda environment
containing the required dependencies:

```
conda activate immune
```

Next, run `snakemake` with the desired options. Setting the `-j` parameter
controls the maximum number of cores used by the workflow:

```
snakemake -j 24
```

See the [snakemake](https://snakemake.readthedocs.io/en/stable/)
documentation for a complete list of available options.

## Datasets

### PBMC

#### pbmc_atac_500

![](./plots/pbmc_atac_500.png)

#### pbmc_atac_1k

![](./plots/pbmc_atac_1k.png)

#### pbmc_atac_5k

![](./plots/pbmc_atac_5k.png)

#### pbmc_atac_10k

![](./plots/pbmc_atac_10k.png)

#### pbmc_atac_10k_chromium

![](./plots/pbmc_atac_10k_chromium.png)

#### pbmc_atac_10k_chromiumX

![](./plots/pbmc_atac_10k_chromiumX.png)

#### pbmc_multiome_3k_sorted

![](./plots/pbmc_multiome_3k_sorted.png)

#### pbmc_multiome_3k_unsorted

![](./plots/pbmc_multiome_3k_unsorted.png)

#### pbmc_multiome_10k_sorted

![](./plots/pbmc_multiome_10k_sorted.png)

#### pbmc_multiome_10k_unsorted

![](./plots/pbmc_multiome_10k_unsorted.png)

#### pbmc_multiome_10k_chromium

![](./plots/pbmc_multiome_10k_chromium.png)

#### pbmc_multiome_10k_chromiumX

![](./plots/pbmc_multiome_10k_chromiumX.png)

#### pbmc_reference

![](./plots/pbmc_reference.png)

### BMMC

#### bmmc_atac

![](./plots/bmmc_atac.png)

#### bmmc_multiome

![](./plots/bmmc_multiome.png)

