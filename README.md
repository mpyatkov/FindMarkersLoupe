# FindMarkersLoupe

<!-- badges: start -->
<!-- badges: end -->

The goal of FindMarkersLoupe is to provide the same output for DE analysis like Loupe Browser.

Based on the following sources:

- [10X cellranger python sources](https://github.com/10XGenomics/cellranger/blob/master/lib/python/cellranger/analysis/diffexp.py)
- [EdgeR, exactTestBetaApprox](https://rdrr.io/bioc/edgeR/src/R/exactTestBetaApprox.R)
- [cellrangerRkit](https://github.com/hb-gitified/cellrangerRkit)

## Installation

You can install the development version of FindMarkersLoupe like so:

``` r
# install.packages("devtools")
devtools::install_github("mpyatkov/FindMarkersLoupe")
```

## Example

This is a basic example which shows you how to compute differential expression using sseq + edger which is the main approach for 10X Loupe Browesr:

``` r
library(FindMarkersLoupe)
## basic example code

library(Seurat)
library(FindMarkersLoupe)

# read Seurat object with precalculated clusters
seurat_obj <- readRDS("/path/to/rds/with/seurat/object")

# choose 2 clusters which are you going to compare (ex. "cluster_id_1" and "cluster_id_2") and run analysis
DE.df <- FindMarkersLoupe(hep01, id.1 = "cluster_id_1", id.2 = "cluster_id_2", formatted = "short" )
```

**DE.df** - is a resulted data.frame which contains gname, log2fc, intesity.1, intensity.2 and adj_pvalue columns. The data.frame will contains intermediate columns with addtional information if "formatted" option set to "full".

