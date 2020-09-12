[![Build Status](https://travis-ci.com/kharchenkolab/sccore.svg?branch=master)](https://travis-ci.com/github/kharchenkolab/sccore)

# sccore
Core utilities for single-cell RNA-seq data analysis. Contained within are utility functions for working with differential expression (DE) matrices and count matrices, a collection of functions for manipulating and plotting data via 'ggplot2', and functions to work with cell graphs and cell embeddings. Graph-based methods include embedding kNN cell graphs into a [UMAP](https://github.com/lmcinnes/umap), collapsing vertices of each cluster in the graph, and propagating graph labels.
 
## Installation

```r
devtools::install_github('kharchenkolab/sccore')
```