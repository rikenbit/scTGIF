# scTGIF
Quality control and cell type annotation for unannotated single-cell RNA-Seq data

Installation of Dependent Packages
======
```r
# CRAN
install.packages("BiocManager", repos="http://cran.r-project.org")

# Bioconductor
library("BiocManager")
BiocManager::install("GSEABase", update=TRUE, ask=FALSE)
BiocManager::install("Biobase", update=TRUE, ask=FALSE)
BiocManager::install("SingleCellExperiment", update=TRUE, ask=FALSE)
BiocManager::install("BiocStyle", update=TRUE, ask=FALSE)
BiocManager::install("plotly", update=TRUE, ask=FALSE)
BiocManager::install("tagcloud", update=TRUE, ask=FALSE)
BiocManager::install("knitr", update=TRUE, ask=FALSE)
BiocManager::install("msigdbr", update=TRUE, ask=FALSE)
BiocManager::install("rmarkdown", update=TRUE, ask=FALSE)
BiocManager::install("Rcpp", update=TRUE, ask=FALSE)
BiocManager::install("S4Vectors", update=TRUE, ask=FALSE)
BiocManager::install("SummarizedExperiment", update=TRUE, ask=FALSE)
BiocManager::install("RColorBrewer", update=TRUE, ask=FALSE)
BiocManager::install("nnTensor", update=TRUE, ask=FALSE)
BiocManager::install("scales", update=TRUE, ask=FALSE)
BiocManager::install("testthat", update=TRUE, ask=FALSE)
```

Installation
======
```r
git clone https://github.com/rikenbit/scTGIF/
R CMD INSTALL scTGIF
```
or type the code below in the R console window
```r
install.packages("devtools", repos="http://cran.r-project.org")
library(devtools)
devtools::install_github("rikenbit/scTGIF")
```

# Prerequisites
- R (3.6.0 or higher)
- pandoc (1.12.3 or higher)

# License
Copyright (c) 2019 Koki Tsuyuzaki and RIKEN Bioinformatics Research Unit Released under the [Artistic License 2.0](http://www.perlfoundation.org/artistic_license_2_0).

# Authors
- Koki Tsuyuzaki
- Itoshi Nikaido