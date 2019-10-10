# scTGIF
Quality control and cell type annotation for unannotated single-cell RNA-Seq data

Installation of Dependent Packages
======
```r
# CRAN
install.packages("BiocManager", repos="http://cran.r-project.org")

# Bioconductor
library("BiocManager")
BiocManager::install(
    c("GSEABase", "Biobase", "SingleCellExperiment", "BiocStyle", "plotly",
    "tagcloud", "knitr", "msigdbr", "rmarkdown", "Rcpp", "S4Vectors",
    "SummarizedExperiment", "RColorBrewer", "nnTensor", "scales", "testthat"),
    update=TRUE, ask=FALSE)
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