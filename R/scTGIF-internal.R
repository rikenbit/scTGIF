.convertGMT <- function (gmt, targetgeneid){
    Y <- matrix(0, nrow = length(targetgeneid), ncol = length(gmt))
    rownames(Y) <- targetgeneid
    colnames(Y) <- names(gmt)
    lapply(seq_len(length(gmt)), function(x) {
        Y[intersect(targetgeneid, gmt[[x]]@geneIds), x] <<- 1
    })
    return(Y)
}

.twoDimBin <- function(sce, common.geneid, assayNames){
    input <- .importAssays(sce, assayNames)
    input <- input[common.geneid, ]
    cID <- sce@metadata$hexbin[[1]]
    out <- t(apply(input, 1, function(x, cID){
        schex:::.make_hexbin_function(x, "mean", cID)
    }, cID=cID))
    colnames(out) <- paste0("Bin", seq(ncol(out)))
    out
}

.importAssays <- function(sce, assayNames){
    if(assayNames %in% names(assays(sce))){
        assays(sce)[[assayNames]]
    }else{
        stop("Please specify the valid assayNames (cf. names(assays(sce)))")
    }
}

.plot_hexbin_pattern <- function(sce, pattern){
    drhex <- data.frame(sce@metadata$hexbin$hexbin.matrix,
        pattern=pattern)
    drhex <- as_tibble(drhex)
    # Plot
    ggplot(drhex, aes_string("x", "y", fill = pattern)) +
    geom_hex(stat = "identity") +
    theme_classic() + theme(legend.position = "bottom") +
    ggtitle("Pattern1") + scale_fill_viridis_c() +
    labs(x = "Dim1", y = "Dim2") + theme(legend.title = element_blank())
}

# Title, Author, Data
.HEADER <- function(author, title){
    HEADER <- paste0("---\ntitle: XXXXX\n",
        "author: YYYYY\ndate:",
        " \"`r Sys.time()`\"\n",
        "output:\n",
        " BiocStyle::html_document:\n",
        "  toc_float: true\n",
        "vignette: >\n ",
        "%\\VignetteIndexEntry{Vignette Title}\n ",
        "%\\VignetteEngine{knitr::rmarkdown}\n ",
        "%\\VignetteEncoding{UTF-8}\n---\n")
    sub("YYYYY", author, sub("XXXXX", title, HEADER))
}

# 1. About Dataset and scTGIF Algorithm
.BODY1 <- paste0("\n\n# About scTGIF Algorithm\n\n",
    "![](Workflow.png)\n",
    "[scTGIF](https://bioconductor.org/packages/release/",
    "bioc/html/scTGIF.html) is the R/Bioconductor package ",
    "for the functional annotation of cells detected ",
    "by single-cell RNA-Seq (scRNA-Seq).\n\n",

    "In the analysis of scRNA-Seq data, ",
    "cell type annotation is one of the most time-consuming ",
    "and subjective step, which requires the domain specific knowledge of ",
    "the cell types.\n\n",

    "Detection of differentially expressed genes (DEGs) and ",
    "the following gene functional analysis such as enrichment analysis ",
    "or geneset enrichment analysis (GSEA), which is used in the bulk-level ",
    "omics analysis might be used for the annotation step ",
    "but this approach heavily depends on the cell type label. ",
    "In the scRNA-Seq analysis, however, the cell type label will updated ",
    "with the update of the knowledge about the cells ",
    "and this approach must be repeatedly performed in each update.\n\n",

    "scTGIF is developed for reducing this try-error cycle; ",
    "using joint non-negative matrix factorization, ",
    "this tool connects the cells and the related gene function ",
    "without the cell type label ",
    "and help users to grasp what kind of cells are included in the dataset. ",
    "scTGIF also contains some quality control metrics ",
    "and provides these as the reasons to remove the low-quality cells or ",
    "artifacts such as doublets.\n\n",

    "For more detail, visit the [vignette](https://bioconductor.org",
    "/packages/devel/bioc/vignettes/scTGIF/inst/doc/scTGIF.html)",
    " of [scTGIF](https://bioconductor.org/",
    "packages/release/bioc/html/scTGIF.html)"
    )

# 2. Global statistics and plots
.BODY2 <- paste0("\n\n# Global statistics and plots\n\n",
    "The result of scTGIF is saved as a R binary file",
    " (reanalysis.RData).\n",
    "```{r, message=F, warning=F}\n", # Top
    "load(\"reanalysis.RData\")\n",
    "# SingleCellExperiment object\n",
    "sce\n",
    "# Reduced data size\n",
    "metadata(sce)$ndim\n",
    "# The result of jNMF\n",
    "str(metadata(sce)$sctgif)\n",
    "# Reconstruction Error of jNMF\n",
    "head(metadata(sce)$recerror)\n",
    "tail(metadata(sce)$recerror)\n",
    "# Relative Change of jNMF\n",
    "head(metadata(sce)$relchange)\n",
    "tail(metadata(sce)$relchange)\n",
    "# Used Gene IDs for jNMF\n",
    "head(metadata(sce)$common.geneid)\n",
    "# Two matrices used for jNMF\n",
    "dim(metadata(sce)$X)\n",
    "dim(metadata(sce)$Y)\n",
    "# The number of bins of schex (Default: 40)\n",
    "metadata(sce)$nbins\n",
    "# GSEABase object for constructing the matrix Y\n",
    "metadata(sce)$gmt\n",
    "# Gene expression matrix\n",
    "is(input)\n",
    "dim(input)\n",
    "input[seq_len(2), seq_len(2)]\n",
    "# The result of 2D dimensional reduction (e.g. t-SNE)\n",
    "is(twoD)\n",
    "dim(twoD)\n",
    "head(twoD)\n",
    "```\n\n", # Bottom

    "\n\n## Two dimensional plot of all cells\n\n",
    "```{r, message=F, warning=F}\n", # Top
    "plot_ly(x=twoD[,1], y=twoD[,2], ",
    "type = \"scatter\", ",
    "text = rownames(twoD), ",
    "mode = \"markers\")\n",
    "```\n\n", # Bottom

    "\n\n## Distribution of core matrix values\n\n",
    "```{r, message=F, warning=F}\n", # Top
    "corenames <- vapply(seq_along(corevalue), ",
    "function(x){paste0(\"Pattern\", x)}, \"\")\n",
    "plot_ly(x=seq_along(corevalue), y=corevalue, ",
    "type=\"bar\", color=names(corevalue), text=corenames, ",
    "colors = c(\"#999999\", \"#E41A1C\"))\n",
    "```\n" # Bottom
    )

# 3. Attension maps and map related gene functions
.BODY3 <- function(ndim, out.dir){
    BODY3 <- paste0("# Attension maps and map related gene functions\n",
        "```{r, message=F, warning=F, message=F, warning=F}\n", # Top
        "library(\"plotly\")\n",
        "```\n" # Bottom
        )

    paste0(BODY3,
        paste(
            vapply(seq_len(ndim), function(x){
                paste(c(
                    paste0("## Pattern ", x, "\n"),
                    ##### H1 #####
                    "![](",
                    paste0(out.dir, "/figures/Hex_", x, ".png"),
                    ")\n",
                    ##### H2 #####
                    "```{r, message=F, warning=F}\n", # Top
                    "d <- ",
                    x, "\n",
                    "value <- H2[,d]\n",
                    "term <- names(H2[,d])[order(value, decreasing=TRUE)]\n",
                    "value <- value[order(value, decreasing=TRUE)]\n",
                    "target <- seq_len(min(20, length(value)))\n",
                    "value <- value[target]\n",
                    "term <- term[target]\n\n",
                    "p <- plot_ly(x=term, y=~value,\n",
                    "type=\"bar\", color=~value, text=term,\n",
                    "colors=c(\"#4b61ba\", \"gray\", ",
                    "\"#a87963\", \"red\"))\n\n",
                    "layout(p, title=paste0(\"Pattern\", d),\n",
                    "margin = list(b = 250),\n",
                    "   xaxis=list(title=\"Term\", type = 'category',\n",
                    "     categoryorder = \"array\", categoryarray = term,\n",
                    "     tickangle = 45),\n",
                    "     yaxis=list(title=\"Value\"))\n",
                    "```\n\n" # Bottom
                ), collapse="")
            }, ""),
        collapse="")
    )
}

# 4. Session Information
.BODY4 <- paste0("\n\n# Session Information\n\n",
    "\n```{r, message=F, warning=F}\n",
    "sessionInfo()\n",
    "```\n")

# 5. License
.BODY5 <- paste0("\n\n# License\n\n",
    "Copyright (c) 2019 Koki Tsuyuzaki and Laboratory for ",
    "Bioinformatics Research, RIKEN Center for Biosystems Dynamics",
    " Reseach Released under the ",
    "[Artistic License 2.0](",
    "http://www.perlfoundation.org/artistic_license_2_0)\n")
