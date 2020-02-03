#
# settingTGIF
#
setGeneric("settingTGIF", function(sce, gmt, reducedDimNames,
    assayNames="counts", nbins=40){
    standardGeneric("settingTGIF")})

setMethod("settingTGIF",
    signature(sce="SingleCellExperiment"),
    function(sce, gmt, reducedDimNames, assayNames="counts", nbins=40){
        userobjects <- deparse(substitute(sce))
        .settingTGIF(userobjects, gmt, reducedDimNames, assayNames,
            nbins, sce)})

.settingTGIF <- function(userobjects, gmt, reducedDimNames, assayNames,
    nbins, ...){
    sce <- list(...)[[1]]
    # class-check
    classCheck <- class(gmt)
    if(classCheck != "GeneSetCollection"){
        stop(paste0("Please specify the gmt as GeneSetCollection object ",
            "defined by GSEABase package"))
    }
    # Rowname-check
    rn <- rownames(assay(sce))
    if(length(rn) != length(unique(rn))){
        stop("Please specify the row names of the input matrix is unique")
    }
    # Only matrix is accepted
    if(!is.matrix(assay(sce))){
        message("The input data is coverted to matrix format by as.matrix")
        assay(sce) <- as.matrix(assay(sce))
    }
    # Import expression matrix
    input <- .importAssays(sce, assayNames)
    # Low dimensional data
    twoD <- reducedDims(sce)[[reducedDimNames]]
    if(ncol(input) != nrow(twoD)){
        stop(paste0("The number of columns in assay(sce) ",
            "and the number of samples in ",
            paste0("reducedDims(sce)[['", reducedDimNames, "']]"),
            "must be same"))
    }
    # two dimensional coordinates
    data.geneid <- unique(rownames(input))
    gmt.geneid <- unique(unlist(lapply(gmt, function(g){
        geneIds(g)
    })))
    common.geneid <- intersect(data.geneid, gmt.geneid)
    if(length(common.geneid) == 0){
        stop(paste0("Please specify the rownames of assay(sce) and ",
            paste0("reducedDims(sce)[['", reducedDimNames, "']]"),
            " must be same and specified as NCBI (Entrez) Gene IDs"))
    }
    # Setting of schex
    sce <- make_hexbin(sce, nbins=nbins,
        dimension_reduction=reducedDimNames)
    # Gene * Bin matrix
    X <- .twoDimBin(sce, common.geneid, assayNames)
    # Gene * Function matrix
    Y <- .convertGMT(gmt, common.geneid)
    X[is.nan(X)] <- 0
    Y[is.nan(Y)] <- 0
    # Overwrite
    metadata(sce)[["gmt"]] <- gmt
    metadata(sce)[["X"]] <- X
    metadata(sce)[["Y"]] <- Y
    metadata(sce)[["common.geneid"]] <- common.geneid
    metadata(sce)[["nbins"]] <- nbins
    metadata(sce)[["reducedDimNames"]] <- reducedDimNames
    metadata(sce)[["assayNames"]] <- assayNames
    assign(userobjects, sce, envir=.GlobalEnv)
}

#
# calcTGIF
#
setGeneric("calcTGIF", function(sce, ndim, verbose=FALSE, droplet=TRUE){
    standardGeneric("calcTGIF")})

setMethod("calcTGIF",
    signature(sce="SingleCellExperiment"),
    function(sce, ndim, verbose=FALSE, droplet=TRUE){
        userobjects <- deparse(substitute(sce))
        .calcTGIF(userobjects, ndim, verbose, droplet, sce)})

.calcTGIF <- function(userobjects, ndim, verbose, droplet, ...){
    sce <- list(...)[[1]]
    # Import expression matrix
    assayNames <- metadata(sce)$assayNames
    input <- .importAssays(sce, assayNames)
    # ここで鬼QC!!!!!
    # QC <- .QCmetrics(sce, input, droplet)
    X <- metadata(sce)$X
    Y <- metadata(sce)$Y
    if(min(ncol(X), ncol(Y)) < ndim){
        stop(paste0("Please specify ndim parameter smaller than ",
            min(ncol(X), ncol(Y))))
    }
    cat(paste0("Gene x Grid matrix (X) has ", nrow(X), " rows and ",
        ncol(X), " columns\n"))
    cat(paste0("Gene x Function matrix (Y) has ", nrow(Y), " rows and ",
        ncol(Y), " columns\n"))
    # Joint NMF
    res.sctgif <- jNMF(list(X=X, Y=Y), J=ndim, algorithm="Frobenius",
        verbose=verbose)
    # Reconstruction error
    recerror <- res.sctgif$RecError
    relchange <- res.sctgif$RelChange
    # Overwrite
    # metadata(sce)[["QC"]] <- QC
    metadata(sce)[["sctgif"]] <- res.sctgif
    metadata(sce)[["ndim"]] <- ndim
    metadata(sce)[["recerror"]] <- recerror
    metadata(sce)[["relchange"]] <- relchange
    assign(userobjects, sce, envir=.GlobalEnv)
}

#
# reportTGIF
#
setGeneric("reportTGIF", function(sce, out.dir=tempdir(), html.open=FALSE,
    title="The result of scTGIF",
    author="The person who runs this script",
    assayNames="counts"){
    standardGeneric("reportTGIF")})

setMethod("reportTGIF",
    signature(sce="SingleCellExperiment"),
    function(sce, out.dir=tempdir(), html.open=FALSE,
    title="The result of scTGIF",
    author="The person who runs this script",
    assayNames="counts"){
        .reportTGIF(out.dir, html.open, title, author, assayNames, sce)})

.reportTGIF <- function(out.dir, html.open, title, author, assayNames, ...){
    sce <- list(...)[[1]]
    if(is.null(metadata(sce)$sctgif)){
        stop("scTGIF did not performed properly")
    }
    # The Directory for saving the analytical result
    dir.create(paste0(out.dir, "/figures"),
        showWarnings = FALSE, recursive = TRUE)
    # File copy
    file.copy(from = system.file("extdata", "Workflow.png",
        package = "scTGIF"), to = paste0(out.dir, "/Workflow.png"),
        overwrite = TRUE)

    # Grid * Dim
    H1 <- metadata(sce)$sctgif$H[[1]]
    norm.H1 <- apply(H1, 2, function(x){norm(as.matrix(x), "F")})
    H1 <- t(t(H1) / norm.H1)
    # Function * Dim
    H2 <- metadata(sce)$sctgif$H[[2]]
    norm.H2 <- apply(H2, 2, function(x){norm(as.matrix(x), "F")})
    H2 <- t(t(H2) / norm.H2)
    # Sort
    corevalue <- norm.H1 * norm.H2
    H1 <- H1[, order(corevalue, decreasing=TRUE)]
    H2 <- H2[, order(corevalue, decreasing=TRUE)]
    corevalue <- corevalue[order(corevalue, decreasing=TRUE)]
    # Import expression matrix
    input <- .importAssays(sce, assayNames)
    # From metadata
    gmt <- metadata(sce)$gmt
    X <- metadata(sce)$X
    Y <- metadata(sce)$Y
    common.geneid <- metadata(sce)$common.geneid
    nbins <- metadata(sce)$nbins
    reducedDimNames <- metadata(sce)$reducedDimNames
    sctgif <- metadata(sce)$sctgif
    ndim <- metadata(sce)$ndim
    recerror <- metadata(sce)$recerror
    relchange <- metadata(sce)$relchange

    # Low dimensional data
    twoD <- reducedDims(sce)[[reducedDimNames]]

    # Plot
    for(i in seq_len(ncol(H1))){
        filename <- paste0(out.dir, "/figures/Hex_", i, ".png")
        g <- .plot_hexbin_pattern(sce, H1[, i])
        ggsave(filename, plot=g, dpi=200, width=6, height=6.5)
        system(paste0("ls ", filename))
    }

    # Save the result of scTGIF
    save(sce, input, twoD, H1, H2, corevalue,
        file=paste0(out.dir, "/reanalysis.RData"))

    # Output index.html
    message("index.Rmd is created...")
    outIdx <- file(paste0(out.dir, "/index.Rmd"), "w")
    writeLines(.HEADER(author, title), outIdx, sep="\n")
    writeLines(.BODY1, outIdx, sep="\n")
    writeLines(.BODY2, outIdx, sep="\n")
    writeLines(.BODY3(ndim, out.dir), outIdx, sep="\n")
    writeLines(.BODY4, outIdx, sep="\n")
    writeLines(.BODY5, outIdx, sep="\n")
    close(outIdx)

    # Rendering
    message("index.Rmd is compiled to index.html...")
    render(paste0(out.dir, "/index.Rmd"), quiet=TRUE)
    # HTML Open
    message(paste0("################################################\n",
        "Data files are saved in\n",
        out.dir, "\n################################################\n"))
    if (html.open) {
        browseURL(paste0(out.dir, "/index.html"))
    }
}
