#
# settingTGIF
#
setGeneric("settingTGIF", function(sce, gmt, reducedDimNames,
    assayNames="counts", grid.size=50){
    standardGeneric("settingTGIF")})

setMethod("settingTGIF",
    signature(sce="SingleCellExperiment"),
    function(sce, gmt, reducedDimNames, assayNames="counts", grid.size=50){
        userobjects <- deparse(substitute(sce))
        .settingTGIF(userobjects, gmt, reducedDimNames, assayNames,
            grid.size, sce)})

.settingTGIF <- function(userobjects, gmt, reducedDimNames, assayNames,
    grid.size, ...){
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
    twoD <- eval(parse(text=paste0("reducedDims(sce)$", reducedDimNames)))
    if(ncol(input) != nrow(twoD)){
        stop(paste0("The number of columns in assay(sce) ",
            "and the number of samples in ",
            paste0("reducedDims(sce)$", reducedDimNames),
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
            paste0("reducedDims(sce)$", reducedDimNames),
            " must be same and specified as NCBI (Entrez) Gene IDs"))
    }
    # Gene * Grid matrix
    X <- .twoDimGrid(input[common.geneid, ], twoD, grid.size)
    # Gene * Function matrix
    Y <- .convertGMT(gmt, common.geneid)
    X[is.nan(X)] <- 0
    Y[is.nan(Y)] <- 0
    # Overwrite
    metadata(sce) <- list(gmt=gmt, X=X, Y=Y,
        common.geneid=common.geneid, grid.size=grid.size,
        reducedDimNames=reducedDimNames)
    assign(userobjects, sce, envir=.GlobalEnv)
}

#
# calcTGIF
#
setGeneric("calcTGIF", function(sce, rank){
    standardGeneric("calcTGIF")})

setMethod("calcTGIF",
    signature(sce="SingleCellExperiment"),
    function(sce, rank){
        userobjects <- deparse(substitute(sce))
        .calcTGIF(userobjects, rank, sce)})

.calcTGIF <- function(userobjects, rank, ...){
    sce <- list(...)[[1]]
    X <- metadata(sce)$X
    Y <- metadata(sce)$Y
    if(min(ncol(X), ncol(Y)) < rank){
        stop(paste0("Please specify rank parameter smaller than ",
            min(ncol(X), ncol(Y))))
    }
    cat(paste0("Gene * Grid matrix X has ", nrow(X), " rows and ",
        ncol(X), " columns\n"))
    cat(paste0("Gene * Function matrix Y has ", nrow(Y), " rows and ",
        ncol(Y), " columns\n"))
    # Joint NMF
    res.sctgif <- jNMF(list(X=X, Y=Y), J=rank, algorithm="Frobenius")
    # Reconstruction error
    recerror <- res.sctgif$RecError
    relchange <- res.sctgif$RelChange
    # Overwrite
    gmt <- metadata(sce)$gmt
    common.geneid <- metadata(sce)$common.geneid
    grid.size <- metadata(sce)$grid.size
    reducedDimNames <- metadata(sce)$reducedDimNames
    metadata(sce) <- list(gmt=gmt, X=X, Y=Y, common.geneid=common.geneid,
        grid.size=grid.size, reducedDimNames=reducedDimNames,
        sctgif=res.sctgif, rank=rank, recerror=recerror,
        relchange=relchange)
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
        userobjects <- deparse(substitute(sce))
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
    grid.size <- metadata(sce)$grid.size
    reducedDimNames <- metadata(sce)$reducedDimNames
    sctgif <- metadata(sce)$sctgif
    rank <- metadata(sce)$rank
    recerror <- metadata(sce)$recerror
    relchange <- metadata(sce)$relchange

    # Low dimensional data
    twoD <- eval(parse(text=paste0("reducedDims(sce)$", reducedDimNames)))

    # Plot
    for(i in seq_len(ncol(H1))){
        par(ask=FALSE)
        png(filename=paste0(out.dir, "/figures/Grid_", i, ".png"),
            width=750, height=750)
        .plot.twoD_grid.SVD(H1, twoD, grid.size, i)
        dev.off()
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
    writeLines(.BODY3(rank, out.dir), outIdx, sep="\n")
    writeLines(.BODY4, outIdx, sep="\n")
    writeLines(.BODY5, outIdx, sep="\n")
    close(outIdx)

    # Rendering
    message("index.Rmd is compiled to index.html...")
    render(paste0(out.dir, "/index.Rmd"), quiet=TRUE)
    if (html.open) {
        browseURL(paste0(out.dir, "/index.html"))
    }
}
