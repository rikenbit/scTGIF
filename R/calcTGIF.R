calcTGIF <-
function (DataMatrix = NA, gmt = NA) 
{
    if (!is.matrix(DataMatrix)) {
        stop("Please specify the data matrix!")
    }
    if (class(gmt) != "GeneSetCollection") {
        stop("Please specify the GeneSet object defined by GSEABase package!\nCheck help(calcTGIF)\n")
    }
    cat(paste0("Data matrix has ", nrow(DataMatrix), " rows and ", 
        ncol(DataMatrix), " columns\n"))
    cat(paste0("GeneSet object has ", length(gmt), " Functional Terms\n"))
    allgeneid <- unique(unlist(lapply(gmt, function(x) {
        x@geneIds
    })))
    rownames(DataMatrix) <- toupper(rownames(DataMatrix))
    targetgeneid <- intersect(rownames(DataMatrix), allgeneid)
    cat(paste0(length(targetgeneid), " / ", nrow(DataMatrix), 
        " common entrez gene ids are found!\nThese are used for following scTGIF calculation\n"))
    if (length(targetgeneid) == 0) {
        stop("None of entrez gene ids are found!\nrownames of data matrix may be not entrez gene ids\n")
    }
    DataMatrix <- DataMatrix[targetgeneid, ]
    DataMatrix <- log10(DataMatrix + 1)
    DataMatrix <- scale(DataMatrix, scale = FALSE)
    DataMatrix[is.nan(DataMatrix)] <- 0
    cat("Gene Set object is converted to matrix\n")
    Y <- .convertGMT(gmt, targetgeneid)
    Y <- scale(Y, scale = FALSE)
    Y[is.nan(Y)] <- 0
    cat("Cell-Gene * Gene-Function => Cell-Function\n")
    R <- t(DataMatrix) %*% Y
    R <- svd(R)
    R$NoCell <- ncol(DataMatrix)
    R$NoGene <- nrow(DataMatrix)
    R$NoTerm <- ncol(Y)
    R$NameCell <- colnames(DataMatrix)
    R$NameGene <- rownames(DataMatrix)
    R$NameTerm <- names(gmt)
    class(R) <- "scTGIF"
    return(R)
}
