.convertGMT <-
function (gmt, targetgeneid) 
{
    Y <- matrix(0, nrow = length(targetgeneid), ncol = length(gmt))
    rownames(Y) <- targetgeneid
    colnames(Y) <- names(gmt)
    lapply(1:length(gmt), function(x) {
        Y[intersect(targetgeneid, gmt[[x]]@geneIds), x] <<- 1
    })
    return(Y)
}
