cellMarkerToGmt <- function(infile, outfile,
  uniq.column=c("tissueType", "cellName"),
  geneid.type=c("geneID", "geneSymbol")){
    # Argument Check
    uniq.column = match.arg(uniq.column)
    geneid.type = match.arg(geneid.type)
    # Load
    cm <- read.delim(infile, stringsAsFactors=FALSE)
    uniq.term <- unique(cm[, uniq.column])
    # Parse
    out <- sapply(uniq.term, function(left, cm){
      right <- cm[which(cm[, uniq.column] == left), geneid.type]
      right <- paste(right, collapse=", ")
      right <- na.omit(strsplit(right, ", ")[[1]])
      right <- right[!right %in% "NA"]
      if(geneid.type == "geneID"){
        right[grep("^[0-9]", right)]
      }else{
        right[grep("^[0-9]", right, invert = TRUE)]
      }
    }, cm=cm)
    # Save
    sink(outfile)
    for(j in seq_len(length(out))){
      cat(paste0(c(names(out)[j], "na", out[[j]]), collapse="\t"))
      cat("\n")
    }
    sink()
}
