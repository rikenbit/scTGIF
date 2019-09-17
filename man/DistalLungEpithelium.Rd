\name{DistalLungEpithelium}
\alias{DistalLungEpithelium}
\docType{data}
\title{
Gene expression matrix of DistalLungEpithelium dataset containing five cluster.
}
\description{
A data frame with 3397 rows (genes) with following 80 columns (cells).

The data is downloaded as supplementary information of the distal lung
epithelium paper (https://www.nature.com/articles/nature13173).

Low-expressed genes are filted.

\emph{All Gene ID is converted to Human Entrez Gene ID for
applying the data to scTGIF.}
}
\usage{data("DistalLungEpithelium")}

\source{
\url{http://www.nature.com/nbt/journal/v33/n2/full/nbt.3102.html}
}
\references{
Treutlein, B. et al. (2014) Reconstructing lineage hierarchies of the
distal lung epithelium using single-cell RNA-seq. \emph{Nature}
\bold{509}, 371-375
}
\examples{
data("DistalLungEpithelium")
}
\keyword{datasets}