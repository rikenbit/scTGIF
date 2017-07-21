\name{eSet.RaceID}
\alias{eSet.RaceID}
\docType{data}
\title{
ExpressionSet of RaceID dataset containing twenty-three cluster.
}
\description{
A data frame with 3724 rows (genes) with following 238 columns (cells).

The data is downloaded as supplementary information of RaceID paper (http://www.nature.com/nbt/journal/v33/n2/full/nbt.3102.html).

Genes having no expression are filted.

The data has twenty-three cluster, which defined by performing RaceID script (http://www.nature.com/nature/journal/v525/n7568/full/nature14966.html).

\emph{All Gene ID is converted to Human Entrez Gene ID for applying the data to scTGIF.}
}
\usage{data("eSet.RaceID")}

\source{
\url{http://www.nature.com/nbt/journal/v33/n2/full/nbt.3102.html}
}
\references{
Dominic Grun, Anna Lyubimova, Lennart Kester, Kay Wiebrands,  Onur Basak, Nobuo Sasaki, Hans Clevers, Alexander van Oudenaarden (2015) Single-cell messenger RNA sequencing reveals rare intestinal cell types. \emph{Nature}, \bold{525}: 251-255
}
\examples{
data(eSet.RaceID)
head(exprs(eSet.RaceID))
head(pData(eSet.RaceID))
pairs(pData(eSet.RaceID)$PCA, col=unlist(pData(eSet.RaceID)$Color), pch=16)
plot(pData(eSet.RaceID)$tSNE.2D, col=unlist(pData(eSet.RaceID)$Color), pch=16)
pairs(pData(eSet.RaceID)$tSNE.3D, col=unlist(pData(eSet.RaceID)$Color), pch=16)
}
\keyword{datasets}
