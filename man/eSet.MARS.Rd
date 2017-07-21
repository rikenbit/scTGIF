\name{eSet.MARS}
\alias{eSet.MARS}
\docType{data}
\title{
ExpressionSet of MARS-Seq dataset containing four different cell types.
}
\description{
A data frame with 353 rows (genes) with following 228 columns (cells).

The data is downloaded from GEO (GSE54006) and genes having no expression are filted.

Cell type information is also downloaded from GSE54006.

The data has four different cell types; B cells (B1-B48), Dendritic cells (DC1-DC89), Monocytes (Mono1-Mono46), and Natural Killer cells (NK1-NK45).

\emph{All Gene ID is converted to Human Entrez Gene ID for applying the data to scTGIF.}
}
\usage{data("eSet.MARS")}
\source{
\url{https://www.sciencemag.org/content/343/6172/776?related-urls=yes&legid=sci;343/6172/776}

\url{http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE54006}
}
\references{
Diego Adhemar Jaitin, Ephraim Kenigsberg, Hadas Keren-Shaul, Naama Elefant, Franziska Paul, Irina Zaretsky, Alexander Mildner, Nadav Cohen, Steffen Jung, Amos Tanay, Ido Amit (2014) Massively Parallel Single-Cell RNA-Seq for Marker-Free Decomposition of Tissues into Cell Types. \emph{Science}, \bold{343 (6172)}: 776-779
}
\examples{
data(eSet.MARS)
head(exprs(eSet.MARS))
head(pData(eSet.MARS))
pairs(pData(eSet.MARS)$PCA, col=unlist(pData(eSet.MARS)$Color), pch=16)
plot(pData(eSet.MARS)$tSNE.2D, col=unlist(pData(eSet.MARS)$Color), pch=16)
pairs(pData(eSet.MARS)$tSNE.3D, col=unlist(pData(eSet.MARS)$Color), pch=16)
}
\keyword{datasets}
