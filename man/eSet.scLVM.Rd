\name{eSet.scLVM}
\alias{eSet.scLVM}
\docType{data}
\title{
ExpressionSet of scLVM dataset containing three different cell states.
}
\description{
A data frame with 7258 rows (genes) with following 270 columns (cells).

We performed sailfish against raw FASTQ data of ERP006670 and retrived gene-cell expression matrix.

Genes having no expression are filted.

The data has three different cell states; G1, S, G2/M of cell cycle stage (https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-2805/samples/).

\emph{All Gene ID is converted to Human Entrez Gene ID for applying the data to scTGIF.}
}
\usage{data("eSet.scLVM")}

\source{
\url{https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-2805}
\url{http://www.ebi.ac.uk/ena/data/view/PRJEB6989}
}
\references{
Florian Buettner, Kedar N Natarajan, F Paolo Casale, Valentina Proserpio, Antonio Scialdone, Fabian J Theis, Sarah A Teichmann, John C Marioni, Oliver Stegle (2015) Computational analysis of cell-to-cell heterogeneity in single-cell RNA-sequencing data reveals hidden subpopulations of cells. \emph{Nature Biotechnology}, \bold{33}, 155-160
}
\examples{
data(eSet.scLVM)
head(exprs(eSet.scLVM))
head(pData(eSet.scLVM))
pairs(pData(eSet.scLVM)$PCA, col=unlist(pData(eSet.scLVM)$Color), pch=16)
plot(pData(eSet.scLVM)$tSNE.2D, col=unlist(pData(eSet.scLVM)$Color), pch=16)
pairs(pData(eSet.scLVM)$tSNE.3D, col=unlist(pData(eSet.scLVM)$Color), pch=16)
}
\keyword{datasets}
