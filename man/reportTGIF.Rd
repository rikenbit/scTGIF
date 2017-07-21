\name{reportTGIF}
\alias{reportTGIF}
\docType{methods}
\title{
Function for reporting the result of \code{\link{calcTGIF}} function
}
\description{
\code{\link{calcTGIF}} function calculates what kind of cellular patterns and functional patterns are contained in single-cell RNA-seq data and \code{\link{reportTGIF}} function generates report of analytic result.
}
\usage{
reportTGIF(result = NA, Comp.Matrix = NA, nDim = 10, label = NA, out.dir = NA, html.open = FALSE, title = "Result of scTGIF", author = "")
}
\arguments{
  \item{result}{
Output object of \code{\link{calcTGIF}}
}
  \item{Comp.Matrix}{
Compressed matrix performed by dimentional reduction such as PCA, Diffusion Map, t-SNE against same data used for \code{\link{calcTGIF}}.
}
  \item{nDim}{
PLS parameter meaning the number of dimention you want to observe (default: 10).
}
  \item{label}{
Cellular label which data is what kind of cell (not necessary)
}
  \item{out.dir}{
Output directory user want to save the report.
}
  \item{html.open}{
Whether html is opened when \code{\link{reportTGIF}} is finished (default: FALSE)
}
  \item{title}{
Title of report (default: "Result of scTGIF")
}
  \item{author}{
The name of user name (default: "")
}
}
\value{
Some file is generated to output directory user specified.
}
\author{
\packageAuthor{scTGIF}
}
\examples{
\dontrun{
    # Load the GSEABase for using getGmt function
    library("scTGIF")
    library("GSEABase")

    # Load the gmt file from MSigDB
    # Only "Entrez Gene ID" can be used in scTGIF
    gmt <- GSEABase::getGmt("/PATH/YOU/SAVED/THE/GMT/h.all.v6.0.entrez.gmt")

    # Test data
    load("MARS.RData")

    # Perform calcTGIF
    tgif.MARS <- calcTGIF(exprs(eSet.MARS), gmt)

    # Report the result of calcTGIF
    reportTGIF(result=tgif.MARS, Comp.Matrix=pData(eSet.MARS)$tSNE.2D, nDim=10, label=pData(eSet.MARS)$Color, out.dir=".", html.open=TRUE, title="scTGIF Report for MARS-Seq Data (t-SNE, 2D)", author="Koki Tsuyuzaki")
    }
}