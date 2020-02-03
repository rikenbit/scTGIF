\name{convertRowID}
\alias{convertRowID}
\docType{methods}
\title{
A function to change the row names of a matrix.
}
\description{
To avoid to specify the duplicated row names against matrix, multiple aggregation rules are implemented.
}
\usage{
convertRowID(input, rowID, LtoR,
    aggr.rule=c("sum", "mean", "large.mean", "large.var", "large.cv2"))
}
\arguments{
\item{input}{
A matrix filled with number (n * m).
}
\item{rowID}{
A vector to specify the row names of input (length: n).
}
\item{LtoR}{
A corresponding table to covert the row names of input as different type of IDs.
(Left: current row names -> Right: new row names)
}
\item{aggr.rule}{
The aggregation rule to change the row names of input and collapse/select the values, if the row names changed by LtoR are duplicated.

sum: Collapses multiple row vectors by summation.
mean: Collapses multiple row vectors by mean.
large.mean: Select a vector having the largest mean in the duplicated vectors.
large.var: Select a vector having the largest variance in the duplicated vectors.
large.cv2: Select a vector having the largest CV2 in the duplicated vectors.
}
}
\value{
\item{output}{
A matrix, in which the row names are changed, according to the aggregation rule user specified.
}
\item{ctable}{
The corresponding table explaining the relationship between previous row names and changed row names of input.
}
}
\author{
\packageAuthor{scTGIF}
}
\examples{
    input <- matrix(1:20, nrow=4, ncol=5)
    rowID <- c("A", "B", "C", "D")
    LtoR <- rbind(
      c("A", "3"),
      c("B", "2"),
      c("C", "4"),
      c("D", "7"))
    (input2 <- convertRowID(input, rowID, LtoR, "sum"))
    (input3 <- convertRowID(input, rowID, LtoR, "mean"))
    (input4 <- convertRowID(input, rowID, LtoR, "large.mean"))
    (input5 <- convertRowID(input, rowID, LtoR, "large.var"))
    (input6 <- convertRowID(input, rowID, LtoR, "large.cv2"))
}

\keyword{methods}