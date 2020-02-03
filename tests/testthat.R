library("testthat")
library("scTGIF")

options(testthat.use_colours = FALSE)

test_file("testthat/test_DistalLungEpithelium.R")
test_file("testthat/test_TGIFFunctions.R")
test_file("testthat/test_convertRowID.R")
test_file("testthat/test_cellMarkerToGmt.R")
