context("TGIFFunctions")

# Package loading
library("SingleCellExperiment")

# Data loading
data("DistalLungEpithelium")
data("pca.DistalLungEpithelium")
data("label.DistalLungEpithelium")

# SingleCellExperiment-class
sce <- SingleCellExperiment(assays = list(counts = DistalLungEpithelium))
reducedDims(sce) <- SimpleList(PCA=pca.DistalLungEpithelium)

# Generic
expect_true("genericFunction" %in% is(settingTGIF))
expect_true("genericFunction" %in% is(calcTGIF))
expect_true("genericFunction" %in% is(reportTGIF))
