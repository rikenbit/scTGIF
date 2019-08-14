context("DistalLungEpithelium")

data("DistalLungEpithelium")
expect_true(nrow(DistalLungEpithelium) == 3397)
expect_true(ncol(DistalLungEpithelium) == 80)

data("pca.DistalLungEpithelium")
expect_true(nrow(pca.DistalLungEpithelium) == 80)
expect_true(ncol(pca.DistalLungEpithelium) == 2)

data("label.DistalLungEpithelium")
expect_true(length(label.DistalLungEpithelium) == 80)
