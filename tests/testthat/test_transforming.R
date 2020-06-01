testthat::context("Testing 'transforming'")

testthat::test_that("transforming-eset", {

  sacurine.eset <- phenomis::reading(system.file("extdata/W4M00001_Sacurine-statistics", package = "phenomis"))
  sacurine.eset <- phenomis::correcting(sacurine.eset, figure.c = "none")
  sacurine.eset <- sacurine.eset[, Biobase::pData(sacurine.eset)[, "sampleType"] != "pool"]
  sacurine.eset <- phenomis::transforming(sacurine.eset)

  testthat::expect_equal(Biobase::exprs(sacurine.eset)["Testosterone glucuronide", "HU_neg_020"],
                         15.4297,
                         tolerance = 1e-6)

})

testthat::test_that("transforming-mset", {

  prometis.mset <- phenomis::reading(system.file("extdata/prometis", package = "phenomis"))

  prometis.mset1 <- phenomis::transforming(prometis.mset)
  testthat::expect_equal(Biobase::assayData(prometis.mset1)[["metabolomics"]][["exprs"]]["v3", "s3"],
                         2.616842,
                         tolerance = 1e-6)

  prometis.mset2 <- phenomis::transforming(prometis.mset, "sqrt")
  testthat::expect_equal(Biobase::assayData(prometis.mset2)[["metabolomics"]][["exprs"]]["v3", "s3"],
                         2.476703,
                         tolerance = 1e-6)

})
