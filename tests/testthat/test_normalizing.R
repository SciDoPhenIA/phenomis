testthat::context("Testing 'normalizing'")

testthat::test_that("normalizing-eset", {
  
  eset <- phenomis::reading(file.path(system.file(package = "phenomis"),
                                      "extdata/W4M00001_Sacurine-statistics"))
  eset <- eset[, Biobase::sampleNames(eset) != 'HU_neg_096_b2']
  norm.eset <- phenomis::normalizing(eset, method.c = "pqn")

  testthat::expect_equal(Biobase::exprs(norm.eset)["Testosterone glucuronide", "HU_neg_020"],
                         148872.7,
                         tolerance = 1e-6)

})

testthat::test_that("normalizing-mset", {

  mset <- phenomis::reading(system.file("extdata/prometis", package = "phenomis"))

  nrom.mset <- phenomis::normalizing(mset, method.c = "pqn")
  
  testthat::expect_equal(Biobase::assayData(nrom.mset)[["metabolomics"]][["exprs"]]["v3", "s3"],
                         6.012876,
                         tolerance = 1e-6)

})
