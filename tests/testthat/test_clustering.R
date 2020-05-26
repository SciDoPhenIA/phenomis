testthat::context("Testing 'clustering'")

testthat::test_that("clustering-eset", {
  
  sacurine.eset <- phenomis::reading(system.file("extdata/sacurine", package = "phenomis"))
  sacurine.eset <- phenomis::correcting(sacurine.eset, figure.c = "none")
  sacurine.eset <- sacurine.eset[, Biobase::pData(sacurine.eset)[, "sampleType"] != "pool"]
  sacurine.eset <- phenomis::transforming(sacurine.eset)
  sacurine.eset <- sacurine.eset[, Biobase::sampleNames(sacurine.eset) != "HU_neg_096_b2"]
  sacurine.eset <- phenomis::clustering(sacurine.eset, clusters.vi = c(10, 10))
  
  testthat::expect_equal(Biobase::pData(sacurine.eset)["HU_neg_021", "hclust"],
                         3,
                         tolerance = 1e-6)
  testthat::expect_equal(Biobase::fData(sacurine.eset)["Testosterone glucuronide", "hclust"],
                         2,
                         tolerance = 1e-6)
  
})

testthat::test_that("clustering-mset", {
  
  prometis.mset <- phenomis::reading(system.file("extdata/prometis", package = "phenomis"))
  
  prometis.mset <- phenomis::clustering(prometis.mset, clusters.vi = c(10, 10))
  
  testthat::expect_equal(Biobase::pData(prometis.mset)[["metabolomics"]]["s6", "hclust"],
                         3,
                         tolerance = 1e-6)
  
})
