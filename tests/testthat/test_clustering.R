testthat::context("Testing 'clustering'")

testthat::test_that("clustering-se", {
  
  sacurine.se <- reading(system.file("extdata/W4M00001_Sacurine-statistics", package = "phenomis"))
  sacurine.se <- correcting(sacurine.se, figure.c = "none")
  sacurine.se <- sacurine.se[, colData(sacurine.se)[, "sampleType"] != "pool"]
  sacurine.se <- transforming(sacurine.se)
  sacurine.se <- sacurine.se[, colnames(sacurine.se) != "HU_neg_096_b2"]
  sacurine.se <- clustering(sacurine.se, clusters.vi = c(10, 10))
  
  testthat::expect_equal(as.numeric(colData(sacurine.se)["HU_neg_021", "hclust"]),
                         3,
                         tolerance = 1e-6)
  testthat::expect_equal(as.numeric(rowData(sacurine.se)["Testosterone glucuronide", "hclust"]),
                         2,
                         tolerance = 1e-6)
  
})

testthat::test_that("clustering-mae", {
  
  prometis.mae <- reading(system.file("extdata/prometis", package = "phenomis"))
  
  prometis.mae <- clustering(prometis.mae, clusters.vi = c(10, 10))
  
  testthat::expect_equal(as.numeric(colData(prometis.mae[["metabolomics"]])["s6", "hclust"]),
                         3,
                         tolerance = 1e-6)
  
})
