testthat::context("Testing 'normalizing'")

testthat::test_that("normalizing-se", {
  
  sacurine.se <- reading(file.path(system.file(package = "phenomis"),
                                   "extdata/W4M00001_Sacurine-statistics"))
  sacurine.se <- sacurine.se[, colnames(sacurine.se) != 'HU_neg_096_b2']
  norm.se <- normalizing(sacurine.se, method.vc = "pqn")
  
  testthat::expect_equivalent(assay(norm.se)["Testosterone glucuronide", "HU_neg_020"],
                              148872.7,
                              tolerance = 1e-6)
  
})

testthat::test_that("normalizing-mae", {
  
  mae <- reading(system.file("extdata/prometis", package = "phenomis"))
  
  nrom.mae <- normalizing(mae, method.vc = "pqn")
  
  testthat::expect_equivalent(assays(nrom.mae)[["metabolomics"]]["v3", "s3"],
                              6.012876,
                              tolerance = 1e-6)
  
})


testthat::test_that("normalizing-mset", {
  
  mset <- reading(system.file("extdata/prometis", package = "phenomis"), output.c = "set")
  
  nrom.mset <- normalizing(mset, method.vc = "pqn")
  
  testthat::expect_equivalent(Biobase::assayData(nrom.mset)[["metabolomics"]][["exprs"]]["v3", "s3"],
                              6.012876,
                              tolerance = 1e-6)
  
})
