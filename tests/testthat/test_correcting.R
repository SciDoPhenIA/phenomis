testthat::context("Testing 'correcting'")

testthat::test_that("correcting-eset", {
  
  sacurine.eset <- phenomis::reading(system.file("extdata/sacurine", package = "phenomis"))
  sacurine.eset <- phenomis::correcting(sacurine.eset)
  
  testthat::expect_equal(Biobase::exprs(sacurine.eset)["Testosterone glucuronide", "HU_neg_020"],
                         44136.83,
                         tolerance = 1e-6)
  
})

testthat::test_that("correcting-mset", {
  
  # to be done
  
})
