testthat::context("Testing 'inspecting'")

testthat::test_that("inspecting-eset", {
  
  sacurine.eset <- phenomis::reading(system.file("extdata/sacurine", package = "phenomis"))
  
  sacurine.eset <- phenomis::inspecting(sacurine.eset,
                                        figure.c = "none",
                                        report.c = "none")
  
  testthat::expect_equal(Biobase::fData(sacurine.eset)["(2-methoxyethoxy)propanoic acid isomer", "pool_CV"],
                         0.3160307,
                         tolerance = 1e-6)
  
  proteo.eset <- phenomis::reading(system.file("extdata/prometis/proteomics", package = "phenomis"))
  
  proteo.eset <- phenomis::inspecting(proteo.eset,
                                      figure.c = "none",
                                      report.c = "none")
  
  testthat::expect_equal(Biobase::pData(proteo.eset)["s1", "deci_pval"],
                         0.106289011,
                         tolerance = 1e-6)
  
})

testthat::test_that("inspecting-mset", {
  
  prometis.mset <- phenomis::reading(system.file("extdata/prometis", package = "phenomis"))
  
  prometis.mset <- phenomis::inspecting(prometis.mset,
                                        figure.c = "none",
                                        report.c = "none")
  
  testthat::expect_equal(Biobase::pData(prometis.mset)[["metabolomics"]]["s1", "hotel_pval"],
                         0.3350008,
                         tolerance = 1e-6)
  
})
