testthat::context("Testing 'inspecting'")

testthat::test_that("inspecting-se", {
  
  sacurine.se <- reading(system.file("extdata/W4M00001_Sacurine-statistics", package = "phenomis"))
  
  sacurine.se <- inspecting(sacurine.se,
                            figure.c = "none",
                            report.c = "none")
  
  testthat::expect_equivalent(rowData(sacurine.se)["(2-methoxyethoxy)propanoic acid isomer", "pool_CV"],
                              0.3160307,
                              tolerance = 1e-6)
  
  proteo.eset <- reading(system.file("extdata/prometis/proteomics", package = "phenomis"))
  
  proteo.eset <- inspecting(proteo.eset,
                            figure.c = "none",
                            report.c = "none")
  
  testthat::expect_equivalent(colData(proteo.eset)["s1", "deci_pval"],
                              0.106289011,
                              tolerance = 1e-6)
  
})

testthat::test_that("inspecting-mae", {
  
  prometis.mae <- reading(system.file("extdata/prometis", package = "phenomis"))
  
  prometis.mae <- inspecting(prometis.mae,
                              figure.c = "none",
                              report.c = "none")
  
  testthat::expect_equivalent(colData(prometis.mae[["metabolomics"]])["s1", "hotel_pval"],
                              0.3350008,
                              tolerance = 1e-6)
  
})

testthat::test_that("inspecting-eset", {
  
  sacurine.eset <- reading(system.file("extdata/W4M00001_Sacurine-statistics", package = "phenomis"),
                           output.c = "set")
  
  sacurine.eset <- inspecting(sacurine.eset,
                            figure.c = "none",
                            report.c = "none")
  
  testthat::expect_equivalent(Biobase::fData(sacurine.eset)["(2-methoxyethoxy)propanoic acid isomer", "pool_CV"],
                              0.3160307,
                              tolerance = 1e-6)
  
  proteo.eset <- reading(system.file("extdata/prometis/proteomics", package = "phenomis"),
                         output.c = "set")
  
  proteo.eset <- inspecting(proteo.eset,
                            figure.c = "none",
                            report.c = "none")
  
  testthat::expect_equivalent(Biobase::pData(proteo.eset)["s1", "deci_pval"],
                              0.106289011,
                              tolerance = 1e-6)
  
})

testthat::test_that("inspecting-mset", {
  
  prometis.mset <- reading(system.file("extdata/prometis", package = "phenomis"), output.c = "set")
  
  prometis.mset <- inspecting(prometis.mset,
                              figure.c = "none",
                              report.c = "none")
  
  testthat::expect_equivalent(Biobase::pData(prometis.mset)[["metabolomics"]]["s1", "hotel_pval"],
                              0.3350008,
                              tolerance = 1e-6)
  
})
