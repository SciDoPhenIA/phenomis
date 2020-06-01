testthat::context("Testing 'reducing'")

testthat::test_that("reducing-eset", {
  
  sacurine.eset <- phenomis::reading(system.file("extdata/W4M00002_Sacurine-comprehensive",
                                                 package = "phenomis"),
                                     report.c = "none")
  sacurine.eset <- phenomis::reducing(sacurine.eset,
                                      rt_tol.n = 6)
  testthat::expect_identical(Biobase::fData(sacurine.eset)["M192.0221T64", "redund_relative"],
                             "[M+13C]")
  testthat::expect_identical(Biobase::fData(sacurine.eset)["M193.0228T64", "redund_relative"],
                             "[M+18O]|[M+13C2]")
  testthat::expect_identical(Biobase::fData(sacurine.eset)["M195.0485T309", "redund_iso_add_frag"],
                             "@M194.0449T309|+13C(1.0033)")
  testthat::expect_identical(as.numeric(table(Biobase::fData(sacurine.eset)[, "redund_group"])),
                             c(3, rep(2, 7)))
  
  
})