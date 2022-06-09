testthat::context("Testing 'reducing'")

testthat::test_that("reducing-se", {
  
  sacurine.se <- reading(system.file("extdata/W4M00002_Sacurine-comprehensive",
                                     package = "phenomis"),
                         report.c = "none")
  sacurine.se <- reducing(sacurine.se)
  testthat::expect_identical(unname(rowData(sacurine.se)["M192.0221T64", "redund_relative"]),
                             "[M+13C]")
  testthat::expect_identical(unname(rowData(sacurine.se)["M193.0228T64", "redund_relative"]),
                             "[M+18O]|[M+13C2]")
  testthat::expect_identical(unname(rowData(sacurine.se)["M195.0485T309", "redund_iso_add_frag"]),
                             "@M194.0449T309|+13C")
  testthat::expect_identical(as.numeric(table(rowData(sacurine.se)[, "redund_group"])),
                             c(3, rep(2, 7)))
  
  
})
