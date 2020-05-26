testthat::context("Testing 'filtering'")

testthat::test_that("filtering-eset", {

  sacurine.eset <- phenomis::reading(system.file("extdata/sacurine", package = "phenomis"))
  
  exprs.mn <- Biobase::exprs(sacurine.eset)
  exprs.mn[exprs.mn < 1e5] <- NA
  Biobase::exprs(sacurine.eset) <- exprs.mn

  sacurine.eset <- phenomis::filtering(sacurine.eset)
  
  expect_dims.mn <- matrix(c(85, 200), nrow = 2, ncol = 1,
                           dimnames = list(c("Features", "Samples"), "exprs"))
  
  testthat::expect_equal(Biobase::dims(sacurine.eset),
                         expect_dims.mn)

})

testthat::test_that("filtering-eset-gender", {
  
  sacurine.eset <- phenomis::reading(system.file("extdata/sacurine", package = "phenomis"))
  
  exprs.mn <- Biobase::exprs(sacurine.eset)
  exprs.mn[exprs.mn < 1e5] <- NA
  Biobase::exprs(sacurine.eset) <- exprs.mn
  
  sacurine.eset <- phenomis::filtering(sacurine.eset, class.c = "gender")
  
  expect_dims.mn <- matrix(c(90, 196), nrow = 2, ncol = 1,
                           dimnames = list(c("Features", "Samples"), "exprs"))
  
  testthat::expect_equal(Biobase::dims(sacurine.eset),
                         expect_dims.mn)
  
})

testthat::test_that("filtering-mset", {

   prometis.mset <- phenomis::reading(system.file("extdata/prometis", package = "phenomis"))
   
   for (set.c in names(prometis.mset)) {
     eset <- prometis.mset[[set.c]]
     exprs.mn <- Biobase::exprs(eset)
     exprs.mn[exprs.mn < quantile(c(exprs.mn), 0.2)] <- NA
     Biobase::exprs(eset) <- exprs.mn
     prometis.mset <- MultiDataSet::add_eset(prometis.mset, eset, dataset.type = set.c,
                                             GRanges = NA, overwrite = TRUE, warnings = FALSE)
   }
   
   prometis.mset <- phenomis::filtering(prometis.mset)
   
   expect_dims.ls <- list(metabolomics = c(Features = 74, Samples = 24),
                          proteomics = c(Features = 77, Samples = 28))
   
   testthat::expect_equal(Biobase::dims(prometis.mset),
                          expect_dims.ls)

})

testthat::test_that("filtering-mset-gene", {
  
  prometis.mset <- phenomis::reading(system.file("extdata/prometis", package = "phenomis"))
  
  for (set.c in names(prometis.mset)) {
    eset <- prometis.mset[[set.c]]
    exprs.mn <- Biobase::exprs(eset)
    exprs.mn[exprs.mn < quantile(c(exprs.mn), 0.2)] <- NA
    Biobase::exprs(eset) <- exprs.mn
    prometis.mset <- MultiDataSet::add_eset(prometis.mset, eset, dataset.type = set.c,
                                            GRanges = NA, overwrite = TRUE, warnings = FALSE)
  }
  
  prometis.mset <- phenomis::filtering(prometis.mset, class.c = "gene")
  
  expect_dims.ls <- list(metabolomics = c(Features = 83, Samples = 23),
                         proteomics = c(Features = 79, Samples = 28))
  
  testthat::expect_equal(Biobase::dims(prometis.mset),
                         expect_dims.ls)
  
})
