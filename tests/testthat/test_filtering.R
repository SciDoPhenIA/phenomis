testthat::context("Testing 'filtering'")

testthat::test_that("filtering-se", {

  sacurine.se <- reading(system.file("extdata/W4M00001_Sacurine-statistics", package = "phenomis"))
  
  assay.mn <- assay(sacurine.se)
  assay.mn[assay.mn < 1e5] <- NA
  assay(sacurine.se) <- assay.mn

  sacurine.se <- filtering(sacurine.se)
  
  expect_dims.vi <- c(85, 200)
  
  testthat::expect_equal(dim(sacurine.se),
                         expect_dims.vi)

})

testthat::test_that("filtering-se-gender", {
  
  sacurine.se <- reading(system.file("extdata/W4M00001_Sacurine-statistics", package = "phenomis"))
  
  assay.mn <- assay(sacurine.se)
  assay.mn[assay.mn < 1e5] <- NA
  assay(sacurine.se) <- assay.mn
  
  sacurine.se <- filtering(sacurine.se, class.c = "gender")
  
  expect_dims.vi <- c(90, 196)
  
  testthat::expect_equal(dim(sacurine.se),
                         expect_dims.vi)
  
})

testthat::test_that("filtering-mae", {
  
  prometis.mae <- reading(system.file("extdata/prometis", package = "phenomis"))
  
  for (set.c in names(prometis.mae)) {
    set.se <- prometis.mae[[set.c]]
    set.mn <- assay(set.se)
    set.mn[set.mn < quantile(c(set.mn), 0.2)] <- NA
    assay(set.se) <- set.mn
    prometis.mae[[set.c]] <- set.se
  }
  
  prometis.mae <- filtering(prometis.mae)
  
  expect_dims.mn <- matrix(c(74, 24, 77, 28),
                           nrow = 2,
                           dimnames = list(NULL, c("metabolomics", "proteomics")))
  
  testthat::expect_equal(sapply(names(prometis.mae), function(set.c) dim(prometis.mae[[set.c]])),
                         expect_dims.mn)
  
})

testthat::test_that("filtering-mset", {

   prometis.mset <- reading(system.file("extdata/prometis", package = "phenomis"), output.c = "set")
   
   for (set.c in names(prometis.mset)) {
     eset <- prometis.mset[[set.c]]
     exprs.mn <- Biobase::exprs(eset)
     exprs.mn[exprs.mn < quantile(c(exprs.mn), 0.2)] <- NA
     Biobase::exprs(eset) <- exprs.mn
     prometis.mset <- MultiDataSet::add_eset(prometis.mset, eset, dataset.type = set.c,
                                             GRanges = NA, overwrite = TRUE, warnings = FALSE)
   }
   
   prometis.mset <- filtering(prometis.mset)
   
   expect_dims.ls <- list(metabolomics = c(Features = 74, Samples = 24),
                          proteomics = c(Features = 77, Samples = 28))
   
   testthat::expect_equal(Biobase::dims(prometis.mset),
                          expect_dims.ls)

})

testthat::test_that("filtering-mset-gene", {
  
  prometis.mset <- reading(system.file("extdata/prometis", package = "phenomis"), output.c = "set")
  
  for (set.c in names(prometis.mset)) {
    eset <- prometis.mset[[set.c]]
    exprs.mn <- Biobase::exprs(eset)
    exprs.mn[exprs.mn < quantile(c(exprs.mn), 0.2)] <- NA
    Biobase::exprs(eset) <- exprs.mn
    prometis.mset <- MultiDataSet::add_eset(prometis.mset, eset, dataset.type = set.c,
                                            GRanges = NA, overwrite = TRUE, warnings = FALSE)
  }
  
  prometis.mset <- filtering(prometis.mset, class.c = "gene")
  
  expect_dims.ls <- list(metabolomics = c(Features = 83, Samples = 23),
                         proteomics = c(Features = 79, Samples = 28))
  
  testthat::expect_equal(Biobase::dims(prometis.mset),
                         expect_dims.ls)
  
})
