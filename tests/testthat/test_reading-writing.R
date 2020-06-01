testthat::context("Testing 'read-write'")

testthat::test_that(".reading", {
  
  sacdir.c <- system.file("extdata/W4M00001_Sacurine-statistics", package = "phenomis")
  
  sacSet1 <- phenomis:::.reading(sacdir.c)
  
  testthat::expect_true(class(sacSet1) == "ExpressionSet")
  
  testthat::expect_equal(Biobase::exprs(sacSet1)[1, 1],
                         477491,
                         tolerance = 1e-3)
  # alternatively
  sacSet2 <- phenomis:::.reading(NA,
                                 file.path(sacdir.c, "Galaxy1_dataMatrix.tabular"),
                                 file.path(sacdir.c, "Galaxy2_sampleMetadata.tabular"),
                                 file.path(sacdir.c, "Galaxy3_variableMetadata.tabular"))
  
  testthat::expect_true(class(sacSet2) == "ExpressionSet")
  
  testthat::expect_equal(Biobase::exprs(sacSet2)[1, 1],
                         477491,
                         tolerance = 1e-3)
  
  
  testthat::expect_error(phenomis:::.reading(NA,
                                             file.path(sacdir.c, "Galaxy1_dataMatrix.tsv"),
                                             file.path(sacdir.c, "Galaxy2_sampleMetadata.tabular"),
                                             file.path(sacdir.c, "Galaxy3_variableMetadata.tabular")))
  
})

testthat::test_that("reading_ExpressionSet", {
  
  sacdir.c <- system.file("extdata/W4M00001_Sacurine-statistics", package = "phenomis")
  
  sacSet2 <- phenomis::reading(sacdir.c)
  
  testthat::expect_true(class(sacSet2) == "ExpressionSet")
  
  testthat::expect_equal(Biobase::exprs(sacSet2)[1, 1],
                         477491,
                         tolerance = 1e-3)
  # alternatively
  sacSet3 <- phenomis::reading(NA,
                               files.ls = list(dataMatrix = file.path(sacdir.c, "Galaxy1_dataMatrix.tabular"),
                                               sampleMetadata = file.path(sacdir.c, "Galaxy2_sampleMetadata.tabular"),
                                               variableMetadata = file.path(sacdir.c, "Galaxy3_variableMetadata.tabular")))
  
  testthat::expect_true(class(sacSet3) == "ExpressionSet")
  
  testthat::expect_equal(Biobase::exprs(sacSet3)[1, 1],
                         477491,
                         tolerance = 1e-3)
  
  
  testthat::expect_error(phenomis::reading(NA,
                                           files.ls = list(dataMatrix = file.path(sacdir.c, "Galaxy1_dataMatrix.tsv"),
                                                           sampleMetadata = file.path(sacdir.c, "Galaxy2_sampleMetadata.tabular"),
                                                           variableMetadata = file.path(sacdir.c, "Galaxy3_variableMetadata.tabular"))))
  
})


testthat::test_that("reading_MultiDataSet", {
  
  prometis_dir.c <- system.file("extdata/prometis", package = "phenomis")
  
  prometis.mset1 <- phenomis::reading(prometis_dir.c)
  
  testthat::expect_true(class(prometis.mset1) == "MultiDataSet")
  
  testthat::expect_identical(names(prometis.mset1), c("metabolomics", "proteomics"))
  
  testthat::expect_equal(Biobase::exprs(prometis.mset1[["metabolomics"]])[1, 1],
                         5.576,
                         tolerance = 1e-3)
  
  
  prometis.mset2 <- phenomis::reading(NA,
                                      files.ls = list(metabolomics = list(dataMatrix = file.path(prometis_dir.c, "metabolomics", "dataMatrix.tsv"),
                                                                          sampleMetadata = file.path(prometis_dir.c, "metabolomics", "sampleMetadata.tsv"),
                                                                          variableMetadata = file.path(prometis_dir.c, "metabolomics", "variableMetadata.tsv")),
                                                      proteomics = list(dataMatrix = file.path(prometis_dir.c, "proteomics", "dataMatrix.tsv"),
                                                                        sampleMetadata = file.path(prometis_dir.c, "proteomics", "sampleMetadata.tsv"),
                                                                        variableMetadata = file.path(prometis_dir.c, "proteomics", "variableMetadata.tsv"))))
  
  testthat::expect_true(class(prometis.mset2) == "MultiDataSet")
  
  testthat::expect_identical(names(prometis.mset2), c("metabolomics", "proteomics"))
  
  testthat::expect_equal(Biobase::exprs(prometis.mset2[["metabolomics"]])[1, 1],
                         5.576,
                         tolerance = 1e-3)
  
  # testthat::expect_identical(colnames(Biobase::pData(prometis.mset2[["metabolomics"]])), c("gene", "id"))
  
  
  
  metMset1 <- phenomis::reading(prometis_dir.c, subsets.vc = "metabolomics")
  
  testthat::expect_true(class(metMset1) == "MultiDataSet")
  
  testthat::expect_identical(names(metMset1), "metabolomics")
  
  testthat::expect_equal(Biobase::exprs(metMset1[["metabolomics"]])[1, 1],
                         5.576,
                         tolerance = 1e-3)
  
  
  metMset2 <- phenomis::reading(NA,
                                files.ls = list(metabolomics = list(dataMatrix = file.path(prometis_dir.c, "metabolomics", "dataMatrix.tsv"),
                                                                    sampleMetadata = file.path(prometis_dir.c, "metabolomics", "sampleMetadata.tsv"),
                                                                    variableMetadata = file.path(prometis_dir.c, "metabolomics", "variableMetadata.tsv")),
                                                proteomics = list(dataMatrix = file.path(prometis_dir.c, "proteomics", "dataMatrix.tsv"),
                                                                  sampleMetadata = file.path(prometis_dir.c, "proteomics", "sampleMetadata.tsv"),
                                                                  variableMetadata = file.path(prometis_dir.c, "proteomics", "variableMetadata.tsv"))),
                                subsets.vc = "metabolomics")
  
  testthat::expect_true(class(metMset2) == "MultiDataSet")
  
  testthat::expect_identical(names(metMset2), "metabolomics")
  
  testthat::expect_equal(Biobase::exprs(metMset2[["metabolomics"]])[1, 1],
                         5.576,
                         tolerance = 1e-3)
  
  
  testthat::expect_error(phenomis::reading(NA,
                                           list(metabolomics = list(dataMatrix = file.path(prometis_dir.c, "metabolomics", "dataMatrix.tsv_XXX"),
                                                                    sampleMetadata = file.path(prometis_dir.c, "metabolomics", "sampleMetadata.tsv"),
                                                                    variableMetadata = file.path(prometis_dir.c, "metabolomics", "variableMetadata.tsv")),
                                                proteomics = list(dataMatrix = file.path(prometis_dir.c, "proteomics", "dataMatrix.tsv_XXX"),
                                                                  sampleMetadata = file.path(prometis_dir.c, "proteomics", "sampleMetadata.tsv"),
                                                                  variableMetadata = file.path(prometis_dir.c, "proteomics", "variableMetadata.tsv")))))
  
  prometis.mset3 <- phenomis::reading(NA,
                                      list(metabolomics = list(dataMatrix = file.path(prometis_dir.c, "metabolomics", "dataMatrix.tsv"),
                                                               sampleMetadata = file.path(prometis_dir.c, "metabolomics", "sampleMetadata.tsv_XXX"),
                                                               variableMetadata = file.path(prometis_dir.c, "metabolomics", "variableMetadata.tsv")),
                                           proteomics = list(dataMatrix = file.path(prometis_dir.c, "proteomics", "dataMatrix.tsv"),
                                                             sampleMetadata = file.path(prometis_dir.c, "proteomics", "sampleMetadata.tsv"),
                                                             variableMetadata = file.path(prometis_dir.c, "proteomics", "variableMetadata.tsv"))))
  
  testthat::expect_true(colnames(Biobase::pData(prometis.mset3[["metabolomics"]])) == "id")
  
  
})

testthat::test_that("writing_MultiDataSet", {
  
  prometis_dir.c <- system.file("extdata/prometis", package = "phenomis")
  
  prometis.mset4 <- reading(prometis_dir.c)
  
  testthat::expect_error(phenomis::writing(prometis.mset4,
                                           dir.c = NA,
                                           files.ls = list(metabolomics = list(dataMatrix = NA,
                                                                               sampleMetadata = file.path(getwd(), "metabolomics_sampleMetadata.tsv"),
                                                                               variableMetadata = file.path(getwd(), "metabolomics_variableMetadata.tsv")),
                                                           proteomics = list(dataMatrix = file.path(getwd(), "proteomics_dataMatrix.tsv"),
                                                                             sampleMetadata = file.path(getwd(), "proteomics_sampleMetadata.tsv"),
                                                                             variableMetadata = file.path(getwd(), "proteomics_variableMetadata.tsv")))))
  
})
