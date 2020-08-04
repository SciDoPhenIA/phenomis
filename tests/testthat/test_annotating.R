testthat::context("Testing 'annotating'")

testthat::test_that("annotating-eset-chebi", {
  
  sacurine.eset <- phenomis::reading(system.file("extdata/W4M00001_Sacurine-statistics", package = "phenomis"))
  sacurine.eset <- sacurine.eset[tail(Biobase::featureNames(sacurine.eset)), ]
  sacurine.eset <- phenomis::annotating(sacurine.eset, database.c = "chebi",
                                        param.ls = list(query.type = "mz",
                                                        query.col = "mass_to_charge",
                                                        ms.mode = "neg",
                                                        mz.tol = 10,
                                                        mz.tol.unit = "ppm"))
  testthat::expect_identical(Biobase::fData(sacurine.eset)["Threonic acid/Erythronic acid", "chebi.id"],
                             "15908")
  
  matchVl <- sapply(1:dim(sacurine.eset)["Features"],
                    function(featI) {
                      dbidC <- Biobase::fData(sacurine.eset)[featI, "database_identifier"]
                      if (dbidC == "") {
                        return(NA)
                      } else {
                        dbidVc <- gsub("CHEBI:", "",
                                       unlist(strsplit(dbidC, split = "|", fixed = TRUE)))
                        chebidC <- Biobase::fData(sacurine.eset)[featI, "chebi.id"]
                        chebidVc <-  unlist(strsplit(chebidC, split = "|", fixed = TRUE))
                        return(any(sapply(dbidVc, function(dbidC) {
                          dbidC %in% chebidVc
                        })))
                      }
                    })
  # table(matchVl)
  # FALSE
  #     3
  testthat::expect_identical(as.numeric(table(matchVl)), c(2, 1))
  # FALSE  TRUE 
  # 2     1
  # for the whole sacurine.eset:
  # FALSE TRUE
  #    68   14  
  
})

testthat::test_that("annotating-eset-localms", {
  
  sacurine.eset <- phenomis::reading(system.file("extdata/W4M00001_Sacurine-statistics", package = "phenomis"))
  msdbDF <- read.table(system.file("extdata/local_ms_db.tsv", package = "phenomis"),
                       header = TRUE,
                       sep = "\t",
                       stringsAsFactors = FALSE)
  sacurine.eset <- phenomis::annotating(sacurine.eset, database.c = "local.ms",
                                        param.ls = list(query.type = "mz",
                                                        query.col = "mass_to_charge",
                                                        ms.mode = "neg",
                                                        mz.tol = 5,
                                                        mz.tol.unit = "ppm",
                                                        local.ms.db = msdbDF))
  testthat::expect_identical(Biobase::fData(sacurine.eset)["Taurine", "local.ms.name"],
                             "Taurine")
  
})

testthat::test_that("annotating-eset-chebiID", {
  
  sacurine.eset <- phenomis::reading(system.file("extdata/W4M00001_Sacurine-statistics", package = "phenomis"))
  sacurine.eset <- sacurine.eset[tail(Biobase::featureNames(sacurine.eset)), ]
  sacurine.eset <- phenomis::annotating(sacurine.eset, database.c = "chebi",
                                        param.ls = list(query.type = "chebi.id", query.col = "database_identifier",
                                                        prefix = "chebiID."))
  
  testthat::expect_identical(Biobase::fData(sacurine.eset)["Tryptophan", "chebiID.formula"],
                             "C11H12N2O2")
  
})



