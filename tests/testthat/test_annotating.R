testthat::context("Testing 'annotating'")

testthat::test_that("annotating-se-chebi", {
  
  sacurine.se <- reading(system.file("extdata/W4M00001_Sacurine-statistics", package = "phenomis"))
  sacurine.se <- sacurine.se[tail(rownames(sacurine.se)), ]
  sacurine.se <- annotating(sacurine.se, database.c = "chebi",
                            param.ls = list(query.type = "mz",
                                            query.col = "mass_to_charge",
                                            ms.mode = "neg",
                                            mz.tol = 10,
                                            mz.tol.unit = "ppm"))
  testthat::expect_identical(rowData(sacurine.se)["Threonic acid/Erythronic acid", "chebi.id"],
                             "15908")
  
  matchVl <- sapply(1:nrow(sacurine.se),
                    function(featI) {
                      dbidC <- rowData(sacurine.se)[featI, "database_identifier"]
                      if (dbidC == "") {
                        return(NA)
                      } else {
                        dbidVc <- gsub("CHEBI:", "",
                                       unlist(strsplit(dbidC, split = "|", fixed = TRUE)))
                        chebidC <- rowData(sacurine.se)[featI, "chebi.id"]
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
  # for the whole sacurine.se:
  # FALSE TRUE
  #    68   14  
  
})

testthat::test_that("annotating-se-localms", {
  
  sacurine.se <- reading(system.file("extdata/W4M00001_Sacurine-statistics", package = "phenomis"))
  msdbDF <- read.table(system.file("extdata/local_ms_db.tsv", package = "phenomis"),
                       header = TRUE,
                       sep = "\t",
                       stringsAsFactors = FALSE)
  sacurine.se <- annotating(sacurine.se, database.c = "local.ms",
                            param.ls = list(query.type = "mz",
                                            query.col = "mass_to_charge",
                                            ms.mode = "neg",
                                            mz.tol = 5,
                                            mz.tol.unit = "ppm",
                                            local.ms.db = msdbDF))
  testthat::expect_identical(rowData(sacurine.se)["Taurine", "local.ms.name"],
                             "Taurine")
  
})

testthat::test_that("annotating-se-chebiID", {
  
  sacurine.se <- reading(system.file("extdata/W4M00001_Sacurine-statistics", package = "phenomis"))
  sacurine.se <- sacurine.se[tail(rownames(sacurine.se)), ]
  sacurine.se <- annotating(sacurine.se, database.c = "chebi",
                            param.ls = list(query.type = "chebi.id", query.col = "database_identifier",
                                            prefix = "chebiID."))
  
  testthat::expect_identical(rowData(sacurine.se)["Tryptophan", "chebiID.formula"],
                             "C11H12N2O2")
  
})

testthat::test_that("annotating-mae", {
  
  sacurine.se <- reading(system.file("extdata/W4M00001_Sacurine-statistics", package = "phenomis"))
  sacurine.se <- sacurine.se[tail(rownames(sacurine.se)), ]
  map.ls <- list(sacurine = data.frame(primary = colnames(sacurine.se),
                                       colname = colnames(sacurine.se)))
  map.df <- MultiAssayExperiment::listToMap(map.ls)
  sac.mae <- MultiAssayExperiment::MultiAssayExperiment(experiments = list(sacurine = sacurine.se),
                                                        colData = colData(sacurine.se),
                                                        sampleMap = map.df)
  sac.mae <- annotating(sac.mae, database.c = "chebi",
                        param.ls = list(query.type = "chebi.id", query.col = "database_identifier",
                                        prefix = "chebiID."))
  
  testthat::expect_identical(rowData(sac.mae[["sacurine"]])["Tryptophan", "chebiID.formula"],
                             "C11H12N2O2")

  
})




