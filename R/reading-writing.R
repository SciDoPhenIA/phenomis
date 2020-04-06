#' reading
#'
#' Reading dataset(s) in the 3 tables 'dataMatrix.tsv', sampleMetadata.tsv' and
#' 'variableMetadata.tsv' format. In case of a single dataset (3 tables in the
#' specified directory), an ExpressionSet instance is returned. In case of a
#' multiple dataset (several subfolders containing 3 tables), a MultiDataSet
#' instance is created.
#'
#' @param dir.c Character: directory containing the 3 .tsv files (single
#' dataset), or containing several subdirectories with 3 .tsv files
#' (multiple datasets)
#' @param subsets.vc Vector of characters: specifying a subset of the
#' subdirectories to be included in the MultiDataSet (by default, all
#' subdirectories containing the 3 tables will be considered as datasets)
#' @param files.ls List: if dir.c is set to NA, the full names of the
#' individual files can be provided; in case of an ExpressionSet, the
#' names of the list must be 'dataMatrix.tsvC', 'sampleMetadata.tsvC',
#' and 'variableMetadata.tsvC' with the corresponding file full names;
#' in case of a MultiDataSet, the list must consists of one such sublist
#' per dataset
#' @param report.c Character: File name for the printed results (call to
#' 'sink()'); if NA (default), messages will be printed on the screen; if NULL,
#' no verbose will be generated
#' @return \code{ExpressionSet} (one dataset) or \code{MultiDataSet} (multiple
#' datasets) instance containing the dataset(s)
#' @examples
#' data_dir.c <- system.file("extdata", package="phenomis")
#' ## 1) Single set
#' sacurine_dir.c <- file.path(data_dir.c, "sacurine")
#' sacurine.eset <- reading(sacurine_dir.c)
#' # or
#' sacurine.eset <- reading(NA,
#'                   files.ls = list(dataMatrix.tsvC = file.path(sacurine_dir.c, "dataMatrix.tsv"),
#'                                  sampleMetadata.tsvC = file.path(sacurine_dir.c, "sampleMetadata.tsv"),
#'                                  variableMetadata.tsvC = file.path(sacurine_dir.c, "variableMetadata.tsv")))
#' ## 2) Multiple sets
#' prometis_dir.c <- file.path(data_dir.c, "prometis")
#' prometis.mset <- reading(prometis_dir.c)
#' metabo.mset <- reading(prometis_dir.c, subsets.vc = "metabolomics")
#' # or
#' prometis.mset <- reading(NA,
#'                        files.ls = list(metabolomics = list(dataMatrix.tsvC = file.path(prometis_dir.c, "metabolomics", "dataMatrix.tsv"),
#'                                                     sampleMetadata.tsvC = file.path(prometis_dir.c, "metabolomics", "sampleMetadata.tsv"),
#'                                                     variableMetadata.tsvC = file.path(prometis_dir.c, "metabolomics", "variableMetadata.tsv")),
#'                                       proteomics = list(dataMatrix.tsvC = file.path(prometis_dir.c, "proteomics", "dataMatrix.tsv"),
#'                                                     sampleMetadata.tsvC = file.path(prometis_dir.c, "proteomics", "sampleMetadata.tsv"),
#'                                                     variableMetadata.tsvC = file.path(prometis_dir.c, "proteomics", "variableMetadata.tsv"))))
#' @rdname reading
#' @export
reading <- function(dir.c,
                    files.ls = NULL,
                    subsets.vc = NA,
                    report.c = c("none", "interactive", "myfile.txt")[2]) {
  
  if (!(report.c %in% c("none", "interactive")))
    sink(report.c, append = TRUE)
  
  x <- NULL
  
  ## Creating the ExpressionSet or building the list for the MultiDataSet
  
  if (!is.na(dir.c)) {
    
    if (!file.exists(dir.c))
      stop("Directory '", dir.c, "' was not found.")
    
    if (!file.info(dir.c)[, "isdir"])
      stop(dir.c, "' is not a directory.")
    
    dirVc <- dir(dir.c, full.names = TRUE)
    
    dirVl <- file.info(dirVc)[, "isdir"]
    
    subDirVc <- dirVc[dirVl]
    
    if (length(subDirVc) == 0) { ## ExpressionSet
      
      x <- .reading(dir.c,
                    dataMatrix.tsvC = NA,
                    sampleMetadata.tsvC = NA,
                    variableMetadata.tsvC = NA,
                    report.c = report.c)
      
    } else {## MultiDataSet
      
      names(subDirVc) <- basename(subDirVc)
      
      subDirVl <- sapply(subDirVc,
                         function(sub_dir.c) {
                           fileC <- list.files(sub_dir.c,
                                               pattern = "dataMatrix.tsv",
                                               full.names = TRUE)
                           length(fileC) == 1 && file.exists(fileC)
                         })
      
      if (sum(subDirVl) == 0) {
        
        stop("All subfolders have none or multiple 'dataMatrix.tsv' file(s):\n",
             paste(subDirVc, collapse = "\n"),
             "\n", call. = FALSE)
        
      } else if (sum(!subDirVl) > 0) {
        
        if (report.c != "none")
          message("No or multiple 'dataMatrix.tsv' file(s) was/were found in the following subfolders:\n",
                  paste(subDirVc[!subDirVl], collapse = "\n"),
                  "\nThe corresponding datasets will be skipped.")
        
        subDirVc <- subDirVc[subDirVl]
        
      }
      
      files.ls <- vector(mode = "list", length = length(subDirVc))
      names(files.ls) <- names(subDirVc)
      
      for (setC in names(files.ls)) {
        
        files.ls[[setC]] <- list(list.files(subDirVc[setC], pattern = "dataMatrix.tsv",
                                            full.names = TRUE),
                                 list.files(subDirVc[setC], pattern = "sampleMetadata.tsv",
                                            full.names = TRUE),
                                 list.files(subDirVc[setC], pattern = "variableMetadata.tsv",
                                            full.names = TRUE))
        
        # files.ls[setC] <- list(file.path(subDirVc[setC],
        #                                 c("dataMatrix.tsv",
        #                                   "sampleMetadata.tsv",
        #                                   "variableMetadata.tsv")))
        names(files.ls[[setC]]) <- c("dataMatrix.tsvC",
                                    "sampleMetadata.tsvC",
                                    "variableMetadata.tsvC")
        
      }
      
    }
    
  } else if (is.na(dir.c)) {
    
    subNamVc <- names(files.ls)
    
    if (sum(sapply(files.ls, is.list)) == 0) { ## ExpressionSet
      
      if (length(subNamVc) == 3 &&
          identical(subNamVc, c("dataMatrix.tsvC",
                                "sampleMetadata.tsvC",
                                "variableMetadata.tsvC"))) {
        
        x <- .reading(NA,
                      dataMatrix.tsvC = files.ls[["dataMatrix.tsvC"]],
                      sampleMetadata.tsvC = files.ls[["sampleMetadata.tsvC"]],
                      variableMetadata.tsvC = files.ls[["variableMetadata.tsvC"]])
        
      } else {
        
        stop("'files.ls does not contain any sublist nor is a list with names 'dataMatrix.tsvC', 'sampleMetadata.tsvC' and 'variableMetadata.tsvC' giving the corresponding file full names.",
             call. = FALSE)
        
      }
      
    } else {## MultiDataSet
      
      for (setC in names(files.ls)) {
        
        setLs <- files.ls[[setC]]
        
        if (!identical(names(setLs),
                       c("dataMatrix.tsvC",
                         "sampleMetadata.tsvC",
                         "variableMetadata.tsvC"))) {
          
          if (report.c != "none")
            message("The names of the following sublist are not 'dataMatrix.tsvC', 'sampleMetadata.tsvC' and 'variableMetadata.tsvC':\n",
                    setC,
                    "\nThe corresponding datasets will be skipped.")
          
          files.ls[[setC]] <- NULL
          
        } else if (!file.exists(setLs[["dataMatrix.tsvC"]])) {
          
          if (report.c != "none")
            message("No 'dataMatrix.tsv' file was found in the following sublist:\n",
                    setC,
                    "\nThe corresponding datasets will be skipped.")
          
          files.ls[[setC]] <- NULL
          
        }
        
      }
      
      if (length(files.ls) == 0)
        stop("None of the provided sublists meets the requirements for the creation of a dataset.",
             call. = FALSE)
      
    }
    
  }
  
  ## Creating the MultiDataSet
  
  if (is.null(x)) {
    
    if (!any(is.na(subsets.vc))) {
      
      misSetVc <- subsets.vc[!(subsets.vc %in% names(files.ls))]
      
      if (length(misSetVc))
        stop("The following selected subsets were not found in the subfolder(s) or sublist(s):\n",
             paste(misSetVc, collapse = "\n"))
      
      files.ls <- files.ls[subsets.vc]
      
    }
    
    # MultiDataSet
    
    x <- MultiDataSet::createMultiDataSet()
    
    for (setC in names(files.ls)) {
      
      setLs <- files.ls[[setC]]
      
      if (report.c != "none")
        message("Reading the '", setC, "' dataset...")
      
      ese <- .reading(NA,
                      setLs[["dataMatrix.tsvC"]],
                      setLs[["sampleMetadata.tsvC"]],
                      setLs[["variableMetadata.tsvC"]])
      
      ese$id <- rownames(Biobase::pData(ese))
      
      x <- MultiDataSet::add_eset(x, ese,
                                  dataset.type = setC,
                                  GRanges = NA)
      
    }
    
  }
  
  methods::validObject(x)
  
  if (report.c != "none")
    print(x)
  
  if (!(report.c %in% c("none", "interactive")))
    sink()
  
  return(invisible(x))
  
}


.reading <- function(dir.c,
                     dataMatrix.tsvC = NA,
                     sampleMetadata.tsvC = NA,
                     variableMetadata.tsvC = NA,
                     report.c = "interactive") {
  
  if (!is.na(dir.c) && !is.na(dataMatrix.tsvC))
    stop("Either 'dir.c' or 'dataMatrix.tsvC' argument must be set to NA",
         call. = FALSE)
  
  if (!is.na(dir.c)) {
    
    dataMatrixFileC <- list.files(dir.c, pattern = "dataMatrix.tsv",
                                  full.names = TRUE)
    
    if (length(dataMatrixFileC) == 0) {
      stop("No 'dataMatrix.tsv' file was found in the directory:\n",
           dir.c, call. = FALSE)
    } else if (length(dataMatrixFileC) > 1) {
      stop("Multiple 'dataMatrix.tsv' files were found in the directory:\n",
           dir.c, call. = FALSE)
    }
    
    sampleMetadataFileC <- list.files(dir.c, pattern = "sampleMetadata.tsv",
                                      full.names = TRUE)
    
    if (length(sampleMetadataFileC) == 0) {
      if (report.c != "none")
        warning("No 'sampleMetadata.tsv' file was found in the directory:\n",
                dir.c,
                "\nThe corresponding 'sampleMetadata' slot will be empty.",
                immediate. = FALSE, call. = FALSE)
      sampleMetadataFileC <- NA_character_
    } else if (length(sampleMetadataFileC) > 1) {
      stop("Multiple 'sampleMetadata.tsv' files were found in the directory:\n",
           dir.c, call. = FALSE)
    }
    
    variableMetadataFileC <- list.files(dir.c, pattern = "variableMetadata.tsv",
                                        full.names = TRUE)
    
    if (length(variableMetadataFileC) == 0) {
      if (report.c != "none")
        warning("No 'variableMetadata.tsv' file was found in the directory:\n",
                dir.c,
                "\nThe corresponding 'variableMetadata' slot will be empty.",
                immediate. = FALSE, call. = FALSE)
      variableMetadataFileC <- NA_character_
    } else if (length(variableMetadataFileC) > 1) {
      stop("Multiple 'variableMetadata.tsv' files were found in the directory:\n",
           dir.c, call. = FALSE)
    }
    
    tabFilVc <- c(dataMatrix = dataMatrixFileC,
                  sampleMetadata = sampleMetadataFileC,
                  variableMetadata = variableMetadataFileC)
    
  } else {
    
    tabFilVc <- c(dataMatrix = dataMatrix.tsvC,
                  sampleMetadata = sampleMetadata.tsvC,
                  variableMetadata = variableMetadata.tsvC)
    
    for (tabC in names(tabFilVc)) {
      
      tabFilC <- tabFilVc[tabC]
      
      if (!file.exists(tabFilC)) {
        
        if (tabC == "dataMatrix") {
          
          stop("The provided dataMatrix file was not found:\n",
               tabFilC, call. = FALSE)
          
        } else if (!file.exists(tabFilC)) {

          if (report.c != "none")
            warning("The following '", tabC,
                    "' file was not found:\n", tabFilC,
                    "\nThe corresponding '",
                    ifelse(tabC == "sampleMetadata", 'phenoData', 'featureData'),
                    "' slot will be empty.",
                    immediate. = TRUE, call. = FALSE)
          
          tabFilVc[tabC] <- NA
          
        }
      }
    }
    
  }
  
  for (tabC in names(tabFilVc)) {
    
    tabFilC <- tabFilVc[tabC]
    
    if (!is.na(tabFilC)) {
      
      ## R standards for row and column names in matrices and data frames
      .checkRformat(tabFilC)
      
      tabDF <- data.frame(data.table::fread(tabFilC,
                                            header = TRUE,
                                            sep = "\t"),
                          check.names = FALSE,
                          row.names = 1,
                          stringsAsFactors = FALSE)
      
      # looking for duplicates in column names
      colname_dup.vi <- which(duplicated(colnames(tabDF)))
      if (length(colname_dup.vi))
        stop("The '", tabC, "' file has several columns named: '",
             paste(colnames(tabDF[colname_dup.vi]), collapse = "', '"), "'",
             call. = FALSE)
      
      switch(tabC,
             dataMatrix = {
               tdatMN <- as.matrix(tabDF)
             },
             sampleMetadata = {
               samDF <- tabDF
             },
             variableMetadata = {
               varDF <- tabDF
             })
      
    } else {
      
      switch(tabC,
             sampleMetadata = {
               samDF <- data.frame(row.names = colnames(tdatMN))
             },
             variableMetadata = {
               varDF <- data.frame(row.names = rownames(tdatMN))
             })
      
    }
  }
  
  chkLs <- .checkW4Mformat(t(tdatMN), samDF, varDF,
                           infCw = report.c)
  
  if (!chkLs[["chkL"]]) {
    stop("Sample and/or variable names do not match between your tables. Use the 'report.c = NA' argument to get more feedback.",
         call. = FALSE)
  } else if (chkLs[["ordL"]]) {
    tdatMN <- t(chkLs[["datMN"]])
  }
  
  eset <- Biobase::ExpressionSet(assayData = tdatMN,
                                 phenoData = Biobase::AnnotatedDataFrame(data = samDF),
                                 featureData = Biobase::AnnotatedDataFrame(data = varDF),
                                 experimentData = Biobase::MIAME(title = ifelse(!is.na(dir.c),
                                                                                basename(dir.c),
                                                                                "")))
  
  methods::validObject(eset)
  
  return(eset)
  
}

#### writing (MultiDataSet) ####

#' @rdname writing
#' @export
setMethod("writing", "MultiDataSet",
          function(x,
                   dir.c,
                   files.ls = NULL,
                   overwrite.l = FALSE,
                   report.c = c("none", "interactive", "myfile.txt")[2]) {
            
            if (!(report.c %in% c("none", "interactive")))
              sink(report.c, append = TRUE)
            
            infTxtC <- report.c
            if (infTxtC != "none")
              infTxtC <- "interactive"
            
            setVc <- names(x)
            
            if (!is.na(dir.c)) {
              
              setDirVc <- file.path(dir.c, setVc)
              names(setDirVc) <- setVc
              
              if (file.exists(dir.c) && file.info(dir.c)[, "isdir"]) {
                
                dirVc <- dir(dir.c, full.names = TRUE)
                
                dirVl <- file.info(dirVc)[, "isdir"]
                
                subDirVc <- dirVc[dirVl]
                
                subDupVc <- intersect(setDirVc,
                                      subDirVc)
                
                if (length(subDupVc) && !overwrite.l)
                  stop("The following subfolder(s) were detected in your directory, please remove them or specify another parent directory to avoid overwriting:\n",
                       paste(subDupVc, collapse = "\n"))
                
              } else {
                
                dir.create(dir.c,
                           showWarnings = report.c != "none")
                
              }
              
              for (setC in names(setDirVc)) {
                
                if (report.c != "none")
                  message("Writing the '", setC, "' dataset...")
                
                set_dir.c <- setDirVc[setC]
                
                if (!(file.exists(set_dir.c) &&
                      file.info(set_dir.c)[, "isdir"]))
                  dir.create(set_dir.c,
                             showWarnings = report.c != "none")
                
                phenomis::writing(x[[setC]],
                                  set_dir.c,
                                  overwrite.l = overwrite.l,
                                  report.c = infTxtC)
                
              }
              
              if (report.c != "none")
                message("The subfolders have been written in the directory:\n", dir.c)
              
            } else if (is.na(dir.c)) {
              
              if (is.null(files.ls))
                stop("'files.ls' must be provided when 'dir.c' is set to NA",
                     call. = FALSE)
              
              if (is.null(names(files.ls)) || any(is.na(names(files.ls))))
                stop("All names of the sublists must be provided (they should match the names of the MultiDataSet datasets)",
                     call. = FALSE)
              
              filLisVl <- sapply(files.ls, is.list)
              
              if (!all(filLisVl))
                stop("The following element(s) of 'files.ls' is/are not sublist(s):\n",
                     paste(names(filLisVl)[!filLisVl], collapse = "\n"),
                     call. = FALSE)
              
              if (!identical(setVc, names(filLisVl)))
                stop("The name(s) of the 'x' MultiDataSet:\n",
                     paste(setVc, collapse = ", "),
                     "\ndo(es) not match the names of the sublists:\n",
                     paste(names(filLisVl), collapse = ", "),
                     call. = FALSE)
              
              for (setC in setVc) {
                
                filLs <- files.ls[[setC]]
                
                if (length(filLs) != 3 ||
                    !identical(names(filLs), c("dataMatrix.tsvC",
                                               "sampleMetadata.tsvC",
                                               "variableMetadata.tsvC")))
                  stop("The names of the '", setC,
                       "' sublist of 'files.ls are not identical to 'dataMatrix.tsvC', 'sampleMetadata.tsvC' and 'variableMetadata.tsvC'.")
                
                if (is.na(filLs[["dataMatrix.tsvC"]]))
                  stop("The 'dataMatrix.tsvC' file name from the '", setC,
                       "' sublist is missing (ie set to NA).")
                
                for (filC in names(filLs)) {
                  
                  filFulNamC <- filLs[[filC]]
                  
                  if (!is.na(filFulNamC)
                      && file.exists(filFulNamC)
                      && !overwrite.l)
                    stop("The following file from the '",
                         setC, "' sublist already exists:\n",
                         filFulNamC,
                         "\nPlease remove it or choose another name to avoid overwriting.")
                  
                }
                
              }
              
              for (setC in setVc) {
                
                filLs <- files.ls[[setC]]
                
                if (report.c != "none")
                  message("Writing the '", setC, "' dataset")
                
                phenomis::writing(x[[setC]],
                                  NA,
                                  files.ls = list(dataMatrix.tsvC = filLs[["dataMatrix.tsvC"]],
                                                 sampleMetadata.tsvC = filLs[["sampleMetadata.tsvC"]],
                                                 variableMetadata.tsvC = filLs[["variableMetadata.tsvC"]]),
                                  overwrite.l = overwrite.l,
                                  report.c = infTxtC)
                
              }
            }
            
            if (!(report.c %in% c("none", "interactive")))
              sink()
            
          })


#### writing (ExpressionSet) ####

#' @rdname writing
#' @export
setMethod("writing", "ExpressionSet",
          function(x,
                   dir.c,
                   files.ls = NULL,
                   overwrite.l = FALSE,
                   report.c = c("none", "interactive", "myfile.txt")[2]){
            
            if (!(report.c %in% c("none", "interactive")))
              sink(report.c, append = TRUE)
            
            if (!is.na(dir.c)) {
              
              if (!(file.exists(dir.c) && file.info(dir.c)[, "isdir"]))
                dir.create(dir.c,
                           showWarnings = report.c != "none")
              
              tabFilVc <- c(dataMatrix = file.path(dir.c, "dataMatrix.tsv"),
                            sampleMetadata = file.path(dir.c, "sampleMetadata.tsv"),
                            variableMetadata = file.path(dir.c, "variableMetadata.tsv"))
              
            } else if (is.na(dir.c)) {
              
              if (is.null(files.ls))
                stop("'files.ls' must be provided when 'dir.c' is set to NA",
                     call. = FALSE)
              
              if (length(files.ls) != 3 ||
                  !identical(names(files.ls), c("dataMatrix.tsvC",
                                               "sampleMetadata.tsvC",
                                               "variableMetadata.tsvC")))
                stop("The names of the 'files.ls' list are not identical to 'dataMatrix.tsvC', 'sampleMetadata.tsvC' and 'variableMetadata.tsvC'.",
                     call. = FALSE)
              
              if (is.na(files.ls[["dataMatrix.tsvC"]]))
                stop("The 'dataMatrix.tsvC' file name from the 'files.ls' list is missing (ie set to NA).")
              
              tabFilVc <- c(dataMatrix = files.ls[["dataMatrix.tsvC"]],
                            sampleMetadata = files.ls[["sampleMetadata.tsvC"]],
                            variableMetadata = files.ls[["variableMetadata.tsvC"]])
              
            }
            
            for (tabC in names(tabFilVc)) {
              
              tabFilC <- tabFilVc[tabC]
              
              if (!is.na(tabFilC) && file.exists(tabFilC) && !overwrite.l)
                stop("The following file already exists:\n", tabFilC,
                     "\nPlease choose another file name.")
              
            }
            
            ## Writing
            
            tdatMN <- Biobase::exprs(x)
            samDF <- Biobase::pData(x)
            varDF <- Biobase::fData(x)
            chkLs <- .checkW4Mformat(t(tdatMN), samDF, varDF, infCw = report.c)
            
            if (!chkLs[["chkL"]]) {
              stop("Sample and/or variable names do not match between your tables. Use the report.c = 'interactive' argument to get more feedback.")
            } else if (chkLs[["ordL"]]) {
              tdatMN <- t(chkLs[["datMN"]])
            }
            
            datDF <- cbind.data.frame(dataMatrix = rownames(tdatMN),
                                      as.data.frame(tdatMN))
            
            utils::write.table(datDF,
                               file = tabFilVc['dataMatrix'],
                               quote = FALSE,
                               row.names = FALSE,
                               sep = "\t")
            
            if (!is.na(dir.c) || !is.na(tabFilVc["sampleMetadata"])) {
              
              samDF <- cbind.data.frame(sampleMetadata = rownames(samDF),
                                        samDF)
              utils::write.table(samDF,
                                 file = tabFilVc['sampleMetadata'],
                                 quote = FALSE,
                                 row.names = FALSE,
                                 sep = "\t")
              
            }
            
            if (!is.na(dir.c) || !is.na(tabFilVc["variableMetadata"])) {
              
              varDF <- cbind.data.frame(variableMetadata = rownames(varDF),
                                        varDF)
              utils::write.table(varDF,
                                 file = tabFilVc['variableMetadata'],
                                 quote = FALSE,
                                 row.names = FALSE,
                                 sep = "\t")
              
            }
            
            if (report.c != "none")
              message("The following file(s) have been written:\n",
                      paste(tabFilVc[!is.na(basename(tabFilVc))], collapse = "\n"))
            
            if (!(report.c %in% c("none", "interactive")))
              sink()
            
          })


.checkRformat <- function(filCa) {
  
  rowVc <- data.table::fread(filCa,
                             header = TRUE,
                             sep = "\t")[[1]]
  
  colVc <- unlist(data.table::fread(filCa,
                                    header = FALSE,
                                    nrows = 1,
                                    sep = "\t")[1])[-1]
  
  if (any(duplicated(rowVc)))
    stop("The following ",
         ifelse(names(filCa) == 'sampleMetadata', 'sample', 'variable'),
         " name(s) is/are duplicated in the ",
         names(filCa),
         ": '",
         paste(rowVc[duplicated(rowVc)], collapse = "', '"), "'",
         call. = FALSE)
  
  if (any(duplicated(colVc)))
    stop("The following ", ifelse(names(filCa) == 'sampleMetadata', 'variable', 'sample'), " name(s) is/are duplicated in the ",
         names(filCa),
         ": '",
         paste(colVc[duplicated(colVc)], collapse = "', '"), "'",
         call. = FALSE)
  
}


.checkW4Mformat <- function(datMNw, samDFw, varDFw,
                            infCw = "interactive") {
  
  chkL <- TRUE
  ordL <- FALSE
  
  if (mode(datMNw) != "numeric") {
    message("The dataMatrix is not of the 'numeric' type.")
    chkL <- FALSE
  }
  
  if (!identical(rownames(datMNw), rownames(samDFw))) {
    ## checking sample names
    
    datSamDifVc <- setdiff(rownames(datMNw), rownames(samDFw))
    
    if (length(datSamDifVc)) {
      if (infCw != "none") {
        message("The following samples were found in the dataMatrix column names but not in the sampleMetadata row names:\n")
        print(cbind.data.frame(col = as.numeric(sapply(datSamDifVc,
                                                       function(samC) which(rownames(datMNw) == samC))),
                               name = datSamDifVc))
      }
      chkL <- FALSE
    }
    
    samDatDifVc <- setdiff(rownames(samDFw), rownames(datMNw))
    
    if (length(samDatDifVc)) {
      if (infCw != "none") {
        message("The following samples were found in the sampleMetadata row names but not in the dataMatrix column names:\n")
        print(cbind.data.frame(row = as.numeric(sapply(samDatDifVc, function(samC) which(rownames(samDFw) == samC))),
                               name = samDatDifVc))
      }
      chkL <- FALSE
    }
    
    if (nrow(datMNw) != nrow(samDFw)) {
      if (infCw != "none")
        message("The dataMatrix has ", nrow(datMNw), " columns (ie samples) whereas the sampleMetadata has ", nrow(samDFw), " rows.")
      chkL <- FALSE
    } else if (identical(gsub("^X", "", rownames(datMNw)), rownames(samDFw))) {
      if (infCw != "none")
        message("The dataMatrix column names start with an 'X' but not the sampleMetadata row names.")
      chkL <- FALSE
    } else if (identical(gsub("^X", "", rownames(samDFw)), rownames(datMNw))) {
      if (infCw != "none")
        message("The sampleMetadata row names start with an 'X' but not the dataMatrix column names.")
      chkL <- FALSE
    } else if (identical(sort(rownames(datMNw)), sort(rownames(samDFw)))) {
      if (infCw != "none")
        message("Re-ordering dataMatrix sample names to match sampleMetadata.")
      datMNw <- datMNw[rownames(samDFw), , drop = FALSE]
      stopifnot(identical(sort(rownames(datMNw)), sort(rownames(samDFw))))
      ordL <- TRUE
    } else {
      if (infCw != "none") {
        message("The dataMatrix column names and the sampleMetadata row names are not identical:\n")
        print(cbind.data.frame(indice = 1:nrow(datMNw),
                               dataMatrix_columnnames = rownames(datMNw),
                               sampleMetadata_rownames = rownames(samDFw))[rownames(datMNw) != rownames(samDFw), , drop = FALSE])
      }
      chkL <- FALSE
    }
    
  }
  
  if (!identical(colnames(datMNw), rownames(varDFw))) {
    ## checking variable names
    
    datVarDifVc <- setdiff(colnames(datMNw), rownames(varDFw))
    
    if (length(datVarDifVc)) {
      if (infCw != "none") {
        message("The following variables were found in the dataMatrix row names but not in the variableMetadata row names:\n")
        print(cbind.data.frame(row = as.numeric(sapply(datVarDifVc, function(varC) which(colnames(datMNw) == varC))),
                               name = datVarDifVc))
      }
      chkL <- FALSE
    }
    
    varDatDifVc <- setdiff(rownames(varDFw), colnames(datMNw))
    
    if (length(varDatDifVc)) {
      if (infCw != "none") {
        message("The following variables were found in the variableMetadata row names but not in the dataMatrix row names:\n")
        print(cbind.data.frame(row = as.numeric(sapply(varDatDifVc, function(varC) which(rownames(varDFw) == varC))),
                               name = varDatDifVc))
      }
      chkL <- FALSE
    }
    
    if (ncol(datMNw) != nrow(varDFw)) {
      if (infCw != "none")
        message("The dataMatrix has ",
                nrow(datMNw),
                " rows (ie variables) whereas the variableMetadata has ",
                nrow(varDFw),
                " rows.")
      chkL <- FALSE
    } else if (identical(sort(colnames(datMNw)), sort(rownames(varDFw)))) {
      if (infCw != "none")
        message("Re-ordering dataMatrix variable names to match variableMetadata.")
      datMNw <- datMNw[, rownames(varDFw), drop = FALSE]
      stopifnot(identical(sort(colnames(datMNw)), sort(rownames(varDFw))))
      ordL <- TRUE
    } else {
      if (infCw != "none") {
        message("The dataMatrix row names and the variableMetadata row names are not identical:\n")
        print(cbind.data.frame(row = 1:ncol(datMNw),
                               dataMatrix_rownames = colnames(datMNw),
                               variableMetadata_rownames = rownames(varDFw))[colnames(datMNw) != rownames(varDFw), , drop = FALSE])
      }
      chkL <- FALSE
    }
  }
  
  return(list(chkL = chkL,
              ordL = ordL,
              datMN = datMNw))
  
}















