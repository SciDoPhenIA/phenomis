#### annotating (MultiAssayExperiment) ####
#' @rdname annotating
#' @export
setMethod("annotating", signature(x = "MultiAssayExperiment"),
          function(x,
                   database.c = c("chebi", "local.ms")[1],
                   param.ls = list(query.type = c("mz", "chebi.id", "kegg.id")[1],
                                   query.col = "mz",
                                   ms.mode = "pos",
                                   mz.tol = 10,
                                   mz.tol.unit = "ppm",
                                   fields = c("chebi.id", "name", "formula", "molecular.mass", "monoisotopic.mass"),
                                   fieldsLimit = 1,
                                   max.results = 3,
                                   local.ms.db = data.frame(),
                                   organism = "hsa",
                                   prefix = paste0(database.c, "."),
                                   sep = "|"),
                   report.c = c("none", "interactive", "myfile.txt")[2]) {
            
            if (!(report.c %in% c("none", "interactive")))
              sink(report.c, append = TRUE)
            
            report_set.c <- report.c
            if (report_set.c != "none")
              report_set.c <- "interactive"
            
            for (set.c in names(x)) {
              
              if (report.c != "none")
                message("Annotating the '", set.c, "' dataset:")
              
              se <- x[[set.c]]
              
              SummarizedExperiment::rowData(se) <- .annotating(SummarizedExperiment::rowData(se),
                                                               database.c = database.c,
                                                               param.ls = param.ls,
                                                               report_set.c != "none")
              
              x[[set.c]] <- se
              
            }
            
            if (!(report.c %in% c("none", "interactive")))
              sink()
            
            methods::validObject(x)
            
            return(invisible(x))
            
          })

#### annotating (SummarizedExperiment) ####
#' @rdname annotating
#' @export
setMethod("annotating", signature(x = "SummarizedExperiment"),
          function(x,
                   database.c = c("chebi", "local.ms")[1],
                   param.ls = list(query.type = c("mz", "chebi.id", "kegg.id")[1],
                                   query.col = "mz",
                                   ms.mode = "pos",
                                   mz.tol = 10,
                                   mz.tol.unit = "ppm",
                                   fields = c("chebi.id", "name", "formula", "molecular.mass", "monoisotopic.mass"),
                                   fieldsLimit = 1,
                                   max.results = 3,
                                   local.ms.db = data.frame(),
                                   organism = "hsa",
                                   prefix = paste0(database.c, "."),
                                   sep = "|"),
                   report.c = c("none", "interactive", "myfile.txt")[2]) {
            
            if (!(report.c %in% c("none", "interactive")))
              sink(report.c, append = TRUE)
            
            if (length(database.c) != 1)
              stop("'database.c' must be a character vector of length 1.",
                   call. = FALSE)
            
            availableDatabasesVc <- c("chebi", "local.ms", "kegg")
            if (!(database.c %in% availableDatabasesVc))
              stop("'database.c' must be in either: '", paste(availableDatabasesVc,
                                                              collapse = "', '"), "'.",
                   call. = FALSE)
            
            SummarizedExperiment::rowData(x) <- .annotating(SummarizedExperiment::rowData(x),
                                                            database.c = database.c,
                                                            param.ls = param.ls,
                                                            report.c != "none")
            
            if (!(report.c %in% c("none", "interactive")))
              sink()
            
            methods::validObject(x)
            
            return(invisible(x))
            
          })

#### annotating (MultiDataSet) ####
#' @rdname annotating
#' @export
setMethod("annotating", signature(x = "MultiDataSet"),
          function(x,
                   database.c = c("chebi", "local.ms")[1],
                   param.ls = list(query.type = c("mz", "chebi.id", "kegg.id")[1],
                                   query.col = "mz",
                                   ms.mode = "pos",
                                   mz.tol = 10,
                                   mz.tol.unit = "ppm",
                                   fields = c("chebi.id", "name", "formula", "molecular.mass", "monoisotopic.mass"),
                                   fieldsLimit = 1,
                                   max.results = 3,
                                   local.ms.db = data.frame(),
                                   organism = "hsa",
                                   prefix = paste0(database.c, "."),
                                   sep = "|"),
                   report.c = c("none", "interactive", "myfile.txt")[2]) {
            
            if (!(report.c %in% c("none", "interactive")))
              sink(report.c, append = TRUE)
            
            report_set.c <- report.c
            if (report_set.c != "none")
              report_set.c <- "interactive"
            
            for (set.c in names(x)) {
              
              if (report.c != "none")
                message("Annotating the '", set.c, "' dataset:")
              
              ese <- x[[set.c]]
              
              Biobase::fData(ese) <- .annotating(Biobase::fData(ese),
                                                 database.c = database.c,
                                                 param.ls = param.ls,
                                                 report_set.c != "none")
              
              x <- MultiDataSet::add_eset(x,
                                          ese,
                                          dataset.type = set.c,
                                          GRanges = NA,
                                          overwrite = TRUE,
                                          warnings = FALSE)
              
            }
            
            if (!(report.c %in% c("none", "interactive")))
              sink()
            
            methods::validObject(x)
            
            return(invisible(x))
            
          })

#### annotating (ExpressionSet) ####
#' @rdname annotating
#' @export
setMethod("annotating", signature(x = "ExpressionSet"),
          function(x,
                   database.c = c("chebi", "local.ms")[1],
                   param.ls = list(query.type = c("mz", "chebi.id", "kegg.id")[1],
                                   query.col = "mz",
                                   ms.mode = "pos",
                                   mz.tol = 10,
                                   mz.tol.unit = "ppm",
                                   fields = c("chebi.id", "name", "formula", "molecular.mass", "monoisotopic.mass"),
                                   fieldsLimit = 1,
                                   max.results = 3,
                                   local.ms.db = data.frame(),
                                   organism = "hsa",
                                   prefix = paste0(database.c, "."),
                                   sep = "|"),
                   report.c = c("none", "interactive", "myfile.txt")[2]) {
            
            if (!(report.c %in% c("none", "interactive")))
              sink(report.c, append = TRUE)
            
            if (length(database.c) != 1)
              stop("'database.c' must be a character vector of length 1.",
                   call. = FALSE)
            
            availableDatabasesVc <- c("chebi", "local.ms", "kegg")
            if (!(database.c %in% availableDatabasesVc))
              stop("'database.c' must be in either: '", paste(availableDatabasesVc,
                                                              collapse = "', '"), "'.",
                   call. = FALSE)
            
            Biobase::fData(x) <- .annotating(Biobase::fData(x),
                                             database.c = database.c,
                                             param.ls = param.ls,
                                             report.c != "none")
            
            if (!(report.c %in% c("none", "interactive")))
              sink()
            
            methods::validObject(x)
            
            return(invisible(x))
            
          })

#' Displays the parameters from the 'annotating' method
#'
#' The parameters and their default values are printed for the selected database
#' 
#' @rdname annotating
#' @export
#' @examples
#' phenomis::annotating_parameters()
#' phenomis::annotating_parameters("chebi")
annotating_parameters <- function(database.c = c("chebi", "local.ms", "kegg")[1]) {
  
  phenomis:::.annot_param_default(database.c = database.c,
                                  printL = TRUE)
  
}


.annotating <- function(feat.df, # variable metadata (dataframe or DataFrame; variables x variable metadata)
                        database.c,
                        param.ls,
                        verboseL) {
  
  paramDefaultLs <- .annot_param_check(database.c = database.c, param.ls = param.ls)
  
  switch(database.c,
         chebi = {feat.df <- .annot_chebi(feat.df = feat.df,
                                          param.ls = param.ls,
                                          paramDefaultLs = paramDefaultLs,
                                          verboseL = verboseL)},
         local.ms = {feat.df <- .annot_local.ms(feat.df = feat.df,
                                                param.ls = param.ls,
                                                paramDefaultLs = paramDefaultLs,
                                                verboseL = verboseL)},
         kegg = {feat.df <- .annot_kegg(feat.df = feat.df,
                                        param.ls = param.ls,
                                        paramDefaultLs = paramDefaultLs,
                                        verboseL = verboseL)})
  
  return(feat.df)
  
}


.annot_param_default <- function(database.c = c("chebi", "local.ms", "kegg")[1],
                                 printL = FALSE) {
  
  paramDefaultDF <- utils::read.table(system.file("extdata/annot_param_default.tsv",
                                                  package = "phenomis"),
                                      header = TRUE,
                                      quote = "",
                                      row.names = 1,
                                      sep = "\t",
                                      stringsAsFactors = FALSE)
  
  paramDefaultDF <- paramDefaultDF[paramDefaultDF[, database.c] == "x", ]
  
  for (dbC in c("chebi", "local.ms", "kegg"))
    paramDefaultDF[, dbC] <- NULL
  
  for (evalParamC in c("prefix"))
    paramDefaultDF[evalParamC, "default"] <- eval(parse(text = paramDefaultDF[evalParamC, "default"]))
  
  if (database.c %in% c("chebi", "kegg")) {
    
    for (fieldColC in c("default", "available")) {
      
      fieldDefaultC <- paramDefaultDF["fields", fieldColC]
      fieldDefaultVc <- unlist(strsplit(fieldDefaultC, split = "|", fixed = TRUE))
      fieldDefaultC <- fieldDefaultVc[grep(paste0(database.c, " = "), fieldDefaultVc)]
      paramDefaultDF["fields", fieldColC] <- gsub(paste0(database.c, " = "), "",
                                                  fieldDefaultC)
      
    }
    
  }
  
  if (printL) {
    
    printParamDF <- paramDefaultDF
    if ("fields" %in% rownames(printParamDF))
      printParamDF["fields", "default"] <- printParamDF["fields", "available"] <- "see below"
    message("Annotating parameters for the query of ", database.c, ":\n")
    print(printParamDF)
    message("'fields' possible values:\n'",
            paste(eval(parse(text = paramDefaultDF["fields", "available"])),
                  collapse = "', '"), "'\n")
    message("'fields' default values:\n'",
            paste(eval(parse(text = paramDefaultDF["fields", "default"])),
                  collapse = "', '"), "'\n")
    
  } else {
    
    return(paramDefaultDF)
    
  }
  
}


.annot_param_check <- function(database.c, param.ls) {
  
  paramDefaultDF <- .annot_param_default(database.c = database.c)
  
  # checking the names of the selected parameters
  
  paramErrorVc <- setdiff(names(param.ls), rownames(paramDefaultDF))
  
  if (length(paramErrorVc) > 0)
    stop("The following parameter name(s) is/are not recognized by the 'annotating' method for the selected '", database.c, "' database:\n'",
         paste(paramErrorVc, collapse = "', '"), "'.\n", call. = FALSE)
  
  # checking the mode of the parameters
  
  for (paramC in names(param.ls)) {
    
    if (paramC == "local.ms.db") {
      
      if (class(param.ls[[paramC]]) != "data.frame")
        stop("Parameter '", paramC, "' is expected to be of 'data.frame' class.",
             call. = FALSE)
      
    } else if (mode(param.ls[[paramC]]) != paramDefaultDF[paramC, "mode"])
      stop("Parameter '", paramC, "' is expected to be of '",
           paramDefaultDF[paramC, "mode"], "' mode.",
           call. = FALSE)
    
  }
  
  # checking the selected values
  
  for (availI in which(paramDefaultDF[, "available"] != "")) {
    
    availParamC <- rownames(paramDefaultDF)[availI]
    availValueVc <- eval(parse(text = paramDefaultDF[availI, "available"]))
    
    if (availParamC %in% names(param.ls) &&
        !(all(param.ls[[availParamC]] %in% availValueVc)))
      stop("Value(s) for the '", availParamC, "' parameter should be in:\n'",
           paste(availValueVc, collapse = "', '"), "'", call. = FALSE)
    
  }
  
  # building the list of default parameters
  
  paramDefaultLs <- lapply(1:nrow(paramDefaultDF),
                           function(parI) {
                             param <- paramDefaultDF[parI, "default"]
                             mode(param) <- paramDefaultDF[parI, "mode"]
                             return(param)
                           })
  names(paramDefaultLs) <- rownames(paramDefaultDF)
  
  if (database.c %in% c("chebi", "kegg")) {
    
    paramDefaultLs[["fields"]] <- eval(parse(text = paramDefaultLs[["fields"]]))
    
  } else if (database.c == "local.ms") {
    
    paramDefaultLs[["local.ms.db"]] <- data.frame(stringsAsFactors = FALSE)
    
  }
  
  return(invisible(paramDefaultLs))
  
}


.annot_chebi <- function(feat.df, param.ls, paramDefaultLs, verboseL) {
  
  param.ls <- c(param.ls,
                paramDefaultLs[setdiff(names(paramDefaultLs),
                                       names(param.ls))])
  
  if (!(param.ls[["query.col"]] %in% colnames(feat.df)))
    stop("The '", param.ls[["query.col"]], "' column could not be found in the fData.",
         call. = FALSE)
  
  queryVcn <- feat.df[, param.ls[["query.col"]]]
  
  if (param.ls[["query.type"]] == "chebi.id")
    queryVcn <- gsub("CHEBI:", "", queryVcn)
  
  optionWarnN <- options()$warn
  options(warn = -1)
  queriableVl <- sapply(queryVcn, function(query) {
    !is.na(query) && query != "" && !is.na(as.numeric(query))
  })
  options(warn = optionWarnN)
  
  if (sum(queriableVl) < 1)
    stop("No numeric were found in the '", param.ls[["query.col"]], "' column for ChEBI query.",
         call. = FALSE)
  
  queryVn <- as.numeric(queryVcn[queriableVl])
  
  fdataTempDF <- cbind.data.frame(.fdatarownames = rownames(feat.df),
                                  feat.df,
                                  stringsAsFactors = FALSE)
  
  fdataTempDF[queriableVl, param.ls[["query.col"]]] <- queryVn
  if (sum(!queriableVl))
    fdataTempDF[!queriableVl, param.ls[["query.col"]]] <- -c(1:sum(!queriableVl))
  
  mybiodb <- biodb::newInst()
  
  if (param.ls[["query.type"]] == "mz") {
    
    chebi <- mybiodb$getFactory()$createConn('chebi')
    
    resultDF <- chebi$annotateMzValues(data.frame(mz = queryVn),
                                       mz.tol = param.ls[["mz.tol"]],
                                       mz.tol.unit = param.ls[["mz.tol.unit"]],
                                       ms.mode = param.ls[["ms.mode"]],
                                       max.results = param.ls[["max.results"]],
                                       fields = param.ls[["fields"]],
                                       fieldsLimit = param.ls[["fieldsLimit"]],
                                       prefix = param.ls[["prefix"]])
    
    resultDF <- stats::aggregate(.~mz, resultDF, paste0, collapse = param.ls[["sep"]])
    
  } else if (param.ls[["query.type"]] == "chebi.id") {
    
    chebi <- mybiodb$getFactory()$createConn('chebi')
    
    entriesLs <- chebi$getEntry(queryVn)
    
    resultDF <- mybiodb$entriesToDataframe(entriesLs,
                                           fields = param.ls[["fields"]],
                                           limit = param.ls[["fieldsLimit"]],
                                           prefix = param.ls[["prefix"]])
    
  } else
    stop("'query.type' must be either 'mz' or 'chebi.id' for query of ChEBI.",
         call. = FALSE)
  
  mybiodb$terminate()
  
  fdataMergeDF <- merge(fdataTempDF, resultDF,
                        by.x = param.ls[["query.col"]],
                        by.y = paste0(ifelse(param.ls[["query.type"]] == "chebi.id",
                                             param.ls[["prefix"]], ""),
                                      param.ls[["query.type"]]),
                        all.x = TRUE)
  rownames(fdataMergeDF) <- fdataMergeDF[, ".fdatarownames"]
  
  fdataMergeDF <- fdataMergeDF[rownames(feat.df), , drop = FALSE]
  fdataMergeDF[, param.ls[["query.col"]]] <- feat.df[, param.ls[["query.col"]]]
  fdataMergeDF[[".fdatarownames"]] <- NULL
  
  return(fdataMergeDF)
  
}


.annot_local.ms <- function(feat.df, param.ls, paramDefaultLs, verboseL) {
  
  param.ls <- c(param.ls,
                paramDefaultLs[setdiff(names(paramDefaultLs),
                                       names(param.ls))])
  
  queryVn <- feat.df[, param.ls[["query.col"]]]
  
  missingVl <- sapply(queryVn, function(mzN) {
    is.na(mzN) || mzN == ""
  })
  
  queryVn <- queryVn[!missingVl]
  
  mybiodb <- biodb::newInst()
  
  conn <- mybiodb$getFactory()$createConn('mass.csv.file')
  conn$setDb(param.ls[["local.ms.db"]])
  
  resultDF <- conn$searchMsPeaks(data.frame(mz = queryVn),
                                 mz.tol = param.ls[["mz.tol"]],
                                 mz.tol.unit = param.ls[["mz.tol.unit"]],
                                 ms.mode = param.ls[["ms.mode"]],
                                 prefix = param.ls[["prefix"]])
  
  mybiodb$terminate()
  
  resultDF <- stats::aggregate(.~mz, resultDF, paste0, collapse = param.ls[["sep"]])
  
  fdataTempDF <- cbind.data.frame(.fdatarownames = rownames(feat.df),
                                  feat.df,
                                  stringsAsFactors = FALSE)
  fdataTempDF[missingVl, param.ls[["query.col"]]] <- -c(1:sum(missingVl))
  
  fdataMergeDF <- merge(fdataTempDF, resultDF, by.x = param.ls[["query.col"]],
                        by.y = "mz", all.x = TRUE)
  rownames(fdataMergeDF) <- fdataMergeDF[, ".fdatarownames"]
  
  fdataMergeDF <- fdataMergeDF[rownames(feat.df), , drop = FALSE]
  fdataMergeDF[, param.ls[["query.col"]]] <- feat.df[, param.ls[["query.col"]]]
  fdataMergeDF[[".fdatarownames"]] <- NULL
  
  return(fdataMergeDF)
  
}


.annot_kegg <- function(feat.df, param.ls, paramDefaultLs, verboseL) {
  
  # param.ls <- list(query.col = "database_identifier", query.type = "chebi")
  # in progress
  
  param.ls <- c(param.ls,
                paramDefaultLs[setdiff(names(paramDefaultLs),
                                       names(param.ls))])
  
  
}

