#### transforming (MultiAssayExperiment) ####

#' @rdname transforming
#' @export
setMethod("transforming", signature(x = "MultiAssayExperiment"),
          function(x,
                   method.vc = c("log2", "log10", "sqrt")[1],
                   report.c = c("none", "interactive", "myfile.txt")[2]) {
            
            if (!(length(method.vc) %in% c(1, length(x)))) {
              stop("'The length of 'method.vc' should either be 1 or equal to the number of datasets")
            } else if (length(method.vc) == 1) {
              method.vc <- rep(method.vc, length(x))
            }
            names(method.vc) <- names(x)
            
            if (!(report.c %in% c("none", "interactive")))
              sink(report.c, append = TRUE)
            
            report_set.c <- report.c
            if (report_set.c != "none")
              report_set.c <- "interactive"
            
            for (set.c in names(x)) {
              
              if (report.c != "none")
                message("Transforming the '", set.c, "' dataset...")
              
              x[[set.c]] <- transforming(x = x[[set.c]],
                                         method.vc = method.vc[set.c],
                                         report.c = report_set.c)
              
            }
            
            if (!(report.c %in% c("none", "interactive")))
              sink()
            
            methods::validObject(x)
            
            return(invisible(x))
            
          })


#### transforming (SummarizedExperiment) ####

#' @rdname transforming
#' @export
setMethod("transforming", signature(x = "SummarizedExperiment"),
          function(x,
                   method.vc = c("log2", "log10", "sqrt")[1],
                   report.c = c("none", "interactive", "myfile.txt")[2]) {
            
            if (length(method.vc) > 1) {
              stop("Only one method should be provided for a 'SummarizedExperiment' object")
            } else {
              method.c <- method.vc
            }
            
            if (!(report.c %in% c("none", "interactive")))
              sink(report.c, append = TRUE)
            
            transf.mn <- .transforming(t(SummarizedExperiment::assay(x)),
                                       method.c = method.c,
                                       report.c != "none")
            
            SummarizedExperiment::assay(x) <- t(transf.mn)
            
            if (!(report.c %in% c("none", "interactive")))
              sink()
            
            methods::validObject(x)
            
            return(invisible(x))
            
          })

#### transforming (MultiDataSet) ####

#' @rdname transforming
#' @export
setMethod("transforming", signature(x = "MultiDataSet"),
          function(x,
                   method.vc = c("log2", "log10", "sqrt")[1],
                   report.c = c("none", "interactive", "myfile.txt")[2]) {
            
            if (!(length(method.vc) %in% c(1, length(x)))) {
              stop("'The length of 'method.vc' should either be 1 or equal to the number of datasets")
            } else if (length(method.vc) == 1) {
              method.vc <- rep(method.vc, length(x))
            }
            names(method.vc) <- names(x)
            
            if (!(report.c %in% c("none", "interactive")))
              sink(report.c, append = TRUE)
            
            report_set.c <- report.c
            if (report_set.c != "none")
              report_set.c <- "interactive"
            
            for (set.c in names(x)) {
              
              if (report.c != "none")
                message("Transforming the '", set.c, "' dataset...")
              
              ese <- transforming(x = x[[set.c]],
                                  method.vc = method.vc[set.c],
                                  report.c = report_set.c)
              
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


#### transforming (ExpressionSet) ####

#' @rdname transforming
#' @export
setMethod("transforming", signature(x = "ExpressionSet"),
          function(x,
                   method.vc = c("log2", "log10", "sqrt")[1],
                   report.c = c("none", "interactive", "myfile.txt")[2]) {
            
            if (length(method.vc) > 1) {
              stop("Only one method should be provided for a 'ExpressionSet' object")
            } else {
              method.c <- method.vc
            }
            
            if (!(report.c %in% c("none", "interactive")))
              sink(report.c, append = TRUE)
            
            transf.mn <- .transforming(t(Biobase::exprs(x)),
                                       method.c = method.c,
                                       report.c != "none")
            
            Biobase::exprs(x) <- t(transf.mn)
            
            if (!(report.c %in% c("none", "interactive")))
              sink()
            
            methods::validObject(x)
            
            return(invisible(x))
            
          })

.transforming <- function(data.mn, ## data (matrix of numerics; samples x variables)
                          method.c = c("log2", "log10", "sqrt")[1],
                          verbose.l = TRUE) {
  
  ## checking
  
  if (length(which(data.mn < 0)))
    stop("The 'dataMatrix' contains negative values")
  
  if (!(method.c %in% c("log2", "log10", "sqrt")))
    stop("The transforming method should be either 'log2', 'log10', or 'sqrt'")
  
  ## Number of missing values
  data_na.ml <- is.na(data.mn)
  data_na.i <- sum(data_na.ml)
  if (data_na.i > 0 && verbose.l)
    message("Missing values in the 'dataMatrix': ", data_na.i,
            " (", round(data_na.i / cumprod(dim(data.mn))[2] * 100), "%)")
  
  ## Number of zero values
  data_zero.ml <- data.mn < .Machine$double.eps # warning: data_zero.ml contains NA values
  data_zero.ml[data_na.ml] <- FALSE # NA values set to FALSE
  data_zero.i <- sum(data_zero.ml)
  if (data_zero.i > 0 && verbose.l)
    message("Zero values in the 'dataMatrix': ", data_zero.i,
            " (", round(data_zero.i / cumprod(dim(data.mn))[2] * 100), "%)")
  
  ## transforming
  
  switch(method.c,
         log2 = {
           
           if (verbose.l)
             message("'log2' transformation")
           
           data.mn[!data_zero.ml] <- log2(data.mn[!data_zero.ml])
           # log applied to non-negative and NA values
           
         },
         log10 = {
           
           if (verbose.l)
             message("'log10' transformation")
           
           data.mn[!data_zero.ml] <- log10(data.mn[!data_zero.ml])
           # log applied to non-negative and NA values
           
         },
         sqrt = {
           
           if (verbose.l)
             message("'Square root' transformation")
           
           data.mn <- sqrt(data.mn)
           
         })
  
  return(data.mn)
  
}
