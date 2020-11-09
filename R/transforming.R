#### transforming (MultiDataSet) ####

#' @rdname transforming
#' @export
setMethod("transforming", signature(x = "MultiDataSet"),
          function(x,
                   method.c = c("log2", "log10", "sqrt")[1],
                   report.c = c("none", "interactive", "myfile.txt")[2]) {
            
            if (!(report.c %in% c("none", "interactive")))
              sink(report.c, append = TRUE)
            
            infTxtC <- report.c
            if (infTxtC != "none")
              infTxtC <- "interactive"
            
            for (setC in names(x)) {
              
              if (report.c != "none")
                message("Transforming the '", setC, "' dataset...")
              
              ese <- phenomis::transforming(x = x[[setC]],
                                            method.c = method.c,
                                            report.c = infTxtC)
              
              x <- MultiDataSet::add_eset(x,
                                          ese,
                                          dataset.type = setC,
                                          GRanges = NA,
                                          overwrite = TRUE,
                                          warnings = FALSE)
              
            }
            
            if (!(report.c %in% c("none", "interactive")))
              sink()
            
            return(invisible(x))
            
          })


#### transforming (ExpressionSet) ####

#' @rdname transforming
#' @export
setMethod("transforming", signature(x = "ExpressionSet"),
          function(x,
                   method.c = c("log2", "log10", "sqrt")[1],
                   report.c = c("none", "interactive", "myfile.txt")[2]) {
            
            if (!(report.c %in% c("none", "interactive")))
              sink(report.c, append = TRUE)
            
            x <- .transforming(x,
                               method.c,
                               report.c != "none")
            
            methods::validObject(x)
            
            if (!(report.c %in% c("none", "interactive")))
              sink()
            
            return(invisible(x))
            
          })

.transforming <- function(eset,
                          methC = c("log2", "log10", "sqrt")[1],
                          verbL = TRUE) {
  
  exprs.mn <- Biobase::exprs(eset)
  
  ## checking
  
  if (length(which(exprs.mn < 0)))
    stop("The 'dataMatrix' contains negative values", call. = FALSE)
  
  ## Number of missing values
  exprs_na.ml <- is.na(exprs.mn)
  exprs_na.i <- sum(exprs_na.ml)
  if (exprs_na.i > 0 && verbL)
    message("Missing values in the 'dataMatrix': ", exprs_na.i,
            " (", round(exprs_na.i / cumprod(dim(exprs.mn))[2] * 100), "%)")
  
  ## Number of zero values
  exprs_zero.ml <- exprs.mn < .Machine$double.eps # warning: exprs_zero.ml contains NA values
  exprs_zero.ml[exprs_na.ml] <- FALSE # NA values set to FALSE
  exprs_zero.i <- sum(exprs_zero.ml)
  if (exprs_zero.i > 0 && verbL)
    message("Zero values in the 'dataMatrix': ", exprs_zero.i,
            " (", round(exprs_zero.i / cumprod(dim(exprs.mn))[2] * 100), "%)")
  
  ## transforming
  
  switch(methC,
         log2 = {
           
           if (verbL)
             message("'log2' transformation")
           
           exprs.mn[!exprs_zero.ml] <- log2(exprs.mn[!exprs_zero.ml])
           # log applied to non-negative and NA values
           
         },
         log10 = {
           
           if (verbL)
             message("'log10' transformation")
           
           exprs.mn[!exprs_zero.ml] <- log10(exprs.mn[!exprs_zero.ml])
           # log applied to non-negative and NA values
           
         },
         sqrt = {
           
           if (verbL)
             message("'Square root' transformation")
           
           exprs.mn <- sqrt(exprs.mn)
           
         })
  
  Biobase::exprs(eset) <- exprs.mn
  
  eset
  
}
