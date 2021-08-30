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
            
            transf.mn <- .transforming(t(Biobase::exprs(x)),
                                       method.c = method.c,
                                       report.c != "none")
            
            Biobase::exprs(x) <- t(transf.mn)
            
            methods::validObject(x)
            
            if (!(report.c %in% c("none", "interactive")))
              sink()
            
            return(invisible(x))
            
          })

.transforming <- function(data.mn,
                          method.c = c("log2", "log10", "sqrt")[1],
                          verbose.l = TRUE) {
  
  ## checking
  
  if (length(which(data.mn < 0)))
    stop("The 'dataMatrix' contains negative values", call. = FALSE)
  
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
