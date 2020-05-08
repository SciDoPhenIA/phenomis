#### filtering (MultiDataSet) ####

#' @rdname filtering
#' @export
setMethod("filtering", signature(x = "MultiDataSet"),
          function(x,
                   class.c = "",
                   max_na_prop.n = 0.2,
                   min_variance.n = .Machine$double.eps,
                   dims.vc = c("features", "samples"),
                   report.c = c("none", "interactive", "myfile.txt")[2]) {
            
            if (!(report.c %in% c("none", "interactive")))
              sink(report.c, append = TRUE)
            
            report_set.c <- report.c
            if (report_set.c != "none")
              report_set.c <- "interactive"
            
            for (setC in names(x)) {
              
              if (report.c != "none")
                message("Filtering the '", setC, "' dataset...")
              
              ese <- phenomis::filtering(x = x[[setC]],
                                         class.c = class.c,
                                         max_na_prop.n = max_na_prop.n,
                                         min_variance.n = min_variance.n,
                                         dims.vc = dims.vc,
                                         report.c = report_set.c)
              
              x <- MultiDataSet::add_eset(x,
                                          ese,
                                          dataset.type = setC,
                                          GRanges = NA,
                                          overwrite = TRUE,
                                          warnings = FALSE)
              
            }
            
            if (!(report.c %in% c("none", "interactive")))
              sink()
            
            invisible(x)
            
          })


#### filtering (ExpressionSet) ####

#' @rdname filtering
#' @export
setMethod("filtering", signature(x = "ExpressionSet"),
          function(x,
                   class.c = "",
                   max_na_prop.n = 0.2,
                   min_variance.n = .Machine$double.eps,
                   dims.vc = c("features", "samples"),
                   report.c = c("none", "interactive", "myfile.txt")[2]) {
            
            if (!(report.c %in% c("none", "interactive")))
              sink(report.c, append = TRUE)
            
            if (class.c != "")
              stopifnot(class.c %in% Biobase::varLabels(x))
            
            x <- .filtering(eset = x,
                            class.c = class.c,
                            max_na_prop.n = max_na_prop.n,
                            min_variance.n = min_variance.n,
                            dims.vc = dims.vc,
                            report.c != "none")
            
            methods::validObject(x)
            
            if (!(report.c %in% c("none", "interactive")))
              sink()
            
            invisible(x)
            
          })


.filtering <- function(eset,
                       class.c,
                       max_na_prop.n,
                       min_variance.n,
                       dims.vc,
                       verbose.l) {
  
  if (!all(dims.vc %in% c("features", "samples")))
    stop("'dims.vc' must be 'features' and/or 'samples'")
  
  for (dim.c in dims.vc)
    eset <- .na_zerovar(eset = eset,
                        dim.c = dim.c,
                        class.c = class.c,
                        max_na_prop.n = max_na_prop.n,
                        min_variance.n = min_variance.n,
                        verbose.l = verbose.l)
  
  eset
  
}

.na_zerovar <- function(eset, dim.c, class.c, max_na_prop.n, min_variance.n, verbose.l) {
  
  filter.vi <- c(nas_and_variance = NA_integer_)
  
  if (class.c != "" && dim.c == "features") {
    class.vc <- as.character(Biobase::pData(eset)[, class.c])
    filter.vl <- apply(Biobase::exprs(eset), 1,
                       function(feat.vn) {
                         all(tapply(feat.vn, class.vc, function(x) sum(is.na(x))/length(x)) > max_na_prop.n) ||
                           any(tapply(feat.vn, class.vc, function(x) stats::var(x, na.rm = TRUE)) < min_variance.n)
                       })
  } else
    filter.vl <- apply(Biobase::exprs(eset), ifelse(dim.c == "features", 1, 2),
                       function(x) {
                         sum(is.na(x))/length(x) > max_na_prop.n ||
                           stats::var(x, na.rm = TRUE) < min_variance.n
                       })
  
  stopifnot(length(filter.vl) == dim(eset)[ifelse(dim.c == "features", "Features", "Samples")] &&
              !any(is.na(filter.vl)))
  
  filter.vi["nas_and_variance"] <- sum(filter.vl)
  
  if (sum(filter.vl)) {
    
    if (all(filter.vl))
      stop("All ", dim.c, " would be discarded (because of too many missing values or too low variances). Please check your dataset or your thresholds.")
    
    if (dim.c == "samples") {
      discard.c <- paste(Biobase::sampleNames(eset)[filter.vl], collapse = ", ")
      eset <- eset[, !filter.vl]
    } else {
      discard.c <- paste(Biobase::featureNames(eset)[filter.vl], collapse = ", ")
      eset <- eset[!filter.vl, ]
    }
    
    if (verbose.l)
      message("Discarded ", sum(filter.vl), " ", dim.c,
              ifelse(Biobase::experimentData(eset)@title != "",
                     paste0(" in '", Biobase::experimentData(eset)@title, "'"),
                     ""),
              ": ",
              discard.c)
    
  }
  
  eset
  
}