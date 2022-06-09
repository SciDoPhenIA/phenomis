#### filtering (MultiAssayExperiment) ####

#' @rdname filtering
#' @export
setMethod("filtering", signature(x = "MultiAssayExperiment"),
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
            
            for (set.c in names(x)) {
              
              if (report.c != "none")
                message("Filtering the '", set.c, "' dataset...")
              
              x[[set.c]] <- filtering(x = x[[set.c]],
                                      class.c = class.c,
                                      max_na_prop.n = max_na_prop.n,
                                      min_variance.n = min_variance.n,
                                      dims.vc = dims.vc,
                                      report.c = report_set.c)
              
            }
            
            if (!(report.c %in% c("none", "interactive")))
              sink()
            
            methods::validObject(x)
            
            return(invisible(x))
            
          })


#### filtering (SummarizedExperiment) ####

#' @rdname filtering
#' @export
setMethod("filtering", signature(x = "SummarizedExperiment"),
          function(x,
                   class.c = "",
                   max_na_prop.n = 0.2,
                   min_variance.n = .Machine$double.eps,
                   dims.vc = c("features", "samples"),
                   report.c = c("none", "interactive", "myfile.txt")[2]) {
            
            if (!(report.c %in% c("none", "interactive")))
              sink(report.c, append = TRUE)
            
            if (class.c != "")
              stopifnot(class.c %in% colnames(SummarizedExperiment::colData(x)))
            
            filt_names.ls <- .filtering(data.mn = t(SummarizedExperiment::assay(x)),
                                        samp.df = SummarizedExperiment::colData(x),
                                        set.c = x@metadata$experimentData@title,
                                        class.c = class.c,
                                        max_na_prop.n = max_na_prop.n,
                                        min_variance.n = min_variance.n,
                                        dims.vc = dims.vc,
                                        report.c != "none")
            
            x <- x[filt_names.ls[[2]], filt_names.ls[[1]]]
            
            if (!(report.c %in% c("none", "interactive")))
              sink()
            
            methods::validObject(x)
            
            return(invisible(x))
            
          })


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
            
            for (set.c in names(x)) {
              
              if (report.c != "none")
                message("Filtering the '", set.c, "' dataset...")
              
              ese <- filtering(x = x[[set.c]],
                               class.c = class.c,
                               max_na_prop.n = max_na_prop.n,
                               min_variance.n = min_variance.n,
                               dims.vc = dims.vc,
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
            
            filt_names.ls <- .filtering(data.mn = t(Biobase::exprs(x)),
                                        samp.df = Biobase::pData(x),
                                        set.c = Biobase::experimentData(x)@title,
                                        class.c = class.c,
                                        max_na_prop.n = max_na_prop.n,
                                        min_variance.n = min_variance.n,
                                        dims.vc = dims.vc,
                                        report.c != "none")
            
            x <- x[filt_names.ls[[2]], filt_names.ls[[1]]]
            
            if (!(report.c %in% c("none", "interactive")))
              sink()
            
            methods::validObject(x)
            
            return(invisible(x))
            
          })


.filtering <- function(data.mn,
                       samp.df,
                       set.c = "",
                       class.c,
                       max_na_prop.n,
                       min_variance.n,
                       dims.vc,
                       verbose.l) {
  
  if (!all(dims.vc %in% c("features", "samples")))
    stop("'dims.vc' must be 'features' and/or 'samples'")
  
  for (dim.c in dims.vc) {
    filt_names.ls <- .na_zerovar(data.mn = data.mn,
                                 samp.df = samp.df,
                                 set.c = set.c,
                                 dim.c = dim.c,
                                 class.c = class.c,
                                 max_na_prop.n = max_na_prop.n,
                                 min_variance.n = min_variance.n,
                                 verbose.l = verbose.l)
    data.mn <- data.mn[filt_names.ls[[1]], filt_names.ls[[2]], drop = FALSE]
    samp.df <- samp.df[filt_names.ls[[1]], , drop = FALSE]
  }
  
  return(dimnames(data.mn))
  
}

.na_zerovar <- function(data.mn,
                        samp.df,
                        set.c,
                        dim.c,
                        class.c,
                        max_na_prop.n,
                        min_variance.n,
                        verbose.l) {
  
  filter.vi <- c(nas_and_variance = NA_integer_)
  
  if (class.c != "" && dim.c == "features") {
    class.vc <- as.character(samp.df[, class.c])
    filter.vl <- apply(data.mn, 2,
                       function(feat.vn) {
                         all(tapply(feat.vn, class.vc, function(x) sum(is.na(x))/length(x)) > max_na_prop.n) ||
                           any(tapply(feat.vn, class.vc, function(x) stats::var(x, na.rm = TRUE)) < min_variance.n)
                       })
  } else
    filter.vl <- apply(data.mn, ifelse(dim.c == "features", 2, 1),
                       function(x) {
                         sum(is.na(x))/length(x) > max_na_prop.n ||
                           stats::var(x, na.rm = TRUE) < min_variance.n
                       })
  
  stopifnot(length(filter.vl) == dim(data.mn)[ifelse(dim.c == "samples", 1, 2)] &&
              !any(is.na(filter.vl)))
  
  filter.vi["nas_and_variance"] <- sum(filter.vl)
  
  dimnames.ls <- dimnames(data.mn)
  
  if (sum(filter.vl)) {
    
    if (all(filter.vl))
      stop("All ", dim.c, " would be discarded (because of too many missing values or too low variances). Please check your dataset or your thresholds.")
    
    if (dim.c == "samples") {
      dimnames.ls[[1]] <- dimnames.ls[[1]][!filter.vl]
      discard.c <- paste(dimnames.ls[[1]][filter.vl], collapse = ", ")
    } else {
      dimnames.ls[[2]] <- dimnames.ls[[2]][!filter.vl]
      discard.c <- paste(dimnames.ls[[2]][filter.vl], collapse = ", ")
    }
    
    if (verbose.l)
      message("Discarded ", sum(filter.vl), " ", dim.c,
              ifelse(set.c != "",
                     paste0(" in '", set.c, "'"),
                     ""),
              ": ",
              discard.c)
    
  }
  
  return(dimnames.ls)
  
}