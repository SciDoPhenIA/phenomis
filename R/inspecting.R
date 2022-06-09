#### inspecting (MultiAssayExperiment) ####

#' @rdname inspecting
#' @export
setMethod("inspecting", signature(x = "MultiAssayExperiment"),
          function(x,
                   pool_as_pool1.l = FALSE,
                   pool_cv.n = 0.3,
                   loess_span.n = 1,
                   sample_intensity.c = c("median", "mean", "sum")[2],
                   title.c = NA,
                   plot_dims.l = TRUE,
                   col_batch.c = "batch",
                   col_injectionOrder.c = "injectionOrder",
                   col_sampleType.c = "sampleType",
                   figure.c = c("none", "interactive", "myfile.pdf")[2],
                   report.c = c("none", "interactive", "myfile.txt")[2]) {
            
            if (!(report.c %in% c("none", "interactive")))
              sink(report.c, append = TRUE)
            
            report_set.c <- report.c
            if (report_set.c != "none")
              report_set.c <- "interactive"
            
            if (!(figure.c %in% c("none", "interactive")))
              grDevices::pdf(figure.c)
            
            figure_set.c <- figure.c
            if (figure_set.c != "none")
              figure_set.c <- "interactive"
            
            if (plot_dims.l)
              .barplot_dims(as.list(MultiAssayExperiment::assays(x)), ifelse(!is.na(title.c), title.c, ""))
            
            for (set.c in names(x)) {
              
              if (report.c != "none")
                message("Inspecting the '", set.c, "' dataset...")
              
              x[[set.c]] <- inspecting(x = x[[set.c]],
                                       loess_span.n = loess_span.n,
                                       pool_as_pool1.l = pool_as_pool1.l,
                                       pool_cv.n = pool_cv.n,
                                       sample_intensity.c = sample_intensity.c,
                                       title.c = set.c,
                                       figure.c = figure_set.c,
                                       report.c = report_set.c)
              
            }
            
            if (!(figure.c %in% c("none", "interactive")))
              grDevices::dev.off()
            
            if (!(report.c %in% c("none", "interactive")))
              sink()
            
            methods::validObject(x)
            
            return(invisible(x))
            
          })


#### inspecting (SummarizedExperiment) ----

#' @rdname inspecting
#' @export
setMethod("inspecting", signature(x = "SummarizedExperiment"),
          function(x,
                   pool_as_pool1.l = FALSE,
                   pool_cv.n = 0.3,
                   loess_span.n = 1,
                   sample_intensity.c = c("median", "mean", "sum")[2],
                   title.c = NA,
                   plot_dims.l = TRUE,
                   col_batch.c = "batch",
                   col_injectionOrder.c = "injectionOrder",
                   col_sampleType.c = "sampleType",
                   figure.c = c("none", "interactive", "myfile.pdf")[2],
                   report.c = c("none", "interactive", "myfile.txt")[2]) {
            
            if (!(report.c %in% c("none", "interactive")))
              sink(report.c, append = TRUE)
            
            if (is.na(title.c))
              title.c <- x@metadata$experimentData@title
            
            
            inspected.ls <- .inspecting(data.mn = t(SummarizedExperiment::assay(x)),
                                        samp.df = SummarizedExperiment::colData(x),
                                        feat.df = SummarizedExperiment::rowData(x),
                                        pool_as_pool1.l = pool_as_pool1.l,
                                        pool_cv.n = pool_cv.n,
                                        loess_span.n = loess_span.n,
                                        sample_intensity.c = sample_intensity.c,
                                        title.c = title.c,
                                        col_batch.c = col_batch.c,
                                        col_injectionOrder.c = col_injectionOrder.c,
                                        col_sampleType.c = col_sampleType.c,
                                        figure.c = figure.c,
                                        report.c = report.c)
            
            SummarizedExperiment::colData(x) <- inspected.ls[["samp.df"]]
            SummarizedExperiment::rowData(x) <- inspected.ls[["feat.df"]]
            
            ## End
            
            if (!(report.c %in% c("none", "interactive")))
              sink()
            
            methods::validObject(x)
            
            return(invisible(x))
            
          })

#### inspecting (MultiDataSet) ####

#' @rdname inspecting
#' @export
setMethod("inspecting", signature(x = "MultiDataSet"),
          function(x,
                   pool_as_pool1.l = FALSE,
                   pool_cv.n = 0.3,
                   loess_span.n = 1,
                   sample_intensity.c = c("median", "mean", "sum")[2],
                   title.c = NA,
                   plot_dims.l = TRUE,
                   col_batch.c = "batch",
                   col_injectionOrder.c = "injectionOrder",
                   col_sampleType.c = "sampleType",
                   figure.c = c("none", "interactive", "myfile.pdf")[2],
                   report.c = c("none", "interactive", "myfile.txt")[2]) {
            
            if (!(report.c %in% c("none", "interactive")))
              sink(report.c, append = TRUE)
            
            report_set.c <- report.c
            if (report_set.c != "none")
              report_set.c <- "interactive"
            
            if (!(figure.c %in% c("none", "interactive")))
              grDevices::pdf(figure.c)
            
            figure_set.c <- figure.c
            if (figure_set.c != "none")
              figure_set.c <- "interactive"
            
            if (plot_dims.l)
              .barplot_dims(MultiDataSet::as.list(x), ifelse(!is.na(title.c), title.c, ""))
            
            for (set.c in names(x)) {
              
              if (report.c != "none")
                message("Inspecting the '", set.c, "' dataset...")
              
              ese <- x[[set.c]]
              
              ese <- inspecting(x = ese,
                                loess_span.n = loess_span.n,
                                pool_as_pool1.l = pool_as_pool1.l,
                                pool_cv.n = pool_cv.n,
                                sample_intensity.c = sample_intensity.c,
                                title.c = set.c,
                                figure.c = figure_set.c,
                                report.c = report_set.c)
              
              x <- MultiDataSet::add_eset(x, ese,
                                          dataset.type = set.c,
                                          GRanges = NA,
                                          overwrite = TRUE,
                                          warnings = FALSE)
              
            }
            
            if (!(figure.c %in% c("none", "interactive")))
              grDevices::dev.off()
            
            if (!(report.c %in% c("none", "interactive")))
              sink()
            
            methods::validObject(x)
            
            return(invisible(x))
            
          })

.barplot_dims <- function(data_mn.ls,
                          title.c = "Dataset dimensions",
                          cex_axis.i = 15,
                          cex_bar.i = 6,
                          bar_just.n = 0.8,
                          cex_title.i = 25) {
  
  set.i <- length(data_mn.ls)

  dims.mn <- sapply(names(data_mn.ls),
                    function(set.c)
                      dim(data_mn.ls[[set.c]]))

  dims.mn <- t(dims.mn)
  colnames(dims.mn) <- c("Features", "Samples")
  gg_barplot(dims.mn, title.c = title.c,
             col_levels.vc = c("Samples", "Features"),
             palette.vc = "Paired",
             cex_axis.i = cex_axis.i,
             cex_bar.i = cex_bar.i,
             cex_title.i = cex_title.i,
             bar_just.n = bar_just.n)
  
}


## inspecting (ExpressionSet) ----

#' @rdname inspecting
#' @export
setMethod("inspecting", signature(x = "ExpressionSet"),
          function(x,
                   pool_as_pool1.l = FALSE,
                   pool_cv.n = 0.3,
                   loess_span.n = 1,
                   sample_intensity.c = c("median", "mean", "sum")[2],
                   title.c = NA,
                   plot_dims.l = TRUE,
                   col_batch.c = "batch",
                   col_injectionOrder.c = "injectionOrder",
                   col_sampleType.c = "sampleType",
                   figure.c = c("none", "interactive", "myfile.pdf")[2],
                   report.c = c("none", "interactive", "myfile.txt")[2]) {
            
            if (!(report.c %in% c("none", "interactive")))
              sink(report.c, append = TRUE)
            
            if (is.na(title.c))
              title.c <- Biobase::experimentData(x)@title
            
            
            inspected.ls <- .inspecting(data.mn = t(Biobase::exprs(x)),
                                        samp.df = Biobase::pData(x),
                                        feat.df = Biobase::fData(x),
                                        pool_as_pool1.l = pool_as_pool1.l,
                                        pool_cv.n = pool_cv.n,
                                        loess_span.n = loess_span.n,
                                        sample_intensity.c = sample_intensity.c,
                                        title.c = title.c,
                                        col_batch.c = col_batch.c,
                                        col_injectionOrder.c = col_injectionOrder.c,
                                        col_sampleType.c = col_sampleType.c,
                                        figure.c = figure.c,
                                        report.c = report.c)
            
            Biobase::pData(x) <- inspected.ls[["samp.df"]]
            Biobase::fData(x) <- inspected.ls[["feat.df"]]
            
            ## End
            
            if (!(report.c %in% c("none", "interactive")))
              sink()
            
            methods::validObject(x)
            
            return(invisible(x))
            
          })


.inspecting <- function(data.mn, ## data (matrix of numerics; samples x variables)
                        samp.df, ## sample metadata (dataframe; samples x metadata)
                        feat.df, ## feature metadata (dataframe; features x metadata)
                        pool_as_pool1.l,
                        pool_cv.n,
                        loess_span.n,
                        sample_intensity.c,
                        title.c = title.c,
                        col_batch.c,
                        col_injectionOrder.c,
                        col_sampleType.c,
                        figure.c,
                        report.c) {
  
  ## Checking sample(s) or feature(s) with NA only or 0 variance
  
  .sampfeat_nas_zerovar(data.mn = data.mn,
                        set.c = title.c)
  
  ## Description
  
  if (report.c != "none") {
    message("Data description:")
    message("observations: ", nrow(data.mn))
    message("variables: ", ncol(data.mn))
    message("missing: ", format(sum(is.na(data.mn)), big.mark = ","),
            " (", round(sum(is.na(data.mn)) / cumprod(dim(data.mn))[2] * 100), "%)")
    message("0 values: ",
            format(sum(abs(data.mn) < .Machine[["double.eps"]], na.rm = TRUE), big.mark = ","),
            " (", round(sum(abs(data.mn) < .Machine[["double.eps"]],
                            na.rm = TRUE) / cumprod(dim(data.mn))[2] * 100), "%)")
    message("min: ", signif(min(data.mn, na.rm = TRUE), 2))
    message("mean: ", signif(mean(data.mn, na.rm = TRUE), 2))
    message("median: ", signif(stats::median(data.mn, na.rm = TRUE), 2))
    message("max: ", signif(max(data.mn, na.rm = TRUE), 2))
    
    if (col_sampleType.c %in% colnames(samp.df)) {
      message("Sample types:")
      print(table(samp.df[, col_sampleType.c]))
    }
  }
  
  ## Sample metrics
  
  sample_metrics.ls <- .sample_metrics(data.mn = data.mn,
                                       samp.df = samp.df,
                                       set.c = title.c)
  samp.df <- sample_metrics.ls[["samp.df"]]
  pca_metrics.ls <- sample_metrics.ls[["pca_metrics.ls"]]
  
  ## Variable metrics
  
  feat.df <- .variableMetrics(data.mn = data.mn,
                              samp.df = samp.df,
                              feat.df = feat.df,
                              pool_as_pool1.l = pool_as_pool1.l,
                              col_sampleType.c = col_sampleType.c)
  
  ## Figure
  
  if (figure.c != "none") {
    
    if (figure.c != "interactive")
      grDevices::pdf(figure.c)
    
    .plotMetrics(data.mn = data.mn,
                 samp.df = samp.df,
                 feat.df = feat.df,
                 pca_metrics.ls = pca_metrics.ls,
                 pool_cv.n = pool_cv.n,
                 loess_span.n = loess_span.n,
                 sample_intensity.c = sample_intensity.c,
                 col_batch.c = col_batch.c,
                 col_injectionOrder.c = col_injectionOrder.c,
                 col_sampleType.c = col_sampleType.c,
                 title.c = title.c)
    
    if (figure.c != "interactive")
      grDevices::dev.off()
    
  }
  
  return(invisible(list(samp.df = samp.df,
                        feat.df = feat.df)))
  
}

.sampfeat_nas_zerovar <- function(data.mn,
                                  set.c) {
  
  check.l <- TRUE
  
  feat_allna.vl <- apply(data.mn, 2,
                         function(feat.vn) all(is.na(feat.vn)))
  if (sum(feat_allna.vl, na.rm = TRUE)) {
    warning("The following feature(s) have NA only",
            ifelse(!is.na(set.c) && set.c != "", paste0(" in the '", set.c, "' dataset"), ""),
            ":\n",
            paste(Biobase::featureNames(x)[feat_allna.vl], collapse = ", "))
    check.l <- FALSE
  }
  
  samp_allna.vl <- apply(data.mn, 1,
                         function(samp.vn) all(is.na(samp.vn)))
  if (sum(samp_allna.vl, na.rm = TRUE)) {
    warning("The following sample(s) have NA only",
            ifelse(!is.na(set.c) && set.c != "", paste0(" in the '", set.c, "' dataset"), ""),
            ":\n",
            paste(Biobase::sampleNames(x)[samp_allna.vl], collapse = ", "))
    check.l <- FALSE
  }
  
  feat_zerovar.vl <- apply(data.mn, 2,
                           function(feat.vn) stats::var(feat.vn, na.rm = TRUE) < .Machine$double.eps)
  if (sum(feat_zerovar.vl, na.rm = TRUE)) {
    warning("The following feature(s) have zero variance",
            ifelse(!is.na(set.c) && set.c != "", paste0(" in the '", set.c, "' dataset"), ""),
            ":\n",
            paste(Biobase::featureNames(x)[feat_zerovar.vl], collapse = ", "))
    check.l <- FALSE
  }
  
  samp_zerovar.vl <- apply(data.mn, 1,
                           function(samp.vn) stats::var(samp.vn, na.rm = TRUE) < .Machine$double.eps)
  if (sum(samp_zerovar.vl, na.rm = TRUE)) {
    warning("The following sample(s) have zero variance",
            ifelse(!is.na(set.c) && set.c != "", paste0(" in the '", set.c, "' dataset"), ""),
            ":\n",
            paste(Biobase::sampleNames(x)[samp_zerovar.vl], collapse = ", "))
    check.l <- FALSE
  }
  
  if (!check.l)
    stop("Please remove the sample(s) and/or feature(s) with NA only or 0 variance by using the 'filtering' method, and apply the 'inspecting' function again on the filtered dataset.")
  
}


.pca_metrics <- function(data.mn,
                         samp.df,
                         pred.i = 2,
                         set.c = "") {
  
  nsamp.i <- nrow(data.mn)
  nfeat.i <- ncol(data.mn)
  
  ## Hotelling: p-value associated to the distance from the center in the first PCA score plane
  
  set.pca <- try(ropls::opls(data.mn, predI = pred.i,
                             crossvalI = min(nsamp.i, 7),
                             fig.pdfC = 'none', info.txtC = 'none'))
  
  if (inherits(set.pca, "try-error")) {
    stop("The PCA could not be computed",
         ifelse(set.c != "",
                paste0(" for set '", set.c, "'"),
                ""),
         ". Check for the presence of features with a high proportion of NA or a low variance and discard them with the 'filtering' method before starting 'inspecting' again.")
  }
  
  relative_var.vn <- ropls::getPcaVarVn(set.pca) / ncol(data.mn) ## for plotting
  
  score_pca.mn <- ropls::getScoreMN(set.pca)
  
  hotelling_df.i <- 2 * (nsamp.i - 1) * (nsamp.i^2 - 1) / (nsamp.i^2 * (nsamp.i - 2))
  
  inverse_covariance.mn <- solve(stats::cov(score_pca.mn[, 1:2]))
  
  hotelling_pval.vn <- apply(score_pca.mn[, 1:2],
                      1,
                      function(x)
                        1 - stats::pf(1 / hotelling_df.i * t(as.matrix(x)) %*% inverse_covariance.mn %*% as.matrix(x), 2, nsamp.i - 2))
  
  list(hotelling_df.i = hotelling_df.i,
       relative_var.vn = relative_var.vn,
       score_pca.mn = score_pca.mn,
       hotelling_pval.vn = hotelling_pval.vn)
  
}


.sample_metrics <- function(data.mn,
                            samp.df,
                            set.c) {
  
  pca_metrics.ls <- .pca_metrics(data.mn = data.mn,
                                 samp.df = samp.df,
                                 pred.i = 2,
                                 set.c = set.c)

  samp.df[, "pca_sco1"] <- pca_metrics.ls[["score_pca.mn"]][, 1]
  samp.df[, "pca_sco2"] <- pca_metrics.ls[["score_pca.mn"]][, 2]
  
  samp.df[, "hotel_pval"] <- pca_metrics.ls[["hotelling_pval.vn"]]
  
  ## p-value associated to number of missing values
  
  missing_zscore.vn <- .zscore(apply(data.mn,
                              1,
                              function(samp.vn) {
                                sum(is.na(samp.vn))
                              }))
  
  samp.df[, "miss_pval"] <- sapply(missing_zscore.vn, function(zsco.n) 2 * (1 - stats::pnorm(abs(zsco.n))))
  
  ## p-value associated to the deciles of the profiles
  
  decile.mn <- t(as.matrix(apply(data.mn,
                              1,
                              function(x) stats::quantile(x, 0.1 * 1:9, na.rm = TRUE))))
  
  decile_zscore.mn <- apply(decile.mn, 2, .zscore)
  
  decile_zscore_max.vn <- apply(decile_zscore.mn, 1, function(samp.vn) samp.vn[which.max(abs(samp.vn))])
  
  samp.df[, "deci_pval"] <- sapply(decile_zscore_max.vn, function(zsco.n) 2 * (1 - stats::pnorm(abs(zsco.n))))
  
  return(list(samp.df = samp.df,
              pca_metrics.ls = pca_metrics.ls))
  
}


.variableMetrics <- function(data.mn,
                             samp.df,
                             feat.df,
                             pool_as_pool1.l,
                             col_sampleType.c) { ## for the call to 'univariate'
  
  temp.df <- feat.df ## some of the intermediate metrics will not be included in feature metadata
  
  ## 'blank' observations
  
  if (col_sampleType.c %in% colnames(samp.df) && "blank" %in% samp.df[, col_sampleType.c]) {
    
    blank.vl <- samp.df[, col_sampleType.c] == "blank"
    
    if (sum(blank.vl) == 1) {
      temp.df[, "blank_mean"] <- data.mn[blank.vl, ]
    } else {
      temp.df[, "blank_mean"] <- apply(data.mn[blank.vl, , drop = FALSE],
                                     2,
                                     function(feat.vn) mean(feat.vn, na.rm = TRUE))
    }
    
    if (sum(blank.vl) == 1) {
      temp.df[, "blank_sd"] <- rep(0, nrow(feat.df))
    } else {
      temp.df[, "blank_sd"] <- apply(data.mn[blank.vl, , drop = FALSE],
                                   2,
                                   function(feat.vn) stats::sd(feat.vn, na.rm = TRUE))
    }
    
    temp.df[, "blank_CV"] <- temp.df[, "blank_sd"] / temp.df[, "blank_mean"]
    
  }
  
  ## 'sample' observations
  
  if (col_sampleType.c %in% colnames(samp.df) &&
      "sample" %in% samp.df[, col_sampleType.c]) {
    
    samp.vl <- samp.df[, col_sampleType.c] == "sample"
    
    if (sum(samp.vl) == 1) {
      temp.df[, "sample_mean"] <- data.mn[samp.vl, ]
    } else {
      temp.df[, "sample_mean"] <- apply(data.mn[samp.vl, , drop = FALSE], 2,
                                      function(feat.vn) mean(feat.vn, na.rm = TRUE))
    }
    
    if (sum(samp.vl) == 1) {
      temp.df[, "sample_sd"] <- rep(0, nrow(feat.df))
    } else {
      temp.df[, "sample_sd"] <- apply(data.mn[samp.vl, , drop = FALSE], 2,
                                    function(feat.vn) stats::sd(feat.vn, na.rm = TRUE))
    }
    
    temp.df[, "sample_CV"] <- temp.df[, "sample_sd"] / temp.df[, "sample_mean"]
    
  }
  
  ## 'blank' mean / 'sample' mean ratio
  
  if (all(c("blank_mean", "sample_mean") %in% colnames(temp.df)))
    feat.df[, "blankMean_over_sampleMean"] <- temp.df[, "blank_mean"] / temp.df[, "sample_mean"]
  
  ## 'pool' observations
  
  if (col_sampleType.c %in% colnames(samp.df) &&
      "pool" %in% samp.df[, col_sampleType.c]) {
    
    pool.vl <- samp.df[, col_sampleType.c] == "pool"
    
    if (sum(pool.vl) == 1) {
      temp.df[, "pool_mean"] <- data.mn[pool.vl, ]
    } else {
      temp.df[, "pool_mean"] <- apply(data.mn[pool.vl, , drop = FALSE], 2,
                                    function(feat.vn) mean(feat.vn, na.rm = TRUE))
    }
    
    if (sum(pool.vl) == 1) {
      temp.df[, "pool_sd"] <- rep(0, nrow(feat.df))
    } else {
      temp.df[, "pool_sd"] <- apply(data.mn[pool.vl, , drop = FALSE], 2,
                                  function(feat.vn) stats::sd(feat.vn, na.rm = TRUE))
    }
    
    feat.df[, "pool_CV"] <- temp.df[, "pool_sd"] / temp.df[, "pool_mean"]
    
  }
  
  ## 'pool' CV / 'sample' CV ratio
  
  if ("pool_CV" %in% colnames(feat.df) && "sample_CV" %in% colnames(temp.df))
    feat.df[, "poolCV_over_sampleCV"] <- feat.df[, "pool_CV"] / temp.df[, "sample_CV"]
  
  ## 'pool' dilutions
  
  if (col_sampleType.c %in% colnames(samp.df) &&
      any(grepl("pool.+", samp.df[, col_sampleType.c]))) {
    
    pool.vi <- grep("pool.*", samp.df[, col_sampleType.c]) ## pool, pool2, pool4, poolInter, ...
    
    pool_name.vc <- samp.df[pool.vi, col_sampleType.c]
    
    if (pool_as_pool1.l) {
      
      pool_name.vc[pool_name.vc == "pool"] <- "pool1" ## 'pool' -> 'pool1'
      
    } else {
      
      pool.vl <- pool_name.vc == "pool"
      pool.vi <- pool.vi[!pool.vl]
      pool_name.vc <- pool_name.vc[!pool.vl]
      
    }
    
    pool_dil.vc <- gsub("pool", "", pool_name.vc)
    
    pool_dil.vl <- sapply(pool_dil.vc, .allDigits)
    
    if (sum(pool_dil.vl)) {
      
      pool_name.vc <- pool_name.vc[pool_dil.vl]
      
      pool.vi <- pool.vi[pool_dil.vl]
      
      dilVn <- 1 / as.numeric(pool_dil.vc[pool_dil.vl])
      
      var.vn <- apply(data.mn[pool.vi, , drop = FALSE], 2,
                      function(feat.vn) stats::var(feat.vn))
      
      var.vl <- !is.na(var.vn) & var.vn > 0
      
      poolDil_pval.vn <- poolDil_cor.vn <- rep(NA, length(var.vn))
      
      poolDil_cor.vn[var.vl] <- apply(data.mn[pool.vi, var.vl, drop = FALSE], 2,
                                       function(feat.vn) stats::cor(dilVn, feat.vn))
      
      feat.df[, "poolDil_cor"] <- poolDil_cor.vn
      
      poolDil_pval.vn[var.vl] <- apply(data.mn[pool.vi, var.vl, drop = FALSE], 2,
                                       function(feat.vn) stats::cor.test(dilVn, feat.vn)[["p.value"]])
      
      feat.df[, "poolDil_pval"] <- poolDil_pval.vn
      
    }
    
  }
  
  return(feat.df)
  
}


.plotMetrics <- function(data.mn,
                         samp.df,
                         feat.df,
                         pca_metrics.ls,
                         pool_cv.n,
                         loess_span.n,
                         sample_intensity.c,
                         col_batch.c,
                         col_injectionOrder.c,
                         col_sampleType.c,
                         title.c) {
  
  ## Constants
  
  marLs <- list(tit = c(0.6, 1.1, 1.1, 0.6),
                sca = c(0.6, 3.1, 4.1, 0.6),
                ima = c(0.6, 2.6, 4.1, 0.9),
                dri = c(3.5, 3.6, 1.1, 0.6),
                pca = c(3.5, 3.6, 1.1, 0.9))
 
  ## Script
  
  opar <- graphics::par(font = 2,
                        font.axis = 2,
                        font.lab = 2,
                        pch = 18)
  
  layVn <- c(1, 2, 3, 3,
             4, 4, 4, 5)
  
  graphics::layout(matrix(layVn,
                          byrow = TRUE,
                          nrow = 2),
                   heights = c(3.5, 3.5),
                   widths = c(1.5, 1, 1, 3.5))
  
  # Colors
  
  sample_color.vc <- .sample_color(samp.df = samp.df,
                                   col_sampleType.c = col_sampleType.c)
 
  ## tit: Title
  
  graphics::par(mar = marLs[["tit"]])
  graphics::plot(0:1, bty = "n", type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
  graphics::text(1, 0.95, adj = 0, cex = 1.2, labels = title.c)
  graphics::text(1, 0.75, adj = 0, labels = paste0("NAs: ",
                                                   round(length(which(is.na(c(data.mn)))) / cumprod(dim(data.mn))[2] * 100), "%"))
  graphics::text(1, 0.68, adj = 0, labels = paste0("0 values: ",
                                                   round(sum(abs(data.mn) < .Machine[["double.eps"]], na.rm = TRUE) / cumprod(dim(data.mn))[2] * 100, 2), "%"))
  graphics::text(1, 0.61, adj = 0, labels = paste0("min: ", signif(min(data.mn, na.rm = TRUE), 2)))
  graphics::text(1, 0.54, adj = 0, labels = paste0("median: ", signif(stats::median(data.mn, na.rm = TRUE), 2)))
  graphics::text(1, 0.47, adj = 0, labels = paste0("mean: ", signif(mean(data.mn, na.rm = TRUE), 2)))
  graphics::text(1, 0.40, adj = 0, labels = paste0("max: ", signif(max(data.mn, na.rm = TRUE), 2)))
  if (col_sampleType.c %in% colnames(samp.df) &&
      "pool" %in% samp.df[, col_sampleType.c])
    graphics::text(1,
                   0.33,
                   adj = 0,
                   labels = paste0("CVpool<",
                                   round(pool_cv.n * 100), "%: ",
                                   round(sum(feat.df[, "pool_CV"] < pool_cv.n, na.rm = TRUE) / ncol(data.mn) * 100),
                                   "%"))
  
  ## sca: Color scale
  
  .plot_color_scale(data.mn = data.mn,
                    mar.vn = marLs[["sca"]])
  
  ## ima: Image
  
  .plot_image(data.mn = data.mn,
              mar.vn = marLs[["ima"]])
  
  ## dri: Analytical drift
  
  .plot_drift(data.mn = data.mn,
              samp.df = samp.df,
              loess_span.n = loess_span.n,
              sample_intensity.c = sample_intensity.c,
              col_batch.c = col_batch.c,
              col_injectionOrder.c = col_injectionOrder.c,
              col_sampleType.c = col_sampleType.c,
              mar.vn = marLs[["dri"]])
  
  ## pca: PCA and Hotelling ellipse
  
  .plot_pca_metrics(data.mn = data.mn,
                    samp.df = samp.df,
                    pca_metrics.ls = pca_metrics.ls,
                    col_sampleType.c = col_sampleType.c,
                    mar.vn = marLs[["pca"]])
  
  ## resetting par
  
  graphics::par(mfrow = c(1, 1))
  graphics::par(opar)
  
}


.palette <- function(palette.c = "heat") {
  if (palette.c == "heat") {
    return(rev(grDevices::rainbow(ceiling(256 * 1.5))[1:256]))
  }
}

.plot_pretty_axis <- function(valVn,
                              lenN) {
  
  if (NA %in% valVn) {
    warning("NA in valVn")
    valVn <- as.vector(stats::na.omit(valVn))
  }
  
  if (lenN < length(valVn))
    stop("The length of in vector must be inferior to the length of the length parameter.")
  
  if (length(valVn) < lenN)
    valVn <- seq(from = min(valVn), to = max(valVn),
                 length.out = lenN)
  
  preValVn <- pretty(valVn)
  
  preLabVn <- preAtVn <- c()
  
  for (n in 1:length(preValVn))
    if (min(valVn) < preValVn[n] && preValVn[n] < max(valVn)) {
      preLabVn <- c(preLabVn, preValVn[n])
      preAtVn <- c(preAtVn, which(abs(valVn - preValVn[n]) == min(abs(valVn - preValVn[n])))[1])
    }
  
  return(list(atVn = preAtVn,
              labVn = preLabVn))
  
}

.plot_color_scale <- function(data.mn,
                              mar.vn = c(0.6, 3.1, 4.1, 0.6)) {
  
  palette.vc <- .palette("heat")
  
  graphics::par(mar = mar.vn)
  
  ylimVn <- c(0, 256)
  ybottomVn <- 0:255
  ytopVn <- 1:256
  
  graphics::plot(x = 0,
                 y = 0,
                 font.axis = 2,
                 font.lab = 2,
                 type = "n",
                 xlim = c(0, 1),
                 ylim = ylimVn,
                 xlab = "",
                 ylab = "",
                 xaxs = "i",
                 yaxs = "i",
                 xaxt = "n",
                 yaxt = "n")
  
  graphics::rect(xleft = 0,
                 ybottom = ybottomVn,
                 xright = 1,
                 ytop = ytopVn,
                 col = palette.vc,
                 border = NA)
  
  eval(parse(text = paste0("axis(at = .plot_pretty_axis(c(ifelse(min(data.mn, na.rm = TRUE) == -Inf, yes = 0, no = min(data.mn, na.rm = TRUE)) , max(data.mn, na.rm = TRUE)), 256)$atVn,
                           font = 2,
                           font.axis = 2,
                           labels = .plot_pretty_axis(c(ifelse(min(data.mn, na.rm = TRUE) == -Inf, yes = 0, no = min(data.mn, na.rm = TRUE)), max(data.mn, na.rm = TRUE)), 256)$labVn,
                           las = 1,
                           lwd = 2,
                           lwd.ticks = 2,
                           side = 2,
                           xpd = TRUE)")))
  
  graphics::arrows(graphics::par("usr")[1],
                   graphics::par("usr")[4],
                   graphics::par("usr")[1],
                   graphics::par("usr")[3],
                   code = 0,
                   lwd = 2,
                   xpd = TRUE)
  
  graphics::box(lwd = 2)
  
}


.plot_image <- function(data.mn,
                        mar.vn = c(0.6, 2.6, 4.1, 0.9)) {
  
  palette.vc <- .palette("heat")
  
  graphics::par(mar = mar.vn)
  
  image.mn <- data.mn[rev(1:nrow(data.mn)), , drop = FALSE]
  
  graphics::image(x = 1:nrow(image.mn),
                  y = 1:ncol(image.mn),
                  z = image.mn,
                  col = palette.vc,
                  font.axis = 2,
                  font.lab = 2,
                  xaxt = "n",
                  yaxt = "n",
                  xlab = "",
                  ylab = "")
  
  if (length(rownames(data.mn)) == 0) {
    rowNamVc <- rep("", times = nrow(data.mn))
  } else
    rowNamVc <- rownames(data.mn)
  
  if (length(colnames(data.mn)) == 0) {
    colNamVc <- rep("", times = ncol(data.mn))
  } else
    colNamVc <- colnames(data.mn)
  
  xlaVc <- paste(paste(rep("[", 2),
                       c(1, nrow(image.mn)),
                       rep("] ", 2),
                       sep = ""),
                 rep("\n", times = 2),
                 c(colNamVc[1], utils::tail(colNamVc, 1)),
                 sep = "")
  
  for (k in 1:2)
    graphics::axis(side = 3,
                   hadj = c(0, 1)[k],
                   at = c(1, nrow(image.mn))[k],
                   cex = 0.8,
                   font = 2,
                   labels = xlaVc[k],
                   line = -0.5,
                   tick = FALSE)
  
  
  ylaVc <- paste(paste(rep("[", times = 2),
                       c(ncol(image.mn), 1),
                       rep("]", times = 2),
                       sep = ""),
                 rep("\n", times = 2),
                 c(utils::tail(rowNamVc, 1), rowNamVc[1]),
                 sep = "")
  
  for (k in 1:2)
    graphics::axis(side = 2,
                   at = c(1, ncol(image.mn))[k],
                   cex = 0.8,
                   font = 2,
                   hadj = c(0, 1)[k],
                   labels = ylaVc[k],
                   las = 0,
                   line = -0.5,
                   lty = "blank",
                   tick = FALSE)
  
  graphics::box(lwd = 2)
  
}


.plot_drift <- function(data.mn,
                        samp.df,
                        loess_span.n = 1,
                        sample_intensity.c = "mean",
                        col_batch.c = "batch",
                        col_injectionOrder.c = "injectionOrder",
                        col_sampleType.c = "sampleType",
                        mar.vn = c(3.5, 3.6, 1.1, 0.6)) {
  
  sample_color.vc <- .sample_color(samp.df = samp.df,
                                   col_sampleType.c = col_sampleType.c)
  
  graphics::par(mar = mar.vn)
  
  ## ordering

  samp.df[, "ordIniVi"] <- 1:nrow(data.mn)
  
  if (col_injectionOrder.c %in% colnames(samp.df)) {
    ordNamC <- "Injection Order"
    if (col_batch.c %in% colnames(samp.df)) {
      ordVi <- order(samp.df[, col_batch.c],
                     samp.df[, col_injectionOrder.c])
    } else
      ordVi <- order(samp.df[, col_injectionOrder.c])
  } else {
    ordNamC <- "Samples"
    ordVi <- 1:nrow(data.mn)
  }
  
  data.mn <- data.mn[ordVi, , drop = FALSE]
  samp.df <- samp.df[ordVi, ]
  sample_color_ordered.vc <- sample_color.vc[ordVi]
  
  if (col_batch.c %in% colnames(samp.df))
    batch.table <- table(samp.df[, col_batch.c])
  
  sample_means.vn <- eval(parse(text = paste0("apply(data.mn, 1, function(obsVn) ",
                                              sample_intensity.c, "(obsVn, na.rm = TRUE))")))
  
  graphics::plot(sample_means.vn,
                 col = sample_color_ordered.vc,
                 pch = 18,
                 # ylim = range(data.mn, na.rm = TRUE),
                 type = "n",
                 xaxs = "i",
                 xlab = "",
                 ylab = "")
  
  # for (samI in 1:nrow(data.mn))
  #   graphics::boxplot(data.mn[samI, ],
  #                     at = samI,
  #                     add = TRUE)
  
  graphics::points(sample_means.vn,
                   col = sample_color_ordered.vc,
                   pch = 18)
  
  graphics::mtext(ordNamC,
                  cex = 0.7,
                  line = 2,
                  side = 1)
  
  graphics::mtext(paste0(toupper(substr(sample_intensity.c, 1, 1)), substr(sample_intensity.c, 2, nchar(sample_intensity.c)), " of variable intensities"),
                  cex = 0.7,
                  line = 2,
                  side = 2)
  
  if (col_batch.c %in% colnames(samp.df) && length(unique(samp.df[, "batch"])) > 1) {
    
    graphics::abline(v = cumsum(batch.table) + 0.5,
                     col = "red")
    
    graphics::mtext(names(batch.table),
                    at = batch.table / 2 + c(0, cumsum(batch.table[-length(batch.table)])),
                    cex = 0.7)
    
    for (batC in names(batch.table)) {
      
      batch_seq.vi <- which(samp.df[, col_batch.c] == batC)
      
      if (col_sampleType.c %in% colnames(samp.df)) {
        batch_sample.vi <- intersect(batch_seq.vi,
                                     grep("sample", samp.df[, col_sampleType.c]))
      } else
        batch_sample.vi <- batch_seq.vi
      
      graphics::lines(batch_seq.vi,
                      .loess(sample_means.vn, batch_sample.vi, batch_seq.vi, loess_span.n),
                      col = .sample_palette("sample"))
      
      if (col_sampleType.c %in% colnames(samp.df) &&
          "pool" %in% samp.df[, col_sampleType.c]) {
        
        batch_pool.vi <- intersect(batch_seq.vi,
                                   grep("^pool$", samp.df[, col_sampleType.c]))
        
        graphics::lines(batch_seq.vi,
                        .loess(sample_means.vn, batch_pool.vi, batch_seq.vi, loess_span.n),
                        col = .sample_palette("pool"))
        
      }
      
    }
    
  } else {
    
    batch_seq.vi <- 1:nrow(samp.df)
    
    if (col_sampleType.c %in% colnames(samp.df)) {
      batch_sample.vi <- intersect(batch_seq.vi,
                                   grep("sample", samp.df[, col_sampleType.c]))
    } else
      batch_sample.vi <- batch_seq.vi
    
    graphics::lines(batch_seq.vi,
                    .loess(sample_means.vn, batch_sample.vi, batch_seq.vi, loess_span.n),
                    col = .sample_palette("sample"))
    
    if (col_sampleType.c %in% colnames(samp.df) &&
        "pool" %in% samp.df[, col_sampleType.c]) {
      
      batch_pool.vi <- intersect(batch_seq.vi,
                                 grep("^pool$", samp.df[, col_sampleType.c]))
      
      graphics::lines(batch_seq.vi,
                      .loess(sample_means.vn, batch_pool.vi, batch_seq.vi, loess_span.n),
                      col = .sample_palette("pool"))
      
      
    }
    
  }
  
  # legend
  
  if (col_sampleType.c %in% colnames(samp.df)) {
    obsColVuc <- sample_color_ordered.vc[sort(unique(names(sample_color_ordered.vc)))]
    legOrdVc <- c("blank", paste0("pool", 8:1), "pool", "other", "sample")
    obsColVuc <- obsColVuc[legOrdVc[legOrdVc %in% names(obsColVuc)]]
    
    graphics::text(rep(graphics::par("usr")[2], times = length(obsColVuc)),
                   graphics::par("usr")[3] + (0.03 + 1:length(obsColVuc) * 0.03) * diff(graphics::par("usr")[3:4]),
                   adj = 1,
                   col = obsColVuc,
                   font = 2,
                   labels = names(obsColVuc),
                   pos = 2)
  }
  
}


.plot_pca_metrics <- function(data.mn,
                              samp.df,
                              pred.i = 2,
                              show_pred.vi = c(1, 2),
                              pca_metrics.ls = NULL,
                              col_sampleType.c = "sampleType",
                              labels.l = TRUE,
                              mar.vn = c(3.5, 3.6, 1.1, 0.9)) {
  
  if (is.null(pca_metrics.ls))
    pca_metrics.ls <- .pca_metrics(data.mn = data.mn,
                                   samp.df = samp.df,
                                   pred.i = pred.i)
  
  score_pca.mn <- pca_metrics.ls[["score_pca.mn"]]
  
  if (ncol(score_pca.mn) < max(show_pred.vi))
    stop("Not enough pca components computed")
  
  sample_color.vc <- .sample_color(samp.df = samp.df, col_sampleType.c = col_sampleType.c)
  if ("injectionOrder" %in% colnames(samp.df) &&
      "sampleType" %in% colnames(samp.df)) {
    palette.vc <- rev(rainbow(100, end = 4/6))
    inj_rank.vi <- rank(samp.df[, "injectionOrder"])
    sample_color.vc <- palette.vc[round((inj_rank.vi - 1) / diff(range(inj_rank.vi)) * 99) + 1]
  }
  
  sample_label.vc <- rownames(data.mn)
  if ("sampleType" %in% colnames(samp.df)) {
    sample_label.vc <- samp.df[, "sampleType"]
    sample_label.vc[sample_label.vc == "sample"] <- "s"
    sample_label.vc[sample_label.vc == "pool"] <- "QC"
  }
  
  
  graphics::par(mar = mar.vn)
  
  graphics::plot(score_pca.mn[, show_pred.vi],
                 type = "n",
                 xlab = "",
                 ylab = "")
  graphics::mtext(paste0("t",
                         show_pred.vi[1],
                         " (",
                         round(pca_metrics.ls[["relative_var.vn"]][show_pred.vi[1]] * 100), "%)"),
                  cex = 0.7,
                  line = 2,
                  side = 1)
  graphics::mtext(paste0("t",
                         show_pred.vi[2],
                         " (",
                         round(pca_metrics.ls[["relative_var.vn"]][show_pred.vi[2]] * 100), "%)"),
                  cex = 0.7,
                  las = 0,
                  line = 2,
                  side = 2)
  graphics::abline(h = 0, lty = "dashed")
  graphics::abline(v = 0, lty = "dashed")
  radVn <- seq(0, 2 * pi, length.out = 100)
  
  hotFisN <- pca_metrics.ls[["hotelling_df.i"]] * stats::qf(0.95, 2, nrow(data.mn) - 2)
  
  graphics::lines(sqrt(stats::var(score_pca.mn[, show_pred.vi[1]]) * hotFisN) * cos(radVn),
                  sqrt(stats::var(score_pca.mn[, show_pred.vi[2]]) * hotFisN) * sin(radVn))
  
  if (identical(show_pred.vi, c(1, 2))) {
    
    hotVn <- pca_metrics.ls[["hotelling_pval.vn"]]
    
  } else if (identical(show_pred.vi, c(3, 4))) {
    
    invCovSco34MN <- solve(stats::cov(score_pca.mn[, 3:4]))
    
    hotVn <- apply(score_pca.mn[, 3:4],
                   1,
                   function(x)
                     1 - stats::pf(1 / pca_metrics.ls[["hotelling_df.i"]] * t(as.matrix(x)) %*% invCovSco34MN %*% as.matrix(x), 2, nrow(data.mn) - 2))
    
  } else
    stop("'show_pred.vi' must be either 'c(1, 2)' or 'c(3, 4)'")
  
  obsHotVi <- which(hotVn < 0.05)
  
  if (labels.l) {
    if (length(obsHotVi))
      graphics::text(score_pca.mn[obsHotVi, show_pred.vi[1]],
                     score_pca.mn[obsHotVi, show_pred.vi[2]],
                     cex = 0.7,
                     col = sample_color.vc[obsHotVi],
                     labels = rownames(data.mn)[obsHotVi])
    graphics::text(score_pca.mn[setdiff(1:nrow(score_pca.mn), obsHotVi), show_pred.vi[1]],
                     score_pca.mn[setdiff(1:nrow(score_pca.mn), obsHotVi), show_pred.vi[2]],
                     col = sample_color.vc[setdiff(1:nrow(score_pca.mn), obsHotVi)],
                     labels = sample_label.vc[setdiff(1:nrow(score_pca.mn), obsHotVi)])
  } else
    graphics::text(score_pca.mn[, show_pred.vi[1]],
                     score_pca.mn[, show_pred.vi[2]],
                     col = sample_color.vc,
                     labels = sample_label.vc)
  
}


.sample_color <- function(samp.df,
                          col_sampleType.c = "sampleType") {
  
  if (col_sampleType.c %in% colnames(samp.df)) {
    
    sample_types.vc <- samp.df[, col_sampleType.c]
    
  } else {
    
    sample_types.vc <- rep("other", nrow(samp.df))
    
  }
  
  .sample_palette(sample_types.vc = sample_types.vc)
  
}

.sample_palette <- function(sample_types.vc) {
  
  type_colors.vc <- c(sample = "green4",
                      pool = RColorBrewer::brewer.pal(9, "Reds")[7],
                      pool1 = RColorBrewer::brewer.pal(9, "Reds")[5],
                      pool2 = RColorBrewer::brewer.pal(9, "Reds")[4],
                      pool4 = RColorBrewer::brewer.pal(9, "Reds")[3],
                      pool8 = RColorBrewer::brewer.pal(9, "Reds")[2],
                      blank = "black",
                      other = RColorBrewer::brewer.pal(9, "Blues")[7])
  
  sample_types.vc[!(sample_types.vc %in% names(type_colors.vc))] <- "other"
  
  type_colors.vc[sample_types.vc]
  
}


.allDigits <- function(string) { ## from the Hmisc package (all.digits)
  k <- length(string)
  result <- logical(k)
  for (i in 1:k) {
    st <- string[i]
    ls <- nchar(st)
    ex <- substring(st, 1:ls, 1:ls)
    result[i] <- all(match(ex, c("0", "1", "2", "3", "4",
                                 "5", "6", "7", "8", "9"), nomatch = 0) > 0)
  }
  result
}


.zscore <- function(x) {
  sdxN <- stats::sd(x, na.rm = TRUE)
  if (sdxN < .Machine[["double.eps"]]) {
    return(rep(0, length(x)))
  } else
    return((x - mean(x, na.rm = TRUE)) / sdxN)
}
