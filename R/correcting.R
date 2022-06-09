## correcting (MultiAssayExperiment) ----

#' @rdname correcting
#' @export
setMethod("correcting", signature(x = "MultiAssayExperiment"),
          function(x,
                   method.vc = c("loess", "serrf")[1],
                   reference.vc = c("pool", "sample")[1],
                   loess_span.vn = 1,
                   serrf_corvar.vi = 10,
                   sample_intensity.c = c("median", "mean", "sum")[2],
                   title.c = NA,
                   col_batch.c = "batch",
                   col_injectionOrder.c = "injectionOrder",
                   col_sampleType.c = "sampleType",
                   figure.c = c("none", "interactive", "myfile.pdf")[2],
                   report.c = c("none", "interactive", "myfile.txt")[2]) {
            
            if (!(length(method.vc) %in% c(1, length(x)))) {
              stop("'The length of 'method.vc' should either be 1 or equal to the number of datasets")
            } else if (length(method.vc) == 1) {
              method.vc <- rep(method.vc, length(x))
            }
            names(method.vc) <- names(x)
            
            if (!(length(reference.vc) %in% c(1, length(x)))) {
              stop("'The length of 'reference.vc' should either be 1 or equal to the number of datasets")
            } else if (length(reference.vc) == 1) {
              reference.vc <- rep(reference.vc, length(x))
            }
            names(reference.vc) <- names(x)
            
            if (!(length(loess_span.vn) %in% c(1, length(x)))) {
              stop("'The length of 'loess_span.vn' should either be 1 or equal to the number of datasets")
            } else if (length(loess_span.vn) == 1) {
              loess_span.vn <- rep(loess_span.vn, length(x))
            }
            names(loess_span.vn) <- names(x)
            
            if (!(length(serrf_corvar.vi) %in% c(1, length(x)))) {
              stop("'The length of 'serrf_corvar.vi' should either be 1 or equal to the number of datasets")
            } else if (length(serrf_corvar.vi) == 1) {
              serrf_corvar.vi <- rep(serrf_corvar.vi, length(x))
            }
            names(serrf_corvar.vi) <- names(x)
            
            if (!(report.c %in% c("none", "interactive")))
              sink(report.c, append = TRUE)
            
            report_set.c <- report.c
            if (report_set.c != "none")
              report_set.c <- "interactive"
            
            if (!(figure.c %in% c("none", "interactive")))
              grDevices::pdf(figure.c, width = 11, height = 7)
            
            figure_set.c <- figure.c
            if (figure_set.c != "none")
              figure_set.c <- "interactive"
            
            for (set.c in names(x)) {
              
              if (report.c != "none")
                message("Correcting the '", set.c, "' dataset:")
              
              x[[set.c]] <- correcting(x = x[[set.c]],
                                       method.vc = method.vc[[set.c]],
                                       reference.vc = reference.vc[[set.c]],
                                       loess_span.vn = loess_span.vn[[set.c]],
                                       serrf_corvar.vi = serrf_corvar.vi[[set.c]],
                                       col_batch.c = col_batch.c,
                                       col_injectionOrder.c = col_injectionOrder.c,
                                       col_sampleType.c = col_sampleType.c,
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



## correcting (SummarizedExperiment) ----

#' @rdname correcting
#' @export
setMethod("correcting", signature(x = "SummarizedExperiment"),
          function(x,
                   method.vc = c("loess", "serrf")[1],
                   reference.vc = c("pool", "sample")[1],
                   loess_span.vn = 1,
                   serrf_corvar.vi = 10,
                   sample_intensity.c = c("median", "mean", "sum")[2],
                   title.c = NA,
                   col_batch.c = "batch",
                   col_injectionOrder.c = "injectionOrder",
                   col_sampleType.c = "sampleType",
                   figure.c = c("none", "interactive", "myfile.pdf")[2],
                   report.c = c("none", "interactive", "myfile.txt")[2]) {
            
            if (length(method.vc) != 1)
              stop("'method.vc' should be of length 1 for an 'SummarizedExperiment'")
            
            if (length(reference.vc) != 1)
              stop("'reference.vc' should be of length 1 for an 'SummarizedExperiment'")
            
            if (length(loess_span.vn) != 1)
              stop("'loess_span.vn' should be of length 1 for an 'SummarizedExperiment'")
            
            if (length(serrf_corvar.vi) != 1)
              stop("'serrf_corvar.vi' should be of length 1 for an 'SummarizedExperiment'")
            
            if (!(report.c %in% c("none", "interactive")))
              sink(report.c, append = TRUE)
            
            if (is.na(title.c))
              title.c <- x@metadata$experimentData@title
            
            norm.mn <- .correcting(data.mn = t(SummarizedExperiment::assay(x)),
                                   samp.df = SummarizedExperiment::colData(x),
                                   method.c = method.vc,
                                   reference.c = reference.vc,
                                   loess_span.n = loess_span.vn,
                                   serrf_corvar.i = serrf_corvar.vi,
                                   col_batch.c = col_batch.c,
                                   col_injectionOrder.c = col_injectionOrder.c,
                                   col_sampleType.c = col_sampleType.c,
                                   sample_intensity.c = sample_intensity.c,
                                   title.c = title.c,
                                   figure.c = figure.c)
            
            SummarizedExperiment::assay(x) <- t(norm.mn)
            
            if (!(report.c %in% c("none", "interactive")))
              sink()
            
            methods::validObject(x)
            
            return(invisible(x))
            
          })


## correcting (MultiDataSet) ----

#' @rdname correcting
#' @export
setMethod("correcting", signature(x = "MultiDataSet"),
          function(x,
                   method.vc = c("loess", "serrf")[1],
                   reference.vc = c("pool", "sample")[1],
                   loess_span.vn = 1,
                   serrf_corvar.vi = 10,
                   sample_intensity.c = c("median", "mean", "sum")[2],
                   title.c = NA,
                   col_batch.c = "batch",
                   col_injectionOrder.c = "injectionOrder",
                   col_sampleType.c = "sampleType",
                   figure.c = c("none", "interactive", "myfile.pdf")[2],
                   report.c = c("none", "interactive", "myfile.txt")[2]) {
            
            if (!(length(method.vc) %in% c(1, length(x)))) {
              stop("'The length of 'method.vc' should either be 1 or equal to the number of datasets")
            } else if (length(method.vc) == 1) {
              method.vc <- rep(method.vc, length(x))
            }
            names(method.vc) <- names(x)
            
            if (!(length(reference.vc) %in% c(1, length(x)))) {
              stop("'The length of 'reference.vc' should either be 1 or equal to the number of datasets")
            } else if (length(reference.vc) == 1) {
              reference.vc <- rep(reference.vc, length(x))
            }
            names(reference.vc) <- names(x)
            
            if (!(length(loess_span.vn) %in% c(1, length(x)))) {
              stop("'The length of 'loess_span.vn' should either be 1 or equal to the number of datasets")
            } else if (length(loess_span.vn) == 1) {
              loess_span.vn <- rep(loess_span.vn, length(x))
            }
            names(loess_span.vn) <- names(x)
            
            if (!(length(serrf_corvar.vi) %in% c(1, length(x)))) {
              stop("'The length of 'serrf_corvar.vi' should either be 1 or equal to the number of datasets")
            } else if (length(serrf_corvar.vi) == 1) {
              serrf_corvar.vi <- rep(serrf_corvar.vi, length(x))
            }
            names(serrf_corvar.vi) <- names(x)
            
            if (!(report.c %in% c("none", "interactive")))
              sink(report.c, append = TRUE)
            
            report_set.c <- report.c
            if (report_set.c != "none")
              report_set.c <- "interactive"
            
            if (!(figure.c %in% c("none", "interactive")))
              grDevices::pdf(figure.c, width = 11, height = 7)
            
            figure_set.c <- figure.c
            if (figure_set.c != "none")
              figure_set.c <- "interactive"
            
            for (set.c in names(x)) {
              
              if (report.c != "none")
                message("Correcting the '", set.c, "' dataset:")
              
              ese <- correcting(x = x[[set.c]],
                                method.vc = method.vc[[set.c]],
                                reference.vc = reference.vc[[set.c]],
                                loess_span.vn = loess_span.vn[[set.c]],
                                serrf_corvar.vi = serrf_corvar.vi[[set.c]],
                                col_batch.c = col_batch.c,
                                col_injectionOrder.c = col_injectionOrder.c,
                                col_sampleType.c = col_sampleType.c,
                                sample_intensity.c = sample_intensity.c,
                                title.c = set.c,
                                figure.c = figure_set.c,
                                report.c = report_set.c)
              
              x <- MultiDataSet::add_eset(x,
                                          ese,
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


## correcting (ExpressionSet) ----

#' @rdname correcting
#' @export
setMethod("correcting", signature(x = "ExpressionSet"),
          function(x,
                   method.vc = c("loess", "serrf")[1],
                   reference.vc = c("pool", "sample")[1],
                   loess_span.vn = 1,
                   serrf_corvar.vi = 10,
                   sample_intensity.c = c("median", "mean", "sum")[2],
                   title.c = NA,
                   col_batch.c = "batch",
                   col_injectionOrder.c = "injectionOrder",
                   col_sampleType.c = "sampleType",
                   figure.c = c("none", "interactive", "myfile.pdf")[2],
                   report.c = c("none", "interactive", "myfile.txt")[2]) {
            
            if (length(method.vc) != 1)
              stop("'method.vc' should be of length 1 for an 'ExpressionSet'")
            
            if (length(reference.vc) != 1)
              stop("'reference.vc' should be of length 1 for an 'ExpressionSet'")
            
            if (length(loess_span.vn) != 1)
              stop("'loess_span.vn' should be of length 1 for an 'ExpressionSet'")
            
            if (length(serrf_corvar.vi) != 1)
              stop("'serrf_corvar.vi' should be of length 1 for an 'ExpressionSet'")
            
            if (!(report.c %in% c("none", "interactive")))
              sink(report.c, append = TRUE)
            
            if (is.na(title.c))
              title.c <- Biobase::experimentData(x)@title
            
            norm.mn <- .correcting(data.mn = t(Biobase::exprs(x)),
                                   samp.df = Biobase::pData(x),
                                   method.c = method.vc,
                                   reference.c = reference.vc,
                                   loess_span.n = loess_span.vn,
                                   serrf_corvar.i = serrf_corvar.vi,
                                   col_batch.c = col_batch.c,
                                   col_injectionOrder.c = col_injectionOrder.c,
                                   col_sampleType.c = col_sampleType.c,
                                   sample_intensity.c = sample_intensity.c,
                                   title.c = title.c,
                                   figure.c = figure.c)
            
            Biobase::exprs(x) <- t(norm.mn)
            
            if (!(report.c %in% c("none", "interactive")))
              sink()
            
            methods::validObject(x)
            
            return(invisible(x))
            
          })

.correcting <- function(data.mn, ## data (matrix of numerics; samples x variables)
                        samp.df, ## sample metadata (data frame; samples x sample metadata)
                        method.c,
                        reference.c,
                        loess_span.n,
                        serrf_corvar.i,
                        title.c,
                        figure.c,
                        col_batch.c,
                        col_injectionOrder.c,
                        col_sampleType.c,
                        sample_intensity.c) {
 
  ## checking
  
  if (method.c == "serrf" && reference.c != "pool")
    stop("'reference' should be set to 'pool' for the serrf method.")
  
  if (sum(grepl(reference.c, samp.df[, col_sampleType.c])) == 0)
    stop("No '", reference.c, "' reference sample type found in the '",
         col_sampleType.c, "' column of the sampleMetadata.")
  
  ref_data.mn <- data.mn[samp.df[, col_sampleType.c] == reference.c, ]
  ref_samp.df <- samp.df[samp.df[, col_sampleType.c] == reference.c, ]
  
  ref_nazeros.vl <- apply(ref_data.mn, 2,
                          function(ref.vn)
                            all(sapply(ref.vn,
                                       function(ref.n) {is.na(ref.n) || ref.n == 0})))
  
  if (sum(ref_nazeros.vl)) {
    message(sum(ref_nazeros.vl), " features with 'NA' or 0 values in all reference samples removed from the data.")
    data.mn <- ref_data.mn[, !ref_nazeros.vl]
  }
  
  ## Computation
  
  ## ordering (batch and injection order)
  
  samp.df[, "initial_order"] <- 1:nrow(samp.df)
  order_batch_inj.vi <- order(samp.df[, col_batch.c],
                              samp.df[, col_injectionOrder.c])
  data.mn <- data.mn[order_batch_inj.vi, ]
  samp.df <- samp.df[order_batch_inj.vi, ]
  
  ## signal drift and batch-effect correction
  
  normalized.mn <- .batch_correct(data.mn = data.mn,
                                  samp.df = samp.df,
                                  method.c = method.c,
                                  reference.c = reference.c,
                                  loess_span.n = loess_span.n,
                                  serrf_corvar.i = serrf_corvar.i,
                                  col_batch.c = col_batch.c,
                                  col_sampleType.c = col_sampleType.c)
  
  ## figure
  
  if (figure.c != "none") {
    
    if (figure.c != "interactive")
      grDevices::pdf(figure.c)
    
    .plot_drift_pca(data.mn = data.mn,
                    samp.df = samp.df,
                    loess_span.n = loess_span.n,
                    sample_intensity.c = sample_intensity.c,
                    title.c = title.c,
                    raw_vs_normalized.c = "raw",
                    col_batch.c = col_batch.c,
                    col_injectionOrder.c = col_injectionOrder.c,
                    col_sampleType.c = col_sampleType.c)
    .plot_drift_pca(data.mn = normalized.mn,
                    samp.df = samp.df,
                    loess_span.n = loess_span.n,
                    sample_intensity.c = sample_intensity.c,
                    title.c = title.c,
                    raw_vs_normalized.c = method.c,
                    col_batch.c = col_batch.c,
                    col_injectionOrder.c = col_injectionOrder.c,
                    col_sampleType.c = col_sampleType.c)
    
    if (figure.c != "interactive")
      grDevices::dev.off()
    
  }
  
  ## returning to initial order
  
  initial_order.vi <- order(samp.df[, "initial_order"])
  normalized.mn <- normalized.mn[samp.df[, "initial_order"], ]
  samp.df <- samp.df[samp.df[, "initial_order"], ] # not used

  ## returning
  
  normalized.mn
  
}


.batch_correct <- function(data.mn,
                           samp.df,
                           method.c,
                           reference.c,
                           loess_span.n,
                           serrf_corvar.i,
                           col_batch.c,
                           col_sampleType.c) {
  
  message("Correction method: ", method.c)
  message("Reference observations: ", reference.c)
  
  ## computing means of all pools (or samples) for each variable (medians used in Fan et al., 2019)
  
  ref_mean.vn <- apply(data.mn[samp.df[, col_sampleType.c] == reference.c, ], 2,
                       function(feat.vn) {
                         mean(feat.vn, na.rm = TRUE)
                       })
  ref_median.vn <- apply(data.mn[samp.df[, col_sampleType.c] == reference.c, ], 2,
                       function(feat.vn) {
                         median(feat.vn, na.rm = TRUE)
                       })
  pred_median.vn <- apply(data.mn[samp.df[, col_sampleType.c] != reference.c, ], 2,
                          function(feat.vn) {
                            median(feat.vn, na.rm = TRUE)
                          })
  
  ## splitting data and sample metadata from each batch
  
  batch_data.ls <- split(as.data.frame(data.mn),
                         f = samp.df[, col_batch.c])
  batch_data.ls <- lapply(batch_data.ls, function(input.df) as.matrix(input.df))
  
  batch_samp.ls <- split(as.data.frame(samp.df),
                         f = samp.df[, col_batch.c])
  
  ## checking extrapolation: are there pools at the first and last observations of each batch
  
  pool_extra.ml <- matrix(FALSE, nrow = 2, ncol = length(batch_data.ls),
                          dimnames = list(c("first", "last"), names(batch_data.ls)))
  
  for (batch.c in names(batch_samp.ls)) {
    batch_sampleType.vc <- batch_samp.ls[[batch.c]][, col_sampleType.c]
    pool_extra.ml["first", batch.c] <- utils::head(batch_sampleType.vc, 1) == reference.c
    pool_extra.ml["last", batch.c] <- utils::tail(batch_sampleType.vc, 1) == reference.c
  }
  
  if (!all(c(pool_extra.ml))) {
    warnings("Reference samples are missing at the first and/or last position of the following batches:\n")
    pool_extra_batch.vi <- which(!apply(pool_extra.ml, 2, all))
    for (i in 1:length(pool_extra_batch.vi))
      message(names(pool_extra_batch.vi)[i], ": ",
              paste(rownames(pool_extra.ml)[!pool_extra.ml[, pool_extra_batch.vi[i]]], collapse = ", "))
    message("Extrapolating loess fits for these batches may result in inaccurate modeling!")
  }
  
  ## normalizing
  
  normalized.mn <- NULL ## normalized data matrix to be computed
  
  message("Processing batch:")
  
  for (batch.c in names(batch_data.ls)) { ## processing each batch individually
    
    message(batch.c)
    
    batch_data.mn <- batch_data.ls[[batch.c]]
    batch_samp.df <- batch_samp.ls[[batch.c]]
    
    batch_all.vi <- 1:nrow(batch_data.mn)
    batch_ref.vi <- which(batch_samp.df[, col_sampleType.c] == reference.c)
    
    
    if (method.c == "loess") {
      
      if (length(batch_ref.vi) < 5)
        message("less than 5 '", reference.c,
                "'; linear regression will be performed instead of loess regression for this batch.")
      
      batch_pred.mn <- .loess_pred(data.mn = batch_data.mn,
                                   ref.vi = batch_ref.vi,
                                   all.vi = batch_all.vi,
                                   span.n = loess_span.n)
      
      ## normalization
      
      batch_pred.mn[batch_pred.mn <= 0] <- NA
      
      batch_normalized.mn <- batch_data.mn / batch_pred.mn
      
      
    } else if (method.c == "serrf") {
      
      batch_normalized.mn <- .serrf_pred(data.mn = batch_data.mn,
                                         ref.vi = batch_ref.vi,
                                         all.vi = batch_all.vi,
                                         corvar.i = serrf_corvar.i,
                                         ref_mean.vn = ref_mean.vn,
                                         ref_median.vn = ref_median.vn,
                                         pred_median.vn = pred_median.vn)
    }  else
      stop("'method.c' must be either 'loess' or 'serrf'")
    
    
    normalized.mn <- rbind(normalized.mn,
                           batch_normalized.mn)
    
  }
  
  
  if (method.c == "loess")
    normalized.mn <- sweep(normalized.mn, MARGIN = 2, STATS = ref_mean.vn, FUN = "*")
  
  
  
  return(normalized.mn)
  
} ## batch_correct


.loess_pred <- function(data.mn,
                        reference.c,
                        ref.vi,
                        all.vi,
                        span.n) {
  
  
  ## prediction of the loess fit
  
  apply(data.mn, 2,
        function(rawVn)
          .loess(raw.vn = rawVn,
                 ref.vi = ref.vi,
                 pred.vi = all.vi,
                 span.n = span.n))
  
}


.loess <- function(raw.vn, ref.vi, pred.vi, span.n) {
  
  if (length(ref.vi) < 5) {
    
    return(stats::predict(stats::lm(raw.vn[ref.vi] ~ ref.vi),
                          newdata = data.frame(ref.vi = pred.vi)))
    
  } else {
    
    return(stats::predict(stats::loess(raw.vn[ref.vi] ~ ref.vi,
                                       control = stats::loess.control(surface = "direct"),
                                       span = span.n),
                          newdata = data.frame(ref.vi = pred.vi)))
    
  }
  
  ## Note:
  ##  the surface = 'direct' argument allows extrapolation
  
}


.serrf_pred <- function(data.mn,
                        ref.vi,
                        all.vi,
                        corvar.i,
                        ref_mean.vn,
                        ref_median.vn,
                        pred_median.vn) {
  # Fan et al. (2019). Systematic Error Removal Using Random Forest for Normalizing Large-Scale Untargeted Lipidomics Data, Anal. Chem., vol. 91, nᵒ 5, p. 3590‑3596, doi: 10.1021/acs.analchem.8b05592.
  
  pred.vi <- setdiff(all.vi, ref.vi)
  
  ref.mn <- data.mn[ref.vi, ]
  pred.mn <- data.mn[pred.vi, ]
  
  ref_scale.mn <- scale(ref.mn)
  pred_scale.mn <- scale(pred.mn)
  
  ref_cor.mn <- cor(ref_scale.mn, method = "spearman")
  pred_cor.mn <- cor(pred_scale.mn, method = "spearman")
  
  ## Replacing 0 and NA values (Fan et al., 2019)
  
  data.mn <- apply(data.mn, 2,
                   function(feat.vn) {
                     zero.vl <- feat.vn == 0
                     na.vl <- is.na(feat.vn)
                     feat.vn[zero.vl] <- rnorm(length(zero.vl),
                                               mean = min(feat.vn[!na.vl]) + 1,
                                               sd = 0.1 * min(feat.vn[!na.vl]) + 0.1)
                     feat.vn[na.vl] <- rnorm(length(na.vl),
                                             mean = 0.5 * min(feat.vn[!na.vl]) + 1,
                                             sd = 0.1 * min(feat.vn[!na.vl]) + 0.1)
                     return(feat.vn)
                   })
  
  serrf.mn <- data.mn
  
  # normalizing each j variable successively
  for (j in 1:ncol(data.mn)) {
    
    # computing the corvar.i closest features to j
    ref_cor_ord.vi <- order(abs(ref_cor.mn[, j]), decreasing = TRUE)
    pred_cor_ord.vi <- order(abs(pred_cor.mn[, j]), decreasing = TRUE)
    
    corvar_j.vi <- integer()
    length.i <- corvar.i
    while (length(corvar_j.vi) < corvar.i) {
      corvar_j.vi = intersect(ref_cor_ord.vi[1:length.i], pred_cor_ord.vi[1:length.i])
      corvar_j.vi = corvar_j.vi[corvar_j.vi != j]
      length.i = length.i + 1
    }
    
    # restricting to these features
    train_x.mn <- ref_scale.mn[, corvar_j.vi]
    train_y.vn <- scale(ref.mn[, j], scale = FALSE)
    
    test_x.mn <- pred_scale.mn[, corvar_j.vi]
    
    # discarding the features with NA only
    nomissing.vl <- apply(rbind(train_x.mn,
                                test_x.mn), 2, function(feat.vn) sum(is.na(feat.vn)) == 0)
    
    train_x.mn <- train_x.mn[, nomissing.vl, drop = FALSE]
    test_x.mn <- test_x.mn[, nomissing.vl, drop = FALSE]
    
    # random forest prediction
    train.df = data.frame(y = train_y.vn, train_x.mn)
    colnames(train.df) = c("y", paste0("V",1:(ncol(train.df)-1)))
    
    model.rf = ranger::ranger(y~., data = train.df,
                              seed = 123)
    
    # predictions
    test.df = data.frame(test_x.mn)
    colnames(test.df) = colnames(train.df)[-1]
    
    serrf_refj.vn <- data.mn[ref.vi, j] / (predict(model.rf, data = train.df)[["predictions"]] + attributes(ref_scale.mn)[["scaled:center"]][j]) * ref_mean.vn[j]
    serrf_prej.vn <- data.mn[pred.vi, j] / (predict(model.rf, data = test.df)[["predictions"]] + attributes(pred_scale.mn)[["scaled:center"]][j]) * pred_median.vn[j]
    serrf_prej.vn[serrf_prej.vn < 0] <- data.mn[pred.vi, j][serrf_prej.vn < 0]
    
    
    serrf_refj.vn <- serrf_refj.vn / median(serrf_refj.vn, na.rm = TRUE) * ref_median.vn[j]
    serrf_prej.vn <- serrf_prej.vn / median(serrf_prej.vn, na.rm = TRUE) * pred_median.vn[j]
    
    serrf.mn[ref.vi, j] <- serrf_refj.vn
    serrf.mn[pred.vi, j] <- serrf_prej.vn
    isfinite.vl <- is.finite(serrf.mn[, j])
    if (sum(!isfinite.vl)) {
      serrf.mn[!isfinite.vl, j] <- rnorm(sum(!isfinite.vl, na.rm = TRUE), sd = sd(serrf.mn[isfinite.vl, j], na.rm = TRUE) * 0.01)
    }
    
    out = boxplot.stats(serrf.mn[, j], coef = 3)$out
    
    serrf.mn[pred.vi, j][serrf.mn[pred.vi, j] %in% out] <- (data.mn[pred.vi, j] - ((predict(model.rf, data = test.df)[["predictions"]] + attributes(pred_scale.mn)[["scaled:center"]][j]) - pred_median.vn[j]))[serrf.mn[pred.vi, j] %in% out]
    serrf.mn[pred.vi, j][serrf.mn[pred.vi, j] < 0] <- data.mn[pred.vi, j][serrf.mn[pred.vi, j] < 0]
    
  }
  
  .fill_naneg <- function(feat.vn) {
    isna.vl <- is.na(feat.vn)
    feat.vn[isna.vl] <- rnorm(sum(isna.vl), mean = min(feat.vn[!isna.vl]), sd = sd(feat.vn[!isna.vl]) * 0.1)
    isneg.vl <- feat.vn < 0
    feat.vn[isneg.vl] <- runif(1) * min(feat.vn[feat.vn > 0])
    return(feat.vn)
  }
  
  serrf.mn[ref.vi, ] <- apply(serrf.mn[ref.vi, ], 2, .fill_naneg)
  serrf.mn[pred.vi, ] <- apply(serrf.mn[pred.vi, ], 2, .fill_naneg)
 
  return(serrf.mn)
  
}
  
  
.plot_drift_pca <- function(data.mn,
                            samp.df,
                            loess_span.n = 1,
                            sample_intensity.c = "mean",
                            title.c = NA,
                            raw_vs_normalized.c = "",
                            col_batch.c = "batch",
                            col_injectionOrder.c = "injectionOrder",
                            col_sampleType.c = "sampleType") {
  
  main.c <- paste0(raw_vs_normalized.c, 
                   ifelse(!is.na(title.c) && title.c != "",
                          paste0(" (", title.c, ")"),
                          ""))
  
  graphics::par(font = 2, font.axis = 2, font.lab = 2, lwd = 2, pch = 18)
  
  graphics::layout(matrix(c(1, 1, 2, 3), nrow = 2),
                   widths = c(0.7, 0.3))
  
  obsNamVc <- rownames(samp.df)
  
  obsColVc <- .sample_color(samp.df = samp.df, col_sampleType.c = col_sampleType.c)
  
  ## Graphic 1: Mean of intensities for each sample
  
  .plot_drift(data.mn = data.mn,
              samp.df = samp.df,
              loess_span.n = loess_span.n,
              sample_intensity.c = sample_intensity.c,
              mar.vn = c(3.6, 3.6, 3.1, 0.6),
              col_batch.c = col_batch.c,
              col_injectionOrder.c = col_injectionOrder.c,
              col_sampleType.c = col_sampleType.c)
  
  title(main.c)
  
  ## Graphics 2 and 3 (right): PCA score plots of components 1-4
  
  pca_metrics.ls <- .pca_metrics(data.mn = data.mn,
                                 samp.df = samp.df,
                                 pred.i = 4)
 
  .plot_pca_metrics(data.mn = data.mn,
                    samp.df = samp.df,
                    pred.i = 4,
                    show_pred.vi = c(1, 2),
                    pca_metrics.ls = pca_metrics.ls,
                    col_sampleType.c = "sampleType",
                    mar.vn = c(3.6, 3.6, 0.6, 1.1))
  
  .plot_pca_metrics(data.mn = data.mn,
                    samp.df = samp.df,
                    pred.i = 4,
                    show_pred.vi = c(3, 4),
                    pca_metrics.ls = pca_metrics.ls,
                    col_sampleType.c = "sampleType",
                    mar.vn = c(3.6, 3.6, 0.6, 1.1))
  
  
} ## plot_batch
