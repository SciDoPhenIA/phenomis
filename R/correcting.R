## correcting (MultiDataSet) ----

#' @rdname correcting
#' @export
setMethod("correcting", signature(x = "MultiDataSet"),
          function(x,
                   reference.c = c("pool", "sample")[1],
                   span.n = 1,
                   sample_intensity.c = c("median", "mean", "sum")[2],
                   title.c = NA,
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
              grDevices::pdf(figure.c, width = 11, height = 7)
            
            figure_set.c <- figure.c
            if (figure_set.c != "none")
              figure_set.c <- "interactive"
            
            for (set.c in names(x)) {
              
              if (report.c != "none")
                message("Correcting the '", set.c, "' dataset:")
              
              ese <- phenomis::correcting(x = x[[set.c]],
                                          reference.c = reference.c,
                                          span.n = span.n,
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
            
            return(invisible(x))
            
          })


## correcting (ExpressionSet) ----

#' @rdname correcting
#' @export
setMethod("correcting", signature(x = "ExpressionSet"),
          function(x,
                   reference.c = c("pool", "sample")[1],
                   span.n = 1,
                   sample_intensity.c = c("median", "mean", "sum")[2],
                   title.c = NA,
                   col_batch.c = "batch",
                   col_injectionOrder.c = "injectionOrder",
                   col_sampleType.c = "sampleType",
                   figure.c = c("none", "interactive", "myfile.pdf")[2],
                   report.c = c("none", "interactive", "myfile.txt")[2]) {
            
            if (!(report.c %in% c("none", "interactive")))
              sink(report.c, append = TRUE)
            
            if (is.na(title.c))
              title.c <- Biobase::experimentData(x)@title
            
            x <- .correcting(eset = x,
                                        reference.c = reference.c,
                                        span.n = span.n,
                                        col_batch.c = col_batch.c,
                                        col_injectionOrder.c = col_injectionOrder.c,
                                        col_sampleType.c = col_sampleType.c,
                                        sample_intensity.c = sample_intensity.c,
                                        title.c = title.c,
                                        figure.c = figure.c)
            
            methods::validObject(x)
            
            if (!(report.c %in% c("none", "interactive")))
              sink()
            
            return(invisible(x))
            
          })

.correcting <- function(eset,
                        reference.c,
                        span.n,
                        title.c,
                        figure.c,
                        col_batch.c,
                        col_injectionOrder.c,
                        col_sampleType.c,
                        sample_intensity.c) {
  
  ## checking
  
  if (sum(grepl(reference.c, Biobase::pData(eset)[, col_sampleType.c])) == 0)
    stop("No '", reference.c, "' reference sample type found in the '",
         col_sampleType.c, "' column of the sampleMetadata.")
  
  ref.eset <- eset[, Biobase::pData(eset)[, col_sampleType.c] == reference.c]
  
  ref_nazeros.vl <- apply(Biobase::exprs(ref.eset), 1,
                       function(refVn)
                         all(sapply(refVn,
                                    function(refN) {is.na(refN) || refN == 0})))
  
  if (sum(ref_nazeros.vl)) {
    message(sum(ref_nazeros.vl), " variables with 'NA' or 0 values in all reference sample will be removed from the data.")
    eset <- eset[!ref_nazeros.vl, ]
  }
  
  ## Computation
  
  ## ordering (batch and injection order)
  
  pdata.df <- Biobase::pData(eset)
  pdata.df[, "initial_order"] <- 1:nrow(pdata.df)
  order_batch_inj.vi <- order(pdata.df[, col_batch.c],
                       pdata.df[, col_injectionOrder.c])
  Biobase::pData(eset) <- pdata.df
  eset <- eset[, order_batch_inj.vi]
  
  ## signal drift and batch-effect correction
  
  normalized.eset <- .batch_correct(eset = eset,
                                               reference.c = reference.c,
                                               span.n = span.n,
                                               col_batch.c = col_batch.c,
                                               col_sampleType.c = col_sampleType.c)
  
  ## figure
  
  if (figure.c != "none") {
    
    if (figure.c != "interactive")
      grDevices::pdf(figure.c)
    
    .plot_batch(eset = eset,
                span.n = span.n,
                sample_intensity.c = sample_intensity.c,
                title.c = title.c,
                raw_vs_normalized.c = "raw",
                col_batch.c = col_batch.c,
                col_injectionOrder.c = col_injectionOrder.c,
                col_sampleType.c = col_sampleType.c)
    .plot_batch(eset = normalized.eset,
                span.n = span.n,
                sample_intensity.c = sample_intensity.c,
                title.c = title.c,
                raw_vs_normalized.c = "normalized",
                col_batch.c = col_batch.c,
                col_injectionOrder.c = col_injectionOrder.c,
                col_sampleType.c = col_sampleType.c)
    
    if (figure.c != "interactive")
      grDevices::dev.off()
    
  }
  
  ## returning to initial order
  
  initial_order.vi <- order(Biobase::pData(normalized.eset)[, "initial_order"])
  normalized.eset <- normalized.eset[, initial_order.vi]
  pdata.df <- Biobase::pData(normalized.eset)
  pdata.df[, "initial_order"] <- NULL
  Biobase::pData(normalized.eset) <- pdata.df

  ## returning
  
  normalized.eset
  
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

.plot_batch <- function(eset,
                        span.n = 1,
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
  
  obsNamVc <- Biobase::sampleNames(eset)
  
  obsColVc <- .sample_color_eset(eset = eset, col_sampleType.c = col_sampleType.c)

  ## Graphic 1: Mean of intensities for each sample
  
  .plot_drift(eset = eset,
              span.n = span.n,
              sample_intensity.c = sample_intensity.c,
              mar.vn = c(3.6, 3.6, 3.1, 0.6),
              col_batch.c = col_batch.c,
              col_injectionOrder.c = col_injectionOrder.c,
              col_sampleType.c = col_sampleType.c)
  
  title(main.c)

  ## Graphics 2 and 3 (right): PCA score plots of components 1-4
  
  pca_metrics.ls <- .pca_metrics(eset = eset, pred.i = 4)
  
  .plot_pca_metrics(eset = eset,
                               pred.i = 4,
                               show_pred.vi = c(1, 2),
                               pca_metrics.ls = pca_metrics.ls,
                               col_sampleType.c = "sampleType",
                               mar.vn = c(3.6, 3.6, 0.6, 1.1))
  
  .plot_pca_metrics(eset = eset,
                               pred.i = 4,
                               show_pred.vi = c(3, 4),
                               pca_metrics.ls = pca_metrics.ls,
                               col_sampleType.c = "sampleType",
                               mar.vn = c(3.6, 3.6, 0.6, 1.1))
 
  
} ## plot_batch

.batch_correct <- function(eset,
                           reference.c,
                           span.n,
                           col_batch.c,
                           col_sampleType.c) {
  
  message("Reference observations are: ", reference.c)
  
  ## computing means of all pools (or samples) for each variable
  
  ref_mean.vn <- rowMeans(Biobase::exprs(eset)[, Biobase::pData(eset)[, col_sampleType.c] == reference.c],
                          na.rm = TRUE)
  
  ## splitting data and sample metadata from each batch
  
  batch_texprs.ls <- split(as.data.frame(t(Biobase::exprs(eset))),
                           f = Biobase::pData(eset)[, col_batch.c])
  batch_texprs.ls <- lapply(batch_texprs.ls, function(input.df) as.matrix(input.df))
  
  batch_pdata.ls <- split(as.data.frame(Biobase::pData(eset)),
                          f = Biobase::pData(eset)[, col_batch.c])
  
  ## checking extrapolation: are there pools at the first and last observations of each batch
  
  pool_extra.ml <- matrix(FALSE, nrow = 2, ncol = length(batch_texprs.ls),
                          dimnames = list(c("first", "last"), names(batch_texprs.ls)))
  
  for (batch.c in names(batch_pdata.ls)) {
    batch_sampleType.vc <- batch_pdata.ls[[batch.c]][, col_sampleType.c]
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
  
  for (batch.c in names(batch_texprs.ls)) { ## processing each batch individually
    
    message(batch.c)
    
    batch_texprs.mn <- batch_texprs.ls[[batch.c]]
    batch_pdata.df <- batch_pdata.ls[[batch.c]]
    
    batch_all.vi <- 1:nrow(batch_texprs.mn)
    
    batch_ref.vi <- which(batch_pdata.df[, col_sampleType.c] == reference.c)
    
    if (length(batch_ref.vi) < 5)
      message("less than 5 '", reference.c,
              "'; linear regression will be performed instead of loess regression for this batch.")
    
    ## prediction of the loess fit
    
    batch_loess.mn <- apply(batch_texprs.mn, 2,
                            function(rawVn)
                              .loess(raw.vn = rawVn,
                                     ref.vi = batch_ref.vi,
                                     pred.vi = batch_all.vi,
                                     span.n = span.n))
    
    ## normalization
    
    batch_loess.mn[batch_loess.mn <= 0] <- NA
    
    batch_normalized.mn <- batch_texprs.mn / batch_loess.mn
    
    normalized.mn <- rbind(normalized.mn,
                           batch_normalized.mn)
    
  }
  
  normalized.mn <- sweep(normalized.mn, MARGIN = 2, STATS = ref_mean.vn, FUN = "*")
  
  normalized.eset <- eset
  stopifnot(identical(rownames(normalized.mn), Biobase::sampleNames(eset)))
  stopifnot(identical(colnames(normalized.mn), Biobase::featureNames(eset)))
  
  Biobase::exprs(normalized.eset) <- t(normalized.mn)
  
  return(normalized.eset)
  
} ## batch_correct
