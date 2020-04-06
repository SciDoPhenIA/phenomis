## correcting (MultiDataSet) ----

#' @rdname correcting
#' @export
setMethod("correcting", signature(x = "MultiDataSet"),
          function(x,
                   reference.c = c("pool", "sample")[1],
                   col_batch.c = "batch",
                   col_injectionOrder.c = "injectionOrder",
                   col_sampleType.c = "sampleType",
                   span.n = 1,
                   title.c = NA,
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
                                          col_batch.c = col_batch.c,
                                          col_injectionOrder.c = col_injectionOrder.c,
                                          col_sampleType.c = col_sampleType.c,
                                          span.n = span.n,
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
            
            return(x)
            
          })


## correcting (ExpressionSet) ----

#' @rdname correcting
#' @export
setMethod("correcting", signature(x = "ExpressionSet"),
          function(x,
                   reference.c = c("pool", "sample")[1],
                   col_batch.c = "batch",
                   col_injectionOrder.c = "injectionOrder",
                   col_sampleType.c = "sampleType",
                   span.n = 1,
                   title.c = NA,
                   figure.c = c("none", "interactive", "myfile.pdf")[2],
                   report.c = c("none", "interactive", "myfile.txt")[2]) {
            
            if (!(report.c %in% c("none", "interactive")))
              sink(report.c, append = TRUE)
            
            if (is.na(title.c))
              title.c <- Biobase::experimentData(x)@title
            
            x <- phenomis:::.correcting(eset = x,
                                        reference.c = reference.c,
                                        col_batch.c = col_batch.c,
                                        col_injectionOrder.c = col_injectionOrder.c,
                                        col_sampleType.c = col_sampleType.c,
                                        span.n = span.n,
                                        title.c = title.c,
                                        figure.c = figure.c)
            
            methods::validObject(x)
            
            if (!(report.c %in% c("none", "interactive")))
              sink()
            
            return(invisible(x))
            
          })

.correcting <- function(eset,
                        reference.c,
                        col_batch.c,
                        col_injectionOrder.c,
                        col_sampleType.c,
                        span.n,
                        title.c,
                        figure.c) {
  
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
  
  normalized.eset <- phenomis:::.batch_correct(eset = eset,
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
                col_batch.c = col_batch.c,
                col_sampleType.c = col_sampleType.c,
                title.c = title.c,
                type.c = "raw")
    .plot_batch(eset = normalized.eset,
                span.n = span.n,
                col_batch.c = col_batch.c,
                col_sampleType.c = col_sampleType.c,
                title.c = title.c,
                type.c = "normalized")
    
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
                        span.n,
                        col_batch.c,
                        col_sampleType.c,
                        title.c = NA,
                        type.c) {
  
  main.c <- paste0(type.c, 
                   ifelse(!is.na(title.c) && title.c != "",
                          paste0(" (", title.c, ")"),
                          ""))
  
  graphics::par(font = 2, font.axis = 2, font.lab = 2, lwd = 2, pch = 18)
  
  graphics::layout(matrix(c(1, 1, 2, 3), nrow = 2),
                   widths = c(0.7, 0.3))
  
  obsNamVc <- Biobase::sampleNames(eset)
  
  obsColVc <- .sample_color(Biobase::pData(eset)[, col_sampleType.c])

  ## Graphic 1: Median of intensities for each sample
  
  graphics::par(mar = c(3.6, 3.6, 3.1, 0.6))
  
  batch.tab <- table(Biobase::pData(eset)[, col_batch.c])
  
  median.vn <- apply(Biobase::exprs(eset), 2,
                     function(samp.vn) median(samp.vn, na.rm = TRUE))
  
  graphics::plot(median.vn,
                 cex = 1.2,
                 col = obsColVc,
                 pch = 18,
                 xaxs = "i",
                 xlab = "",
                 ylab = "")
  
  graphics::mtext("Injection order",
                  line = 2.2,
                  side = 1)
  graphics::mtext("Median of variable intensities",
                  line = 2.2,
                  side = 2)
  
  graphics::mtext(main.c, cex = 1.2, line = 1.5, side = 3)
  
  graphics::abline(v = cumsum(batch.tab) + 0.5,
                   col = "black")
  
  graphics::mtext(names(batch.tab),
                  at = batch.tab / 2 + c(0, cumsum(batch.tab[-length(batch.tab)])))
  
  obsColVuc <- obsColVc[sort(unique(names(obsColVc)))]
  
  graphics::text(rep(batch.tab[1], times = length(obsColVuc)),
                 graphics::par("usr")[3] + (0.97 - length(obsColVuc) * 0.03 + 1:length(obsColVuc) * 0.03) * diff(graphics::par("usr")[3:4]),
                 col = obsColVuc,
                 font = 2,
                 labels = names(obsColVuc),
                 pos = 2)
  
  for (batch.c in names(batch.tab)) {
    
    batch_seq.vi <- which(Biobase::pData(eset)[, col_batch.c] == batch.c)
    batch_pool.vi <- intersect(batch_seq.vi,
                               which(Biobase::pData(eset)[, col_sampleType.c] == "pool"))
    batch_samp.vi <- intersect(batch_seq.vi,
                               which(Biobase::pData(eset)[, col_sampleType.c] == "sample"))
    if (length(batch_pool.vi))
      graphics::lines(batch_seq.vi,
                      .loess(median.vn, batch_pool.vi, batch_seq.vi, span.n),
                      col = .sample_color("pool"))
    if (length(batch_samp.vi))
      graphics::lines(batch_seq.vi,
                      .loess(median.vn, batch_samp.vi, batch_seq.vi, span.n),
                      col = .sample_color("sample"))
    
  }
  
  ## Graphics 2 and 3 (right): PCA score plots of components 1-4
  
  rad.vn <- seq(0, 2 * pi, length.out = 100)
  
  pca.mn <- t(Biobase::exprs(eset))
  
  if (any(is.na(pca.mn))) {
    minN <- min(pca.mn, na.rm = TRUE)
    pca.mn[is.na(pca.mn)] <- minN
  }
  
  pca_model <- ropls::opls(pca.mn, predI = 4, algoC = "svd",
                        fig.pdfC = "none", info.txtC = "none")
  score.mn <- ropls::getScoreMN(pca_model)
  vRelVn <- pca_model@modelDF[, "R2X"]
  
  n <- nrow(score.mn)
  hotN <- 2 * (n - 1) * (n^2 - 1) / (n^2 * (n - 2))
  
  hotFisN <- hotN * stats::qf(0.95, 2, n - 2)
  
  pcsLs <- list(c(1, 2), c(3, 4))
  
  graphics::par(mar = c(3.6, 3.6, 0.6, 1.1))
  
  for (pcsN in 1:length(pcsLs)) {
    
    pcsVn <- pcsLs[[pcsN]]
    
    tcsMN <- score.mn[, pcsVn]
    
    micMN <- solve(stats::cov(tcsMN))
    
    n <- nrow(score.mn)
    hotN <- 2 * (n - 1) * (n^2 - 1) / (n^2 * (n - 2))
    
    hotFisN <- hotN * stats::qf(0.95, 2, n - 2)
    
    hotVn <- apply(tcsMN,
                   1,
                   function(x) 1 - stats::pf(1 / hotN * t(as.matrix(x)) %*% micMN %*% as.matrix(x), 2, n - 2))
    
    obsHotVi <- which(hotVn < 0.05)
    
    xLabC <- paste("t",
                   pcsVn[1],
                   "(",
                   round(vRelVn[pcsVn[1]] * 100),
                   "%)",
                   sep = "")
    
    yLabC <- paste("t",
                   pcsVn[2],
                   "(",
                   round(vRelVn[pcsVn[2]] * 100),
                   "%)",
                   sep = "")
    
    xLimVn <- c(-1, 1) * max(sqrt(stats::var(tcsMN[, 1]) * hotFisN), max(abs(tcsMN[, 1])))
    yLimVn <- c(-1, 1) * max(sqrt(stats::var(tcsMN[, 2]) * hotFisN), max(abs(tcsMN[, 2])))
    
    graphics::plot(tcsMN,
                   main = "",
                   type = "n",
                   xlab = "",
                   ylab = "",
                   xlim = xLimVn,
                   ylim = yLimVn)
    
    graphics::mtext(xLabC,
                    line = 2.2,
                    side = 1)
    graphics::mtext(yLabC,
                    line = 2.2,
                    side = 2)
    
    graphics::par(lwd = 1)
    
    graphics::abline(v = graphics::axTicks(1),
                     col = "grey")
    
    graphics::abline(h = graphics::axTicks(2),
                     col = "grey")
    
    graphics::abline(v = 0)
    graphics::abline(h = 0)
    
    graphics::lines(sqrt(stats::var(tcsMN[, 1]) * hotFisN) * cos(rad.vn),
                    sqrt(stats::var(tcsMN[, 2]) * hotFisN) * sin(rad.vn))
    
    graphics::points(tcsMN,
                     col = obsColVc,
                     pch = 18)
    
    if (length(obsHotVi))
      graphics::text(tcsMN[obsHotVi, 1],
                     tcsMN[obsHotVi, 2],
                     col = obsColVc[obsHotVi],
                     labels = obsNamVc[obsHotVi],
                     pos = 3)
    
  } ## for(pcsN in 1:length(pcsLs)) {
  
  invisible(list(median.vn = median.vn,
                 tcsMN = tcsMN))
  
} ## plotBatchF

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
