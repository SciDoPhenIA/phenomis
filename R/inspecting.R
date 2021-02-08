#### inspecting (MultiDataSet) ####

#' @rdname inspecting
#' @export
setMethod("inspecting", signature(x = "MultiDataSet"),
          function(x,
                   pool_as_pool1.l = FALSE,
                   pool_cv.n = 0.3,
                   span.n = 1,
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
            
            figPdfC <- figure.c
            if (figPdfC != "none")
              figPdfC <- "interactive"
            
            if (plot_dims.l)
              .barplot_dims(x, ifelse(!is.na(title.c), title.c, ""))
            
            for (set.c in names(x)) {
              
              if (report.c != "none")
                message("Inspecting the '", set.c, "' dataset...")
              
              ese <- x[[set.c]]
              
              ese <- inspecting(x = ese,
                                span.n = span.n,
                                pool_as_pool1.l = pool_as_pool1.l,
                                pool_cv.n = pool_cv.n,
                                sample_intensity.c = sample_intensity.c,
                                title.c = set.c,
                                figure.c = figPdfC,
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
            
            return(invisible(x))
            
          })

.barplot_dims <- function(mset,
                          title.c = "Dataset dimensions",
                          cex_axis.i = 15,
                          cex_bar.i = 6,
                          bar_just.n = 0.8,
                          cex_title.i = 25,
                          set_names.vc = NA) {
  
  set.i <- length(mset)

  dims.mn <- sapply(names(mset),
                    function(set.c)
                      Biobase::dims(mset[[set.c]]))
  
  if (length(set_names.vc) > 1 || !is.na(set_names.vc)) {
    stopifnot(length(set_names.vc) == set.i)
    colnames(dims.mn) <- set_names.vc
  }
    
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
                   span.n = 1,
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
            
            ## Checking
            
            check.l <- TRUE
            
            feat_allna.vl <- apply(Biobase::exprs(x), 1,
                                    function(feat.vn) all(is.na(feat.vn)))
            if (sum(feat_allna.vl, na.rm = TRUE)) {
              warning("The following feature(s) have NA only",
                      ifelse(!is.na(title.c) && title.c != "", paste0(" in the '", title.c, "' dataset"), ""),
                      ":\n",
                      paste(Biobase::featureNames(x)[feat_allna.vl], collapse = ", "))
              check.l <- FALSE
            }
            
            samp_allna.vl <- apply(Biobase::exprs(x), 2,
                                    function(samp.vn) all(is.na(samp.vn)))
            if (sum(samp_allna.vl, na.rm = TRUE)) {
              warning("The following sample(s) have NA only",
                      ifelse(!is.na(title.c) && title.c != "", paste0(" in the '", title.c, "' dataset"), ""),
                      ":\n",
                      paste(Biobase::sampleNames(x)[samp_allna.vl], collapse = ", "))
              check.l <- FALSE
            }
            
            feat_zerovar.vl <- apply(Biobase::exprs(x), 1,
                                    function(feat.vn) stats::var(feat.vn, na.rm = TRUE) < .Machine$double.eps)
            if (sum(feat_zerovar.vl, na.rm = TRUE)) {
              warning("The following feature(s) have zero variance",
                      ifelse(!is.na(title.c) && title.c != "", paste0(" in the '", title.c, "' dataset"), ""),
                      ":\n",
                      paste(Biobase::featureNames(x)[feat_zerovar.vl], collapse = ", "))
              check.l <- FALSE
            }
            
            samp_zerovar.vl <- apply(Biobase::exprs(x), 2,
                                     function(samp.vn) stats::var(samp.vn, na.rm = TRUE) < .Machine$double.eps)
            if (sum(samp_zerovar.vl, na.rm = TRUE)) {
              warning("The following sample(s) have zero variance",
                      ifelse(!is.na(title.c) && title.c != "", paste0(" in the '", title.c, "' dataset"), ""),
                      ":\n",
                      paste(Biobase::sampleNames(x)[samp_zerovar.vl], collapse = ", "))
              check.l <- FALSE
            }
            
            if (!check.l)
              stop("Please remove the sample(s) and/or feature(s) with NA only or 0 variance by using the 'filtering' method, and apply the 'inspecting' function again on the filtered dataset.")
            
            ## Description
            
            if (report.c != "none") {
              message("Data description:")
              message("observations: ", ncol(Biobase::exprs(x)))
              message("variables: ", nrow(Biobase::exprs(x)))
              message("missing: ", sum(is.na(Biobase::exprs(x))))
              message("0 values (%): ",
                      signif(sum(abs(Biobase::exprs(x)) < .Machine[["double.eps"]],
                                 na.rm = TRUE) / cumprod(dim(Biobase::exprs(x)))[2] * 100, 1), 2)
              message("min: ", signif(min(Biobase::exprs(x), na.rm = TRUE), 2))
              message("mean: ", signif(mean(Biobase::exprs(x), na.rm = TRUE), 2))
              message("median: ", signif(stats::median(Biobase::exprs(x), na.rm = TRUE), 2))
              message("max: ", signif(max(Biobase::exprs(x), na.rm = TRUE), 2))
              
              if (col_sampleType.c %in% colnames(Biobase::pData(x))) {
                message("Sample types:")
                print(table(Biobase::pData(x)[, col_sampleType.c]))
              }
            }
            
            ## Sample metrics
            
            sample_metrics.ls <- .sampleMetrics(eset = x)
            x <- sample_metrics.ls[["eset"]]
            pca_metrics.ls <- sample_metrics.ls[["pca_metrics.ls"]]
            
            ## Variable metrics
            
            x <- .variableMetrics(eset = x,
                                  pool_as_pool1.l = pool_as_pool1.l,
                                  col_sampleType.c = col_sampleType.c)
            
            ## Figure
            
            if (figure.c != "none") {
              
              if (figure.c != "interactive")
                grDevices::pdf(figure.c)
              
              .plotMetrics(eset = x,
                           pca_metrics.ls = pca_metrics.ls,
                           pool_cv.n = pool_cv.n,
                           span.n = span.n,
                           sample_intensity.c = sample_intensity.c,
                           col_batch.c = col_batch.c,
                           col_injectionOrder.c = col_injectionOrder.c,
                           col_sampleType.c = col_sampleType.c,
                           title.c = title.c)
              
              if (figure.c != "interactive")
                grDevices::dev.off()
              
            }
            
            ## End
            
            if (!(report.c %in% c("none", "interactive")))
              sink()
            
            return(invisible(x))
            
          })


.pca_metrics <- function(eset, pred.i = 2) {
  
  ## Hotelling: p-value associated to the distance from the center in the first PCA score plane
  
  pcaMod <- try(ropls::opls(t(Biobase::exprs(eset)), predI = pred.i,
                            crossvalI = min(ncol(Biobase::exprs(eset)), 7),
                            fig.pdfC = 'none', info.txtC = 'none'))
  
  if (inherits(pcaMod, "try-error")) {
    stop("The PCA could not be computed",
         ifelse(Biobase::experimentData(eset)@title != "",
                paste0(" for set '", Biobase::experimentData(eset)@title, "'"),
                ""),
         ". Check for the presence of features with a high proportion of NA or a low variance and discard them with the 'filtering' method before starting 'inspecting' again.")
  }
  
  varRelVn <- ropls::getPcaVarVn(pcaMod) / nrow(Biobase::exprs(eset)) ## for plotting
  
  pcaScoreMN <- ropls::getScoreMN(pcaMod)
  
  n <- ncol(Biobase::exprs(eset))
  hotN <- 2 * (n - 1) * (n^2 - 1) / (n^2 * (n - 2))
  
  invCovSco12MN <- solve(stats::cov(pcaScoreMN[, 1:2]))
  
  hotPva12Vn <- apply(pcaScoreMN[, 1:2],
                    1,
                    function(x)
                      1 - stats::pf(1 / hotN * t(as.matrix(x)) %*% invCovSco12MN %*% as.matrix(x), 2, n - 2))
  
  list(hotN = hotN,
       varRelVn = varRelVn,
       pcaScoreMN = pcaScoreMN,
       hotPva12Vn = hotPva12Vn)
  
}


.sampleMetrics <- function(eset) {
  
  pca_metrics.ls <- .pca_metrics(eset = eset, pred.i = 2)
  
  pdaDF <- Biobase::pData(eset)

  pdaDF[, "pca_sco1"] <- pca_metrics.ls[["pcaScoreMN"]][, 1]
  pdaDF[, "pca_sco2"] <- pca_metrics.ls[["pcaScoreMN"]][, 2]
  
  pdaDF[, "hotel_pval"] <- pca_metrics.ls[["hotPva12Vn"]]
  
  ## p-value associated to number of missing values
  
  missZscoVn <- .zscore(apply(t(Biobase::exprs(eset)),
                              1,
                              function(rowVn) {
                                sum(is.na(rowVn))
                              }))
  
  pdaDF[, "miss_pval"] <- sapply(missZscoVn, function(zscoN) 2 * (1 - stats::pnorm(abs(zscoN))))
  
  ## p-value associated to the deciles of the profiles
  
  deciMN <- t(as.matrix(apply(t(Biobase::exprs(eset)),
                              1,
                              function(x) stats::quantile(x, 0.1 * 1:9, na.rm = TRUE))))
  
  deciZscoMN <- apply(deciMN, 2, .zscore)
  
  deciZscoMaxVn <- apply(deciZscoMN, 1, function(rowVn) rowVn[which.max(abs(rowVn))])
  
  pdaDF[, "deci_pval"] <- sapply(deciZscoMaxVn, function(zscoN) 2 * (1 - stats::pnorm(abs(zscoN))))
  
  Biobase::pData(eset) <- pdaDF
  
  return(list(eset = eset,
              pca_metrics.ls = pca_metrics.ls))
  
}


.variableMetrics <- function(eset,
                             pool_as_pool1.l,
                             col_sampleType.c) { ## for the call to 'univariate'
  
  fdaDF <- Biobase::fData(eset)
  
  tmpDF <- fdaDF ## some of the intermediate metrics will not be included in fData
  
  ## 'blank' observations
  
  if (col_sampleType.c %in% Biobase::varLabels(eset) && "blank" %in% Biobase::pData(eset)[, col_sampleType.c]) {
    
    blkVl <- Biobase::pData(eset)[, col_sampleType.c] == "blank"
    
    if (sum(blkVl) == 1) {
      tmpDF[, "blank_mean"] <- t(Biobase::exprs(eset))[blkVl, ]
    } else {
      tmpDF[, "blank_mean"] <- apply(t(Biobase::exprs(eset))[blkVl, , drop = FALSE],
                                     2,
                                     function(varVn) mean(varVn, na.rm = TRUE))
    }
    
    if (sum(blkVl) == 1) {
      tmpDF[, "blank_sd"] <- rep(0, nrow(fdaDF))
    } else {
      tmpDF[, "blank_sd"] <- apply(t(Biobase::exprs(eset))[blkVl, , drop = FALSE],
                                   2,
                                   function(varVn) stats::sd(varVn, na.rm = TRUE))
    }
    
    tmpDF[, "blank_CV"] <- tmpDF[, "blank_sd"] / tmpDF[, "blank_mean"]
    
  }
  
  ## 'sample' observations
  
  if (col_sampleType.c %in% Biobase::varLabels(eset) &&
      "sample" %in% Biobase::pData(eset)[, col_sampleType.c]) {
    
    samVl <- Biobase::pData(eset)[, col_sampleType.c] == "sample"
    
    if (sum(samVl) == 1) {
      tmpDF[, "sample_mean"] <- t(Biobase::exprs(eset))[samVl, ]
    } else {
      tmpDF[, "sample_mean"] <- apply(t(Biobase::exprs(eset))[samVl, , drop = FALSE], 2,
                                      function(varVn) mean(varVn, na.rm = TRUE))
    }
    
    if (sum(samVl) == 1) {
      tmpDF[, "sample_sd"] <- rep(0, nrow(fdaDF))
    } else {
      tmpDF[, "sample_sd"] <- apply(t(Biobase::exprs(eset))[samVl, , drop = FALSE], 2,
                                    function(varVn) stats::sd(varVn, na.rm = TRUE))
    }
    
    tmpDF[, "sample_CV"] <- tmpDF[, "sample_sd"] / tmpDF[, "sample_mean"]
    
  }
  
  ## 'blank' mean / 'sample' mean ratio
  
  if (all(c("blank_mean", "sample_mean") %in% colnames(tmpDF)))
    fdaDF[, "blankMean_over_sampleMean"] <- tmpDF[, "blank_mean"] / tmpDF[, "sample_mean"]
  
  ## 'pool' observations
  
  if (col_sampleType.c %in% Biobase::varLabels(eset) &&
      "pool" %in% Biobase::pData(eset)[, col_sampleType.c]) {
    
    pooVl <- Biobase::pData(eset)[, col_sampleType.c] == "pool"
    
    if (sum(pooVl) == 1) {
      tmpDF[, "pool_mean"] <- t(Biobase::exprs(eset))[pooVl, ]
    } else {
      tmpDF[, "pool_mean"] <- apply(t(Biobase::exprs(eset))[pooVl, , drop = FALSE], 2,
                                    function(varVn) mean(varVn, na.rm = TRUE))
    }
    
    if (sum(pooVl) == 1) {
      tmpDF[, "pool_sd"] <- rep(0, nrow(fdaDF))
    } else {
      tmpDF[, "pool_sd"] <- apply(t(Biobase::exprs(eset))[pooVl, , drop = FALSE], 2,
                                  function(varVn) stats::sd(varVn, na.rm = TRUE))
    }
    
    fdaDF[, "pool_CV"] <- tmpDF[, "pool_sd"] / tmpDF[, "pool_mean"]
    
  }
  
  ## 'pool' CV / 'sample' CV ratio
  
  if ("pool_CV" %in% colnames(fdaDF) && "sample_CV" %in% colnames(tmpDF))
    fdaDF[, "poolCV_over_sampleCV"] <- fdaDF[, "pool_CV"] / tmpDF[, "sample_CV"]
  
  ## 'pool' dilutions
  
  if (col_sampleType.c %in% Biobase::varLabels(eset) &&
      any(grepl("pool.+", Biobase::pData(eset)[, col_sampleType.c]))) {
    
    pooVi <- grep("pool.*", Biobase::pData(eset)[, col_sampleType.c]) ## pool, pool2, pool4, poolInter, ...
    
    pooNamVc <- Biobase::pData(eset)[pooVi, col_sampleType.c]
    
    if (pool_as_pool1.l) {
      
      pooNamVc[pooNamVc == "pool"] <- "pool1" ## 'pool' -> 'pool1'
      
    } else {
      
      pooVl <- pooNamVc == "pool"
      pooVi <- pooVi[!pooVl]
      pooNamVc <- pooNamVc[!pooVl]
      
    }
    
    pooDilVc <- gsub("pool", "", pooNamVc)
    
    pooDilVl <- sapply(pooDilVc, .allDigits)
    
    if (sum(pooDilVl)) {
      
      pooNamVc <- pooNamVc[pooDilVl]
      
      pooVi <- pooVi[pooDilVl]
      
      dilVn <- 1 / as.numeric(pooDilVc[pooDilVl])
      
      var.vn <- apply(t(Biobase::exprs(eset))[pooVi, , drop = FALSE], 2,
                      function(varVn) stats::var(varVn))
      
      var.vl <- !is.na(var.vn) & var.vn > 0
      
      poolDil_pval.vn <- poolDil_cor.vn <- rep(NA, length(var.vn))
      
      poolDil_cor.vn[var.vl] <- apply(t(Biobase::exprs(eset))[pooVi, var.vl, drop = FALSE], 2,
                                       function(varVn) stats::cor(dilVn, varVn))
      
      fdaDF[, "poolDil_cor"] <- poolDil_cor.vn
      
      # fdaDF[, "poolDil_cor"] <- apply(t(Biobase::exprs(eset))[pooVi, , drop = FALSE], 2,
      #                                 function(varVn) stats::cor(dilVn, varVn))
      
      poolDil_pval.vn[var.vl] <- apply(t(Biobase::exprs(eset))[pooVi, var.vl, drop = FALSE], 2,
                                       function(varVn) stats::cor.test(dilVn, varVn)[["p.value"]])
      
      fdaDF[, "poolDil_pval"] <- poolDil_pval.vn
      
      # fdaDF[, "poolDil_pval"] <- apply(t(Biobase::exprs(eset))[pooVi, , drop = FALSE], 2,
      #                                  function(varVn) {
      #                                    pval.n <- try(stats::cor.test(dilVn, varVn)[["p.value"]], silent = TRUE)
      #                                    if (!inherits(pval.n, "try-error")) {
      #                                      return(pval.n)
      #                                    } else {
      #                                      return(NA)
      #                                    }
      #                                  })
      
    }
    
  }
  
  Biobase::fData(eset) <- fdaDF
  
  return(eset)
  
}


.na_zerovar_filter <- function(eset,
                               max_na_prop.n,
                               report.c) {
  
  feat_init.vc <- Biobase::featureNames(eset)
  
  na_zerovar <- function(eset)
    apply(Biobase::exprs(eset), 1,
          function(var.vn)
            sum(is.na(var.vn))/length(var.vn) > max_na_prop.n ||
            stats::var(var.vn, na.rm = TRUE) < .Machine$double.eps)
  
  na_zero.vl <- na_zerovar(eset)
  
  while (sum(na_zero.vl)) {
    
    if (sum(!na_zero.vl) == 0)
      stop("No more features left in the dataset.", call. = FALSE)
    
    eset <- eset[!na_zero.vl, ]
    
    na_zero.vl <- na_zerovar(eset)
    
  }
  
  feat_final.vc <- Biobase::featureNames(eset)
  
  feat_diff.vc <- setdiff(feat_init.vc, feat_final.vc)
  
  if (length(feat_diff.vc) && report.c != "none")
    message(length(feat_diff.vc),
            " feature(s) was/were discarded because of a proportion of NAs > ",
            max_na_prop.n, " or a variance of 0.")
  
  eset
  
}


.plotMetrics <- function(eset,
                         pca_metrics.ls,
                         pool_cv.n,
                         span.n,
                         sample_intensity.c,
                         col_batch.c,
                         col_injectionOrder.c,
                         col_sampleType.c,
                         title.c) {
  
  ## Constants
  
  marLs <- list(tit = c(0.6, 1.1, 1.1, 0.6),
                drivol = c(3.5, 3.6, 4.1, 0.6),
                sca = c(0.6, 3.1, 4.1, 0.6),
                scavol = c(3.5, 3.1, 4.1, 0.6),
                ima = c(0.6, 2.6, 4.1, 0.9),
                imavol = c(3.5, 2.6, 4.1, 0.9),
                dri = c(3.5, 3.6, 1.1, 0.6),
                vol = c(3.5, 3.6, 1.1, 0.9),
                pca = c(3.5, 3.6, 1.1, 0.9))
  
  heat_palette.vc <- rev(grDevices::rainbow(ceiling(256 * 1.5))[1:256])
 
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
  
  sample_color.vc <- .sample_color_eset(eset = eset,
                                        col_sampleType.c = col_sampleType.c)
 
  ## tit: Title
  
  graphics::par(mar = marLs[["tit"]])
  graphics::plot(0:1, bty = "n", type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
  graphics::text(1, 0.95, adj = 0, cex = 1.2, labels = title.c)
  graphics::text(1, 0.75, adj = 0, labels = paste0("NAs: ",
                                                   round(length(which(is.na(c(t(Biobase::exprs(eset)))))) / cumprod(dim(t(Biobase::exprs(eset))))[2] * 100), "%"))
  graphics::text(1, 0.68, adj = 0, labels = paste0("0 values: ",
                                                   round(sum(abs(t(Biobase::exprs(eset))) < .Machine[["double.eps"]], na.rm = TRUE) / cumprod(dim(t(Biobase::exprs(eset))))[2] * 100, 2), "%"))
  graphics::text(1, 0.61, adj = 0, labels = paste0("min: ", signif(min(t(Biobase::exprs(eset)), na.rm = TRUE), 2)))
  graphics::text(1, 0.54, adj = 0, labels = paste0("median: ", signif(stats::median(t(Biobase::exprs(eset)), na.rm = TRUE), 2)))
  graphics::text(1, 0.47, adj = 0, labels = paste0("mean: ", signif(mean(t(Biobase::exprs(eset)), na.rm = TRUE), 2)))
  graphics::text(1, 0.40, adj = 0, labels = paste0("max: ", signif(max(t(Biobase::exprs(eset)), na.rm = TRUE), 2)))
  if (col_sampleType.c %in% Biobase::varLabels(eset) &&
      "pool" %in% Biobase::pData(eset)[, col_sampleType.c])
    graphics::text(1,
                   0.33,
                   adj = 0,
                   labels = paste0("CVpool<",
                                   round(pool_cv.n * 100), "%: ",
                                   round(sum(Biobase::fData(eset)[, "pool_CV"] < pool_cv.n, na.rm = TRUE) / nrow(Biobase::exprs(eset)) * 100),
                                   "%"))
  
  ## sca: Color scale
  
  .plot_color_scale(eset = eset,
                    palette.vc = heat_palette.vc,
                    mar.vn = marLs[["sca"]])
  
  ## ima: Image
  
  .plot_image(eset = eset,
              palette.vc = heat_palette.vc,
              mar.vn = marLs[["ima"]])
  
  ## dri: Analytical drift
  
  .plot_drift(eset = eset,
              span.n = span.n,
              sample_intensity.c = sample_intensity.c,
              mar.vn = marLs[["dri"]],
              col_batch.c,
              col_injectionOrder.c,
              col_sampleType.c)
  
  ## pca: PCA and Hotelling ellipse
  
  .plot_pca_metrics(eset = eset,
                    col_sampleType.c = col_sampleType.c,
                    labels.l = TRUE,
                    mar.vn = marLs[["pca"]],
                    pca_metrics.ls = pca_metrics.ls)
  
  ## resetting par
  
  graphics::par(mfrow = c(1, 1))
  graphics::par(opar)
  
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

.plot_color_scale <- function(eset,
                              palette.vc,
                              mar.vn) {
  
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
  
  eval(parse(text = paste0("axis(at = .plot_pretty_axis(c(ifelse(min(t(Biobase::exprs(eset)), na.rm = TRUE) == -Inf, yes = 0, no = min(t(Biobase::exprs(eset)), na.rm = TRUE)) , max(t(Biobase::exprs(eset)), na.rm = TRUE)), 256)$atVn,
                           font = 2,
                           font.axis = 2,
                           labels = .plot_pretty_axis(c(ifelse(min(t(Biobase::exprs(eset)), na.rm = TRUE) == -Inf, yes = 0, no = min(t(Biobase::exprs(eset)), na.rm = TRUE)), max(t(Biobase::exprs(eset)), na.rm = TRUE)), 256)$labVn,
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


.plot_image <- function(eset,
                        palette.vc,
                        mar.vn) {
  
  graphics::par(mar = mar.vn)
  
  image.mn <- Biobase::exprs(eset)[, rev(1:ncol(Biobase::exprs(eset))), drop = FALSE]
  
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
  
  if (length(rownames(t(Biobase::exprs(eset)))) == 0) {
    rowNamVc <- rep("", times = ncol(Biobase::exprs(eset)))
  } else
    rowNamVc <- rownames(t(Biobase::exprs(eset)))
  
  if (length(colnames(t(Biobase::exprs(eset)))) == 0) {
    colNamVc <- rep("", times = nrow(Biobase::exprs(eset)))
  } else
    colNamVc <- colnames(t(Biobase::exprs(eset)))
  
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


.plot_drift <- function(eset,
                        span.n,
                        sample_intensity.c,
                        mar.vn,
                        col_batch.c,
                        col_injectionOrder.c,
                        col_sampleType.c) {
  
  sample_color.vc <- .sample_color_eset(eset = eset,
                                        col_sampleType.c = col_sampleType.c)
  
  graphics::par(mar = mar.vn)
  
  ## ordering
  
  texprs.mn <- t(Biobase::exprs(eset))
  pdata.df <- Biobase::pData(eset)
  
  pdata.df[, "ordIniVi"] <- 1:nrow(texprs.mn)
  
  if (col_injectionOrder.c %in% colnames(pdata.df)) {
    ordNamC <- "Injection Order"
    if (col_batch.c %in% colnames(pdata.df)) {
      ordVi <- order(pdata.df[, col_batch.c],
                     pdata.df[, col_injectionOrder.c])
    } else
      ordVi <- order(pdata.df[, col_injectionOrder.c])
  } else {
    ordNamC <- "Samples"
    ordVi <- 1:nrow(texprs.mn)
  }
  
  texprs.mn <- texprs.mn[ordVi, ]
  pdata.df <- pdata.df[ordVi, ]
  sample_color_ordered.vc <- sample_color.vc[ordVi]
  
  if (col_batch.c %in% colnames(pdata.df))
    batch.table <- table(pdata.df[, col_batch.c])
  
  sample_means.vn <- eval(parse(text = paste0("apply(texprs.mn, 1, function(obsVn) ",
                                              sample_intensity.c, "(obsVn, na.rm = TRUE))")))
  
  graphics::plot(sample_means.vn,
                 col = sample_color_ordered.vc,
                 pch = 18,
                 # ylim = range(texprs.mn, na.rm = TRUE),
                 type = "n",
                 xaxs = "i",
                 xlab = "",
                 ylab = "")
  
  # for (samI in 1:nrow(texprs.mn))
  #   graphics::boxplot(texprs.mn[samI, ],
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
  
  if (col_batch.c %in% colnames(pdata.df) && length(unique(pdata.df[, "batch"])) > 1) {
    
    graphics::abline(v = cumsum(batch.table) + 0.5,
                     col = "red")
    
    graphics::mtext(names(batch.table),
                    at = batch.table / 2 + c(0, cumsum(batch.table[-length(batch.table)])),
                    cex = 0.7)
    
    for (batC in names(batch.table)) {
      
      batch_seq.vi <- which(pdata.df[, col_batch.c] == batC)
      
      if (col_sampleType.c %in% colnames(pdata.df)) {
        batch_sample.vi <- intersect(batch_seq.vi,
                              grep("sample", pdata.df[, col_sampleType.c]))
      } else
        batch_sample.vi <- batch_seq.vi
      
      graphics::lines(batch_seq.vi,
                      .loess(sample_means.vn, batch_sample.vi, batch_seq.vi, span.n),
                      col = .sample_color_vector("sample"))
      
      if (col_sampleType.c %in% colnames(pdata.df) &&
          "pool" %in% pdata.df[, col_sampleType.c]) {
        
        batch_pool.vi <- intersect(batch_seq.vi,
                              grep("^pool$", pdata.df[, col_sampleType.c]))
        
        graphics::lines(batch_seq.vi,
                        .loess(sample_means.vn, batch_pool.vi, batch_seq.vi, span.n),
                        col = .sample_color_vector("pool"))
        
      }
      
    }
    
  } else {
    
    batch_seq.vi <- 1:nrow(pdata.df)
    
    if (col_sampleType.c %in% colnames(pdata.df)) {
      batch_sample.vi <- intersect(batch_seq.vi,
                            grep("sample", pdata.df[, col_sampleType.c]))
    } else
      batch_sample.vi <- batch_seq.vi
    
    graphics::lines(batch_seq.vi,
                    .loess(sample_means.vn, batch_sample.vi, batch_seq.vi, span.n),
                    col = .sample_color_vector("sample"))
    
    if (col_sampleType.c %in% colnames(pdata.df) &&
        "pool" %in% pdata.df[, col_sampleType.c]) {
      
      batch_pool.vi <- intersect(batch_seq.vi,
                            grep("^pool$", pdata.df[, col_sampleType.c]))
      
      graphics::lines(batch_seq.vi,
                      .loess(sample_means.vn, batch_pool.vi, batch_seq.vi, span.n),
                      col = .sample_color_vector("pool"))
      
      
    }
    
  }
  
}


.plot_pca_metrics <- function(eset,
                              pred.i = 2,
                              show_pred.vi = c(1, 2),
                              col_sampleType.c = "sampleType",
                              labels.l = TRUE,
                              mar.vn,
                              pca_metrics.ls = NULL) {
  
  if (is.null(pca_metrics.ls))
    pca_metrics.ls <- .pca_metrics(eset = eset, pred.i = pred.i)
  
  pcaScoreMN <- pca_metrics.ls[["pcaScoreMN"]]

  if (ncol(pcaScoreMN) < max(show_pred.vi))
    stop("Not enough pca components computed")
  
  sample_color.vc <- .sample_color_eset(eset = eset, col_sampleType.c = col_sampleType.c)
  
  graphics::par(mar = mar.vn)
  
  graphics::plot(pcaScoreMN[, show_pred.vi],
                 type = "n",
                 xlab = "",
                 ylab = "")
  graphics::mtext(paste("t1 (", round(pca_metrics.ls[["varRelVn"]][show_pred.vi[1]] * 100), "%)", sep = ""),
                  cex = 0.7,
                  line = 2,
                  side = 1)
  graphics::mtext(paste("t2 (", round(pca_metrics.ls[["varRelVn"]][show_pred.vi[2]] * 100), "%)", sep = ""),
                  cex = 0.7,
                  las = 0,
                  line = 2,
                  side = 2)
  graphics::abline(h = 0, lty = "dashed")
  graphics::abline(v = 0, lty = "dashed")
  radVn <- seq(0, 2 * pi, length.out = 100)
  
  hotFisN <- pca_metrics.ls[["hotN"]] * stats::qf(0.95, 2, Biobase::dims(eset)["Samples", 1] - 2)
  
  graphics::lines(sqrt(stats::var(pcaScoreMN[, show_pred.vi[1]]) * hotFisN) * cos(radVn),
                  sqrt(stats::var(pcaScoreMN[, show_pred.vi[2]]) * hotFisN) * sin(radVn))
  
  if (identical(show_pred.vi, c(1, 2))) {
    
    hotVn <- pca_metrics.ls[["hotPva12Vn"]]
    
  } else if (identical(show_pred.vi, c(3, 4))) {
    
    invCovSco34MN <- solve(stats::cov(pcaScoreMN[, 3:4]))
    
    hotVn <- apply(pcaScoreMN[, 3:4],
                   1,
                   function(x)
                     1 - stats::pf(1 / pca_metrics.ls[["hotN"]] * t(as.matrix(x)) %*% invCovSco34MN %*% as.matrix(x), 2, Biobase::dims(eset)["Samples", 1] - 2))
    
  } else
    stop("'show_pred.vi' must be either 'c(1, 2)' or 'c(3, 4)'")

  obsHotVi <- which(hotVn < 0.05)

  if (labels.l) {
    graphics::text(pcaScoreMN[obsHotVi, show_pred.vi[1]],
                   pcaScoreMN[obsHotVi, show_pred.vi[2]],
                   cex = 0.7,
                   col = sample_color.vc[obsHotVi],
                   labels = Biobase::sampleNames(eset)[obsHotVi])
    graphics::points(pcaScoreMN[setdiff(1:nrow(pcaScoreMN), obsHotVi), show_pred.vi[1]],
                     pcaScoreMN[setdiff(1:nrow(pcaScoreMN), obsHotVi), show_pred.vi[2]],
                     col = sample_color.vc[setdiff(1:nrow(pcaScoreMN), obsHotVi)],
                     pch = 16)
  } else
    graphics::points(pcaScoreMN[, show_pred.vi[1]],
                     pcaScoreMN[, show_pred.vi[2]],
                     col = sample_color.vc,
                     pch = 16)
  
  if (col_sampleType.c %in% Biobase::varLabels(eset)) {
    obsColVuc <- sample_color.vc[sort(unique(names(sample_color.vc)))]
    legOrdVc <- c("blank", paste0("pool", 8:1), "pool", "other", "sample")
    obsColVuc <- obsColVuc[legOrdVc[legOrdVc %in% names(obsColVuc)]]
    
    graphics::text(rep(graphics::par("usr")[1], times = length(obsColVuc)),
                   graphics::par("usr")[3] + (0.97 - length(obsColVuc) * 0.03 + 1:length(obsColVuc) * 0.03) * diff(graphics::par("usr")[3:4]),
                   col = obsColVuc,
                   font = 2,
                   labels = names(obsColVuc),
                   pos = 4)
  }
  
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
