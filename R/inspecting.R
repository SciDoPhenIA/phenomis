#### inspecting (MultiDataSet) ####

#' @rdname inspecting
#' @export
setMethod("inspecting", signature(x = "MultiDataSet"),
          function(x,
                   pool_as_pool1.l = FALSE,
                   pool_cv.n = 0.3,
                   span.n = 1,
                   sample_intensity.c = c("median", "mean")[1],
                   title.c = NA,
                   plot_dims.l = TRUE,
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
                                pool_as_pool1.l = pool_as_pool1.l,
                                pool_cv.n = pool_cv.n,
                                span.n = span.n,
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
            
            return(x)
            
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
                   sample_intensity.c = c("median", "mean")[1],
                   title.c = NA,
                   plot_dims.l = TRUE,
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
              
              if ("sampleType" %in% colnames(Biobase::pData(x))) {
                message("Sample types:")
                print(table(Biobase::pData(x)[, "sampleType"]))
              }
            }
            
            ## Sample metrics
            
            sampMetLs <- .sampleMetrics(x = x)
            
            x <- sampMetLs[["x"]]
            pcaLs <- sampMetLs[["pcaLs"]]
            
            ## Variable metrics
            
            x <- .variableMetrics(x = x, pool_as_pool1.l = pool_as_pool1.l)
            
            ## Figure
            
            if (figure.c != "none") {
              
              if (figure.c != "interactive")
                grDevices::pdf(figure.c)
              
              .plotMetrics(x = x,
                           pcaLs = pcaLs,
                           pool_cv.n = pool_cv.n,
                           title.c = title.c,
                           span.n = span.n,
                           sample_intensity.c = sample_intensity.c)
              
              if (figure.c != "interactive")
                grDevices::dev.off()
              
            }
            
            ## End
            
            if (!(report.c %in% c("none", "interactive")))
              sink()
            
            return(invisible(x))
            
          })


.sampleMetrics <- function(x) {
  
  ## Hotelling: p-value associated to the distance from the center in the first PCA score plane
  
  pcaMod <- try(ropls::opls(t(Biobase::exprs(x)), predI = 2,
                            crossvalI = min(ncol(Biobase::exprs(x)), 7),
                            fig.pdfC = 'none', info.txtC = 'none'))
  
  if (inherits(pcaMod, "try-error")) {
    stop("The PCA could not be computed",
             ifelse(Biobase::experimentData(x)@title != "",
                    paste0(" for set '", Biobase::experimentData(x)@title, "'"),
                    ""),
             ". Check for the presence of features with a high proportion of NA or a low variance and discard them with the 'filtering' method before starting 'inspecting' again.")
  }
  
  varRelVn <- ropls::getPcaVarVn(pcaMod) / nrow(Biobase::exprs(x)) ## for plotting
  
  pcaScoreMN <- ropls::getScoreMN(pcaMod)
  
  pdaDF <- Biobase::pData(x)
  
  pdaDF[, "pca_sco1"] <- pcaScoreMN[, 1]
  pdaDF[, "pca_sco2"] <- pcaScoreMN[, 2]
  
  invCovScoMN <- solve(stats::cov(pcaScoreMN))
  
  n <- ncol(Biobase::exprs(x))
  hotN <- 2 * (n - 1) * (n^2 - 1) / (n^2 * (n - 2))
  
  hotPvaVn <- apply(pcaScoreMN,
                    1,
                    function(x)
                      1 - stats::pf(1 / hotN * t(as.matrix(x)) %*% invCovScoMN %*% as.matrix(x), 2, n - 2))
  
  pdaDF[, "hotel_pval"] <- hotPvaVn
  
  pcaLs <- list(hotN = hotN,
                varRelVn = varRelVn)
  
  
  ## p-value associated to number of missing values
  
  missZscoVn <- .zscore(apply(t(Biobase::exprs(x)),
                              1,
                              function(rowVn) {
                                sum(is.na(rowVn))
                              }))
  
  pdaDF[, "miss_pval"] <- sapply(missZscoVn, function(zscoN) 2 * (1 - stats::pnorm(abs(zscoN))))
  
  ## p-value associated to the deciles of the profiles
  
  deciMN <- t(as.matrix(apply(t(Biobase::exprs(x)),
                              1,
                              function(x) stats::quantile(x, 0.1 * 1:9, na.rm = TRUE))))
  
  deciZscoMN <- apply(deciMN, 2, .zscore)
  
  deciZscoMaxVn <- apply(deciZscoMN, 1, function(rowVn) rowVn[which.max(abs(rowVn))])
  
  pdaDF[, "deci_pval"] <- sapply(deciZscoMaxVn, function(zscoN) 2 * (1 - stats::pnorm(abs(zscoN))))
  
  Biobase::pData(x) <- pdaDF
  
  resLs <- list(x = x,
                pcaLs = pcaLs)
  
  return(resLs)
  
}


.variableMetrics <- function(x,
                             pool_as_pool1.l) { ## for the call to 'univariate'
  
  fdaDF <- Biobase::fData(x)
  
  tmpDF <- fdaDF ## some of the intermediate metrics will not be included in fData
  
  ## 'blank' observations
  
  if ("sampleType" %in% colnames(Biobase::pData(x)) && "blank" %in% Biobase::pData(x)[, "sampleType"]) {
    
    blkVl <- Biobase::pData(x)[, "sampleType"] == "blank"
    
    if (sum(blkVl) == 1) {
      tmpDF[, "blank_mean"] <- t(Biobase::exprs(x))[blkVl, ]
    } else {
      tmpDF[, "blank_mean"] <- apply(t(Biobase::exprs(x))[blkVl, , drop = FALSE],
                                     2,
                                     function(varVn) mean(varVn, na.rm = TRUE))
    }
    
    if (sum(blkVl) == 1) {
      tmpDF[, "blank_sd"] <- rep(0, nrow(fdaDF))
    } else {
      tmpDF[, "blank_sd"] <- apply(t(Biobase::exprs(x))[blkVl, , drop = FALSE],
                                   2,
                                   function(varVn) stats::sd(varVn, na.rm = TRUE))
    }
    
    tmpDF[, "blank_CV"] <- tmpDF[, "blank_sd"] / tmpDF[, "blank_mean"]
    
  }
  
  ## 'sample' observations
  
  if ("sampleType" %in% colnames(Biobase::pData(x)) &&
      "sample" %in% Biobase::pData(x)[, "sampleType"]) {
    
    samVl <- Biobase::pData(x)[, "sampleType"] == "sample"
    
    if (sum(samVl) == 1) {
      tmpDF[, "sample_mean"] <- t(Biobase::exprs(x))[samVl, ]
    } else {
      tmpDF[, "sample_mean"] <- apply(t(Biobase::exprs(x))[samVl, , drop = FALSE], 2,
                                      function(varVn) mean(varVn, na.rm = TRUE))
    }
    
    if (sum(samVl) == 1) {
      tmpDF[, "sample_sd"] <- rep(0, nrow(fdaDF))
    } else {
      tmpDF[, "sample_sd"] <- apply(t(Biobase::exprs(x))[samVl, , drop = FALSE], 2,
                                    function(varVn) stats::sd(varVn, na.rm = TRUE))
    }
    
    tmpDF[, "sample_CV"] <- tmpDF[, "sample_sd"] / tmpDF[, "sample_mean"]
    
  }
  
  ## 'blank' mean / 'sample' mean ratio
  
  if (all(c("blank_mean", "sample_mean") %in% colnames(tmpDF)))
    fdaDF[, "blankMean_over_sampleMean"] <- tmpDF[, "blank_mean"] / tmpDF[, "sample_mean"]
  
  ## 'pool' observations
  
  if ("sampleType" %in% colnames(Biobase::pData(x)) &&
      "pool" %in% Biobase::pData(x)[, "sampleType"]) {
    
    pooVl <- Biobase::pData(x)[, "sampleType"] == "pool"
    
    if (sum(pooVl) == 1) {
      tmpDF[, "pool_mean"] <- t(Biobase::exprs(x))[pooVl, ]
    } else {
      tmpDF[, "pool_mean"] <- apply(t(Biobase::exprs(x))[pooVl, , drop = FALSE], 2,
                                    function(varVn) mean(varVn, na.rm = TRUE))
    }
    
    if (sum(pooVl) == 1) {
      tmpDF[, "pool_sd"] <- rep(0, nrow(fdaDF))
    } else {
      tmpDF[, "pool_sd"] <- apply(t(Biobase::exprs(x))[pooVl, , drop = FALSE], 2,
                                  function(varVn) stats::sd(varVn, na.rm = TRUE))
    }
    
    fdaDF[, "pool_CV"] <- tmpDF[, "pool_sd"] / tmpDF[, "pool_mean"]
    
  }
  
  ## 'pool' CV / 'sample' CV ratio
  
  if ("pool_CV" %in% colnames(fdaDF) && "sample_CV" %in% colnames(tmpDF))
    fdaDF[, "poolCV_over_sampleCV"] <- fdaDF[, "pool_CV"] / tmpDF[, "sample_CV"]
  
  ## 'pool' dilutions
  
  if ("sampleType" %in% colnames(Biobase::pData(x)) &&
      any(grepl("pool.+", Biobase::pData(x)[, "sampleType"]))) {
    
    pooVi <- grep("pool.*", Biobase::pData(x)[, "sampleType"]) ## pool, pool2, pool4, poolInter, ...
    
    pooNamVc <- Biobase::pData(x)[pooVi, "sampleType"]
    
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
      
      var.vn <- apply(t(Biobase::exprs(x))[pooVi, , drop = FALSE], 2,
                      function(varVn) stats::var(varVn))
      
      var.vl <- !is.na(var.vn) & var.vn > 0
      
      poolDil_pval.vn <- poolDil_cor.vn <- rep(NA, length(var.vn))
      
      poolDil_cor.vn[var.vl] <- apply(t(Biobase::exprs(x))[pooVi, var.vl, drop = FALSE], 2,
                                       function(varVn) stats::cor(dilVn, varVn))
      
      fdaDF[, "poolDil_cor"] <- poolDil_cor.vn
      
      # fdaDF[, "poolDil_cor"] <- apply(t(Biobase::exprs(x))[pooVi, , drop = FALSE], 2,
      #                                 function(varVn) stats::cor(dilVn, varVn))
      
      poolDil_pval.vn[var.vl] <- apply(t(Biobase::exprs(x))[pooVi, var.vl, drop = FALSE], 2,
                                       function(varVn) stats::cor.test(dilVn, varVn)[["p.value"]])
      
      fdaDF[, "poolDil_pval"] <- poolDil_pval.vn
      
      # fdaDF[, "poolDil_pval"] <- apply(t(Biobase::exprs(x))[pooVi, , drop = FALSE], 2,
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
  
  Biobase::fData(x) <- fdaDF
  
  return(x)
  
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


.plotMetrics <- function(x,
                         pcaLs,
                         pool_cv.n,
                         title.c,
                         span.n,
                         sample_intensity.c) {
  
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
  
  palHeaVc <- rev(grDevices::rainbow(ceiling(256 * 1.5))[1:256])
  
  ## Functions
  
  .prettyAxis <- function(valVn,
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
  
  if ("sampleType" %in% colnames(Biobase::pData(x))) {
    obsColVc <- .sample_color(Biobase::pData(x)[, "sampleType"])
  } else
    obsColVc <- rep("black", nrow(Biobase::pData(x)))
  
  ## tit: Title
  
  graphics::par(mar = marLs[["tit"]])
  graphics::plot(0:1, bty = "n", type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
  graphics::text(1, 0.95, adj = 0, cex = 1.2, labels = title.c)
  graphics::text(1, 0.75, adj = 0, labels = paste0("NAs: ",
                                                   round(length(which(is.na(c(t(Biobase::exprs(x)))))) / cumprod(dim(t(Biobase::exprs(x))))[2] * 100), "%"))
  graphics::text(1, 0.68, adj = 0, labels = paste0("0 values: ",
                                                   round(sum(abs(t(Biobase::exprs(x))) < .Machine[["double.eps"]], na.rm = TRUE) / cumprod(dim(t(Biobase::exprs(x))))[2] * 100, 2), "%"))
  graphics::text(1, 0.61, adj = 0, labels = paste0("min: ", signif(min(t(Biobase::exprs(x)), na.rm = TRUE), 2)))
  graphics::text(1, 0.54, adj = 0, labels = paste0("median: ", signif(stats::median(t(Biobase::exprs(x)), na.rm = TRUE), 2)))
  graphics::text(1, 0.47, adj = 0, labels = paste0("mean: ", signif(mean(t(Biobase::exprs(x)), na.rm = TRUE), 2)))
  graphics::text(1, 0.40, adj = 0, labels = paste0("max: ", signif(max(t(Biobase::exprs(x)), na.rm = TRUE), 2)))
  if ("sampleType" %in% colnames(Biobase::pData(x)) &&
      "pool" %in% Biobase::pData(x)[, "sampleType"])
    graphics::text(1,
                   0.33,
                   adj = 0,
                   labels = paste0("CVpool<",
                                   round(pool_cv.n * 100), "%: ",
                                   round(sum(Biobase::fData(x)[, "pool_CV"] < pool_cv.n, na.rm = TRUE) / nrow(Biobase::exprs(x)) * 100),
                                   "%"))
  
  ## sca: Color scale
  
  graphics::par(mar = marLs[["sca"]])
  
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
                 col = palHeaVc,
                 border = NA)
  
  eval(parse(text = paste0("axis(at = .prettyAxis(c(ifelse(min(t(Biobase::exprs(x)), na.rm = TRUE) == -Inf, yes = 0, no = min(t(Biobase::exprs(x)), na.rm = TRUE)) , max(t(Biobase::exprs(x)), na.rm = TRUE)), 256)$atVn,
                           font = 2,
                           font.axis = 2,
                           labels = .prettyAxis(c(ifelse(min(t(Biobase::exprs(x)), na.rm = TRUE) == -Inf, yes = 0, no = min(t(Biobase::exprs(x)), na.rm = TRUE)), max(t(Biobase::exprs(x)), na.rm = TRUE)), 256)$labVn,
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
  
  ## ima: Image
  
  graphics::par(mar = marLs[["ima"]])
  
  imaMN <- Biobase::exprs(x)[, rev(1:ncol(Biobase::exprs(x))), drop = FALSE]
  
  graphics::image(x = 1:nrow(imaMN),
                  y = 1:ncol(imaMN),
                  z = imaMN,
                  col = palHeaVc,
                  font.axis = 2,
                  font.lab = 2,
                  xaxt = "n",
                  yaxt = "n",
                  xlab = "",
                  ylab = "")
  
  if (length(rownames(t(Biobase::exprs(x)))) == 0) {
    rowNamVc <- rep("", times = ncol(Biobase::exprs(x)))
  } else
    rowNamVc <- rownames(t(Biobase::exprs(x)))
  
  if (length(colnames(t(Biobase::exprs(x)))) == 0) {
    colNamVc <- rep("", times = nrow(Biobase::exprs(x)))
  } else
    colNamVc <- colnames(t(Biobase::exprs(x)))
  
  xlaVc <- paste(paste(rep("[", 2),
                       c(1, nrow(imaMN)),
                       rep("] ", 2),
                       sep = ""),
                 rep("\n", times = 2),
                 c(colNamVc[1], utils::tail(colNamVc, 1)),
                 sep = "")
  
  for (k in 1:2)
    graphics::axis(side = 3,
                   hadj = c(0, 1)[k],
                   at = c(1, nrow(imaMN))[k],
                   cex = 0.8,
                   font = 2,
                   labels = xlaVc[k],
                   line = -0.5,
                   tick = FALSE)
  
  
  ylaVc <- paste(paste(rep("[", times = 2),
                       c(ncol(imaMN), 1),
                       rep("]", times = 2),
                       sep = ""),
                 rep("\n", times = 2),
                 c(utils::tail(rowNamVc, 1), rowNamVc[1]),
                 sep = "")
  
  for (k in 1:2)
    graphics::axis(side = 2,
                   at = c(1, ncol(imaMN))[k],
                   cex = 0.8,
                   font = 2,
                   hadj = c(0, 1)[k],
                   labels = ylaVc[k],
                   las = 0,
                   line = -0.5,
                   lty = "blank",
                   tick = FALSE)
  
  graphics::box(lwd = 2)
  
  ## dri: Analytical drift
  
  graphics::par(mar = marLs[["dri"]])
  
  ## ordering
  
  driProMN <- t(Biobase::exprs(x))
  driObsDF <- Biobase::pData(x)
  
  driObsDF[, "ordIniVi"] <- 1:nrow(driProMN)
  
  if ("injectionOrder" %in% colnames(driObsDF)) {
    ordNamC <- "Injection Order"
    if ("batch" %in% colnames(driObsDF)) {
      ordVi <- order(driObsDF[, "batch"],
                     driObsDF[, "injectionOrder"])
    } else
      ordVi <- order(driObsDF[, "injectionOrder"])
  } else {
    ordNamC <- "Samples"
    ordVi <- 1:nrow(driProMN)
  }
  
  driProMN <- driProMN[ordVi, ]
  driObsDF <- driObsDF[ordVi, ]
  driColVc <- obsColVc[ordVi]
  
  if ("batch" %in% colnames(driObsDF))
    batTab <- table(driObsDF[, "batch"])
  
  driSumVn <- eval(parse(text = paste0("apply(driProMN, 1, function(obsVn) ", sample_intensity.c, "(obsVn, na.rm = TRUE))")))
  
  graphics::plot(driSumVn,
                 col = driColVc,
                 pch = 18,
                 ylim = range(driProMN, na.rm = TRUE),
                 type = "n",
                 xaxs = "i",
                 xlab = "",
                 ylab = "")
  
  for (samI in 1:nrow(driProMN))
    graphics::boxplot(driProMN[samI, ],
                      at = samI,
                      add = TRUE)
  
  graphics::points(driSumVn,
                   col = driColVc,
                   pch = 18)
  
  graphics::mtext(ordNamC,
                  cex = 0.7,
                  line = 2,
                  side = 1)
  
  graphics::mtext(paste0(toupper(substr(sample_intensity.c, 1, 1)), substr(sample_intensity.c, 2, nchar(sample_intensity.c)), " of variable intensities"),
                  cex = 0.7,
                  line = 2,
                  side = 2)
  
  if ("batch" %in% colnames(driObsDF)) {
    
    graphics::abline(v = cumsum(batTab) + 0.5,
                     col = "red")
    
    graphics::mtext(names(batTab),
                    at = batTab / 2 + c(0, cumsum(batTab[-length(batTab)])),
                    cex = 0.7)
    
    for (batC in names(batTab)) {
      
      batSeqVi <- which(driObsDF[, "batch"] == batC)
      
      if ("sampleType" %in% colnames(driObsDF)) {
        batSamVi <- intersect(batSeqVi,
                              grep("sample", driObsDF[, "sampleType"]))
      } else
        batSamVi <- batSeqVi
      
      graphics::lines(batSeqVi,
                      .loess(driSumVn, batSamVi, batSeqVi, span.n),
                      col = .sample_color("sample"))
      
      if ("sampleType" %in% colnames(driObsDF) &&
          "pool" %in% driObsDF[, "sampleType"]) {
        
        batPooVi <- intersect(batSeqVi,
                              grep("^pool$", driObsDF[, "sampleType"]))
        
        graphics::lines(batSeqVi,
                        .loess(driSumVn, batPooVi, batSeqVi, span.n),
                        col = .sample_color("pool"))
        
      }
      
    }
    
  } else {
    
    batSeqVi <- 1:nrow(driObsDF)
    
    if ("sampleType" %in% colnames(driObsDF)) {
      batSamVi <- intersect(batSeqVi,
                            grep("sample", driObsDF[, "sampleType"]))
    } else
      batSamVi <- batSeqVi
    
    graphics::lines(batSeqVi,
                    .loess(driSumVn, batSamVi, batSeqVi, span.n),
                    col = .sample_color("sample"))
    
    if ("sampleType" %in% colnames(driObsDF) &&
        "pool" %in% driObsDF[, "sampleType"]) {
      
      batPooVi <- intersect(batSeqVi,
                            grep("^pool$", driObsDF[, "sampleType"]))
      
      graphics::lines(batSeqVi,
                      .loess(driSumVn, batPooVi, batSeqVi, span.n),
                      col = .sample_color("pool"))
      
      
    }
    
  }
  
  ## pca: PCA and Hotelling ellipse
  
  graphics::par(mar = marLs[["pca"]])
  
  graphics::plot(Biobase::pData(x)[, c("pca_sco1", "pca_sco2")],
                 type = "n",
                 xlab = "",
                 ylab = "",
                 xlim = range(Biobase::pData(x)[, "pca_sco1"]) * 1.1)
  graphics::mtext(paste("t1 (", round(pcaLs[["varRelVn"]][1] * 100), "%)", sep = ""),
                  cex = 0.7,
                  line = 2,
                  side = 1)
  graphics::mtext(paste("t2 (", round(pcaLs[["varRelVn"]][2] * 100), "%)", sep = ""),
                  cex = 0.7,
                  las = 0,
                  line = 2,
                  side = 2)
  graphics::abline(h = 0, lty = "dashed")
  graphics::abline(v = 0, lty = "dashed")
  radVn <- seq(0, 2 * pi, length.out = 100)
  
  hotFisN <- pcaLs[["hotN"]] * stats::qf(1 - 0.05, 2, ncol(Biobase::exprs(x)) - 2)
  
  graphics::lines(sqrt(stats::var(Biobase::pData(x)[, "pca_sco1"]) * hotFisN) * cos(radVn),
                  sqrt(stats::var(Biobase::pData(x)[, "pca_sco2"]) * hotFisN) * sin(radVn))
  
  graphics::text(Biobase::pData(x)[, "pca_sco1"],
                 Biobase::pData(x)[, "pca_sco2"],
                 cex = 0.7,
                 col = obsColVc,
                 labels = Biobase::sampleNames(x))
  
  if ("sampleType" %in% colnames(Biobase::pData(x))) {
    obsColVuc <- obsColVc[sort(unique(names(obsColVc)))]
    legOrdVc <- c("blank", paste0("pool", 8:1), "pool", "other", "sample")
    obsColVuc <- obsColVuc[legOrdVc[legOrdVc %in% names(obsColVuc)]]
    
    graphics::text(rep(graphics::par("usr")[1], times = length(obsColVuc)),
                   graphics::par("usr")[3] + (0.97 - length(obsColVuc) * 0.03 + 1:length(obsColVuc) * 0.03) * diff(graphics::par("usr")[3:4]),
                   col = obsColVuc,
                   font = 2,
                   labels = names(obsColVuc),
                   pos = 4)
  }
  
  graphics::par(mfrow = c(1, 1))
  graphics::par(opar)
  
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
