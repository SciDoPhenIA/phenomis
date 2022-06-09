#### hypotesting (MultiAssayExperiment) ####

#' @rdname hypotesting
#' @export
setMethod("hypotesting", signature(x = "MultiAssayExperiment"),
          function(x,
                   test.c = c("ttest", "limma", "wilcoxon",
                              "anova", "kruskal",
                              "pearson", "spearman",
                              "limma2ways", "limma2waysInter",
                              "anova2ways", "anova2waysInter")[2],
                   factor_names.vc,
                   factor_levels.ls = list(factor1.vc = "default", factor2.vc = "default"),
                   adjust.c = c("holm", "hochberg", "hommel", "bonferroni", "BH",
                                "BY", "fdr", "none")[5],
                   adjust_thresh.n = 0.05,
                   signif_maxprint.i = NA,
                   title.c = NA,
                   display_signif.l = FALSE,
                   prefix.c = "",
                   figure.c = c("none", "interactive", "interactive_plotly", "myfile.pdf")[2],
                   report.c = c("none", "interactive", "myfile.txt")[2]) {
            
            if (!(report.c %in% c("none", "interactive")))
              sink(report.c, append = TRUE)
            
            report_set.c <- report.c
            if (report_set.c != "none")
              report_set.c <- "interactive"
            
            if (!(figure.c %in% c("none", "interactive", "interactive_plotly"))) {
              grDevices::pdf(figure.c)
              figure_set.c <- "interactive"
            } else
              figure_set.c <- figure.c
            
            for (set.c in names(x)) {
              
              main.c <- paste0(ifelse(!is.na(title.c), paste0(title.c, ", "), ""),
                               set.c)
              
              if (report.c != "none")
                message("Hypothesis testing for the '", set.c, "' dataset:")
              
              x[[set.c]] <- hypotesting(x = x[[set.c]],
                                        test.c = test.c,
                                        factor_names.vc = factor_names.vc,
                                        factor_levels.ls = factor_levels.ls,
                                        adjust.c = adjust.c,
                                        adjust_thresh.n = adjust_thresh.n,
                                        signif_maxprint.i = signif_maxprint.i,
                                        title.c = main.c,
                                        display_signif.l = display_signif.l,
                                        prefix.c = prefix.c,
                                        figure.c = figure_set.c,
                                        report.c = report_set.c)
              
            }
            
            if (!(figure.c %in% c("none", "interactive", "interactive_plotly")))
              grDevices::dev.off()
            
            if (!(report.c %in% c("none", "interactive")))
              sink()
            
            methods::validObject(x)
            
            return(invisible(x))
            
          })


#### hypotesting (SummarizedExperiment) ####

#' @rdname hypotesting
#' @export
setMethod("hypotesting", signature(x = "SummarizedExperiment"),
          function(x,
                   test.c = c("ttest", "limma", "wilcoxon",
                              "anova", "kruskal",
                              "pearson", "spearman",
                              "limma2ways", "limma2waysInter",
                              "anova2ways", "anova2waysInter")[2],
                   factor_names.vc,
                   factor_levels.ls = list(factor1.vc = "default", factor2.vc = "default"),
                   adjust.c = c("holm", "hochberg", "hommel", "bonferroni", "BH",
                                "BY", "fdr", "none")[5],
                   adjust_thresh.n = 0.05,
                   signif_maxprint.i = NA,
                   title.c = NA,
                   display_signif.l = FALSE,
                   prefix.c = "",
                   figure.c = c("none", "interactive", "interactive_plotly", "myfile.pdf")[2],
                   report.c = c("none", "interactive", "myfile.txt")[2]) {
            
            if (!(report.c %in% c("none", "interactive")))
              sink(report.c, append = TRUE)
            
            # Checks
            
            factorVl <- factor_names.vc %in% colnames(colData(x))
            
            if (sum(!factorVl))
              stop("The following factor(s) '", paste(factor_names.vc[!factorVl], collapse = "', '"),
                   "' was/were not found in the column names of colData(x): '",
                   paste(colnames(colData(x)), collapse = "', '"), "'")
            
            testVc <- c("ttest", "limma", "wilcoxon",
                        "anova", "kruskal",
                        "pearson", "spearman",
                        "limma2ways", "limma2waysInter",
                        "anova2ways", "anova2waysInter")
            
            if (!test.c %in% testVc)
              stop("'test.c = ", test.c, ": 'test.c' must be in '", paste(testVc, collapse = "', '"), "'")
            
            # Hypothesis testing
            
            if (report.c != "none")
              message("Performing '", test.c, "'...")
            
            if (test.c %in% c("ttest", "limma", "wilcoxon",
                              "pearson", "spearman")) {
              
              testLs <- .twoSampCorTests(data.mn = t(assay(x)),
                                         samp.df = colData(x),
                                         feat.df = rowData(x),
                                         test.c = test.c,
                                         factorNameC = factor_names.vc,
                                         factorLevelsVc = factor_levels.ls[[1]],
                                         adjust.c = adjust.c,
                                         adjust_thresh.n = adjust_thresh.n,
                                         title.c = title.c,
                                         display_signif.l = display_signif.l,
                                         prefix.c = prefix.c,
                                         figure.c = figure.c)
              
            } else if (test.c %in% c("anova", "kruskal")) {
              
              testLs <- .anovas(data.mn = t(assay(x)),
                                samp.df = colData(x),
                                feat.df = rowData(x),
                                test.c = test.c,
                                factorNameC = factor_names.vc,
                                factorLevelsVc = factor_levels.ls[[1]],
                                adjust.c = adjust.c,
                                adjust_thresh.n = adjust_thresh.n,
                                title.c = title.c,
                                prefix.c = prefix.c,
                                figure.c = figure.c)
              
            } else if (test.c %in% c("limma2ways", "limma2waysInter",
                                     "anova2ways", "anova2waysInter")) {
              
              testLs <- .anovas2ways(data.mn = t(assay(x)),
                                     samp.df = colData(x),
                                     feat.df = rowData(x),
                                     test.c = test.c,
                                     factor_names.vc = factor_names.vc,
                                     factor_levels.ls = factor_levels.ls,
                                     adjust.c = adjust.c,
                                     adjust_thresh.n = adjust_thresh.n,
                                     title.c = title.c,
                                     prefix.c = prefix.c,
                                     figure.c = figure.c)
              
            }
            
            rowData(x) <- testLs[["feat.df"]]
            metric.mn <- testLs[["metric.mn"]]
            
            ## Printing significant features
            
            signiColVl <- grepl("_signif", colnames(metric.mn))
            signiVl <- rowSums(metric.mn[, signiColVl, drop = FALSE], na.rm = TRUE) > 0
            
            if (report.c != "none") {
              if (sum(signiVl) > 0) {
                
                metric.mn <- metric.mn[, !signiColVl, drop = FALSE]
                
                .signifPrint(metric.mn = metric.mn,
                             signiVl = signiVl,
                             signif_maxprint.i,
                             adjust.c = adjust.c,
                             adjust_thresh.n = adjust_thresh.n)
                
              } else
                message("\nNo significant variable found at the selected ", adjust_thresh.n, " level.")
              
            }
            
            if (!(report.c %in% c("none", "interactive")))
              sink()
            
            methods::validObject(x)
            
            return(invisible(x))
            
          })




#### hypotesting (MultiDataSet) ####

#' @rdname hypotesting
#' @export
setMethod("hypotesting", signature(x = "MultiDataSet"),
          function(x,
                   test.c = c("ttest", "limma", "wilcoxon",
                             "anova", "kruskal",
                             "pearson", "spearman",
                             "limma2ways", "limma2waysInter",
                             "anova2ways", "anova2waysInter")[2],
                   factor_names.vc,
                   factor_levels.ls = list(factor1.vc = "default", factor2.vc = "default"),
                   adjust.c = c("holm", "hochberg", "hommel", "bonferroni", "BH",
                               "BY", "fdr", "none")[5],
                   adjust_thresh.n = 0.05,
                   signif_maxprint.i = NA,
                   title.c = NA,
                   display_signif.l = FALSE,
                   prefix.c = "",
                   figure.c = c("none", "interactive", "interactive_plotly", "myfile.pdf")[2],
                   report.c = c("none", "interactive", "myfile.txt")[2]) {
            
            if (!(report.c %in% c("none", "interactive")))
              sink(report.c, append = TRUE)
            
            report_set.c <- report.c
            if (report_set.c != "none")
              report_set.c <- "interactive"
            
            if (!(figure.c %in% c("none", "interactive", "interactive_plotly"))) {
              grDevices::pdf(figure.c)
              figure_set.c <- "interactive"
            } else
              figure_set.c <- figure.c
            
            for (set.c in names(x)) {
              
              main.c <- paste0(ifelse(!is.na(title.c), paste0(title.c, ", "), ""),
                               set.c)
              
              if (report.c != "none")
                message("Hypothesis testing for the '", set.c, "' dataset:")
              
              ese <- x[[set.c]]
              
              ese <- hypotesting(x = ese,
                                 test.c = test.c,
                                 factor_names.vc = factor_names.vc,
                                 factor_levels.ls = factor_levels.ls,
                                 adjust.c = adjust.c,
                                 adjust_thresh.n = adjust_thresh.n,
                                 signif_maxprint.i = signif_maxprint.i,
                                 title.c = main.c,
                                 display_signif.l = display_signif.l,
                                 prefix.c = prefix.c,
                                 figure.c = figure_set.c,
                                 report.c = report_set.c)
              
              x <- MultiDataSet::add_eset(x,
                                          ese,
                                          dataset.type = set.c,
                                          GRanges = NA,
                                          overwrite = TRUE,
                                          warnings = FALSE)
              
            }
            
            if (!(figure.c %in% c("none", "interactive", "interactive_plotly")))
              grDevices::dev.off()
            
            if (!(report.c %in% c("none", "interactive")))
              sink()
            
            methods::validObject(x)
            
            return(invisible(x))
            
          })


#### hypotesting (ExpressionSet) ####

#' @rdname hypotesting
#' @export
setMethod("hypotesting", signature(x = "ExpressionSet"),
          function(x,
                   test.c = c("ttest", "limma", "wilcoxon",
                             "anova", "kruskal",
                             "pearson", "spearman",
                             "limma2ways", "limma2waysInter",
                             "anova2ways", "anova2waysInter")[2],
                   factor_names.vc,
                   factor_levels.ls = list(factor1.vc = "default", factor2.vc = "default"),
                   adjust.c = c("holm", "hochberg", "hommel", "bonferroni", "BH",
                               "BY", "fdr", "none")[5],
                   adjust_thresh.n = 0.05,
                   signif_maxprint.i = NA,
                   title.c = NA,
                   display_signif.l = FALSE,
                   prefix.c = "",
                   figure.c = c("none", "interactive", "interactive_plotly", "myfile.pdf")[2],
                   report.c = c("none", "interactive", "myfile.txt")[2]) {
            
            if (!(report.c %in% c("none", "interactive")))
              sink(report.c, append = TRUE)
            
            # Checks
            
            factorVl <- factor_names.vc %in% colnames(Biobase::pData(x))
            
            if (sum(!factorVl))
              stop("The following factor(s) '", paste(factor_names.vc[!factorVl], collapse = "', '"),
                   "' was/were not found in the column names of pData(x): '",
                   paste(colnames(Biobase::pData(x)), collapse = "', '"), "'")
            
            testVc <- c("ttest", "limma", "wilcoxon",
                        "anova", "kruskal",
                        "pearson", "spearman",
                        "limma2ways", "limma2waysInter",
                        "anova2ways", "anova2waysInter")
            
            if (!test.c %in% testVc)
              stop("'test.c = ", test.c, ": 'test.c' must be in '", paste(testVc, collapse = "', '"), "'")
            
            # Hypothesis testing
            
            if (report.c != "none")
              message("Performing '", test.c, "'...")
            
            if (test.c %in% c("ttest", "limma", "wilcoxon",
                             "pearson", "spearman")) {
              
              testLs <- .twoSampCorTests(data.mn = t(Biobase::exprs(x)),
                                         samp.df = Biobase::pData(x),
                                         feat.df = Biobase::fData(x),
                                         test.c = test.c,
                                         factorNameC = factor_names.vc,
                                         factorLevelsVc = factor_levels.ls[[1]],
                                         adjust.c = adjust.c,
                                         adjust_thresh.n = adjust_thresh.n,
                                         title.c = title.c,
                                         display_signif.l = display_signif.l,
                                         prefix.c = prefix.c,
                                         figure.c = figure.c)
              
            } else if (test.c %in% c("anova", "kruskal")) {
              
              testLs <- .anovas(data.mn = t(Biobase::exprs(x)),
                                samp.df = Biobase::pData(x),
                                feat.df = Biobase::fData(x),
                                test.c = test.c,
                                factorNameC = factor_names.vc,
                                factorLevelsVc = factor_levels.ls[[1]],
                                adjust.c = adjust.c,
                                adjust_thresh.n = adjust_thresh.n,
                                title.c = title.c,
                                prefix.c = prefix.c,
                                figure.c = figure.c)
              
            } else if (test.c %in% c("limma2ways", "limma2waysInter",
                                    "anova2ways", "anova2waysInter")) {
              
              testLs <- .anovas2ways(data.mn = t(Biobase::exprs(x)),
                                     samp.df = Biobase::pData(x),
                                     feat.df = Biobase::fData(x),
                                     test.c = test.c,
                                     factor_names.vc = factor_names.vc,
                                     factor_levels.ls = factor_levels.ls,
                                     adjust.c = adjust.c,
                                     adjust_thresh.n = adjust_thresh.n,
                                     title.c = title.c,
                                     prefix.c = prefix.c,
                                     figure.c = figure.c)
              
            }
            
            Biobase::fData(x) <- testLs[["feat.df"]]
            metric.mn <- testLs[["metric.mn"]]
            
            ## Printing significant features
            
            signiColVl <- grepl("_signif", colnames(metric.mn))
            signiVl <- rowSums(metric.mn[, signiColVl, drop = FALSE], na.rm = TRUE) > 0
            
            if (report.c != "none") {
              if (sum(signiVl) > 0) {
                
                metric.mn <- metric.mn[, !signiColVl, drop = FALSE]
                
                .signifPrint(metric.mn = metric.mn,
                             signiVl = signiVl,
                             signif_maxprint.i,
                             adjust.c = adjust.c,
                             adjust_thresh.n = adjust_thresh.n)
                
              } else
                message("\nNo significant variable found at the selected ", adjust_thresh.n, " level.")
              
            }
            
            if (!(report.c %in% c("none", "interactive")))
              sink()
            
            methods::validObject(x)
            
            return(invisible(x))
            
          })


# ttest, wilcoxon, limma, pearson, spearman: test
.twoSampCorTests <- function(data.mn, ## data (matrix of numerics; samples x variables)
                             samp.df, ## sample metadata (dataframe; samples x metadata)
                             feat.df, ## feature metadata (dataframe; features x metadata)
                             test.c,
                             factorNameC,
                             factorLevelsVc,
                             adjust.c,
                             adjust_thresh.n,
                             title.c,
                             display_signif.l,
                             prefix.c,
                             figure.c) {
  
  checkLs <- .oneFactorCheck(samp.df = samp.df,
                             test.c = test.c,
                             factorNameC = factorNameC,
                             factorLevelsVc = factorLevelsVc)
  
  factorVc <- checkLs[["factorVc"]]
  factorFc <- checkLs[["factorFc"]]
  factorVn <- checkLs[["factorVn"]]
  
  statVn <- .twoSampCorStat(data.mn = data.mn,
                            test.c = test.c,
                            factorFc = factorFc,
                            factorVn = factorVn)
  
  pvalVn <- .twoSampCorPval(data.mn = data.mn,
                            test.c = test.c,
                            factorFc = factorFc,
                            factorVn = factorVn)
  
  pvalAdjustVn <- stats::p.adjust(pvalVn, adjust.c)
  
  adjustSigniVn <- as.numeric(pvalAdjustVn <= adjust_thresh.n)
  
  metric.mn <- cbind(statVn,
                     pvalAdjustVn,
                     adjustSigniVn)
  
  rownames(metric.mn) <- colnames(data.mn)
  colnames(metric.mn) <- .twoSampCorNames(test.c,
                                          factorNameC,
                                          factorFc,
                                          adjust.c,
                                          prefix.c)
  
  for (colC in colnames(metric.mn))
    feat.df[, colC] <- metric.mn[, colC]
  
  metric.mn <- metric.mn[order(pvalVn), ,
                         drop = FALSE]
  
  if (figure.c != "none") {
    
    if (!grepl("interactive", figure.c))
      grDevices::pdf(figure.c)
    
    .twoSampCorPlot(data.mn = data.mn,
                    test.c = test.c,
                    factorNameC = factorNameC,
                    factorFc = factorFc,
                    factorVn = factorVn,
                    metric.mn = metric.mn,
                    adjust.c = adjust.c,
                    adjust_thresh.n = adjust_thresh.n,
                    title.c = title.c,
                    display_signif.l = display_signif.l,
                    figure.c = figure.c)
    
    if (!grepl("interactive", figure.c))
      grDevices::dev.off()
    
  }
  
  return(list(feat.df = feat.df,
              metric.mn = metric.mn))
  
}

# ttest, wilcoxon, limma, pearson, spearman, anova, kruskal: checking factor name and levels
.oneFactorCheck <- function(samp.df,
                            test.c,
                            factorNameC,
                            factorLevelsVc) {
  
  factorLevelsSelectedVc <- factorLevelsVc
  
  factorNameLengthN <- length(factorNameC)
  if (factorNameLengthN != 1)
    stop("'length(factorNameC) = ", factorNameLengthN, "': A single 'factorName' is required")
  
  factorVcn <- samp.df[, factorNameC]
  # either of 'character' or 'numeric' mode at this stage
  
  if (test.c %in% c("ttest", "limma", "wilcoxon", "anova", "kruskal",
                   "anova2ways",
                   "anova2waysInter",
                   "limma2ways",
                   "limma2waysInter")) {
    
    factorVc <- factorVcn
    factorVn <- NULL
    
    if (is.factor(factorVc)) {
      factorFc <- factorVc
      factorVc <- as.character(factorVc)
    } else {
      factorModeC <- mode(factorVc)
      if (factorModeC != "character")
        stop("'mode(factorVc) = ", factorModeC, "': Factor must be of 'character' mode.")
      factorFc <- factor(factorVc)
    }
    
    factorLevelsDefaultVc <- levels(factorFc)
    levelDotDashVl <- grepl(".", factorLevelsDefaultVc, fixed = TRUE) |
      grepl("-", factorLevelsDefaultVc, fixed = TRUE)
    if (sum(levelDotDashVl))
      stop("'Factor levels must not contain '.' or '-' characters.")
    
    factorLevelsDefaultN <- nlevels(factorFc)
    
    if (test.c %in% c("ttest", "limma", "wilcoxon") && factorLevelsDefaultN != 2) {
      stop("'nlevels(factorFc) = ", factorLevelsDefaultN, "': The factor must have 2 levels.")
    } else if (test.c %in% "anova" && factorLevelsDefaultN < 3) {
      stop("'nlevels(factorFc) = ", factorLevelsDefaultN, "': The factor must have 3 levels or more.")
    }
    
    if (length(factorLevelsSelectedVc) > 1 || factorLevelsSelectedVc != "default") {
      
      factorLevelsSelectedN <- length(factorLevelsSelectedVc)
      if (factorLevelsDefaultN != factorLevelsSelectedN)
        stop("'nlevels(factorFc) = ", factorLevelsDefaultN, "'and 'length(factorLevelsVc) = ", factorLevelN,
             "': 'factorLevelsVc' must have the same number of levels as the factor")
      
      if (!identical(sort(factorLevelsVc), sort(factorLevelsDefaultVc)))
        stop("'levels(factorFc) = ", paste(factorLevelsDefaultVc, collapse = ', '),
             "' and 'factorLevelsVc = ", paste(factorLevelsSelectedVc, collapse = ', '),
             "': refLeveslVc' must match all level names of the factor.")
      
      factorFc <- factor(factorVc, factorLevelsSelectedVc)
      
    }
    
  } else if (test.c %in% c("pearson", "spearman")) {
    
    factorVn <- factorVcn
    factorFc <- factorVc <- NULL
    
    factorModeC <- mode(factorVn)
    if (factorModeC != "numeric")
      stop("'mode(factorVn) = ", factorModeC, "': Factor must be of 'numeric' mode.")
    
  }
  
  return(list(factorVc = factorVc,
              factorFc = factorFc,
              factorVn = factorVn))
  
}

.twoSampDiffMean <- function(data.mn,
                             factorFc) {
  apply(data.mn, 2, function(y)
    diff(tapply(y, factorFc, function(x) mean(x, na.rm = TRUE))))
}

.twoSampDiffMedian <- function(data.mn,
                               factorFc) {
  apply(data.mn, 2, function(y)
    diff(tapply(y, factorFc, function(x) stats::median(x, na.rm = TRUE))))
}

.cor <- function(data.mn,
                 factorVn,
                 methodC)
  apply(data.mn, 2, function(y)
    stats::cor(factorVn, y, method = methodC, use = "pairwise.complete.obs"))

.twoSampCorStat <- function(data.mn,
                            test.c,
                            factorFc,
                            factorVn) {
  
  if (test.c %in% c("ttest", "limma")) {
    return(.twoSampDiffMean(data.mn = data.mn, factorFc = factorFc))
  } else if (test.c == "wilcoxon") {
    return(.twoSampDiffMedian(data.mn = data.mn, factorFc = factorFc))
  } else if (test.c %in% c("pearson", "spearman")) {
    return(.cor(data.mn = data.mn, factorVn = factorVn, methodC = test.c))
  }
  
}

.twoSampCorPval <- function(data.mn,
                            test.c,
                            factorFc,
                            factorVn) {
  
  if (test.c %in% c("ttest", "wilcoxon", "pearson", "spearman")) {
    
    if (test.c == "ttest") {
      hypotest <- function(y) stats::t.test(y ~ factorFc)[["p.value"]]
    } else if (test.c == "wilcoxon") {
      hypotest <- function(y) stats::wilcox.test(y ~ factorFc, exact = FALSE)[["p.value"]]
    } else if (test.c %in% c("pearson", "spearman")) {
      hypotest <- function(y) stats::cor.test(factorVn, y, method = test.c,
                                              use = "pairwise.complete.obs", exact = FALSE)[["p.value"]]
    }
    
    return(apply(data.mn, 2, hypotest))
    
  } else if (test.c == "limma") {
    designMN <- cbind(level2 = rep(1, nrow(data.mn)),
                      level1.level2 = as.numeric(as.numeric(factorFc) == 1))
    rownames(designMN) <- rownames(data.mn)
    limmaFitLs <- limma::lmFit(t(data.mn), designMN)
    limmaBayesLs <- limma::eBayes(limmaFitLs)
    
    return(limmaBayesLs[["p.value"]][, "level1.level2"])
    
  }
  
}


.twoSampCorNames <- function(test.c,
                             factorNameC,
                             factorFc,
                             adjust.c,
                             prefix.c) {
  
  prefix.c <- paste0(prefix.c, test.c, "_", make.names(factorNameC), "_")
  
  if (test.c %in% c("ttest", "limma", "wilcoxon",
                    "limma2ways", "limma2waysInter",
                    "anova2ways", "anova2waysInter")) {
    
    prefix.c <- paste0(prefix.c,
                      paste(levels(factorFc), collapse = "."), "_")
    statC <- "diff"
    
  } else if (test.c %in% c("pearson", "spearman")) {
    
    statC <- "cor"
    
  } else
    stop()
  
  paste0(prefix.c, c(statC, adjust.c, "signif"))
  
}


# ttest, wilcoxon, limma: volcano and boxplots
.twoSampCorPlot <- function(data.mn,
                            test.c,
                            factorNameC,
                            factorFc,
                            factorVn,
                            metric.mn,
                            adjust.c,
                            adjust_thresh.n,
                            title.c,
                            display_signif.l,
                            figure.c) {
  
  statC <- ifelse(test.c %in% c("pearson", "spearman"), "cor", "diff")
  
  adjPvalVn <- metric.mn[, grep(adjust.c, colnames(metric.mn))]
  signiVl <- adjPvalVn <= adjust_thresh.n
  signiVl[is.na(signiVl)] <- FALSE
  
  mainC <- paste0(ifelse(!is.na(title.c), paste0(title.c, "\n"), ""),
                  test.c, " '", factorNameC, "', ",
                  adjust.c, " (signif.: ",
                  format(sum(signiVl),
                         big.mark = ","),
                  ")")
  
  if (!is.null(factorFc)) {
    group.vc <- levels(factorFc)
  } else
    group.vc = ""
  
  label.vc <- rownames(metric.mn)
  label.vc[!signiVl] <- ""
  
  if (sum(signiVl) > 30)
    label.vc[adjPvalVn > max(head(sort(adjPvalVn), 30))] <- ""
  
  gg_volcanoplot(fold_change.vn = metric.mn[, grep(paste0("_", statC), colnames(metric.mn))],
                 adjusted_pvalue.vn = adjPvalVn,
                 adjust_method.c = adjust.c,
                 adjust_thresh.n = adjust_thresh.n,
                 label.vc = label.vc,
                 title.c = mainC,
                 xlab.c = switch(statC,
                                 diff = "Fold Change",
                                 cor = "Correlation"),
                 class_name.vc = group.vc,
                 figure.c = figure.c)
  
  
  if (display_signif.l) {
    
    varSigniVi <- which(metric.mn[, grep(paste0("_signif$"), colnames(metric.mn))] > 0)
    
    for (varI in varSigniVi) {
      
      varC <- rownames(metric.mn)[varI]
      
      if (test.c %in% c("pearson", "spearman")) {
        
        mod <- stats::lm(data.mn[, varC] ~  factorVn)
        
        mainC <- paste0(varC, "\n(", statC, " = ",
                        signif(metric.mn[varI,
                                        grep(paste0("_", statC, "$"), colnames(metric.mn))], 2),
                        ", ",
                        adjust.c, " = ",
                        signif(metric.mn[varI,
                                        grep(paste0("_", adjust.c, "$"), colnames(metric.mn))], 2),
                        ", R2 = ", signif(summary(mod)$r.squared, 2), ")")
        if (!is.na(title.c))
          mainC <- paste0(title.c, ", ", mainC)
        
        graphics::plot(factorVn, data.mn[, varC],
                       xlab = factorNameC,
                       ylab = "",
                       pch = 18,
                       main = mainC)
        
        graphics::abline(mod, col = "red")
        
      } else {
        mainC <- paste0(varC, "\n(", statC, " = ",
                        signif(metric.mn[varI,
                                        grep(paste0("_", statC, "$"),
                                             colnames(metric.mn))], 2),
                        ", ",
                        adjust.c, " = ",
                        signif(metric.mn[varI,
                                        grep(paste0("_", adjust.c, "$"),
                                             colnames(metric.mn))], 2), ")")
        if (!is.na(title.c))
          mainC <- paste0(title.c, ", ", mainC)
        
        gg_boxplot(data.frame(factor = factorFc,
                              response = data.mn[, varC]),
                   x.c = "factor",
                   y.c = "response",
                   color.c = "factor",
                   title.c = mainC,
                   xlab.c = factorNameC)
        
        
      }
      
    }
    
  }
  
}

.anovas <- function(data.mn, ## data (matrix of numerics; samples x variables)
                    samp.df, ## sample metadata (dataframe; samples x metadata)
                    feat.df, ## feature metadata (dataframe; features x metadata)
                    test.c,
                    factorNameC,
                    factorLevelsVc,
                    adjust.c,
                    adjust_thresh.n,
                    title.c,
                    prefix.c,
                    figure.c) {
  
  checkLs <- .oneFactorCheck(samp.df = samp.df,
                             test.c = test.c,
                             factorNameC = factorNameC,
                             factorLevelsVc = factorLevelsVc)
  
  factorVc <- checkLs[["factorVc"]]
  factorFc <- checkLs[["factorFc"]]
  
  # metric and pairwise names
  namesLs <- .anovasNames(test.c = test.c,
                          factorNameC = factorNameC,
                          factorFc = factorFc,
                          adjust.c = adjust.c,
                          prefix.c = prefix.c)
  pairNamesVc <- namesLs[["pairNamesVc"]]
  metricNamesVc <- namesLs[["aovNamesVc"]]
  
  ## omnibus and post-hoc tests
  
  diffPvalLs <- .anovasDiffPval(data.mn = data.mn,
                                test.c = test.c,
                                factorFc = factorFc,
                                pairNamesVc = pairNamesVc)
  
  pvalVn <- diffPvalLs[["pvalVn"]]
  pairDiffMN <- diffPvalLs[["pairDiffMN"]]
  pairPvalMN <- diffPvalLs[["pairPvalMN"]]
  
  # main ANOVA test
  pvalAdjustVn <- stats::p.adjust(pvalVn, method = adjust.c)
  pvalSigniVi <- as.integer(pvalAdjustVn <= adjust_thresh.n)
  
  # pairwise test
  pairAdjustMN <- apply(pairPvalMN, 2,
                        function(pairPvalVn) stats::p.adjust(pairPvalVn, method = adjust.c))
  pairSigniMI <- pairAdjustMN <= adjust_thresh.n
  mode(pairSigniMI) <- "integer"
  
  metric.mn <- cbind(pvalAdjustVn,
                    pvalSigniVi)
  for (pairI in 1:ncol(pairAdjustMN)) {
    metric.mn <- cbind(metric.mn,
                      pairDiffMN[, pairI],
                      pairAdjustMN[, pairI],
                      pairSigniMI[, pairI])
  }
  rownames(metric.mn) <- colnames(data.mn)
  colnames(metric.mn) <- metricNamesVc
  
  for (colC in colnames(metric.mn))
    feat.df[, colC] <- metric.mn[, colC]
  
  metric.mn <- metric.mn[order(pvalVn), , drop = FALSE]
  
  ## graphic
  
  if (figure.c != "none") {
    
    if (!grepl("interactive", figure.c))
      grDevices::pdf(figure.c)
    
    .anovasPlot(data.mn = data.mn,
                factorNameC = factorNameC,
                factorFc = factorFc,
                metric.mn = metric.mn,
                pairNamesVc = pairNamesVc,
                adjust.c = adjust.c,
                adjust_thresh.n = adjust_thresh.n,
                title.c = title.c,
                test.c = test.c,
                figure.c = figure.c)
    
    if (!grepl("interactive", figure.c))
      grDevices::dev.off()
    
  }
  
  return(list(feat.df = feat.df,
              metric.mn = metric.mn))
  
}

.anovasDiffPval <- function(data.mn,
                            test.c,
                            factorFc,
                            pairNamesVc) {
  
  if (test.c == "anova") {
    
    aovMN <- t(apply(data.mn, 2,
                     function(varVn) {
                       
                       aovModel <- stats::aov(varVn ~ factorFc)
                       aovPvalN <- summary(aovModel)[[1]][1, "Pr(>F)"]
                       tukeyHsdMN <- stats::TukeyHSD(aovModel)[["factorFc"]]
                       
                       padjVn <- diffVn <- rep(NA_real_, length(pairNamesVc))
                       names(padjVn) <- names(diffVn) <- sapply(pairNamesVc, function(pairC) {
                         paste(rev(unlist(strsplit(pairC, ".", fixed = TRUE))), collapse = "-")
                       })
                       
                       diffVn[rownames(tukeyHsdMN)] <- tukeyHsdMN[, "diff"]
                       padjVn[rownames(tukeyHsdMN)] <- tukeyHsdMN[, "p adj"]
                       c(aovPval = aovPvalN,
                         diffVn,
                         padjVn)
                       
                     }))
    
    # main ANOVA test
    pvalVn <- aovMN[, "aovPval"]
    
    # pairwise differences
    pairI <- as.integer(nlevels(factorFc) * (nlevels(factorFc) - 1) / 2)
    pairDiffMN <- aovMN[, 1 + 1:pairI, drop = FALSE]
    
    # pairwise test
    pairPvalMN <- aovMN[, (2 + pairI):ncol(aovMN), drop = FALSE]
    
  } else if (test.c == "kruskal") {
    
    kruskalMN <- t(apply(data.mn, 2, function(varVn) {
      
      kruskalPvalN <- stats::kruskal.test(varVn ~ factorFc)[["p.value"]]
      nemenyiPvalMN <- PMCMRplus::kwAllPairsNemenyiTest(varVn, factorFc, "Tukey")[["p.value"]]
      
      nemenyiPvalVn <- c(nemenyiPvalMN[lower.tri(nemenyiPvalMN, diag = TRUE)])
      names(nemenyiPvalVn) <- paste0(rep(rownames(nemenyiPvalMN), ncol(nemenyiPvalMN)),
                                     "-", rep(colnames(nemenyiPvalMN), each = nrow(nemenyiPvalMN)))[c(lower.tri(nemenyiPvalMN, diag = T))]
      padjVn <- rep(NA_real_, length(pairNamesVc))
      names(padjVn) <- sapply(pairNamesVc, function(pairC) {
        paste(rev(unlist(strsplit(pairC, ".", fixed = TRUE))), collapse = "-")
      })
      
      padjVn[names(nemenyiPvalVn)] <- nemenyiPvalVn
      
      c(kruskalPval = kruskalPvalN,
        padjVn)
      
    }))
    
    # main test
    pvalVn <- kruskalMN[, "kruskalPval"]
    
    ## difference of the medians for each pairwise comparison
    
    pairDiffMN <- vapply(pairNamesVc, function(pairC) {
      pairVc <- unlist(strsplit(pairC, ".", fixed = TRUE))
      pairSamplesVi <- which(factorFc %in% pairVc)
      pairFactorFc <- factor(as.character(factorFc)[pairSamplesVi], levels = pairVc)
      apply(data.mn[pairSamplesVi, ], 2,
            function(varVn)
              diff(as.numeric(tapply(varVn, pairFactorFc,
                                     function(x)
                                       stats::median(x, na.rm = TRUE)))))
    }, numeric(length(pvalVn)))
    
    # pairwise test
    pairPvalMN <- kruskalMN[, -1, drop = FALSE]
    
  }
  
  return(list(pvalVn = pvalVn,
              pairDiffMN = pairDiffMN,
              pairPvalMN = pairPvalMN))
  
}

.anovasNames <- function(test.c,
                         factorNameC,
                         factorFc,
                         adjust.c,
                         prefix.c) {
  
  prefix.c <- paste0(prefix.c, test.c, "_", make.names(factorNameC), "_")
  
  # getting the names of the pairwise comparisons
  factLevVc <- levels(factorFc)
  pairMC <- matrix("", nrow = length(factLevVc), ncol = length(factLevVc))
  for (i in 1:length(factLevVc)) {
    for (j in i:length(factLevVc))
      pairMC[i, j] <- paste0(factLevVc[i], ".", factLevVc[j])
  }
  pairMC <- t(pairMC)
  pairNamesVc <- pairMC[lower.tri(pairMC)]
  
  aovNamesVc <- paste0(prefix.c,
                       c(adjust.c, "signif"))
  
  for (pairC in pairNamesVc)
    aovNamesVc <- c(aovNamesVc,
                    paste0(prefix.c, pairC, c("_diff", paste0("_", adjust.c), "_signif")))
  
  return(list(pairNamesVc = pairNamesVc,
              aovNamesVc = aovNamesVc))
  
}

.anovasPlot <- function(data.mn,
                        factorNameC,
                        factorFc,
                        metric.mn,
                        pairNamesVc,
                        adjust.c,
                        adjust_thresh.n,
                        title.c,
                        test.c,
                        figure.c) {
  
  for (pairC in pairNamesVc) {
    
    pairmetric.mn <- metric.mn[, grep(pairC, colnames(metric.mn), fixed = TRUE)]
    
    pairFactorFc <- factorFc[factorFc %in% unlist(strsplit(pairC, ".", fixed = TRUE)), drop = TRUE]
    
    pairAdjPvalVn <- pairmetric.mn[, grep(adjust.c, colnames(pairmetric.mn))]
    pairSigniVl <- pairAdjPvalVn <= adjust_thresh.n
    pairSigniVl[is.na(pairSigniVl)] <- FALSE
    
    mainC <- paste0(ifelse(!is.na(title.c), paste0(title.c, "\n"), ""),
                    test.c, " '", factorNameC, "' (", pairC, ") ",
                    adjust.c, " (signif.: ",
                    format(sum(pairSigniVl),
                           big.mark = ","),
                    ")")
    
    if (!is.null(pairFactorFc)) {
      group.vc <- levels(pairFactorFc)
    } else
      group.vc = ""
    
    label.vc <- rownames(pairmetric.mn)
    label.vc[!pairSigniVl] <- ""
    
    gg_volcanoplot(fold_change.vn = pairmetric.mn[, 1],
                   adjusted_pvalue.vn = pairmetric.mn[, 2],
                   adjust_method.c = adjust.c,
                   adjust_thresh.n = adjust_thresh.n,
                   label.vc = label.vc,
                   title.c = mainC,
                   xlab.c = "Fold Change",
                   class_name.vc = group.vc,
                   figure.c = figure.c)
    
  }
  
  pvalAdjVn <- metric.mn[, paste0(test.c, "_", factorNameC, "_", adjust.c)]
  pvalSigniVi <- which(metric.mn[, paste0(test.c, "_", factorNameC, "_signif")] > 0)
  
  if (sum(pvalSigniVi)) {
    
    for (varI in pvalSigniVi) {
      
      varNameC <- rownames(metric.mn)[varI]
      
      mainC <- paste0(varNameC, " (", adjust.c, " = ",
                      signif(pvalAdjVn[varI], 2), ")")
      if (!is.na(title.c))
        mainC <- paste0(title.c, ", ", mainC)
      
      graphics::boxplot(data.mn[, varNameC] ~ factorFc,
                        main = mainC,
                        xlab = factorNameC, ylab = "")
      
    }
  }
}

# ANOVA 2 ways without interaction (2 levels for each factor)
.anovas2ways <- function(data.mn, ## data (matrix of numerics; samples x variables)
                         samp.df, ## sample metadata (dataframe; samples x metadata)
                         feat.df, ## feature metadata (dataframe; features x metadata)
                         test.c,
                         factor_names.vc,
                         factor_levels.ls,
                         adjust.c,
                         adjust_thresh.n,
                         title.c,
                         prefix.c,
                         figure.c) {
  
  checkLs <- .twoFactorsCheck(samp.df,
                              test.c,
                              factor_names.vc,
                              factor_levels.ls,
                              adjust.c,
                              adjust_thresh.n)
  
  factor1.vc <- checkLs[["factor1.vc"]]
  factor1.fc <- checkLs[["factor1.fc"]]
  factor2.vc <- checkLs[["factor2.vc"]]
  factor2.fc <- checkLs[["factor2.fc"]]
  
  diffPvalLs <- .anovas2waysDiffPval(data.mn = data.mn,
                                     test.c = test.c,
                                     factor1.fc = factor1.fc,
                                     factor2.fc = factor2.fc)
  
  diffMN <- diffPvalLs[["diffMN"]]
  pvalMN <- diffPvalLs[["pvalMN"]]
  
  pvalAdjustMN <- apply(pvalMN, 2,
                        function(pvalVn)
                          stats::p.adjust(pvalVn, adjust.c))
  
  adjustSigniMN <- apply(pvalAdjustMN, 2,
                         function(adjustVn)
                           as.numeric(adjustVn <= adjust_thresh.n))
  
  metric.mn <- NULL
  
  for (factorI in 1:2) {
    metric.mn <- cbind(metric.mn,
                      diffMN[, factorI],
                      pvalAdjustMN[, factorI],
                      adjustSigniMN[, factorI])
  }
  
  if (grepl("Inter$", test.c))
    metric.mn <- cbind(metric.mn, pvalAdjustMN[, 3], adjustSigniMN[, 3])
  
  rownames(metric.mn) <- colnames(data.mn)
  colnames(metric.mn) <- .anovas2waysNames(test.c = test.c,
                                           factor_names.vc = factor_names.vc,
                                           factor1.fc = factor1.fc,
                                           factor2.fc = factor2.fc,
                                           adjust.c = adjust.c,
                                           prefix.c = prefix.c)
  
  for (colC in colnames(metric.mn))
    feat.df[, colC] <- metric.mn[, colC]
  
  metric.mn <- metric.mn[do.call("order", as.list(as.data.frame(pvalAdjustMN))), ,
                       drop = FALSE]
  
  ## graphic
  
  if (figure.c != "none")
    .anovas2waysPlot(data.mn = data.mn,
                     test.c = test.c,
                     factor_names.vc = factor_names.vc,
                     adjust.c = adjust.c,
                     adjust_thresh.n = adjust_thresh.n,
                     metric.mn = metric.mn,
                     title.c = title.c,
                     figure.c = figure.c)
  
  return(list(feat.df = feat.df,
              metric.mn = metric.mn))
  
}

# anova2ways: checking factor names and levels
.twoFactorsCheck <- function(samp.df,
                             test.c,
                             factor_names.vc,
                             factor_levels.ls,
                             adjust.c,
                             adjust_thresh.n) {
  
  if (length(factor_names.vc) != 2)
    stop("Two factors are required.")
  
  if (!is.list(factor_levels.ls) || length(factor_levels.ls) != 2)
    stop("'factor_levels.ls' must be a list of length 2.", call. = FALSE)
  
  factorLs <- vector(mode = "list", length = 2)
  
  for (i in 1:length(factor_names.vc)) {
    
    factorLs[[i]] <- .oneFactorCheck(samp.df = samp.df,
                                     test.c = test.c,
                                     factorNameC = factor_names.vc[i],
                                     factorLevelsVc = factor_levels.ls[[i]])
    
  }
  
  factor1.vc <- factorLs[[1]][["factorVc"]]
  factor1.fc <- factorLs[[1]][["factorFc"]]
  factor2.vc <- factorLs[[2]][["factorVc"]]
  factor2.fc <- factorLs[[2]][["factorFc"]]
  
  if (nlevels(factor1.fc) != 2 || nlevels(factor2.fc) != 2)
    stop("Only factors with 2 levels are handled by the current implementation
         of '2ways' tests")
  
  return(list(factor1.vc = factor1.vc,
              factor1.fc = factor1.fc,
              factor2.vc = factor2.vc,
              factor2.fc = factor2.fc))
  
}

.anovas2waysDiffPval <- function(data.mn,
                                 test.c,
                                 factor1.fc,
                                 factor2.fc) {
  
  # Difference of means
  
  anovas2waysDiffMN <- cbind(.twoSampDiffMean(data.mn = data.mn,
                                              factorFc = factor1.fc),
                             .twoSampDiffMean(data.mn = data.mn,
                                              factorFc = factor2.fc))
  
  # Test statistic
  
  if (test.c == "anova2ways") {
    
    anova2waysPvalMN <- t(apply(data.mn, 2,
                                function(varVn) {
                                  aovModel <- stats::lm(varVn ~ factor1.fc + factor2.fc)
                                  stats::anova(aovModel)[c("factor1.fc",
                                                           "factor2.fc"), "Pr(>F)"]
                                }))
    
  } else if (test.c == "anova2waysInter") {
    
    anova2waysPvalMN <- t(apply(data.mn, 2,
                                function(varVn) {
                                  aovModel <- stats::lm(varVn ~ factor1.fc * factor2.fc)
                                  stats::anova(aovModel)[c("factor1.fc",
                                                           "factor2.fc",
                                                           "factor1.fc:factor2.fc"), "Pr(>F)"]
                                }))
    
  } else if (test.c == "limma2ways") {
    
    # create design
    designMN <- stats::model.matrix(~ factor1.fc + factor2.fc)
    rownames(designMN) <- rownames(data.mn)
    
    # apply limma test
    limmaFitLs <- limma::lmFit(t(data.mn), designMN)
    limmaBayesLs <- limma::eBayes(limmaFitLs)
    
    anova2waysPvalMN <- limmaBayesLs[["p.value"]][, 2:3]
    
  } else if (test.c == "limma2waysInter") {
    
    # create the design
    designMN <- stats::model.matrix(~ 0 + factor1.fc:factor2.fc)
    colnames(designMN) <- make.names(gsub("factor1.fc", "",
                                          gsub("factor2.fc", "",
                                               colnames(designMN))))
    rownames(designMN) <- rownames(data.mn)
    
    # apply linear model
    limmaFit <- limma::lmFit(t(data.mn), designMN)
    
    # create contrast matrix
    factInterVc <- colnames(designMN)
    fact1AllC <- paste(paste(factInterVc[c(2, 4)], collapse = "+"),
                       paste(factInterVc[c(1, 3)], collapse = "-"),
                       sep = "-")
    fact2AllC <- paste(paste(factInterVc[3:4], collapse = "+"),
                       paste(factInterVc[1:2], collapse = "-"),
                       sep = "-")
    factInterC <- paste(paste(factInterVc[c(1, 4)], collapse = "+"),
                        paste(factInterVc[c(2, 3)], collapse = "-"),
                        sep = "-")
    
    contrC <- c(fact1AllC, fact2AllC, factInterC)
    
    contrastMN <- limma::makeContrasts(contrasts = contrC,
                                       levels = factInterVc)
    
    # apply the contrast in linear model
    contrastFitLs <- limma::contrasts.fit(limmaFit, contrastMN)
    limmaBayesLs <- limma::eBayes(contrastFitLs)
    
    # recover signifcant variables
    anova2waysPvalMN <- limmaBayesLs[["p.value"]]
    
  }
  
  return(list(diffMN = anovas2waysDiffMN,
              pvalMN = anova2waysPvalMN))
  
}

.anovas2waysNames <- function(test.c,
                              factor_names.vc,
                              factor1.fc,
                              factor2.fc,
                              adjust.c,
                              prefix.c) {
  
  anova2waysNamesC <- c(.twoSampCorNames(test.c = test.c,
                                         factorNameC = factor_names.vc[1],
                                         factorFc = factor1.fc,
                                         adjust.c = adjust.c,
                                         prefix.c = prefix.c),
                        .twoSampCorNames(test.c = test.c,
                                         factorNameC = factor_names.vc[2],
                                         factorFc = factor2.fc,
                                         adjust.c = adjust.c,
                                         prefix.c = prefix.c))
  
  if (grepl("Inter$", test.c)) {
    
    anova2waysNamesC <- c(anova2waysNamesC,
                          paste0(paste0(test.c, "_",
                                        factor_names.vc[1], ":", factor_names.vc[2],
                                        "_"),
                                 c(adjust.c, "signif")))
    
  }
  
  anova2waysNamesC
  
}

.anovas2waysPlot <- function(data.mn,
                             test.c,
                             factor_names.vc,
                             adjust.c,
                             adjust_thresh.n,
                             metric.mn,
                             title.c,
                             figure.c) {
  
  metricSigniMN <- metric.mn[, grep("_signif$", colnames(metric.mn))]
  metricSigniVi <- which(rowSums(metricSigniMN) > 0)
  
  if (sum(metricSigniVi)) {
    
    metricSigniLs <- apply(metricSigniMN, 2, function(signiVn) which(signiVn > 0))
    
    metricSigniNamesVc <- factor_names.vc
    if (grepl("Inter$", test.c))
      metricSigniNamesVc <- c(metricSigniNamesVc, paste(factor_names.vc, collapse = ":"))
    names(metricSigniLs) <- metricSigniNamesVc
    
    mainC <- test.c
    if (!is.na(title.c))
      mainC <- paste0(title.c, "\n", mainC)
    
    if (!grepl("interactive", figure.c)) {
      file.tiffC <- gsub(".pdf", ".tiff", figure.c, fixed = TRUE)
    } else
      file.tiffC <- "none"
    
    vennplot(input.ls = metricSigniLs,
             title.c = mainC,
             figure.c = file.tiffC)
    
  }
  
}

# Print significant features
.signifPrint <- function(metric.mn,
                         signiVl,
                         signif_maxprint.i,
                         adjust.c,
                         adjust_thresh.n) {
  
  metricSigniMN <- metric.mn[signiVl, , drop = FALSE]
  
  if (is.na(signif_maxprint.i))
    signif_maxprint.i <- nrow(metricSigniMN)
  
  message("\n", nrow(metricSigniMN), " variable",
          ifelse(nrow(metricSigniMN) > 1, "s", ""),
          " (", round(nrow(metricSigniMN) / nrow(metric.mn) * 100), "%) ",
          ifelse(nrow(metricSigniMN) > 1, "were", "was"),
          " found significant at the ", adjust_thresh.n,
          " level (after the '", adjust.c, "' correction).\n",
          ifelse(nrow(metricSigniMN) == 1,
                 "It is",
                 ifelse(nrow(metricSigniMN) <= signif_maxprint.i,
                        "They are",
                        paste0("The first ", signif_maxprint.i, " are"))),
          " displayed below",
          ifelse(nrow(metricSigniMN) > 1,
                 " (sorted by increasing corrected p-values)",
                 ""),
          ":\n")
  print(metricSigniMN[1:min(signif_maxprint.i, nrow(metricSigniMN)), , drop = FALSE])
  
}

