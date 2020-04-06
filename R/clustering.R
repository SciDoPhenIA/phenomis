#### clustering (MultiDataSet) ####

#' @rdname clustering
#' @export
setMethod("clustering", signature(x = "MultiDataSet"),
          function(x,
                   dissym.c = c("euclidean",
                                "maximum",
                                "manhattan",
                                "canberra",
                                "binary",
                                "minkowski",
                                "1-cor",
                                "1-abs(cor)")[7],
                   correl.c = c("pearson",
                                "kendall",
                                "spearman")[1],
                   agglo.c = c("ward.D",
                               "ward.D2",
                               "single",
                               "complete",
                               "average",
                               "mcquitty",
                               "median",
                               "centroid")[2],
                   clusters.vi = c(2, 2),
                   cex.vn = c(1, 1),
                   palette.c = c("blueOrangeRed",
                                 "redBlackGreen")[1],
                   scale_plot.l = TRUE,
                   title.c = NA,
                   figure.c = c("none", "interactive", "myfile.pdf")[2],
                   report.c = c("none", "interactive", "myfile.txt")[2]) {
            
            if (!(report.c %in% c("none", "interactive")))
              sink(report.c, append = TRUE)
            
            infTxtC <- report.c
            if (infTxtC != "none")
              infTxtC <- "interactive"
            
            if (!(figure.c %in% c("none", "interactive")))
              grDevices::pdf(figure.c)
            
            figPdfC <- figure.c
            if (figPdfC != "none")
              figPdfC <- "interactive"
            
            for (setC in names(x)) {
              
              if (report.c != "none")
                message("Performing hierarchical clustering on the '", setC, "' dataset...")
              
              ese <- x[[setC]]
              
              ese <- clustering(x = ese,
                                dissym.c = dissym.c,
                                correl.c = correl.c,
                                agglo.c = agglo.c,
                                clusters.vi = clusters.vi,
                                cex.vn = cex.vn,
                                palette.c = palette.c,
                                scale_plot.l = scale_plot.l,
                                title.c = paste0("[", setC, "]"),
                                figure.c = figPdfC,
                                report.c = infTxtC)
              
              x <- MultiDataSet::add_eset(x,
                                          ese,
                                          dataset.type = setC,
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


#### clustering (ExpressionSet) ####

#' @rdname clustering
#' @export
setMethod("clustering", signature(x = "ExpressionSet"),
          function(x,
                   dissym.c = c("euclidean",
                                "maximum",
                                "manhattan",
                                "canberra",
                                "binary",
                                "minkowski",
                                "1-cor",
                                "1-abs(cor)")[7],
                   correl.c = c("pearson",
                                "kendall",
                                "spearman")[1],
                   agglo.c = c("ward.D",
                               "ward.D2",
                               "single",
                               "complete",
                               "average",
                               "mcquitty",
                               "median",
                               "centroid")[2],
                   clusters.vi = c(2, 2),
                   cex.vn = c(1, 1),
                   palette.c = c("blueOrangeRed",
                                 "redBlackGreen")[1],
                   scale_plot.l = TRUE,
                   title.c = NA,
                   figure.c = c("none", "interactive", "myfile.pdf")[2],
                   report.c = c("none", "interactive", "myfile.txt")[2]) {
            
            if (!(report.c %in% c("none", "interactive")))
              sink(report.c, append = TRUE)
            
            x <- .clustering(eset = x,
                             dissym.c = dissym.c,
                             correl.c = correl.c,
                             agglo.c = agglo.c,
                             clusters.vi = clusters.vi,
                             cex.vn = cex.vn,
                             palette.c = palette.c,
                             scale_plot.l = scale_plot.l,
                             title.c = title.c,
                             figure.c = figure.c)
            
            methods::validObject(x)
            
            if (!(report.c %in% c("none", "interactive")))
              sink()
            
            return(invisible(x))
            
          })

.clustering <- function(eset,
                        dissym.c,    ## dissimilarity
                        correl.c, ## correlation method
                        agglo.c, ## agglomeration method
                        clusters.vi,
                        cex.vn,
                        palette.c,    ## color scale
                        scale_plot.l,
                        title.c,
                        figure.c) {
  
  proMN <- t(Biobase::exprs(eset))
  
  ncaN <- 14 ## Sample and variable name truncature for display
  
  if (dissym.c %in% c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski")) {
    
    obsHcl <- stats::hclust(stats::dist(proMN, method = dissym.c),
                            method = agglo.c)
    
    feaHcl <- stats::hclust(stats::dist(t(proMN), method = dissym.c),
                            method = agglo.c)
    
  } else if (dissym.c == "1-cor") {
    
    obsHcl <- stats::hclust(stats::as.dist(1 - stats::cor(t(proMN),
                                                          method = correl.c,
                                                          use = "pairwise.complete.obs")),
                            method = agglo.c)
    
    feaHcl <- stats::hclust(stats::as.dist(1 - stats::cor(proMN,
                                                          method = correl.c,
                                                          use = "pairwise.complete.obs")),
                            method = agglo.c)
    
  } else if (dissym.c == "1-abs(cor)") {
    
    obsHcl <- stats::hclust(stats::as.dist(1 - abs(stats::cor(t(proMN),
                                                              method = correl.c,
                                                              use = "pairwise.complete.obs"))),
                            method = agglo.c)
    
    feaHcl <- stats::hclust(stats::as.dist(1 - abs(stats::cor(proMN,
                                                              method = correl.c,
                                                              use = "pairwise.complete.obs"))),
                            method = agglo.c)
    
  }
  
  heaMN <- proMN <- proMN[obsHcl[["order"]], feaHcl[["order"]]]
  
  if (scale_plot.l)
    heaMN <- scale(heaMN)
  
  heaMN <- heaMN[, rev(1:ncol(heaMN)), drop = FALSE]
  
  switch(palette.c,
         blueOrangeRed = {
           imaPalVn <- grDevices::colorRampPalette(c("blue", "orange", "red"),
                                                   space = "rgb")(5)[1:5]
         },
         redBlackGreen = {
           imaPalVn <- grDevices::colorRampPalette(c("red", "black", "green"),
                                                   space = "rgb")(5)[1:5]
         })
  
  ## figure
  
  if (figure.c != "none") {
    
    if (figure.c != "interactive")
      grDevices::pdf(figure.c)
    
    graphics::layout(matrix(1:4, nrow = 2),
                     widths = c(1, 4), heights = c(1, 4))
    
    ## Color scale
    
    scaN <- length(imaPalVn)
    
    graphics::par(mar = c(0.6, 0.6, 0.6, 4.1))
    
    ylimVn <- c(0, scaN)
    ybottomVn <- 0:(scaN - 1)
    ytopVn <- 1:scaN
    
    graphics::plot(x = 0,
                   y = 0,
                   bty = "n",
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
    
    graphics::rect(xleft = 0.8,
                   ybottom = ybottomVn,
                   xright = 1,
                   ytop = ytopVn,
                   col = imaPalVn,
                   border = NA)
    
    prtVn <- pretty(range(heaMN, na.rm = TRUE))
    graphics::axis(at = scaN / diff(range(prtVn)) * (prtVn - min(prtVn)),
                   font = 2,
                   font.axis = 2,
                   labels = prtVn,
                   las = 1,
                   lwd = 2,
                   lwd.ticks = 2,
                   side = 4,
                   xpd = TRUE)
    
    graphics::arrows(graphics::par("usr")[2],
                     graphics::par("usr")[4],
                     graphics::par("usr")[2],
                     graphics::par("usr")[3],
                     code = 0,
                     lwd = 2,
                     xpd = TRUE)
    
    ## Feature dendrogram
    
    graphics::par(mar = c(4.1, 0.6, 0, 0.1),
                  lwd = 2)
    
    graphics::plot(rev(stats::as.dendrogram(feaHcl)), horiz = TRUE,
                   leaflab = "none",
                   main = "", xaxs = "i", yaxs = "i",
                   xaxt = "n", yaxt = "n", xlab = "", ylab = "")
    
    revFeaHcl <- list(merge = cbind(feaHcl[["merge"]][, 2], feaHcl[["merge"]][, 1]),
                      height = feaHcl[["height"]],
                      order = rev(feaHcl[["order"]]),
                      labels = feaHcl[["labels"]])
    
    if (clusters.vi[2] > 1) {
      cluFeaVn <- stats::cutree(revFeaHcl, k = clusters.vi[2])[revFeaHcl[["order"]]]
      cutFeaVn <- which(abs(diff(cluFeaVn)) > 0)
      cutFeaTxtVn <- c(cutFeaVn[1] / 2, cutFeaVn + diff(c(cutFeaVn, length(cluFeaVn))) / 2) + 0.5
      cutFeaLinVn <- cutFeaVn + 0.5
      graphics::text(graphics::par("usr")[1] + 0.2 * diff(graphics::par("usr")[1:2]),
                     cutFeaTxtVn,
                     labels = unique(cluFeaVn),
                     cex = 2,
                     font = 2,
                     las = 2)
    }
    
    ## Observation dendrogram
    
    graphics::par(mar = c(0.1, 0, 0.6, 4.1),
                  lwd = 2)
    
    graphics::plot(stats::as.dendrogram(obsHcl), leaflab = "none",
                   main = "", xaxs = "i", yaxs = "i",
                   yaxt = "n", xlab = "", ylab = "")
    
    if (clusters.vi[1] > 1) {
      cluObsVn <- stats::cutree(obsHcl, k = clusters.vi[1])[obsHcl[["order"]]]
      cutObsVn <- which(abs(diff(cluObsVn)) > 0)
      cutObsTxtVn <- c(cutObsVn[1] / 2, cutObsVn + diff(c(cutObsVn, length(cluObsVn))) / 2) + 0.5
      cutObsLinVn <- cutObsVn + 0.5
      graphics::text(cutObsTxtVn,
                     0.8 * graphics::par("usr")[4],
                     labels =  unique(cluObsVn),
                     cex = 2,
                     font = 2)
    }
    
    ## Heatmap
    
    graphics::par(mar = c(4.1, 0, 0, 4.1))
    
    graphics::image(x = 1:nrow(heaMN),
                    y = 1:ncol(heaMN),
                    z = round(heaMN),
                    col = imaPalVn,
                    font.axis = 2,
                    font.lab = 2,
                    xaxt = "n",
                    yaxt = "n",
                    xlab = "",
                    ylab = "")
    
    obsOrdVc <- obsHcl[["labels"]][obsHcl[["order"]]]
    obsOrdLenVn <- sapply(obsOrdVc, nchar)
    obsOrdVc <- substr(obsOrdVc, 1, ncaN)
    obsOrdVc <- paste0(obsOrdVc, ifelse(obsOrdLenVn > ncaN, ".", ""), " ")
    
    graphics::mtext(obsOrdVc,
                    at = 1:nrow(heaMN),
                    cex = cex.vn[1],
                    las = 2,
                    side = 1)
    
    feaOrdVc <- feaHcl[["labels"]][feaHcl[["order"]]]
    feaOrdLenVn <- sapply(feaOrdVc, nchar)
    feaOrdVc <- substr(feaOrdVc, 1, ncaN)
    feaOrdVc <- paste0(" ", feaOrdVc, ifelse(feaOrdLenVn > ncaN, ".", ""))
    
    graphics::mtext(feaOrdVc,
                    at = ncol(heaMN):1,
                    cex = cex.vn[2],
                    las = 2,
                    side = 4)
    
    if (clusters.vi[2] > 1)
      graphics::abline(h = cutFeaLinVn)
    if (clusters.vi[1] > 1)
      graphics::abline(v = cutObsLinVn)
    
    graphics::box()
    
    if (!is.na(title.c))
      graphics::title(paste0(title.c, " "), line = -1, adj = 1, outer = TRUE)
    
    if (figure.c != "interactive")
      grDevices::dev.off()
    
  }
  
  ## Returning
  
  if (clusters.vi[1] > 1) ## number of sample clusters
    Biobase::pData(eset)[, "hclust"] <- stats::cutree(obsHcl, k = clusters.vi[1])
  
  if (clusters.vi[2] > 1) ## number of variable clusters
    Biobase::fData(eset)[, "hclust"] <- stats::cutree(feaHcl, k = clusters.vi[2])
  
  return(invisible(eset))
  
}
