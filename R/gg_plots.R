#' Barplot with ggplot2
#'
#' Barplot with ggplot2
#' 
#' @param data.mn Matrix of numerics: values to be barplotted
#' @param row_levels.vc Vector of characters: levels of rownames (default: NA: alphabetical order will be used)
#' @param col_levels.vc Vector of characters: levels of colnames (default: NA: alphabetical order will be used)
#' @param title.c Character: plot title
#' @param xlab.c Character: x label
#' @param ylab.c Character: y label
#' @param palette.vc Character: either the name of an RColorBrewer palette
#' (default: 'Set1'; 'Paired' can be useful for parallel plotting) or a vector
#' manually defining the colors
#' @param cex_axis.i Integer: size of axis text (default: 18)
#' @param cex_bar.i Integer: size of bar value text (default: 10)
#' @param cex_title.i Integer: size of title text (default: 28)
#' @param bar_just.n Numeric: adjustment of bar value text (default : 0.9)
#' @param figure.c Character: either 'interactive' for interactive display,
#' 'my_barplot.pdf' for figure saving (only the extension matters), or 'none' to prevent plotting
#' @return invisible ggplot2 object
#' @export
#' @examples
#' prometis.mset <- phenomis::reading(system.file("extdata/prometis", package = "phenomis"))
#' dims.mn <- sapply(names(prometis.mset), function(set.c) { Biobase::dims(prometis.mset[[set.c]])})
#' dims.mn <- t(dims.mn)
#' colnames(dims.mn) <- c("Features", "Samples")
#' gg_barplot(dims.mn, title.c = "ProMetIS data")
#' prometis_dims.ls <- Biobase::dims(prometis.mset)
gg_barplot <- function(data.mn,
                       row_levels.vc = NA,
                       col_levels.vc = NA,
                       title.c = "",
                       xlab.c = "",
                       ylab.c = "",
                       palette.vc = "Set1",
                       cex_axis.i = 18,
                       cex_bar.i = 10,
                       cex_title.i = 28,
                       bar_just.n = 0.9,
                       figure.c = c("interactive",
                                    "my_barplot.pdf",
                                    "none")[1]) {
  
  data.df <- as.data.frame(data.mn)
  data.df[, "rownames"] <- rownames(data.df)
  
  longer.tib <- tidyr::pivot_longer(data.df, colnames(data.df)[colnames(data.df) != "rownames"])
  
  if (length(row_levels.vc) > 1)
    longer.tib$rownames <- factor(longer.tib$rownames, levels = row_levels.vc)
  
  if (length(col_levels.vc) > 1)
    longer.tib$name <- factor(longer.tib$name, levels = col_levels.vc)
  
  p <- .gg_barplot(longer.tib,
                   x.c = "rownames",
                   y.c = "value",
                   color.c = "name",
                   title.c = title.c,
                   xlab.c = xlab.c,
                   ylab.c = ylab.c,
                   palette.vc = palette.vc,
                   geom_text.ls = list(axis.i = cex_axis.i, bar.i = cex_bar.i,
                                       bar_just.n = bar_just.n, title.i = cex_title.i),
                   figure.c = figure.c)
  
  return(invisible(p))
  
}


.gg_barplot <- function(data.tb,
                       x.c = "",
                       y.c = "",
                       color.c = "",
                       title.c = "",
                       xlab.c = "",
                       ylab.c = "",
                       palette.vc = "Set1",
                       legend_position.c = c("none", "top")[2],
                       geom_text.ls = list(axis.i = 18, bar.i = 10,
                                           bar_just.n = 1.5, title.i = 28),
                       flip.l = TRUE,
                       position_dodge.l = TRUE,
                       figure.c = c("interactive",
                                    "my_barplot.pdf",
                                    "none")[1]) {
  # http://www.sthda.com/french/wiki/ggplot2-barplots-guide-de-demarrage-rapide-logiciel-r-et-visualisation-de-donnees
  
  if (!tibble::is_tibble(data.tb))
    data.tb <- tibble::as_tibble(data.tb)
 
  geom_text_default.vn <- c(axis.i = 18, bar.i = 10,
                            bar_just.n = 1.5, title.i = 28)
  for (geom_text.c in names(geom_text_default.vn)) {
    if (!(geom_text.c %in% names(geom_text.ls)))
      geom_text.ls[[geom_text.c]] <- geom_text_default.vn[geom_text.c]
  }
  
  filename_ext.c <-  utils::tail(unlist(strsplit(basename(figure.c),
                                                 ".", fixed = TRUE)), 1)
  
  if (!is.factor(data.tb[[x.c]]))
    data.tb[[x.c]] <- factor(data.tb[[x.c]], levels = unique(data.tb[[x.c]]))
  
  if (color.c != "" && !is.factor(data.tb[[color.c]]))
    data.tb[[color.c]] <- factor(data.tb[[color.c]],
                                 levels = unique(data.tb[[color.c]]))
  
  if (flip.l) {

    data.tb[[x.c]] <- factor(data.tb[[x.c]],
                             levels = rev(levels(data.tb[[x.c]])))
    
    if (color.c != "")
      data.tb[[color.c]] <- factor(data.tb[[color.c]],
                                   levels = rev(levels(data.tb[[color.c]])))

  }
  
  aes.c <- paste0("ggplot2::ggplot(data.tb, ggplot2::aes(x = ",
                  x.c, ", y = ", y.c, ", fill = ", color.c, "))")
  
  p <- eval(parse(text = aes.c))
  
  # color palette
  if (length(palette.vc) == 1 &&
      palette.vc %in% rownames(RColorBrewer::brewer.pal.info)) {
    p <- p + ggplot2::scale_fill_brewer(palette = palette.vc)
  } else
    p <- p + ggplot2::scale_fill_manual(values = palette.vc)
  
  # parallel bars in case of paired statistics
  if (position_dodge.l) {
    p <- p + ggplot2::geom_bar(stat = "identity",
                               position = ggplot2::position_dodge())
  } else
    p <- p + ggplot2::geom_bar(stat = "identity")
  
  # numbers super-imposed on the bars
  p <- p +
    eval(parse(text = paste0("ggplot2::geom_text(ggplot2::aes(label = ",
                             y.c,
                             "), ",
                             ifelse(flip.l, "h", "v"),
                             "just = ", geom_text.ls[['bar_just.n']], ", fontface = 'bold',
                               position = ggplot2::position_dodge(0.9),
                               size = ", geom_text.ls[['bar.i']], ")")))
  
  p <- p +
    ggplot2::labs(title = title.c, x = xlab.c, y = ylab.c) +
    ggplot2::theme(legend.position = legend_position.c,
                   legend.text = ggplot2::element_text(size = geom_text.ls[['axis.i']],
                                                       face = "bold"),
                   legend.title = ggplot2::element_blank(),
                   axis.text = ggplot2::element_text(size = geom_text.ls[['axis.i']],
                                                     face = "bold"),
                   plot.title = ggplot2::element_text(size = geom_text.ls[['title.i']],
                                                      face = "bold"))
  
  if (flip.l) {
    p <- p + ggplot2::coord_flip()
    if (color.c != "")
      p <- p + ggplot2::guides(fill = ggplot2::guide_legend(reverse = TRUE))
  }
  
  if (filename_ext.c != "none") {
    
    if (filename_ext.c == "pdf")
      grDevices::pdf(figure.c)
    
    print(p)
    
    if (filename_ext.c == "pdf")
      grDevices::dev.off()
    
  }
  
  return(invisible(p))
  
}

#' Boxplot with ggplot2
#' 
#' Boxplot with ggplot2
#' 
#' @param data.tb Data frame (or tibble) containing the information
#' @param x.c Character: name of the column with qualitative levels
#' @param y.c Character: name of the column with quantitative values
#' @param color.c Character: optional name of the column for color information
#' @param title.c Character: plot title
#' @param xlab.c Character: x label
#' @param ylab.c Character: y label
#' @param label.vc Character (vector): either the name of a character column
#' from the data or a character vector of the same length as the rown number of
#' the data, containing the feature labeling for outlier display
#' @param palette.vc Character: either the name of an RColorBrewer palette
#' (default: 'Set1'; 'Paired' can be useful for parallel plotting) or a vector
#' manually defining the colors
#' @param size.ls List of sizes for dots (default is 0.7), labels (default is 16),
#' ticks (14) and title (20)
#' @param figure.c Character: either 'interactive' for interactive display or
#' 'my_barplot.pdf' for figure saving (only the extension matters)
#' @return character vector of outlier labels (same dimension as the number of rows from data.tb)
#' @export
#' @examples
#' sacurine.eset <- phenomis::reading(system.file("extdata/W4M00001_Sacurine-statistics", package = "phenomis"))
#' sacurine_pda.df <- Biobase::pData(sacurine.eset)
#' sacurine_pda.df <- sacurine_pda.df[!grepl("QC", rownames(sacurine_pda.df)), ]
#' phenomis::gg_boxplot(sacurine_pda.df, y.c = "age")
#' phenomis::gg_boxplot(sacurine_pda.df, x.c = "gender", y.c = "bmi", color.c = "gender")
#' phenomis::gg_boxplot(sacurine_pda.df, x.c = "gender", y.c = "bmi", color.c = "gender", label.vc = rownames(sacurine_pda.df))
gg_boxplot <- function(data.tb,
                       x.c = "",
                       y.c = "",
                       color.c = "",
                       title.c = NA,
                       xlab.c = NA,
                       ylab.c = "",
                       label.vc = "",
                       palette.vc = "Set1",
                       size.ls = list(dot.n = 0.7,
                                      lab.i = 20,
                                      tick.i = 20,
                                      title.i = 20),
                       figure.c = c("interactive",
                                    "my_boxplot.pdf")[1]) {
  # http://www.sthda.com/french/wiki/ggplot2-box-plot-guide-de-demarrage-rapide-logiciel-r-et-visualisation-de-donnees
  
  stopifnot(is.numeric(data.tb[[y.c]]))
  
  size_default.vi <- c(class.i = 5, dot.n = 0.7, lab.i = 16, point.i = 3,
                       tick.i = 14, title.i = 20)
  
  for (size.c in names(size_default.vi)) {
    if (!(size.c %in% names(size.ls)))
      size.ls[[size.c]] <- size_default.vi[size.c]
  }
  
  filename_ext.c <-  utils::tail(unlist(strsplit(basename(figure.c), ".", fixed = TRUE)), 1)
  
  if (is.na(xlab.c))
    xlab.c <- ifelse(x.c == "", y.c, "")
  
  if (is.na(title.c))
    title.c <- ifelse(x.c == "", "", y.c)
  
  # labels (for outliers)
  
  if (length(label.vc) == 1) {
    if (label.vc != "") {
      stopifnot(label.vc %in% colnames(data.tb))
      label.vc <- as.character(data.tb[[label.vc]])
    } else
      label.vc <- as.character(data.tb[[y.c]])
  } else {
    stopifnot(length(label.vc) == nrow(data.tb))
  }
  
  # color palette
  if (color.c != "")
    aes.c <- paste0("ggplot2::ggplot(data.tb, ggplot2::aes(x = ",
                    ifelse(x.c == "", "''", x.c), ", y = ", y.c, ", color = ", color.c, "))")
  else
    aes.c <- paste0("ggplot2::ggplot(data.tb, ggplot2::aes(x = ",
                    ifelse(x.c == "", "''", x.c), ", y = ", y.c, "))")
  
  p <- eval(parse(text = aes.c))
  
  p <- p +
    ggplot2::geom_boxplot(lwd = 2)
  

  is_outlier <- function(values.vn, index.vc) {
    out.vl <- values.vn < quantile(values.vn, 0.25, na.rm = TRUE) - 1.5 * IQR(values.vn, na.rm = TRUE) |
      values.vn > quantile(values.vn, 0.75, na.rm = TRUE) + 1.5 * IQR(values.vn, na.rm = TRUE)
    out.vi <- as.numeric(index.vc[out.vl])
    out.vi[!is.na(out.vi)]
  }

  if (x.c == "") {
    outlier.vi <- is_outlier(data.tb[[y.c]], as.character(1:nrow(data.tb)))
  } else {
    outlier.vi <- integer()
    x.vc <- as.character(data.tb[[x.c]])
    x.vuc <- unique(x.vc)
    for (x.uc in x.vuc) {
      index.vi <- which(x.vc == x.uc)
      outlier.vi <- c(outlier.vi, is_outlier(data.tb[[y.c]][index.vi], as.character(index.vi)))
    }
  }

  outlier.vc <- character(nrow(data.tb))
  outlier.vc[outlier.vi] <- label.vc[outlier.vi]

  p <- p + eval(parse(text = paste0("ggrepel::geom_text_repel(ggplot2::aes(x = ",
                                    ifelse(x.c == "", "''", x.c), ", y = ", y.c, ", label = outlier.vc))")))

  if (length(palette.vc) == 1 &&
      palette.vc %in% rownames(RColorBrewer::brewer.pal.info))
    p <- p + ggplot2::scale_colour_brewer(palette = palette.vc)
  else {
    if (is.factor(data.tb[[color.c]]))
      names(palette.vc) <- levels(data.tb[[color.c]])
    p <- p + ggplot2::scale_color_manual(values = palette.vc)
  }
  
  p <- p +
    ggplot2::labs(title = title.c,
                  x = xlab.c, y = ylab.c) +
    ggplot2::geom_dotplot(binaxis = "y", stackdir = "center",
                          dotsize = size.ls[["dot.n"]]) +
    ggplot2::theme(plot.title = ggplot2::element_text(size = size.ls[["title.i"]], face = "bold"),
                   axis.title.x = ggplot2::element_text(size = size.ls[["lab.i"]], face = "bold"),
                   axis.title.y = ggplot2::element_text(size = size.ls[["lab.i"]], face = "bold"),
                   axis.text = ggplot2::element_text(size = size.ls[["tick.i"]]),
                   legend.title = ggplot2::element_blank(),
                   legend.position = "none")
  
  if (filename_ext.c == "pdf")
    grDevices::pdf(figure.c)
  
  print(p)
  
  if (filename_ext.c == "pdf")
    grDevices::dev.off()
  
  return(invisible(outlier.vc))
  
}


#' Pie with ggplot2
#'
#' Pie with ggplot2
#' 
#' @param data.tb Tibble (or data frame) containing the information
#' @param y.c Character: name of the column with the factor to be displayed; alternatively,
#' name of the column with the counts (in this case set the name of the column with
#' the names of the factor levels with the 'color.c' argument)
#' @param color.c Character: optional name of the column with the names of the factor levels
#' @param title.c Character: plot title
#' @param palette.vc Character: either the name of an RColorBrewer palette
#' (default: 'Set1'; 'Paired' can be useful for parallel plotting) or a vector
#' manually defining the colors
#' @param label.c Character: (relative) counts to be displayed on the pie; either
#' 'none' (default), 'value' or 'percent'
#' @param geom_text.ls List of sizes for lab.i (default 7), legend_title.i (16),
#' legend_text.i (14), and title.i (16)
#' @param figure.c Character: either 'interactive' for interactive display,
#' 'my_pie.pdf' for figure saving (only the extension matters), or 'none' to prevent plotting
#' @return invisible ggplot2 object
#' @export
#' @examples
#' sacurine.eset <- phenomis::reading(system.file("extdata/W4M00001_Sacurine-statistics", package = "phenomis"))
#' sacurine_pda.df <- Biobase::pData(sacurine.eset)
#' sacurine_pda.df <- sacurine_pda.df[!grepl("QC", rownames(sacurine_pda.df)), ]
#' gg_pie(sacurine_pda.df, y.c = "gender", label.c = "value")
gg_pie <- function(data.tb,
                   y.c = "",
                   color.c = "",
                   title.c = "",
                   palette.vc = "Set1",
                   label.c = c("none", "value", "percent")[1],
                   geom_text.ls = list(lab.i = 7, legend_title.i = 16, legend_text.i = 14, title.i = 16),
                   figure.c = c("interactive",
                                "my_pie.pdf",
                                "none")[1]) {
  
  if (!tibble::is_tibble(data.tb))
    data.tb <- tibble::as_tibble(data.tb)
  
  geom_text_default.vn <- c(lab.i = 7, legend_title.i = 16, legend_text.i = 14, title.i = 16)
  for (geom_text.c in names(geom_text_default.vn)) {
    if (!(geom_text.c %in% names(geom_text.ls)))
      geom_text.ls[[geom_text.c]] <- geom_text_default.vn[geom_text.c]
  }
  
  filename_ext.c <-  utils::tail(unlist(strsplit(basename(figure.c),
                                                 ".", fixed = TRUE)), 1)
  
  if (color.c == "") {
    if (y.c == "")
      stop("When color.c is '', y.c must be specified.", call. = FALSE)
    y.fc <- data.tb[[y.c]]
    if (!is.factor(y.fc)) {
      if (is.character(y.fc)) {
        y.fc <- factor(y.fc)
      } else
        stop("When color.c is '', the y.c column of data.frame must be a factor or character vector.",
             call. = FALSE)
    }
    data.tb <- eval(parse(text = paste0("dplyr::summarize(dplyr::group_by(data.tb, ",
                                        y.c, "), n = dplyr::n())")))
    color.c <- y.c
    y.c <- "n"
  }
  
  aes.c <- paste0("ggplot2::ggplot(data.tb, ggplot2::aes(x = '', y = ", y.c, ", fill = ", color.c, "))")
  p <- eval(parse(text = aes.c)) + ggplot2::geom_bar(width = 1, stat = "identity")
  
  # color palette
  if (length(palette.vc) == 1 &&
      palette.vc %in% rownames(RColorBrewer::brewer.pal.info)) {
    
    palette_col.i <- RColorBrewer::brewer.pal.info[palette.vc, "maxcolors"]
    
    if (y.c == "") {
      data_col.i <- nlevels(data.tb[[color.c]])
    } else
      data_col.i <- nrow(data.tb)
    
    if (data_col.i <= palette_col.i) {
      p <- p + ggplot2::scale_fill_brewer(palette = palette.vc)
    } else {
      fill_values.vc <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(palette_col.i, palette.vc))(data_col.i)
      p <- p + ggplot2::scale_fill_manual(values = fill_values.vc)
    }
    
  } else
    p <- p + ggplot2::scale_fill_manual(values = palette.vc)
  
  p <- p + ggplot2::coord_polar("y", start = 0, direction = -1) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(size = geom_text.ls[["title.i"]], face = "bold"),
      legend.title = ggplot2::element_text(size = geom_text.ls[["legend_title.i"]]), 
      legend.text = ggplot2::element_text(size = geom_text.ls[["legend_text.i"]])
    ) +
    ggplot2::labs(title = title.c)
  
  if (label.c != "none") {
    if (label.c == "value") {
      p <- p + eval(parse(text = paste0("ggplot2::geom_text(ggplot2::aes(label = ",
                                        y.c, "), position = ggplot2::position_stack(0.5), size = ",
                                        geom_text.ls[["lab.i"]], ")")))
      # p <- p + eval(parse(text = paste0("ggplot2::geom_text(ggplot2::aes(y = ",
      #                                   y.c, "/3 + c(0, cumsum(",
      #                                   y.c, ")[-length(", y.c, ")]), label = ",
      #                                   y.c, "), position = ggplot2::position_stack(0.5), size = ",
      #                                   geom_text.ls[["lab.i"]], ")")))
    } else if (label.c == "percent") {
      p <- p + eval(parse(text = paste0("ggplot2::geom_text(ggplot2::aes(label = scales::percent(",
                                        y.c, "/100)), position = ggplot2::position_stack(0.5)size = ",
                                        geom_text.ls[["lab.i"]], ")")))
    } else
      stop("'label.c' must be either 'none', 'value', or 'percent'", call. = FALSE)
  } 
  
  if (filename_ext.c != "none") {
    
    if (filename_ext.c == "pdf")
      grDevices::pdf(figure.c)
    
    print(p)
    
    if (filename_ext.c == "pdf")
      grDevices::dev.off()
    
  }
  
  return(invisible(p))
  
}

#' Volcano plot with ggplot2
#'
#' Volcano plot with ggplot2
#' 
#' @param fold_change.vn Numeric vector: fold changes
#' @param adjusted_pvalue.vn Numeric vector: (adjusted) p-values
#' @param adjust_method.c Character: method for multiple testing correction
#' @param adjust_thresh.n Numeric: significance threshold
#' @param label.vc Character (vector): either the name of a character column
#' from the data or a character vector of the same length as the row number of
#' the data, containing the feature labeling
#' @param title.c Character: plot title
#' @param xlab.c Character: x label (default: "Fold Change")
#' @param signif_palette.vc Character vector: color palette (default 'green4' for
#' significant features and 'gray' otherwise
#' @param signif_shape.vi Integer vector: shapes for significant (respectively,
#' non significant) features; default is 16 (respectively, 1)
#' @param class_name.vc Character vector: names of the two compared class labels
#' @param class_color.vc Character vector: colors of the two compared class labels
#' @param size.ls List of sizes for classes (default: 5), xy labels (default: 16),
#' points (default: 3), ticks (default: 14) and title (default: 20)
#' @param figure.c Character: either 'interactive' (respectively,
#' 'interactive_plotly') for interactive display with ggplot2 (respectively,
#' with plotly::ggplotly [default]), or 'my_volcanoplot.pdf' (respectively
#' 'my_volcanoplot.html') for figure saving (only the extension matters) with
#' ggplot2 (respectively, with plotly::ggplotly)
#' @return invisible ggplot2 object
#' @export
#' @examples
#' sacurine.eset <- phenomis::reading(system.file("extdata/W4M00001_Sacurine-statistics", package = "phenomis"))
#' sacurine.eset <- phenomis::correcting(sacurine.eset, figure.c = "none")
#' sacurine.eset <- sacurine.eset[, Biobase::pData(sacurine.eset)[, "sampleType"] != "pool"]
#' sacurine.eset <- phenomis::transforming(sacurine.eset)
#' sacurine.eset <- phenomis::hypotesting(sacurine.eset, test.c = "wilcoxon",
#'                                       factor_names.vc = "gender",
#'                                       figure.c = "none", report.c = "none")
#' fold.vn <- Biobase::fData(sacurine.eset)[, "wilcoxon_gender_Female.Male_diff"]
#' fdr.vn <- Biobase::fData(sacurine.eset)[, "wilcoxon_gender_Female.Male_BH"]
#' feat.vc <- Biobase::featureNames(sacurine.eset)
#' phenomis::gg_volcanoplot(fold.vn,
#'                          fdr.vn,
#'                          label.vc = feat.vc,
#'                          adjust_method.c = "BH")
#' feat_signif.vc <-  sapply(seq_along(feat.vc),
#'                           function(feat.i)
#'                            ifelse(fdr.vn[feat.i] <= 0.05, feat.vc[feat.i], ""))
#' phenomis::gg_volcanoplot(fold.vn,
#'                          fdr.vn,
#'                          label.vc = feat_signif.vc,
#'                          adjust_method.c = "BH",
#'                          figure.c = "interactive")
gg_volcanoplot <- function(fold_change.vn,
                           adjusted_pvalue.vn,
                           adjust_method.c = "",
                           adjust_thresh.n = 0.05,
                           label.vc = "",
                           title.c = "",
                           xlab.c = "Fold Change",
                           signif_palette.vc = c(yes = RColorBrewer::brewer.pal(9, "Greens")[8],
                                                 no = RColorBrewer::brewer.pal(9, "Greys")[7]),
                           signif_shape.vi = c(yes = 16,
                                               no = 1),
                           class_name.vc = "",
                           class_color.vc = "",
                           size.ls = list(class.i = 5,
                                          lab.i = 16,
                                          point.i = 3,
                                          tick.i = 14,
                                          title.i = 20),
                           figure.c = c("interactive",
                                      "interactive_plotly",
                                      "my_volcanoplot.pdf",
                                      "my_volcanoplot.html")[2]) {
  
  size_default.vi <- c(class.i = 5, lab.i = 16, point.i = 3, tick.i = 14,
                       title.i = 20)
  
  for (size.c in names(size_default.vi)) {
    if (!(size.c %in% names(size.ls)))
      size.ls[[size.c]] <- size_default.vi[size.c]
  }
  
  filename_ext.c <-  utils::tail(unlist(strsplit(basename(figure.c), ".", fixed = TRUE)), 1)
  
  stopifnot(identical(length(fold_change.vn), length(adjusted_pvalue.vn)))
  
  volcano.df <- data.frame(fold_change = fold_change.vn,
                           log_pval = -log10(adjusted_pvalue.vn))
  
  if (length(label.vc) == 1) {
    if (label.vc != "") {
      stopifnot(label.vc %in% colnames(volcano.df))
      label.vc <- as.character(volcano.df[[label.vc]])
    } else
      label.vc <- rep("", nrow(volcano.df))
  } else {
    stopifnot(length(label.vc) == nrow(volcano.df))
  }
  
  stopifnot(length(signif_palette.vc) == 2)
  stopifnot(length(signif_shape.vi) == 2)
  
  volcano.df[, "shape"] <- volcano.df[, "color"] <- ifelse(adjusted_pvalue.vn <= adjust_thresh.n,
                                                           "yes", "no")
  
  p <- ggplot2::ggplot(volcano.df,
                       ggplot2::aes(x = fold_change,
                                    y = log_pval,
                                    color = color,
                                    shape = shape,
                                    text = label.vc))
  
  if (figure.c == "interactive" || filename_ext.c == "pdf")
    p <- p + ggrepel::geom_text_repel(ggplot2::aes(x = fold_change,
                                                   y = log_pval,
                                                   label = label.vc))
  
  p <- p + ggplot2::scale_color_manual(values = signif_palette.vc)
  
  p <- p + ggplot2::scale_shape_manual(values = signif_shape.vi)
  
  p <- p + ggplot2::geom_hline(ggplot2::aes(yintercept = -log10(adjust_thresh.n))) +
    ggplot2::geom_vline(ggplot2::aes(xintercept = 0)) +
    ggplot2::geom_point(size = size.ls[["point.i"]])
  
  # title and axis labels
  p <- p + ggplot2::labs(title = title.c,
                         x = xlab.c,
                         y = paste0("-log10(", ifelse(adjust_method.c != "",
                                                      adjust_method.c, "p-value"), ")")) +
    ggplot2::theme(plot.title = ggplot2::element_text(size = size.ls[["title.i"]], face = "bold"),
                   axis.title.x = ggplot2::element_text(size = size.ls[["lab.i"]], face = "bold"),
                   axis.title.y = ggplot2::element_text(size = size.ls[["lab.i"]], face = "bold"),
                   axis.text = ggplot2::element_text(size = size.ls[["tick.i"]]),
                   legend.title = ggplot2::element_blank(),
                   legend.position = "none")
  
  
  if (length(class_name.vc) == 2) {
    
    if (length(class_color.vc) == 1 || any(class_color.vc == ""))
      class_color.vc <- rep("black", 2)
    
    p <- p + ggplot2::annotate(geom = "text",
                               x = 0.95 * min(volcano.df[, "fold_change"], na.rm = TRUE),
                               y = 0,
                               label = class_name.vc[1],
                               color = class_color.vc[1],
                               size = size.ls[["class.i"]]) +
      ggplot2::annotate(geom = "text",
                        x = 0.95 * max(volcano.df[, "fold_change"], na.rm = TRUE),
                        y = 0,
                        label = class_name.vc[2],
                        color = class_color.vc[2],
                        size = size.ls[["class.i"]])
    
  }
  
  if (figure.c == "interactive_plotly" || filename_ext.c == "html") {
    
    p <- plotly::ggplotly(p, tooltip = "text")
    
    p <- plotly::layout(p, hoverlabel = list(font = list(size = 20)))
    
    if (filename_ext.c == "html") {
      
      htmlwidgets::saveWidget(plotly::as_widget(p), figure.c)
      
      return(invisible(p))
      
      
    } else {
      
      print(p)
      
      return(invisible(p))
      
    }
    
  } else {
    
    if (filename_ext.c == "pdf")
      grDevices::pdf(figure.c)
    
    print(p)
    
    if (filename_ext.c == "pdf")
      grDevices::dev.off()
    
    return(invisible(p))
    
  }
  
}

#' Venn diagram with VennDiagram
#'
#' Venn diagram with VennDiagram
#' 
#' @param input.ls Named list of vectors to be compared
#' @param palette.vc Character vector: Color palette
#' @param title.c Character: Plot title
#' @param sub.c Character: Plot subtitle
#' @param cat_pos.vi Integer vector giving the position (in degrees) of each
#' category name along the circle, with 0 at 12 o'clock; if NA, (-50, 50),
#' (-40, 40, 180), (-15, 15, 0, 0), and (0, 287.5, 215, 145, 70) values are used
#' @param label_col.c Character: Label color
#' @param lwd.i Integer: Width of the circle's circumference
#' @param inverted.l Logical: Should the Venn diagram be flipped along its
#' vertical axis (pairwise venn only)
#' @param figure.c Character: Filename for image output (with either .tiff, .png,
#' or .svg extensions); if 'none' (default) the grid object is displayed interactively
#' @return invisible grid object
#' @export
#' @examples
#' sacurine.eset <- phenomis::reading(system.file("extdata/W4M00001_Sacurine-statistics", package = "phenomis"))
#' sacurine.eset <- phenomis::correcting(sacurine.eset, figure.c = 'none')
#' sacurine.eset <- sacurine.eset[, Biobase::pData(sacurine.eset)[, "sampleType"] != "pool"]
#' sacurine.eset <- phenomis::transforming(sacurine.eset)
#' sacurine.eset <- sacurine.eset[, Biobase::sampleNames(sacurine.eset) != "HU_neg_096_b2"]
#' # Student's T test
#' sacurine.eset <- phenomis::hypotesting(sacurine.eset, "ttest", "gender")
#' # Wilcoxon T test
#' sacurine.eset <- phenomis::hypotesting(sacurine.eset, "wilcoxon", "gender")
#' signif.ls <- list(ttest = which(Biobase::fData(sacurine.eset)[, "ttest_gender_Female.Male_signif"] > 0),
#' wilcoxon =  which(Biobase::fData(sacurine.eset)[, "wilcoxon_gender_Female.Male_signif"] > 0))
#' vennplot(signif.ls, label_col.c = "black", title.c = "Signif. features\nwith Student or Wilcoxon tests")
vennplot <- function(input.ls,
                     palette.vc = RColorBrewer::brewer.pal(9, "Set1")[1:5],
                     title.c = NA,
                     sub.c = "",
                     cat_pos.vi = NA,
                     label_col.c = "black",
                     lwd.i = 2,
                     inverted.l = FALSE,
                     figure.c = "none") {
  
  if (is.na(title.c))
    title.c <- gsub("MN", "",
                    gsub("Ls", "",
                         deparse(substitute(input.ls))))
  
  if (class(input.ls) != "list")
    stop("'input.ls' must be a list for Venn plot", call. = FALSE)
  
  if (length(input.ls) > 5)
    stop("'input.ls' list must be of maximum length 5 for Venn plot", call. = FALSE)
  
  cat.i <- length(input.ls)
  
  if (length(palette.vc) < cat.i)
    stop("'palette.vc' must contain at least 'length(input.ls)' colors", call. = FALSE)
  
  if (figure.c == "none") {
    filename.c <- NULL
  } else {
    filename.c <- figure.c
    ext.c <- tail(unlist(strsplit(basename(filename.c), ".", fixed = TRUE)), 1)
    if (!(ext.c %in% c("tiff", "png", "svg")))
      stop("Filename extension must be either 'tiff', 'png', or 'svg'",
           call. = FALSE)
  }
  
  futile.logger::flog.threshold(futile.logger::ERROR,
                                name = "VennDiagramLogger")
  
  if (any(is.na(cat_pos.vi)))
    if (cat.i == 2) {
      cat_pos.vi <- c(-40, 40)
      # cat_pos.vi <- c(-50, 50)
    } else if (cat.i == 3) {
      cat_pos.vi <- c(-40, 40, 180)
    } else if (cat.i == 4) {
      cat_pos.vi <- c(-15, 15, 0, 0)
    }
  
  if (inverted.l) {
    if (cat.i != 2)
      stop("'inverted.l' option is only available for pairwise venn",
           call. = FALSE)
    cat_pos.vi <- -cat_pos.vi
  }
  
  if (cat.i <= 3) {
    ven <- VennDiagram::venn.diagram(x = input.ls,
                                     
                                     alpha = 0.8,
                                     
                                     col = palette.vc[1:cat.i],
                                     
                                     cat.cex = 1.7, ## 2.5
                                     cat.col = palette.vc[1:cat.i],
                                     
                                     ## additional argument for cat.i <= 3
                                     cat.dist = c(rep(ifelse(cat.i == 2,
                                                             0.03, 0.05), 2),
                                                  ifelse(cat.i == 3, 0.04, 0.05),
                                                  rep(0.04, 2))[1:cat.i],
                                     
                                     cat.pos = cat_pos.vi,
                                     
                                     cex = 1.7, ## 3
                                     
                                     fill = palette.vc[1:cat.i],
                                     
                                     label.col = label_col.c,
                                     lwd = lwd.i, ## 4
                                     
                                     main = title.c,
                                     main.cex = 1.7,
                                     main.pos = c(0.5, 1.05),
                                     
                                     margin = 0.01,
                                     
                                     filename = filename.c,
                                     cat.fontfamily = "sans",
                                     cat.fontface = "bold",
                                     fontfamily = "sans",
                                     fontface = "bold",
                                     main.fontfamily = "sans",
                                     main.fontface = "bold",
                                     sub = sub.c,
                                     sub.cex = 1,
                                     sub.fontfamily = "sans",
                                     sub.pos = c(0.5, 1.0),
                                     inverted = inverted.l)
  } else {
    ven <- VennDiagram::venn.diagram(x = input.ls,
                                     
                                     alpha = 0.8,
                                     
                                     col = palette.vc[1:cat.i],
                                     
                                     cat.cex = 1.7, ## 2.5
                                     cat.col = palette.vc[1:cat.i],
                                     
                                     cex = 1.7, ## 3
                                     
                                     fill = palette.vc[1:cat.i],
                                     
                                     label.col = label_col.c,
                                     lwd = lwd.i, ## 4
                                     
                                     main = title.c,
                                     main.cex = 1.7,
                                     main.pos = c(0.5, 1.05),
                                     
                                     margin = 0.01,
                                     
                                     filename = filename.c,
                                     cat.fontfamily = "sans",
                                     cat.fontface = "bold",
                                     fontfamily = "sans",
                                     fontface = "bold",
                                     main.fontfamily = "sans",
                                     main.fontface = "bold",
                                     sub = sub.c,
                                     sub.cex = 1,
                                     sub.fontfamily = "sans",
                                     sub.pos = c(0.5, 1.0))
  }
  
  if (is.null(filename.c)) {
    grid::grid.newpage()
    grid::grid.draw(ven)
  }
  
  return(invisible(ven))
  
}
