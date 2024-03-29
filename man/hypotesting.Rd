% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/allGenerics.R, R/hypotesting.R
\name{hypotesting}
\alias{hypotesting}
\alias{hypotesting,MultiAssayExperiment-method}
\alias{hypotesting,SummarizedExperiment-method}
\alias{hypotesting,MultiDataSet-method}
\alias{hypotesting,ExpressionSet-method}
\title{Univariate hypothesis testing}
\usage{
hypotesting(
  x,
  test.c = c("ttest", "limma", "wilcoxon", "anova", "kruskal", "pearson", "spearman",
    "limma2ways", "limma2waysInter", "anova2ways", "anova2waysInter")[2],
  factor_names.vc,
  factor_levels.ls = list(factor1Vc = "default", factor2Vc = "default"),
  adjust.c = c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr",
    "none")[5],
  adjust_thresh.n = 0.05,
  signif_maxprint.i = NA,
  title.c = NA,
  display_signif.l = FALSE,
  prefix.c = "",
  figure.c = c("none", "interactive", "interactive_plotly", "myfile.pdf")[2],
  report.c = c("none", "interactive", "myfile.txt")[2]
)

\S4method{hypotesting}{MultiAssayExperiment}(
  x,
  test.c = c("ttest", "limma", "wilcoxon", "anova", "kruskal", "pearson", "spearman",
    "limma2ways", "limma2waysInter", "anova2ways", "anova2waysInter")[2],
  factor_names.vc,
  factor_levels.ls = list(factor1.vc = "default", factor2.vc = "default"),
  adjust.c = c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr",
    "none")[5],
  adjust_thresh.n = 0.05,
  signif_maxprint.i = NA,
  title.c = NA,
  display_signif.l = FALSE,
  prefix.c = "",
  figure.c = c("none", "interactive", "interactive_plotly", "myfile.pdf")[2],
  report.c = c("none", "interactive", "myfile.txt")[2]
)

\S4method{hypotesting}{SummarizedExperiment}(
  x,
  test.c = c("ttest", "limma", "wilcoxon", "anova", "kruskal", "pearson", "spearman",
    "limma2ways", "limma2waysInter", "anova2ways", "anova2waysInter")[2],
  factor_names.vc,
  factor_levels.ls = list(factor1.vc = "default", factor2.vc = "default"),
  adjust.c = c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr",
    "none")[5],
  adjust_thresh.n = 0.05,
  signif_maxprint.i = NA,
  title.c = NA,
  display_signif.l = FALSE,
  prefix.c = "",
  figure.c = c("none", "interactive", "interactive_plotly", "myfile.pdf")[2],
  report.c = c("none", "interactive", "myfile.txt")[2]
)

\S4method{hypotesting}{MultiDataSet}(
  x,
  test.c = c("ttest", "limma", "wilcoxon", "anova", "kruskal", "pearson", "spearman",
    "limma2ways", "limma2waysInter", "anova2ways", "anova2waysInter")[2],
  factor_names.vc,
  factor_levels.ls = list(factor1.vc = "default", factor2.vc = "default"),
  adjust.c = c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr",
    "none")[5],
  adjust_thresh.n = 0.05,
  signif_maxprint.i = NA,
  title.c = NA,
  display_signif.l = FALSE,
  prefix.c = "",
  figure.c = c("none", "interactive", "interactive_plotly", "myfile.pdf")[2],
  report.c = c("none", "interactive", "myfile.txt")[2]
)

\S4method{hypotesting}{ExpressionSet}(
  x,
  test.c = c("ttest", "limma", "wilcoxon", "anova", "kruskal", "pearson", "spearman",
    "limma2ways", "limma2waysInter", "anova2ways", "anova2waysInter")[2],
  factor_names.vc,
  factor_levels.ls = list(factor1.vc = "default", factor2.vc = "default"),
  adjust.c = c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr",
    "none")[5],
  adjust_thresh.n = 0.05,
  signif_maxprint.i = NA,
  title.c = NA,
  display_signif.l = FALSE,
  prefix.c = "",
  figure.c = c("none", "interactive", "interactive_plotly", "myfile.pdf")[2],
  report.c = c("none", "interactive", "myfile.txt")[2]
)
}
\arguments{
\item{x}{An S4 object of class \code{SummarizedExperiment} or \code{MultiAssayExperiment}}

\item{test.c}{Character: One of the 9 available hypothesis tests can be selected
(either 'ttest', 'limma', 'wilcoxon', 'anova', 'kruskal', 'pearson', 'spearman',
'limma2ways', 'limma2waysInter', 'anova2ways', 'anova2waysInter')}

\item{factor_names.vc}{(Vector of) character(s): Factor(s) of interest (up to two),
i.e. name(s) of a column from the pData(x)}

\item{factor_levels.ls}{List: for each factor of interest (up to two), the levels
of the factor can be specified (i.e. re-ordered) by including a character vector
with those levels in the list; by default (no specification), the two vectors
are set to "default".}

\item{adjust.c}{Character: Name of the method for correction of multiple testing
(the p.adjust function is used)}

\item{adjust_thresh.n}{Numeric: Threshold for (corrected) p-values}

\item{signif_maxprint.i}{Interger: Maximum number of significant feature to display
on the screen (by default, 'NA', all significant features are displayed)}

\item{title.c}{Character: Title of the graphics}

\item{display_signif.l}{Character: In case of two sample tests (or
correlation test), should individual boxplots (or scatterplots) of significant
features be shown?}

\item{prefix.c}{Character: prefix to be added to the supplementary columns from
the variableMetadata to prevent overwriting of pre-existing columns with identical
names [default: ""]}

\item{figure.c}{Character: File name with '.pdf' extension for the figure (for
venn diagrams, e.g. in the 'anova2ways' test, the extension will be internally
changed to '.tiff' for compatibility with the VennDiagram package); if
'interactive' (default), figures will be displayed interactively; if 'none',
no figure will be generated}

\item{report.c}{Character: File name with '.txt' extension for the printed
results (call to sink()'); if 'interactive' (default), messages will be
printed on the screen; if 'none', no verbose will be generated}
}
\value{
\code{SummarizedExperiment} or \code{MultiAssayExperiment} including the difference in
means/medians or correlations and the adjusted p-values in feature metadata
}
\description{
Univariate hypothesis testing
}
\examples{
sacurine.se <- reading(system.file("extdata/W4M00001_Sacurine-statistics", package = "phenomis"))
sacurine.se <- correcting(sacurine.se, figure.c = 'none')
sacurine.se <- sacurine.se[, colData(sacurine.se)[, "sampleType"] != "pool"]
sacurine.se <- transforming(sacurine.se)
sacurine.se <- sacurine.se[, colnames(sacurine.se) != "HU_neg_096_b2"]
# Student's T test
sacurine.se <- hypotesting(sacurine.se, "ttest", "gender")
# Pearson correlation test
sacurine.se <- hypotesting(sacurine.se, "pearson", "age")
# ANOVA
colData(sacurine.se)[, "ageGroup"] <- vapply(colData(sacurine.se)[, "age"],
                                                 function(x) {
                                                   if (x < 35) {
                                                     return("thirty")
                                                   } else if (x < 50) {
                                                     return("fourty")
                                                   } else {
                                                     return("fifty")}},
                                                 FUN.VALUE = character(1))
sacurine.se <- hypotesting(sacurine.se, "anova", "ageGroup")
prometis.mset <- reading(system.file("extdata/prometis", package="phenomis"))
prometis.mset <- hypotesting(prometis.mset, "limma", "gene")
}
