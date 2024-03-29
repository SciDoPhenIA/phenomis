% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/allGenerics.R, R/inspecting.R
\name{inspecting}
\alias{inspecting}
\alias{inspecting,MultiAssayExperiment-method}
\alias{inspecting,SummarizedExperiment-method}
\alias{inspecting,MultiDataSet-method}
\alias{inspecting,ExpressionSet-method}
\title{Inspecting}
\usage{
inspecting(
  x,
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
  report.c = c("none", "interactive", "myfile.txt")[2]
)

\S4method{inspecting}{MultiAssayExperiment}(
  x,
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
  report.c = c("none", "interactive", "myfile.txt")[2]
)

\S4method{inspecting}{SummarizedExperiment}(
  x,
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
  report.c = c("none", "interactive", "myfile.txt")[2]
)

\S4method{inspecting}{MultiDataSet}(
  x,
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
  report.c = c("none", "interactive", "myfile.txt")[2]
)

\S4method{inspecting}{ExpressionSet}(
  x,
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
  report.c = c("none", "interactive", "myfile.txt")[2]
)
}
\arguments{
\item{x}{An S4 object of class \code{SummarizedExperiment} or \code{MultiAssayExperiment}}

\item{pool_as_pool1.l}{Logical: should pool be included (as pool1) in the correlation
with the dilution factor?}

\item{pool_cv.n}{Numeric: threshold for the coefficient of variation of the pools}

\item{loess_span.n}{Numeric: span parameter used in the loess trend estimation}

\item{sample_intensity.c}{Character: function to be used to display the global
sample intensity; default: 'mean'}

\item{title.c}{Character: MultiAssayExperiment: title of the barplot showing the number
of samples and variables in each dataset; ExpressionSet: title of the multipanel
graphic displaying the metrics (if NA -default- the title slot from the experimentData
will be used)}

\item{plot_dims.l}{(MultiAssayExperiment) Logical: should an overview of the number of samples and
variables in all datasets be barplotted?}

\item{col_batch.c}{Character: name of the column from colData(x) containing
the batch information (encoded as characters)}

\item{col_injectionOrder.c}{Character: name of the column from colData(x)
containing the injection order information (encoded as numerics)}

\item{col_sampleType.c}{Character:  name of the column from colData(x)
containing the sample type information (encoded as characters)}

\item{figure.c}{Character: File name with '.pdf' extension for the figure;
if 'interactive' (default), figures will be displayed interactively; if 'none',
no figure will be generated}

\item{report.c}{Character: File name with '.txt' extension for the printed
results (call to sink()'); if 'interactive' (default), messages will be
printed on the screen; if 'none', no verbose will be generated}
}
\value{
\code{SummarizedExperiment} or \code{MultiAssayExperiment} including the computed
in rowData and colData sample and variable metrics
}
\description{
Provides numerical metrics and graphical overview of an ExpressionSet or MultiDataSet instance
Please note that all variables with a proportion of missing values > 'max_na_prop.n'
or a variance of 0 will be filtered out at the beginning of the method and
therefore in the output SummarizedExperiments(s)
}
\examples{
sacurine.se <- reading(system.file("extdata/W4M00001_Sacurine-statistics", package = "phenomis"))
sacurine.se <- inspecting(sacurine.se)
sacurine.se <- correcting(sacurine.se)
sacurine.se <- inspecting(sacurine.se)
sacurine.se <- transforming(sacurine.se)
sacurine.se <- inspecting(sacurine.se)
# MultiAssayExperiment
prometis.mae <- reading(system.file("extdata/prometis", package = "phenomis"))
\dontrun{
prometis.mae <- inspecting(prometis.mae)
}
}
