#### annotating ####

#' MS annotation
#'
#' Annotation with chemical and biological databases
#' 
#' @param x An S4 object of class \code{ExpressionSet} or \code{MultiDataSet}
#' @param database.c character: database to be used for annotation
#' @param param.ls list: parameters for database query; see the example below for
#' the name and related database of each of them
#' @param report.c Character: File name with '.txt' extension for the printed
#' results (call to sink()'); if 'interactive' (default), messages will be
#' printed on the screen; if 'none', no verbose will be generated
#' @return \code{ExpressionSet} or \code{MultiDataSet} including the appended fData
#' data frame(s)
#' @rdname annotating
#' @export
#' @examples
#' sacurine.eset <- phenomis::reading(system.file("extdata/sacurine", package = "phenomis"))
#' # see the (default) parameters (e.g. for ChEBI query)
#' phenomis::annotating_parameters("chebi")
#' # mz annotation with ChEBI
#' sacurine.eset <- phenomis::annotating(sacurine.eset, database.c = "chebi",
#' param.ls = list(query.type = "mz", query.col = "mass_to_charge", ms.mode = "neg",
#' prefix = "chebiMZ."))
#' # mz annotation with local database
#' msdbDF <- read.table(system.file("extdata/local_ms_db.tsv", package = "phenomis"),
#' header = TRUE, sep = "\t", stringsAsFactors = FALSE)
#' sacurine.eset <- phenomis::annotating(sacurine.eset, database.c = "local.ms",
#' param.ls = list(query.type = "mz", query.col = "mass_to_charge", ms.mode = "neg",
#' mz.tol = 5, mz.tol.unit = "ppm", local.ms.db = msdbDF, prefix = "localMS."))
#' Biobase::fData(sacurine.eset)[!is.na(Biobase::fData(sacurine.eset)[, "localMS.accession"]), ]
#' # annotation from ChEBI identifiers
#' sacurine.eset <- phenomis::annotating(sacurine.eset, database.c = "chebi",
#' param.ls = list(query.type = "chebi.id", query.col = "database_identifier",
#' prefix = "chebiID."))
setGeneric("annotating",
           function(x,
                    database.c = c("chebi", "local.ms")[1],
                    param.ls = list(query.type = c("mz", "chebi.id", "kegg.id")[1],
                                    query.col = "mz",
                                    ms.mode = "pos",
                                    mz.tol = 10,
                                    mz.tol.unit = "ppm",
                                    fields = c("chebi.id", "name", "formula", "molecular.mass", "monoisotopic.mass"),
                                    fieldsLimit = 1,
                                    max.results = 3,
                                    local.ms.db = data.frame(),
                                    organism = "hsa",
                                    prefix = paste0(database.c, "."),
                                    sep = "|"),
                    report.c = c("none", "interactive", "myfile.txt")[2])
           standardGeneric("annotating"))


#### clustering ####

#' clustering
#'
#' Hierarchical clustering of both samples and variables
#'
#' @param x An S4 object of class \code{ExpressionSet} or \code{MultiDataSet}
#' @param dissym.c character: [hclust]
#' @param correl.c character: correlation coefficient (in case
#' '1-cor' or '1-abs(cor)' are selected as dissymilarity)
#' @param agglo.c character: agglomeration method
#' @param clusters.vi tupple of integers: number of sample and variable clusters, respectively
#' @param cex.vn tupple of numerics [Plot parameter]; size of the sample and variable labels
#' @param palette.c character [Plot parameter]: color palette
#' @param scale_plot.l logical [Plot parameter]: scaling (mean-centering and unit
#' variance scaling) to enhance contrast (for plotting only)
#' @param title.c character [Plot parameter]: Graphic the subtitle
#' @param figure.c character: File name with '.pdf' extension for the figure;
#' if 'interactive' (default), figures will be displayed interactively; if 'none',
#' no figure will be generated
#' @param report.c character: File name with '.txt' extension for the printed
#' results (call to sink()'); if 'interactive' (default), messages will be
#' printed on the screen; if 'none', no verbose will be generated
#' @return \code{ExpressionSet} or \code{MultiDataSet} including columns indicating
#' the clusters in pData and fData if 'clusters.vi' has been specified
#' @rdname clustering
#' @export
#' @examples
#' sacurine.eset <- phenomis::reading(system.file("extdata/sacurine", package = "phenomis"))
#' sacurine.eset <- phenomis::correcting(sacurine.eset)
#' sacurine.eset <- sacurine.eset[, Biobase::pData(sacurine.eset)[, "sampleType"] != "pool"]
#' sacurine.eset <- phenomis::transforming(sacurine.eset)
#' sacurine.eset <- sacurine.eset[, Biobase::sampleNames(sacurine.eset) != "HU_neg_096_b2"]
#' sacurine.eset <- phenomis::clustering(sacurine.eset)
#' utils::head(Biobase::fData(sacurine.eset))
#' # MultiDataSet
#' prometis.mset <- phenomis::reading(system.file("extdata/prometis/", package="phenomis"))
#' prometis.mset <- phenomis::clustering(prometis.mset, clusters.vi = c(3, 3))
setGeneric("clustering",
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
                    report.c = c("none", "interactive", "myfile.txt")[2])
           standardGeneric("clustering"))


#### correcting ####

#' correcting
#'
#' Signal drift and batch effect correction
#'
#' @param x An S4 object of class \code{ExpressionSet} or \code{MultiDataSet}
#' @param reference.c Character: sample type to be used as reference for
#' the correction (as indicated in the 'colnameSampleType' column from the
#' pData(x); e.g. 'pool')
#' @param col_batch.c Character: name of the column from pData(x) containing
#' the batch information (encoded as characters)
#' @param col_injectionOrder.c Character: name of the column from pData(x)
#' containing the injection order information (encoded as numerics)
#' @param col_sampleType.c Character:  name of the column from pData(x)
#' containing the sample type information (encoded as characters)
#' @param span.n Numeric: smoothing parameter for the loess regression;
#' between 0 and 1; (default set to 1)
#' @param title.c Character (ExpressionSet): Graphic title: if NA [default] the
#' title slot from the experimentData will be used
#' @param figure.c Character: File name with '.pdf' extension for the figure;
#' if 'interactive' (default), figures will be displayed interactively; if 'none',
#' no figure will be generated
#' @param report.c Character: File name with '.txt' extension for the printed
#' results (call to sink()'); if 'interactive' (default), messages will be
#' printed on the screen; if 'none', no verbose will be generated
#' @return \code{ExpressionSet} or \code{MultiDataSet} including the corrected
#' intensities in the expr matrix (matrices)
#' @rdname correcting
#' @export
#' @examples
#' sacurine.eset <- phenomis::reading(system.file("extdata/sacurine", package = "phenomis"))
#' sacurine.eset <- phenomis::correcting(sacurine.eset)
#' # MultiDataSet (to be done)
setGeneric("correcting",
           function(x,
                    reference.c = c("pool", "sample")[1],
                    col_batch.c = "batch",
                    col_injectionOrder.c = "injectionOrder",
                    col_sampleType.c = "sampleType",
                    span.n = 1,
                    title.c = NA,
                    figure.c = c("none", "interactive", "myfile.pdf")[2],
                    report.c = c("none", "interactive", "myfile.txt")[2])
             standardGeneric("correcting"))

#### filtering ####

#' Filtering of the features (or samples) with a high proportion of NAs or a low variance
#'
#' Filtering of the features (or samples) with a high proportion of NAs or a low variance
#'
#' @param x An S4 object of class \code{ExpressionSet} or \code{MultiDataSet}
#' @param class.c character(1): name of the column of the sample metadata giving
#' the classification groups: the filtering will be applied on each class
#' (default: "" meaning that there are no specific classes to consider)
#' @param max_na_prop.n numeric(1): maximum proportion of NAs for a feature (or
#' sample) to be kept [default: 0.2] (in case 'class.c' is provided, the maximum
#' proportion of NAs for a feature must be achieved in at least one sample class)
#' @param min_variance.n numeric(1): minimum variance for a feature (or sample) to be kept
#' [default: .Machine$double.eps] (in case 'class.c' is provided, the minimum variance
#' for a feature must be achieved in all sample classes)
#' @param dims.vc Vector of one or two characters: dimension(s) to which the
#' filtering should be applied; either 'features', 'samples', or c('features', 'samples')
#' @param report.c character(1): File name with '.txt' extension for the printed
#' results (call to sink()'); if 'interactive' (default), messages will be
#' printed on the screen; if 'none', no verbose will be generated
#' @return \code{ExpressionSet} or \code{MultiDataSet} including the exprs matrix
#' (list of matrices) with transformed intensities
#' @rdname filtering
#' @export
#' @examples
#' sacurine.eset <- phenomis::reading(system.file("extdata/sacurine", package = "phenomis"))
#' exprs.mn <- Biobase::exprs(sacurine.eset)
#' ropls::view(exprs.mn)
#' phenomis::filtering(sacurine.eset)
#' exprs.mn[exprs.mn < 1e5] <- NA
#' ropls::view(exprs.mn)
#' Biobase::exprs(sacurine.eset) <- exprs.mn
#' phenomis::filtering(sacurine.eset)
#' phenomis::filtering(sacurine.eset, class.c = "gender")
#' phenomis::filtering(sacurine.eset, class.c = "sampleType")
#' # MultiDataSet
#' prometis.mset <- phenomis::reading(system.file("extdata/prometis", package="phenomis"))
#' phenomis::filtering(prometis.mset)
#' for (set.c in names(prometis.mset)) {
#' eset <- prometis.mset[[set.c]]
#' exprs.mn <- Biobase::exprs(eset)
#' exprs.mn[exprs.mn < quantile(c(exprs.mn), 0.2)] <- NA
#' Biobase::exprs(eset) <- exprs.mn
#' prometis.mset <- MultiDataSet::add_eset(prometis.mset, eset, dataset.type = set.c,
#'                                         GRanges = NA, overwrite = TRUE, warnings = FALSE)
#' }
#' phenomis::filtering(prometis.mset)
setGeneric("filtering",
           function(x,
                    class.c = "",
                    max_na_prop.n = 0.2,
                    min_variance.n = .Machine$double.eps,
                    dims.vc = c("features", "samples"),
                    report.c = c("none", "interactive", "myfile.txt")[2])
             standardGeneric("filtering"))


#### hypotesting ####

#' Univariate hypothesis testing
#'
#' Provides numerical metrics and graphical overview of an ExpressionSet or MultiDataSet instance
#'
#' @param x An S4 object of class \code{ExpressionSet} or \code{MultiDataSet}
#' @param test.c Character: One of the 9 available hypothesis tests can be selected
#' (either 'ttest', 'limma', 'wilcoxon', 'anova', 'kruskal', 'pearson', 'spearman',
#' 'limma2ways', 'limma2waysInter', 'anova2ways', 'anova2waysInter')
#' @param factor_names.vc (Vector of) character(s): Factor(s) of interest (up to two),
#' i.e. name(s) of a column from the pData(x)
#' @param factor_levels.ls List: for each factor of interest (up to two), the levels
#' of the factor can be specified (i.e. re-ordered) by including a character vector
#' with those levels in the list; by default (no specification), the two vectors
#' are set to "default".
#' @param adjust.c Character: Name of the method for correction of multiple testing
#' (the p.adjust function is used)
#' @param adjust_thresh.n Numeric: Threshold for (corrected) p-values
#' @param signif_maxprint.i Interger: Maximum number of significant feature to display
#' on the screen (by default, 'NA', all significant features are displayed)
#' @param title.c Character: Title of the graphics
#' @param display_signif.l Character: In case of two sample tests (or
#' correlation test), should individual boxplots (or scatterplots) of significant
#' features be shown?
#' @param prefix.c Character: prefix to be added to the supplementary columns from
#' the variableMetadata to prevent overwriting of pre-existing columns with identical
#' names [default: ""]
#' @param figure.c Character: File name with '.pdf' extension for the figure (for
#' venn diagrams, e.g. in the 'anova2ways' test, the extension will be internally
#' changed to '.tiff' for compatibility with the VennDiagram package); if
#' 'interactive' (default), figures will be displayed interactively; if 'none',
#' no figure will be generated
#' @param report.c Character: File name with '.txt' extension for the printed
#' results (call to sink()'); if 'interactive' (default), messages will be
#' printed on the screen; if 'none', no verbose will be generated
#' @return \code{ExpressionSet} or \code{MultiDataSet} including the difference in
#' means/medians or correlations and the adjusted p-values in fData
#' @rdname hypotesting
#' @export
#' @examples
#' sacurine.eset <- phenomis::reading(system.file("extdata/sacurine", package = "phenomis"))
#' sacurine.eset <- phenomis::correcting(sacurine.eset, figure.c = 'none')
#' sacurine.eset <- sacurine.eset[, Biobase::pData(sacurine.eset)[, "sampleType"] != "pool"]
#' sacurine.eset <- phenomis::transforming(sacurine.eset)
#' sacurine.eset <- sacurine.eset[, Biobase::sampleNames(sacurine.eset) != "HU_neg_096_b2"]
#' # Student's T test
#' sacurine.eset <- hypotesting(sacurine.eset, "ttest", "gender")
#' # Pearson correlation test
#' sacurine.eset <- hypotesting(sacurine.eset, "pearson", "age")
#' # ANOVA
#' Biobase::pData(sacurine.eset)[, "ageGroup"] <- vapply(Biobase::pData(sacurine.eset)[, "age"],
#'                                                  function(x) {
#'                                                    if (x < 35) {
#'                                                      return("thirty")
#'                                                    } else if (x < 50) {
#'                                                      return("fourty")
#'                                                    } else {
#'                                                      return("fifty")}},
#'                                                  FUN.VALUE = character(1))
#' sacurine.eset <- phenomis::hypotesting(sacurine.eset, "anova", "ageGroup")
#' prometis.mset <- phenomis::reading(system.file("extdata/prometis", package="phenomis"))
#' prometis.mset <- phenomis::hypotesting(prometis.mset, "limma", "gene")
setGeneric("hypotesting",
           function(x,
                    test.c = c("ttest", "limma", "wilcoxon",
                               "anova", "kruskal",
                               "pearson", "spearman",
                               "limma2ways", "limma2waysInter",
                               "anova2ways", "anova2waysInter")[2],
                    factor_names.vc,
                    factor_levels.ls = list(factor1Vc = "default", factor2Vc = "default"),
                    adjust.c = c("holm", "hochberg", "hommel", "bonferroni", "BH",
                                 "BY", "fdr", "none")[5],
                    adjust_thresh.n = 0.05,
                    signif_maxprint.i = NA,
                    title.c = NA,
                    display_signif.l = FALSE,
                    prefix.c = "",
                    figure.c = c("none", "interactive", "interactive_plotly", "myfile.pdf")[2],
                    report.c = c("none", "interactive", "myfile.txt")[2])
           standardGeneric("hypotesting"))


#### inspecting ####

#' Inspecting
#'
#' Provides numerical metrics and graphical overview of an ExpressionSet or MultiDataSet instance
#' Please note that all variables with a proportion of missing values > 'max_na_prop.n'
#' or a variance of 0 will be filtered out at the beginning of the method and
#' therefore in the output ExpressionSet(s)
#'
#' @param x An S4 object of class \code{ExpressionSet} or \code{MultiDataSet}
#' @param pool_as_pool1.l Logical: should pool be included (as pool1) in the correlation
#' with the dilution factor?
#' @param pool_cv.n Numeric: threshold for the coefficient of variation of the pools
#' @param span.n Numeric: span parameter used in the loess trend estimation
#' @param sample_intensity.c Character: function to be used to display the global
#' sample intensity
#' @param title.c Character: MultiDataSet: title of the barplot showing the number
#' of samples and variables in each dataset; ExpressionSet: title of the multipanel
#' graphic displaying the metrics (if NA -default- the title slot from the experimentData
#' will be used)
#' @param plot_dims.l (MultiDataSet) Logical: should an overview of the number of samples and
#' variables in all datasets be barplotted? 
#' @param figure.c Character: File name with '.pdf' extension for the figure;
#' if 'interactive' (default), figures will be displayed interactively; if 'none',
#' no figure will be generated
#' @param report.c Character: File name with '.txt' extension for the printed
#' results (call to sink()'); if 'interactive' (default), messages will be
#' printed on the screen; if 'none', no verbose will be generated
#' @return \code{ExpressionSet} or \code{MultiDataSet} including the computed
#' in pData and fData sample and variable metrics
#' @rdname inspecting
#' @examples
#' sacurine.eset <- reading(system.file("extdata/sacurine", package = "phenomis"))
#' sacurine.eset <- inspecting(sacurine.eset)
#' sacurine.eset <- correcting(sacurine.eset)
#' sacurine.eset <- inspecting(sacurine.eset)
#' sacurine.eset <- transforming(sacurine.eset)
#' sacurine.eset <- inspecting(sacurine.eset)
#' # MultiDataSet
#' prometis.mset <- reading(system.file("extdata/prometis", package = "phenomis"))
#'\dontrun{
#' prometis.mset <- inspecting(prometis.mset)
#'}
setGeneric("inspecting",
           function(x,
                    pool_as_pool1.l = FALSE,
                    pool_cv.n = 0.3,
                    span.n = 1,
                    sample_intensity.c = c("median", "mean")[1],
                    title.c = NA,
                    plot_dims.l = TRUE,
                    figure.c = c("none", "interactive", "myfile.pdf")[2],
                    report.c = c("none", "interactive", "myfile.txt")[2])
             standardGeneric("inspecting"))


#### transforming ####

#' Transformation of the dataMatrix
#'
#' Transformation of the dataMatrix intensities
#'
#' @param x An S4 object of class \code{ExpressionSet} or \code{MultiDataSet}
#' @param method.c Character: Factor of interest (name of a column from the
#' pData(x))
#' @param report.c Character: File name with '.txt' extension for the printed
#' results (call to sink()'); if 'interactive' (default), messages will be
#' printed on the screen; if 'none', no verbose will be generated
#' @return \code{ExpressionSet} or \code{MultiDataSet} including the exprs matrix
#' (list of matrices) with transformed intensities
#' @rdname transforming
#' @export
#' @examples
#' sacurine.eset <- phenomis::reading(system.file("extdata/sacurine", package = "phenomis"))
#' sacurine.eset <- phenomis::correcting(sacurine.eset)
#' sacurine.eset <- sacurine.eset[, Biobase::pData(sacurine.eset)[, "sampleType"] != "pool"]
#' sacurine.eset <- phenomis::transforming(sacurine.eset)
#' # MultiDataSet
#' prometis.mset <- phenomis::reading(system.file("extdata/prometis", package="phenomis"))
#' prometis.mset <- phenomis::transforming(prometis.mset)
setGeneric("transforming",
           function(x,
                    method.c = c("log2", "log10", "sqrt")[1],
                    report.c = c("none", "interactive", "myfile.txt")[2])
             standardGeneric("transforming"))


#### writing ####

#' Exporting an ExpressionSet (or MultiDataSet) instance into (subfolders with)
#' the 3 tabulated files 'dataMatrix.tsv', sampleMetadata.tsv', 'variableMetadata.tsv'
#'
#' Note that the \code{dataMatrix} is transposed before
#' export (e.g., the samples are written column wise in the 'dataMatrix.tsv'
#' exported file).
#'
#' @param x An S4 object of class \code{ExpressionSet} or \code{MultiDataSet}
#' @param dir.c Character: directory where each dataset should be written
#' @param files.ls List: alternatively to the dir.c argument, the full names of the files can be provided as a list
#' @param overwrite.l Logical: should existing files be overwritten?
#' @param report.c Character: File name with '.txt' extension for the printed
#' results (call to sink()'); if 'interactive' (default), messages will be
#' printed on the screen; if 'none', no verbose will be generated
#' @return No object returned.
#' @rdname writing
#' @export
#' @examples
#' metabo.eset <- phenomis::reading(system.file("extdata/prometis/metabolomics", package="phenomis"))
#'\dontrun{
#' writing(metabo.eset, dir.c = file.path(getwd(), "metabolomics"))
#'}
#'# MultiDataSet
#' prometis.mset <- reading(system.file("extdata/prometis",package="phenomis"))
#'\dontrun{
#' writing(prometis.mset, dir.c = file.path(getwd(), "prometis"))
#' # alternatively
#' writing(prometis.mset,
#'          dir.c = NA,
#'          files.ls = list(metabolomics = list(dataMatrix.tsvC = file.path(getwd(), "met_dataMatrix.tsv"),
#'                                       sampleMetadata.tsvC = file.path(getwd(), "met_sampleMetadata.tsv"),
#'                                       variableMetadata.tsvC = file.path(getwd(), "met_variableMetadata.tsv")),
#'                         proteomics = list(dataMatrix.tsvC = file.path(getwd(), "pro_dataMatrix.tsv"),
#'                                       sampleMetadata.tsvC = file.path(getwd(), "pro_sampleMetadata.tsv"),
#'                                       variableMetadata.tsvC = file.path(getwd(), "pro_variableMetadata.tsv"))))
#'}
setGeneric("writing",
           function(x,
                    dir.c,
                    files.ls = NULL,
                    overwrite.l = FALSE,
                    report.c = c("none", "interactive", "myfile.txt")[2])
             standardGeneric("writing"))
