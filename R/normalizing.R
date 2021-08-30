#### normalizing (MultiDataSet) ####

#' @rdname normalizing
#' @export
setMethod("normalizing", signature(x = "MultiDataSet"),
          function(x,
                   method.c = "pqn",
                   report.c = c("none", "interactive", "myfile.txt")[2]) {
            
            if (!(report.c %in% c("none", "interactive")))
              sink(report.c, append = TRUE)
            
            inf_txt.c <- report.c
            if (inf_txt.c != "none")
              inf_txt.c <- "interactive"
            
            for (setC in names(x)) {
              
              if (report.c != "none")
                message("normalizing the '", setC, "' dataset...")
              
              ese <- phenomis::normalizing(x = x[[setC]],
                                           method.c = method.c,
                                           report.c = inf_txt.c)
              
              x <- MultiDataSet::add_eset(x,
                                          ese,
                                          dataset.type = setC,
                                          GRanges = NA,
                                          overwrite = TRUE,
                                          warnings = FALSE)
              
            }
            
            if (!(report.c %in% c("none", "interactive")))
              sink()
            
            return(invisible(x))
            
          })


#### normalizing (ExpressionSet) ####

#' @rdname normalizing
#' @export
setMethod("normalizing", signature(x = "ExpressionSet"),
          function(x,
                   method.c = "pqn",
                   report.c = c("none", "interactive", "myfile.txt")[2]) {
            
            if (!(report.c %in% c("none", "interactive")))
              sink(report.c, append = TRUE)
            
            if (method.c == "pqn") {
              if ("sampleType" %in% Biobase::varLabels(x) &&
                  "pool" %in% Biobase::pData(x)[, "sampleType"]) {
                reference.vl <- Biobase::pData(x)[, "sampleType"] == "pool"
                message("PQN normalization on pools")
              } else {
                reference.vl <- rep(TRUE, Biobase::dims(x)["Samples", "exprs"])
                message("PQN normalization on all samples")
              }
            } else
              reference.vl <- FALSE
            
            norm.mn <- .normalizing(t(Biobase::exprs(x)),
                                    method.c = method.c,
                                    reference.vl = reference.vl)
            
            Biobase::exprs(x) <- t(norm.mn)
            
            methods::validObject(x)
            
            if (!(report.c %in% c("none", "interactive")))
              sink()
            
            return(invisible(x))
            
          })


.normalizing <- function(data.mn,
                         method.c = "pqn",
                         reference.vl) {
  
  stopifnot(is.matrix(data.mn))
  stopifnot(mode(data.mn) == "numeric")
  stopifnot(is.logical(reference.vl))
  stopifnot(identical(nrow(data.mn), length(reference.vl)))
  stopifnot(method.c %in% c("pqn"))
  
  # Dieterle, F., Ross, A., Schlotterbeck, G., & Senn, H. (2006).
  # Probabilistic Quotient Normalization as Robust Method to Account for Dilution
  # of Complex Biological Mixtures. Application in1H NMR Metabonomics.
  # Analytical Chemistry, 78(13), 4281–4290. https://doi.org/10.1021/ac051632c
  
  # 1. Perform an integral normalization (typically a constant integral of 100 is used).
  
  # data_norm.mn <- sweep(data.mn,
  #                       MARGIN = 1,
  #                       STATS = rowSums(data.mn, na.rm = TRUE) / 100,
  #                       FUN = "/")
 
  data_norm.mn <- data.mn
  
  # 2. Choose/calculate the reference spectrum (the best approach is the 
  # calculation of the median spectrum of control samples).
  
  ref_spec.vn <- apply(data_norm.mn[reference.vl, ], 2,
                       function(feat.vn) median(feat.vn, na.rm = TRUE))
  
  # 3. Calculate the quotients of all variables of interest of the test
  # spectrum with those of the reference spectrum.
  
  data_quot.mn <- sweep(data_norm.mn,
                        MARGIN = 2,
                        STATS = ref_spec.vn,
                        FUN = "/")
  
  # 4. Calculate the median of these quotients.
  
  data_quot_med.vn <- apply(data_quot.mn, 1,
                            function(quot.vn) median(quot.vn, na.rm = TRUE))
  
  # 5. Divide all variables of the test spectrum by this median.
  
  data_pqn.mn <- sweep(data_norm.mn,
                      MARGIN = 1,
                      STATS = data_quot_med.vn,
                      FUN = "/")
  
  return(data_pqn.mn)

}