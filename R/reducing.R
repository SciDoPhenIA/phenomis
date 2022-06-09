#### reducing (MultiAssayExperiment) ####

#' @rdname reducing
#' @export
setMethod("reducing", signature(x = "MultiAssayExperiment"),
          function(x,
                   cor_method.c = "pearson",
                   cor_threshold.n = 0.9,
                   rt_tol.n = 6,
                   rt_colname.c = "rt",
                   mzdiff_tol.n = 0.005,
                   mz_colname.c = "mz",
                   return_adjacency.l = FALSE,
                   report.c = c("none", "interactive", "myfile.txt")[2]) {
            
            if (!(report.c %in% c("none", "interactive")))
              sink(report.c, append = TRUE)
            
            report_set.c <- report.c
            if (report_set.c != "none")
              report_set.c <- "interactive"
            
            if (return_adjacency.l)
              adjacency.ls <- list()
            
            for (set.c in names(x)) {
              
              if (report.c != "none")
                message("Reducing the '", set.c, "' dataset...")
              
              reduce_res <- reducing(x = x[[set.c]],
                                     cor_method.c = cor_method.c,
                                     cor_threshold.n = cor_threshold.n,
                                     rt_tol.n = rt_tol.n,
                                     rt_colname.c = rt_colname.c,
                                     mzdiff_tol.n = mzdiff_tol.n,
                                     mz_colname.c = mz_colname.c,
                                     return_adjacency.l = return_adjacency.l,
                                     report.c = report_set.c)
              
              if (!return_adjacency.l) {
                
                x[[set.c]] <- reduce_res
                
              } else {
                
                x[[set.c]] <- reduce_res[["se"]]
                
                adjacency.ls[[set.c]] <- reduce_res[["adjacency.mi"]]
                
              }
              
            }
            
            if (!(report.c %in% c("none", "interactive")))
              sink()
            
            methods::validObject(x)
            
            if (return_adjacency.l) {
              
              return(invisible(list(mae = x,
                                    adjacency.ls = adjacency.ls)))
              
            } else {
              
              return(invisible(x))
              
            }
            
          })

#### reducing (SummarizedExperiment) ####

#' @rdname reducing
#' @export
setMethod("reducing", signature(x = "SummarizedExperiment"),
          function(x,
                   cor_method.c = "pearson",
                   cor_threshold.n = 0.9,
                   rt_tol.n = 6,
                   rt_colname.c = "rt",
                   mzdiff_tol.n = 0.005,
                   mz_colname.c = "mz",
                   return_adjacency.l = FALSE,
                   report.c = c("none", "interactive", "myfile.txt")[2]) {
            
            if (!(report.c %in% c("none", "interactive")))
              sink(report.c, append = TRUE)
            
            redund.ls <- .reducing(data.mn = t(assay(x)),
                                   feat.df = rowData(x),
                                   cor_method.c = cor_method.c,
                                   cor_threshold.n = cor_threshold.n,
                                   rt_tol.n = rt_tol.n,
                                   rt_colname.c = rt_colname.c,
                                   mzdiff_tol.n = mzdiff_tol.n,
                                   mz_colname.c = mz_colname.c)
            
            # re-ordering features
            x <- x[colnames(redund.ls[["data.mn"]]), ]
            assay(x) <- t(redund.ls[["data.mn"]])
            rowData(x) <- redund.ls[["feat.df"]]
            adjacency.mi <- redund.ls[["corrtmz.mi"]]
 
            if (!(report.c %in% c("none", "interactive")))
              sink()
            
            methods::validObject(x)
            
            if (return_adjacency.l) {
              
              return(invisible(list(se = x,
                                    adjacency.mi = adjacency.mi)))
              
            } else {
              
              return(invisible(x))
              
            }
            
          })


#### reducing (MultiDataSet) ####

#' @rdname reducing
#' @export
setMethod("reducing", signature(x = "MultiDataSet"),
          function(x,
                   cor_method.c = "pearson",
                   cor_threshold.n = 0.9,
                   rt_tol.n = 6,
                   rt_colname.c = "rt",
                   mzdiff_tol.n = 0.005,
                   mz_colname.c = "mz",
                   return_adjacency.l = FALSE,
                   report.c = c("none", "interactive", "myfile.txt")[2]) {
            
            if (!(report.c %in% c("none", "interactive")))
              sink(report.c, append = TRUE)
            
            report_set.c <- report.c
            if (report_set.c != "none")
              report_set.c <- "interactive"
            
            if (return_adjacency.l)
              adjacency.ls <- list()
            
            for (set.c in names(x)) {
              
              if (report.c != "none")
                message("Reducing the '", set.c, "' dataset...")
              
              reduce_res <- phenomis::reducing(x = x[[set.c]],
                                               cor_method.c = cor_method.c,
                                               cor_threshold.n = cor_threshold.n,
                                               rt_tol.n = rt_tol.n,
                                               rt_colname.c = rt_colname.c,
                                               mzdiff_tol.n = mzdiff_tol.n,
                                               mz_colname.c = mz_colname.c,
                                               return_adjacency.l = return_adjacency.l,
                                               report.c = report_set.c)
              
              if (!return_adjacency.l) {
                
                x <- MultiDataSet::add_eset(x,
                                            reduce_res,
                                            dataset.type = set.c,
                                            GRanges = NA,
                                            overwrite = TRUE,
                                            warnings = FALSE)
                
              } else {
                
                x <- MultiDataSet::add_eset(x,
                                            reduce_res[["eset"]],
                                            dataset.type = set.c,
                                            GRanges = NA,
                                            overwrite = TRUE,
                                            warnings = FALSE)
                
                adjacency.ls[[set.c]] <- reduce_res[["adjacency.mi"]]
                
              }
              
            }
            
            if (!(report.c %in% c("none", "interactive")))
              sink()
            
            methods::validObject(x)
            
            if (return_adjacency.l) {
              
              return(list(mset = x,
                          adjacency.ls = adjacency.ls))
              
            } else {
              
              return(invisible(x))
              
            }
            
          })


#### reducing (ExpressionSet) ####

#' @rdname reducing
#' @export
setMethod("reducing", signature(x = "ExpressionSet"),
          function(x,
                   cor_method.c = "pearson",
                   cor_threshold.n = 0.9,
                   rt_tol.n = 6,
                   rt_colname.c = "rt",
                   mzdiff_tol.n = 0.005,
                   mz_colname.c = "mz",
                   return_adjacency.l = FALSE,
                   report.c = c("none", "interactive", "myfile.txt")[2]) {
            
            if (!(report.c %in% c("none", "interactive")))
              sink(report.c, append = TRUE)
            
            redund.ls <- .reducing(data.mn = t(Biobase::exprs(x)),
                                   feat.df = Biobase::fData(x),
                                   cor_method.c = cor_method.c,
                                   cor_threshold.n = cor_threshold.n,
                                   rt_tol.n = rt_tol.n,
                                   rt_colname.c = rt_colname.c,
                                   mzdiff_tol.n = mzdiff_tol.n,
                                   mz_colname.c = mz_colname.c)
            
            # re-ordering features
            x <- x[colnames(redund.ls[["data.mn"]]), ]
            Biobase::exprs(x) <- t(redund.ls[["data.mn"]])
            Biobase::fData(x) <- redund.ls[["feat.df"]]
            adjacency.mi <- redund.ls[["corrtmz.mi"]]
            
            if (!(report.c %in% c("none", "interactive")))
              sink()
            
            methods::validObject(x)
            
            if (return_adjacency.l) {
              
              return(list(eset = x,
                          adjacency.mi = adjacency.mi))
              
            } else {
              
              return(invisible(x))
              
            }
            
          })

.reducing <- function(data.mn, ## data (matrix of numerics; samples x variables)
                      feat.df, ## feature metadata (dataframe; features x metadata)
                      cor_method.c = "pearson",
                      cor_threshold.n = 0.9,
                      rt_tol.n = 6,
                      rt_colname.c = "rt",
                      mzdiff_tol.n = 0.005,
                      mz_colname.c = "mz") {
  
  # table of referenced losses (fragments, adducts, isotopes)
  
  mzdiff_db.df <- .mzdiff_db()
  
  mzdiff_db.df[, "delta_mass"] <- abs(mzdiff_db.df[, "delta_mass"])
  
  mzdiff_db.df <- mzdiff_db.df[mzdiff_db.df[, "add_and_frag_dup"] < 1, ]
  
  # mean intensity for each feature
  
  feat.df[, "redund_samp_mean"] <- colMeans(data.mn, na.rm = TRUE)
  
  # correlation between features
  
  if (any(is.na(data.mn))) {
    cor.mn <- stats::cor(data.mn,
                         method = cor_method.c,
                         use = "pairwise.complete.obs")
  } else {
    cor.mn <- stats::cor(data.mn,
                         method = cor_method.c) 
  }
  
  if (any(is.na(cor.mn)))
    warning("Correlation matrix contains missing values")
  
  cor.mi <- cor.mn
  diag(cor.mi) <- 0
  cor.mi[cor.mi < cor_threshold.n] <- 0
  cor.mi[cor.mi > 0] <- 1
  mode(cor.mi) <- "integer"
  cor_sel.vi <- rowSums(cor.mi, na.rm = TRUE) > 0
  cor.mi <- cor.mi[cor_sel.vi, cor_sel.vi]
  
  # RT difference
  
  if (!(rt_colname.c %in% colnames(feat.df)))
    stop("The 'rt_colname.c' value '", rt_colname.c,
         "' was not found in the columns from the variable metadata:\n",
         paste(colnames(feat.df), collapse = ", "))
  
  rt.vn <- feat.df[rownames(cor.mi), rt_colname.c]
  if (any(is.na(rt.vn)))
    stop("The '", rt_colname.c, "' column contains missing values")
  
  rt_tol_over_max.n <- -log10(rt_tol.n / max(rt.vn))
  
  if (rt_tol_over_max.n < 1 || rt_tol_over_max.n > 3)
    warning("The 'rt_tol.n' value (", rt_tol.n,
            ") corresponds to ", signif(10^(-rt_tol_over_max.n) * 100, 2),
            "% of the maximum value (", max(rt.vn),
            ") from the '", rt_colname.c,
            "' values; you might check your 'rt_tol.n' unit")
  
  rt_diff.mn <- abs(outer(rt.vn, rt.vn, "-"))
  dimnames(rt_diff.mn) <- dimnames(cor.mi)
  
  rt.mi <- rt_diff.mn <= rt_tol.n
  mode(rt.mi) <- "integer"
  
  corrt.mi <- cor.mi * rt.mi
  mode(corrt.mi) <- "integer"
  
  corrt_sel.vi <- rowSums(corrt.mi, na.rm = TRUE) > 0
  corrt.mi <- corrt.mi[corrt_sel.vi, corrt_sel.vi]
  
  # known m/z loss
  
  if (!(mz_colname.c %in% colnames(feat.df)))
    stop("The 'mz_colname.c' value '", mz_colname.c,
         "' was not found in the columns from the variable metadata:\n",
         paste(colnames(feat.df), collapse = ", "))
  
  ## annotation of all pairwise differences of corrt.mi features
  mz.vn <- feat.df[rownames(corrt.mi), mz_colname.c]
  if (any(is.na(mz.vn)))
    stop("The '", mz_colname.c, "' column contains missing values")
  
  mzdiff.mn <- -outer(mz.vn, mz.vn, "-")
  dimnames(mzdiff.mn) <- dimnames(corrt.mi)
  
  mzdiff_upper.mi <- which(upper.tri(mzdiff.mn) & corrt.mi > 0, arr.ind = TRUE)
  
  mzdiff.vn <- mzdiff.mn[mzdiff_upper.mi]
  
  mzmin.vn <- feat.df[rownames(mzdiff.mn), mz_colname.c]
  
  mzannot_upper.mc <- vapply(seq_along(mzdiff.vn),
                             function(mzdiff.i) {
                               mzdiff.n <- mzdiff.vn[mzdiff.i]
                               mzmin.n <- mzmin.vn[mzdiff.i]
                               paste(which(abs(c(mzdiff_db.df[, "delta_mass"], mzmin.n) - mzdiff.n) <= mzdiff_tol.n),
                                     collapse = "_")
                             }, character(1))
  
  corrtmz.mc <- matrix("", nrow = nrow(corrt.mi), ncol = ncol(corrt.mi),
                       dimnames = dimnames(corrt.mi))
  corrtmz.mc[mzdiff_upper.mi] <- mzannot_upper.mc
  corrtmz.mc <- paste0(t(corrtmz.mc), corrtmz.mc)
  dim(corrtmz.mc) <- dim(corrt.mi)
  dimnames(corrtmz.mc) <- dimnames(corrt.mi)
  
  ## adjacency matrix of pairs of features with correlation, RT match and loss annotation
  corrtmz.mi <- corrtmz.mc != ""
  mode(corrtmz.mi) <- "integer"
  
  corrtmz_sel.vi <- which(rowSums(corrtmz.mi, na.rm = TRUE) > 0)
  corrtmz.mi <- corrtmz.mi[corrtmz_sel.vi, corrtmz_sel.vi]
  corrtmz.mc <- corrtmz.mc[corrtmz_sel.vi, corrtmz_sel.vi]
  
  # connex components of the corr + RT adjacency matrix
  
  corrtmz.igraph <- igraph::graph_from_adjacency_matrix(corrtmz.mi,
                                                        mode = "undirected")
  
  components.ls <- igraph::components(corrtmz.igraph)
  
  group.vi <- rep(NA_integer_, nrow(feat.df))
  names(group.vi) <- rownames(feat.df)
  group.vi[names(components.ls[["membership"]])] <- components.ls[["membership"]]
  
  repres.vc <- character(nrow(feat.df))
  names(repres.vc) <- rownames(feat.df)
  
  isoadfra.vc <- relat.vc <- repres.vc
  
  for (group.i in seq_along(table(group.vi))) {
    
    # group.i <- 1
    # group.i <- 62
    
    ## features from the component
    group_feat.vc <- names(group.vi)[which(group.vi == group.i)]
    
    ion_sel.c <- group_feat.vc[which.max(feat.df[group_feat.vc, "redund_samp_mean"])]
    
    repres.vc[group_feat.vc] <- ion_sel.c
    
    ## annotation of pairwise m/z differences within a group
    
    for (group_feat.c in group_feat.vc) {
      
      group_feat_isoadfra.vc <- character()
      
      for (linked.c in setdiff(group_feat.vc, group_feat.c)) {
        
        sign.i <- sign(mzdiff.mn[linked.c, group_feat.c])
        
        annot.c <- corrtmz.mc[linked.c, group_feat.c]
        
        if (annot.c != "") {
          
          annot_split.vc <- unlist(strsplit(annot.c, split = "_"))
         
          if (length(annot_split.vc)) {
            
            annot.c <- paste(vapply(annot_split.vc,
                                    function(annot_split.c) {
                                      
                                      if (as.numeric(annot_split.c) > nrow(mzdiff_db.df)) {
                                        ## The 2M annotation is coded as nrow(mzdiff_db.df) + 1
                                        
                                        if (sign.i > 0) {
                                          
                                          code.c <- paste0("+(", floor(feat.df[linked.c, mz_colname.c] * 1e4) / 1e4, ")")
                                          
                                        } else {
                                          
                                          code.c <- paste0("-(", floor(feat.df[group_feat.c, mz_colname.c] * 1e4) / 1e4, ")")
                                          
                                        }
                                        
                                        return(code.c)
                                                         
                                        
                                      } else {
                                        
                                        code.c <- mzdiff_db.df[as.numeric(annot_split.c), "losses_or_gains"]
                                        
                                        if (sign.i > 0) {
                                          
                                          if (!(substr(code.c, 1, 1) %in% c("+", "-"))) {
                                            code.c <- paste0("+", code.c)
                                          } else if (substr(code.c, 1, 1) == "-") {
                                            code.c <- gsub("@", "-",
                                                           gsub("-", "+",
                                                                gsub("+", "@", code.c, fixed = TRUE), fixed = TRUE),
                                                           fixed = TRUE)
                                          }
                                          
                                          # if (!(substr(code.c, 1, 1) %in% c("+", "-")))
                                          #   code.c <- paste0("+", code.c)
                                          
                                        } else {
                                          
                                          if (!(substr(code.c, 1, 1) %in% c("+", "-"))) {
                                            if (code.c != "(H2O-CO)")
                                              code.c <- gsub("@", "-",
                                                             gsub("-", "+",
                                                                  gsub("+", "@", code.c,
                                                                       fixed = TRUE),
                                                                  fixed = TRUE),
                                                             fixed = TRUE)
                                            code.c <- paste0("-", code.c)
                                          } else if (substr(code.c, 1, 1) == "+") {
                                            code.c <- gsub("@", "-",
                                                           gsub("-", "+",
                                                                gsub("+", "@", code.c,
                                                                     fixed = TRUE),
                                                                fixed = TRUE),
                                                           fixed = TRUE)
                                            
                                          }
                                          
                                          # code.c <- gsub("@", "-",
                                          #                gsub("-", "+",
                                          #                     gsub("+", "@", code.c, fixed = TRUE), fixed = TRUE),
                                          #                fixed = TRUE)
                                          # 
                                          # if (!(substr(code.c, 1, 1) %in% c("+", "-")))
                                          #   code.c <- paste0("-", code.c)
                                          
                                        }
                                        
                                        return(code.c)
                                        
                                        # return(paste0(code.c, "(",
                                        #               floor(mzdiff_db.df[as.numeric(annot_split.c), "delta_mass"] * 1e4) / 1e4,
                                        #               ")")) 
                                        
                                      }
                                      
                                    }, character(1)), collapse = "|")
            
            
            
          } else
            annot.c <- ""
          
        }
        
        if (annot.c != "")
          group_feat_isoadfra.vc <- c(group_feat_isoadfra.vc, paste0("@", linked.c, "|", annot.c))
        
        
      }
      
      isoadfra.vc[group_feat.c] <- paste(group_feat_isoadfra.vc, collapse = "|")
    }
    
    ## annotation relative to the representative ion within a group
    
    relat.vc[ion_sel.c] <- "M"
    
    for (linked.c in setdiff(group_feat.vc, ion_sel.c)) {
      
      sign.i <- sign(mzdiff.mn[ion_sel.c, linked.c])
      
      annot.c <- corrtmz.mc[ion_sel.c, linked.c]
      
      if (annot.c != "") {
        
        annot_split.vc <- unlist(strsplit(annot.c, split = "_"))
        
        if (sign.i < 0)
          annot_split.vc <- annot_split.vc[mzdiff_db.df[as.numeric(annot_split.vc), "type"] != "isotope"]
        
        if (length(annot_split.vc)) {
          
          annot.c <- paste(vapply(annot_split.vc,
                                  function(annot_split.c) {
                                    
                                    if (as.numeric(annot_split.c) > nrow(mzdiff_db.df)) {
                                      ## The 2M annotation is coded as nrow(mzdiff_db.df) + 1
                                      
                                      if (sign.i > 0) {
                                        
                                        code.c <- "[2M]"
                                        
                                      } else {
                                        
                                        code.c <- "[0.5M]"
                                        
                                      }
                                      
                                      return(code.c)
                                      
                                      
                                    } else {
                                      
                                      code.c <- mzdiff_db.df[as.numeric(annot_split.c), "losses_or_gains"]
                                      
                                      if (sign.i > 0) {
                                        
                                        if (!(substr(code.c, 1, 1) %in% c("+", "-"))) {
                                          code.c <- paste0("+", code.c)
                                        } else if (substr(code.c, 1, 1) == "-") {
                                          code.c <- gsub("@", "-",
                                                         gsub("-", "+",
                                                              gsub("+", "@", code.c, fixed = TRUE), fixed = TRUE),
                                                         fixed = TRUE)
                                        }
                                        
                                      } else {
                                        
                                        if (!(substr(code.c, 1, 1) %in% c("+", "-"))) {
                                          if (code.c != "(H2O-CO)")
                                            code.c <- gsub("@", "-",
                                                           gsub("-", "+",
                                                                gsub("+", "@", code.c,
                                                                     fixed = TRUE),
                                                                fixed = TRUE),
                                                           fixed = TRUE)
                                          code.c <- paste0("-", code.c)
                                        } else if (substr(code.c, 1, 1) == "+") {
                                          code.c <- gsub("@", "-",
                                                         gsub("-", "+",
                                                              gsub("+", "@", code.c,
                                                                   fixed = TRUE),
                                                              fixed = TRUE),
                                                         fixed = TRUE)
                                          
                                        }
                                        
                                      }
                                      
                                      return(paste0("[M", code.c, "]"))
                                      
                                    }
                                    
                                  }, character(1)), collapse = "|")
          
          relat.vc[linked.c] <- annot.c
          
          
        }
        
      }
      
    }
    
  }
  
  feat.df[, "redund_is"] <- as.numeric(!(relat.vc %in% c("", "M")))
  feat.df[, "redund_group"] <- group.vi
  feat.df[, "redund_iso_add_frag"] <- isoadfra.vc
  feat.df[, "redund_repres"] <- repres.vc
  feat.df[, "redund_relative"] <- relat.vc
  
  # ordering the features according to the component
  
  data.mn <- data.mn[, order(group.vi), drop = FALSE]
  feat.df <- feat.df[order(group.vi), , drop = FALSE]
  
  group.vi <- feat.df[, "redund_group"]
  group_na.vi <- which(is.na(group.vi))
  
  if (length(group_na.vi) > 1) {
    
    group_na_min.i <- group_na.vi[1]
    group_order.vi <- c(seq_len(group_na_min.i - 1),
                        rep(group_na_min.i, length(group.vi) - group_na_min.i + 1))
    feat_order.vi <- order(group_order.vi,
                           feat.df[, mz_colname.c],
                           feat.df[, rt_colname.c])
    data.mn <- data.mn[, feat_order.vi, drop = FALSE]
    feat.df <- feat.df[feat_order.vi, , drop = FALSE]
    
  }
  
  message(length(table(feat.df[, "redund_group"])), " groups")
  
  message(sum(feat.df[, "redund_is"]), " chemically redundant features (",
          round(sum(feat.df[, "redund_is"]) / nrow(feat.df) * 100), "%)")
  
  list(data.mn = data.mn,
       feat.df = feat.df,
       adjacency.mi = corrtmz.mi)
  
}


.mzdiff_db <- function() {
  # table of referenced losses (fragments, adducts, isotopes)

  mzdiff_db.df <- read.table(system.file("extdata/mzdiff_db.tsv", package = "phenomis"),
                             header = TRUE,
                             quote = "",
                             sep = "\t",
                             stringsAsFactors = FALSE)
  
  mzdiff_db.df[ , "losses_or_gains"] <- gsub(" ", "", mzdiff_db.df[ , "losses_or_gains"])
  

  
  add_and_frag.vi <- integer(nrow(mzdiff_db.df))
  
  loss_or_gain.vc <- sapply(mzdiff_db.df[ , "losses_or_gains"],
                            function(loss_or_gain.c) {
                              if (substr(loss_or_gain.c, 1, 1) %in% c("+", "-"))
                                loss_or_gain.c <- substr(loss_or_gain.c, 2, nchar(loss_or_gain.c))
                              loss_or_gain.c
                            })
  
  for (add_and_frag.c in c("(CH3OH)",
                           "(H2O)",
                           "(HCOOH)",
                           "(NaCl)",
                           "2(H2O)",
                           "2(HCOOH)",
                           "(NaCl)")) {
    
    dup.vi <- which(loss_or_gain.vc == add_and_frag.c)
    stopifnot(length(dup.vi) == 2)
    
    add_and_frag.vi[dup.vi[2]] <- 1
    
  }

  mzdiff_db.df[, "add_and_frag_dup"] <- add_and_frag.vi
  
  mzdiff_db.df
}