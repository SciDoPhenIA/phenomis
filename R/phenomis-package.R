#' Phenomics data analysis and integration
#'
#' Package description
#'
#' @name phenomis-package
#' @aliases phenomis phenomis-package
#' @docType package
#' @author E. A. Thevenot (CEA, LIST, MetaboHUB)
#'
#' Maintainer: Etienne Thevenot <etienne.thevenot@@cea.fr>
#' @keywords package
#' @examples
#' # See the package vignette
#'
NULL

#' Searching patterns in code
#'
#' This function is useful for code inspection, development, and debugging
#' 
#' @param pattern.c Character: pattern to be searched in the files
#' @param dir.c Character: parent directory of the files
#' @param fixed.l Logical: should the pattern be searched 'as is'
#' @param filename_pattern.c Character: optional filter of file names to be searched
#' @return nothing is returned
#' @export
#' @examples
#' \dontrun{
#' phenomis::search_code("aes", filename_pattern.c = ".R")
#' }
search_code <- function(pattern.c,
                        dir.c = getwd(),
                        fixed.l = TRUE,
                        filename_pattern.c = NULL) {
  
  warn_option.n <- getOption("warn")
  options(warn = -1)
  
  if (is.na(file.info(dir.c)[["isdir"]]))
    stop(dir.c, " is not a readable directory")
  
  file.vc <- list.files(dir.c,
                        full.names = TRUE,
                        pattern = filename_pattern.c,
                        recursive = TRUE)
  
  found.l <- FALSE
  
  for (file.c in file.vc) {
    
    lines.vc <- readLines(file.c)
    
    pattern.vi <- grep(pattern.c,
                       x = lines.vc,
                       fixed = fixed.l)
    
    if (length(pattern.vi) > 0) {
      
      message("\nPattern found in '",
              gsub(dir.c, "", file.c, fixed = TRUE),
              "' in the following lines:")
 
        for (pattern.i in pattern.vi)
          message("line ", pattern.i, ": ", lines.vc[pattern.i])
      
      found.l <- TRUE
      
    }
  }
  
  if (!found.l) {
    message("The pattern was not found in any of the searched files:")
    message(paste(file.vc, collapse = "\n"))
  }
  
  options(warn = warn_option.n)
  
  invisible(NA)
  
}

.sample_color_eset <- function(eset,
                               col_sampleType.c = "sampleType") {
  
  if (col_sampleType.c %in% Biobase::varLabels(eset)) {
    
    sample_types.vc <- Biobase::pData(eset)[, col_sampleType.c]
    
  } else {
    
    sample_types.vc <- rep("other", Biobase::dims(eset)["Samples", 1])
    
  }
  
  .sample_color_vector(sample_types.vc = sample_types.vc)
  
}

.sample_color_vector <- function(sample_types.vc) {
  
  type_colors.vc <- c(sample = "green4",
                      pool = RColorBrewer::brewer.pal(9, "Reds")[7],
                      pool1 = RColorBrewer::brewer.pal(9, "Reds")[5],
                      pool2 = RColorBrewer::brewer.pal(9, "Reds")[4],
                      pool4 = RColorBrewer::brewer.pal(9, "Reds")[3],
                      pool8 = RColorBrewer::brewer.pal(9, "Reds")[2],
                      blank = "black",
                      other = RColorBrewer::brewer.pal(9, "Blues")[7])
  
  sample_types.vc[!(sample_types.vc %in% names(type_colors.vc))] <- "other"
  
  type_colors.vc[sample_types.vc]
  
}



















