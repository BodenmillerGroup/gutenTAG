#' Non-negative matrix factorisation
#'
#' @param x A targeted intensity dataframe (rows: pixels, cols: features)
#' @param comp The number of principal components to compute.
#' @param cntr Center the data or not (bool)
#'
#' @return ?
#' @export
#'
#' @examples
#' nmf(x = IntensityDF, comp = 5)

nmf <- function(x, comp){

  library(RcppML)

  nmf_temp <- RcppML::nmf(scale(x, center = FALSE), k = comp)

}


