#' Non-negative matrix factorisation
#'
#' @param x A targeted intensity dataframe (rows: pixels, cols: features)
#' @param comp The number of principal components to compute.
#' @param cntr Center the data or not (bool)
#'
#' @return Non-ordered factors of the input dataset
#'
#' @importFrom RcppML nmf
#' @export
#'
#' @examples
#' nmf(x = IntensityDF, comp = 5)

nmf <- function(x, comp){

  nmf_temp <- RcppML::nmf(scale(x, center = FALSE), k = comp)
  return(nmf_temp)

}


