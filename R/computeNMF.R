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
#' computeNMF(x = IntensityDF, comp = 5)

computeNMF <- function(x, comp, seed = 234, cntr = FALSE){

  set.seed(seed)
  .valid.computeNMF(x, comp, cntr)

  nmf_temp <- RcppML::nmf(A = scale(x, center = cntr), k = comp)
  return(nmf_temp)

}


