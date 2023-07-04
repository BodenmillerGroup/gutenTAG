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

nmf <- function(x, comp, cntr = FALSE){

  nmf_temp = nmf(scale(Final_intensity_matrix_targeted, center = cntr), k = comp)

}
