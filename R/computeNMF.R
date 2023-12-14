#' Non-negative matrix factorisation
#'
#' @param x A targeted intensity dataframe (rows: pixels, cols: features)
#' @param comp The number of principal components to compute.
#' @param seed A number to generate a random seed
#' @param cntr Center the data or not (bool)
#'
#' @return Non-ordered factors of the input dataset
#'
#' @importFrom RcppML nmf
#' @export
#'
#' @examples
#' path <- system.file("extdata/Example_data.imzML", package = "maldipackage")
#' panel_path <- system.file("extdata/ref_list.csv", package = "maldipackage")
#' panel <- readPanel(path = panel_path)
#' raw <- readMSIData(path)
#' pre <- preProcess(raw, cores = 2)
#' peaks <- peakDetection(pre, core = 2)
#' metapeaks <- generateMetapeaks(peaks)
#' processed <- assignMetapeaks(x = metapeaks, pre, panel)
#' computeNMF(x = processed$IntensityDF, comp = 5)

computeNMF <- function(x, comp, seed = 234, cntr = FALSE){

  set.seed(seed)
  .valid.computeNMF(x, comp, cntr)

  nmf_temp <- RcppML::nmf(A = scale(x, center = cntr), k = comp)
  return(nmf_temp)

}


