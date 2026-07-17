
#' Title asAnnData
#'
#' @param x the output from `gutenTAG` function `assignMetapeaks()`
#'
#' @return a Python \code{anndata.AnnData} object (a \pkg{reticulate} reference)
#'
#' @export
#'
#' @examples
#' path <- system.file("extdata/Example_data/Example_data.imzML", package = "gutenTAG")
#' panel_path <- system.file("extdata/ref_list.csv", package = "gutenTAG")
#' panel <- readPanel(path = panel_path)
#' raw <- readMSIData(path)
#' pre <- preProcess(raw)
#' peaks <- peakDetection(pre)
#' metapeaks <- generateMetapeaks(peaks, hist_smooth_factor = 1)
#' final <- assignMetapeaks(x = metapeaks, pre, panel)
#'

asAnnData <- function(x){

  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("Package 'reticulate' is required for asAnnData(). Please install it ",
         "and ensure the Python 'anndata' module is available.")
  }

  # extract data from object
  exp_data <- as.matrix(x$IntensityDF)
  var_data <- x$CorrespondenceMatrix
  rownames(var_data) <- rep(0:(nrow(var_data) - 1))

  # create obsdata
  obs_data <- x$SpatialCoords

  # Build the AnnData object directly via reticulate. The R `anndata` wrapper
  # injects a `dtype` argument that Python anndata >= 0.11 removed, so we call
  # the Python class ourselves and omit it.
  if (utils::packageVersion("reticulate") >= "1.41.0") {
    reticulate::py_require("anndata")
  }
  ad <- reticulate::import("anndata")
  output <- ad$AnnData(X = exp_data,
                       obs = data.frame(obs_data),
                       var = data.frame(var_data))

  return(output)
}

