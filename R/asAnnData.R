#' Title asAnnData
#'
#' @param x the output from `gutenTAG` function `assignMetapeaks()`
#'
#' @return an object of class AnnData
#'
#' @importFrom anndata AnnData
#' @export
#'
#' @examples
#' path <- system.file("extdata/Example_data.imzML", package = "gutenTAG")
#' panel_path <- system.file("extdata/ref_list.csv", package = "gutenTAG")
#' panel <- readPanel(path = panel_path)
#' raw <- readMSIData(path)
#' pre <- preProcess(raw, cores = 2)
#' peaks <- peakDetection(pre, core = 2)
#' metapeaks <- generateMetapeaks(peaks)
#' final <- assignMetapeaks(x = metapeaks, pre, panel)
#'
asAnnData <- function(x) {
  # extract data from object
  exp_data <- as.matrix(x$IntensityDF)
  coords <- x$SpatialCoords
  var_data <- x$CorrespondenceMatrix

  rownames(var_data) <- rep(0:(nrow(var_data) - 1))

  # create obsdata
  obs_data <- coords

  # create anndata object
  output <- anndata::AnnData(
    X = exp_data,
    obs = data.frame(obs_data),
    var = data.frame(var_data)
  )

  return(output)
}
