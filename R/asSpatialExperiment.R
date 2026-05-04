#' Make SpatialExperiment object
#'
#' @param x The output from 'assignMetapeaks' function
#'
#' @return An object of the class SpatialExperiment from Cardinal.
#'
#' @importFrom SpatialExperiment SpatialExperiment
#' @importFrom SpatialExperiment spatialCoordsNames
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
#' processed <- assignMetapeaks(metapeaks, pre = pre, refList = panel)
#'
#' spe <- asSpatialExperiment(processed)
#'
asSpatialExperiment <- function(x){

  # validity checks
  .valid.asSpatialExperiment(x)

  # put intensity df and coordinates in correct form
  mat <- t(matter::as.matrix(x$IntensityDF))
  crds <- as.matrix(x$SpatialCoords)

  # put into SpatialExperiment
  spe <- SpatialExperiment(
    assay = list(intensity = mat),
    colData = crds,
    spatialCoordsNames = c("x", "y"))

  return(spe)

}

