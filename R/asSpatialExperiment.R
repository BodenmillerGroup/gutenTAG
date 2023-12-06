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
#' spe <- asSpatialExperiment(final)
#'
asSpatialExperiment <- function(x){

  # validity checks
  .valid.asSpatialExperiment(x)

  # put intensity df and coordinates in correct form
  mat <- t(matter::as.matrix(x$IntensityDF))
  crds <- as.matrix(x$SpatialCoords)

  # put into SpatialExperiment
  spe <- SpatialExperiment(
    assay = list(intenity = mat),
    colData = crds,
    spatialCoordsNames = c("x", "y"))

  return(spe)

}

