#' Make SpatialExperiment object
#'
#' @param final The output from assignMetapeaks() function
#'
#' @return
#'
#' @importFrom SpatialExperiment SpatialExperiment
#' @importFrom SpatialExperiment spatialCoordsNames
#'
#' @export
#'
#' @examples
#' spe <- asSpatialExperiment(final)
#'
asSpatialExperiment <- function(final){

  # put intensity df and coordinates in correct form
  mat <- t(matter::as.matrix(final$IntensityDF))
  crds <- as.matrix(final$SpatialCoords)

  # put into SpatialExperiment
  spe <- SpatialExperiment(
    assay = list(intenity = mat),
    colData = crds,
    spatialCoordsNames = c("x", "y"))

  return(spe)

}

