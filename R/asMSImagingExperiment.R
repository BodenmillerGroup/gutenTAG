#' Make MSImagingExperiment object
#'
#' @param x The output from 'assignMetapeaks' function
#'
#' @return An object of the class MSImagingExperiment from Cardinal.
#'
#' @import Cardinal
#' @importFrom matter as.matrix
#' @export
#'
#' @examples
#' asMSImagingExperiment(x)

asMSImagingExperiment <- function(x){

  # validity checks
  .valid.asMSImagingExperiment(x)

  # pdata: where spatial coordinates go
  coord <- x$SpatialCoords
  run <- factor(rep("run0", nrow(coord)))
  pdata <- PositionDataFrame(run = run, coord = coord)

  # idata: where intensity df goes
  idata <- t(as.matrix(x$IntensityDF)) # intensity dataframe must first be a matrix, then transposed

  # ordered vector of mass locations. metapeak mz stored as mz, expected mz stored as expected_mz
  fdata <- MassDataFrame(mz = sort(x$CorrespondenceMatrix$mz_location),
                         expected_mz = sort(x$CorrespondenceMatrix$expected_mz_location))

  # put them all together
  out <- MSImagingExperiment(imageData=idata,
                             featureData=fdata,
                             pixelData=pdata)

  return(out)


}
