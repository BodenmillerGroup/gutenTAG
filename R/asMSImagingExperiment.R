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
#' path <- system.file("extdata/Example_data.imzML", package = "maldipackage")
#' panel_path <- system.file("extdata/ref_list.csv", package = "maldipackage")
#' panel <- readPanel(path = panel_path)
#' raw <- readMSIData(path)
#' pre <- preProcess(raw, cores = 2)
#' peaks <- peakDetection(pre, core = 2)
#' metapeaks <- generateMetapeaks(peaks)
#' processed <- assignMetapeaks(metapeaks, pre = pre, refList = panel)
#'
#' asMSImagingExperiment(processed)

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
  fdata <- MassDataFrame(mz = sort(x$CorrespondenceMatrix$expected_mz_location),
                         observed_mz = x$CorrespondenceMatrix$mz_location,
                         expected_mz = sort(x$CorrespondenceMatrix$expected_mz_location))

  # put them all together
  out <- MSImagingExperiment(imageData = idata,
                             featureData = fdata,
                             pixelData = pdata)

  return(out)


}
