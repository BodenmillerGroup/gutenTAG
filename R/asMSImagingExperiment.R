#' Make MSImagingExperiment object
#'
#' @param x The output from 'assignMetapeaks' function
#' @param remove.na Removes NA values from the experiment. NA values are generated when a marker in the panel is not assigned a metapeak. The default is FALSE.
#'
#' @return An object of the class MSImagingExperiment from Cardinal.
#'
#' @import Cardinal
#' @importFrom matter as.matrix
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
#' processed <- assignMetapeaks(metapeaks, pre = pre, refList = panel)
#'

asMSImagingExperiment <- function(x, remove.na = FALSE){

  # validity checks
  .valid.asMSImagingExperiment(x)

  # pdata: where spatial coordinates go
  coord <- x$SpatialCoords
  run <- factor(rep("run0", nrow(coord)))
  pdata <- PositionDataFrame(run = run, coord = coord)

  # ordered vector of mass locations. metapeak mz stored as mz, expected mz stored as expected_mz
  fdata <- MassDataFrame(mz = sort(x$CorrespondenceMatrix$expected_mz_location),
                         observed_mz = x$CorrespondenceMatrix$mz_location,
                         expected_mz = sort(x$CorrespondenceMatrix$expected_mz_location))


  # idata: where intensity df goes
  idata <- t(as.matrix(x$IntensityDF)) # intensity dataframe must first be a matrix, then transposed
  # get m/z order of marker names
  correct_order <- x$CorrespondenceMatrix$marker
  # reorder rows (markers) of idata according to m/z order
  idata <- idata[correct_order, ]

  # remove NA clause
  if (remove.na == TRUE){
    keep <- !is.na(x$CorrespondenceMatrix$mz_location)

    # remove from idata
    idata <- idata[keep, ]
    # remove from fdata
    fdata<- fdata[keep, ]
  }

  # put them all together
  out <- MSImagingExperiment(spectraData = idata,
                             featureData = fdata,
                             pixelData = pdata)

  return(out)


}
