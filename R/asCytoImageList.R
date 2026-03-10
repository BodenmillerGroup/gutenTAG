#' asCytoImageList
#'
#' @param x Output from assignMetapeaks function
#' @param name The name of the image
#'
#' @return A CytoImageList object
#'
#' @importFrom abind abind
#' @importFrom cytomapper CytoImageList
#' @importFrom cytomapper channelNames
#' @importFrom EBImage Image
#'
#' @export
#' @examples
#' path <- system.file("extdata/Example_data.imzML", package = "gutenTAG")
#' panel_path <- system.file("extdata/ref_list.csv", package = "gutenTAG")
#' panel <- readPanel(path = panel_path)
#' raw <- readMSIData(path)
#' pre <- preProcess(raw, cores = 2)
#' peaks <- peakDetection(pre, core = 2)
#' metapeaks <- generateMetapeaks(peaks, hist_smooth_factor = 1)
#' processed <- assignMetapeaks(metapeaks, pre = pre, refList = panel)
#'
#' cil <- asCytoImageList(processed, name = "test")

asCytoImageList <- function(x, name = "test"){

  # run validity checks
  .valid.asCytoImageList(x, name)

  dataframe <- x$IntensityDF
  coords <- x$SpatialCoords

  matrix_list <- lapply(colnames(dataframe), .curateMatrix, coords = coords, dataframe = dataframe)
  # Create multidimensional array
  my_array <- do.call(abind, list(matrix_list, along = 3))
  # Create multichannel image from array
  my_image <- Image(my_array)
  # Create CytoImageList object from image
  my_image <- CytoImageList(my_image)
  # Assign image and channel names
  names(my_image) <- name
  cytomapper::channelNames(my_image) <- colnames(dataframe)

  return(my_image)

}
