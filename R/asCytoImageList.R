#' asCytoImageList
#'
#' @param final Output from assignMetapeaks function
#' @param name The name of the image
#'
#' @return A CytoImageList object
#'
#' @importFrom cytomapper CytoImageList
#'
#' @export
#'

asCytoImageList <- function(x, name = "test"){

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
