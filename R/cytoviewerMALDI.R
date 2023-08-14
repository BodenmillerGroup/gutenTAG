#' cytoviewerMALDI
#'
#' @param final Output from getIntensityDF function
#' @param name The name of the image
#'
#' @return Cytoviewer shiny app
#'
#' @importFrom cytomapper CytoImageList
#' @import cytoviewer
#'
#' @export
#'

cytoviewerMALDI <- function(final, name = "test"){

  dataframe <- final$IntensityDF
  coords <- final$SpatialCoords

  matrix_list <- lapply(colnames(dataframe), .curateMatrix, coords = coords, dataframe = dataframe)
  # Create multidimensional array
  my_array <- do.call(abind, list(matrix_list, along = 3))
  # Create multichannel image from array
  my_image <- Image(my_array)
  # Create CytoImageList object from image
  my_image <- CytoImageList(my_image)
  # Assign image and channel names
  names(my_image) <- name
  channelNames(my_image) <- colnames(dataframe)

  cytoviewer(image = my_image)


}
