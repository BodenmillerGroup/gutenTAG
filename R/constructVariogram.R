#' Construct a variogram for each dataframe channel.
#'
#' @param x A targeted intensity dataframe.
#' @param coords A spatial coordinates dataframe.
#'
#' @return A list of variograms for each marker.
#'
#' @import gstat
#' @export
#'
#' @examples
#' constructVariogram(IntensityDF, coords)
constructVariogram <- function(x, coords){

  # Include make.names to remove problematic characters, such as spaces, +, etc.
  names(df) <- make.names(names(df))
  names <- colnames(df)

  # Create list of variograms for each marker
  variograms <- list()
  for (i in 1:ncol(df)){

    marker <- data.frame(df[i])
    colnames(marker) <- names[i]

    # Assign coordinates
    marker$x <- coords$x
    marker$y <- coords$y
    coordinates(marker) <- ~x+y

    # Construct a formula from the current column name
    formula_string <- paste0(names[i], "~1")
    formula_object <- as.formula(formula_string)

    # Calculate variogram for the marker
    marker_variogram <- gstat::variogram(formula_object, marker)

    # Store the result in the list
    variograms[[names[i]]] <- marker_variogram

  }

  return(variograms)

}

