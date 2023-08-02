#' Read .csv file panel
#'
#' @param path The path to the panel.
#'
#' @return The cleaned and ordered panel.
#'
#' @importFrom dplyr arrange
#' @importFrom dplyr relocate
#'
#' @export
#'
#' @examples
#' readPanel("/path/to/panel.csv")
readPanel <- function(path){

  # read in .csv
  panel <- read.delim(path, sep = ",", header = TRUE, col.names = c("Name", "FeatureMass"))

  # sanity check: columns are correctly named
  if (is.numeric(panel$Name) == TRUE){

    colnames(panel) = c("FeatureMass","Name")

  }

  # arrange feature mass in increasing order
  panel <- dplyr::arrange(panel, panel$FeatureMass)

  panel <- .cleanPanel(panel)

  return(panel)
}
