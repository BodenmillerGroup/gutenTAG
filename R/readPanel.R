#' Read .csv file panel
#'
#' @param path The path to the panel.
#'
#' @return A clean and ordered data frame with three columns:
#'   \code{Name} (cleaned marker name), \code{OriginalName} (name as
#'   supplied in the CSV before cleaning), and \code{FeatureMass}
#'   (target m/z value).
#'
#' @examples
#' panel_path <- system.file("extdata/ref_list.csv", package = "gutenTAG")
#' panel <- readPanel(path = panel_path)
#'
#' @importFrom dplyr arrange
#' @importFrom dplyr relocate
#' @importFrom utils read.delim
#'
#' @export


readPanel <- function(path){

  # read in .csv
  panel <- read.delim(path, sep = ",", header = TRUE, col.names = c("Name", "FeatureMass"))

  # sanity check: columns are correctly named
  if (is.numeric(panel$Name)){

    colnames(panel) = c("FeatureMass","Name")

  }

  # arrange feature mass in increasing order
  panel <- dplyr::arrange(panel, panel$FeatureMass)

  panel <- .cleanPanel(panel)

  return(panel)
}
