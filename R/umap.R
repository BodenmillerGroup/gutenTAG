#' Uniform Manifold Approximate Projection (UMAP)
#'
#' @param x A targeted intensity dataframe (rows: pixels, cols: features)
#'
#' @return object of class umap, containing at least a component with an embedding and a component with configuration settings
#' @export
#'
#' @examples
#' umap(IntensityDF)

umap <- function(x){

  umap_representation <- umap(x)

  return(umap_representation)

}


