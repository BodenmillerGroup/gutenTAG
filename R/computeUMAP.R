#' Uniform Manifold Approximate Projection (UMAP)
#'
#' @param x A targeted intensity dataframe (rows: pixels, cols: features)
#'
#' @return a dataframe containing UMAP coordinatees
#'
#' @importFrom uwot umap
#' @export
#'

computeUMAP <- function(x){

  umap_representation <- uwot::umap(x)

  df <- data.frame(UMAP1 = umap_representation[,1],
                   UMAP2 = umap_representation[,2])

  return(df)
}


