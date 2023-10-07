#' Uniform Manifold Approximate Projection (UMAP)
#'
#' @param x A targeted intensity dataframe (rows: pixels, cols: features)
#' @param seed A single number containing to the random number generator state for controlled randomness.
#'
#' @return a dataframe containing UMAP coordinates
#'
#' @importFrom uwot umap
#' @export
#'

computeUMAP <- function(x, seed = NULL){

  # validity checks
  .valid.computeUMAP(x, seed)

  # compute UMAP
  umap_representation <- umap(x, seed = seed)

  # store as data frame
  df <- data.frame(UMAP1 = umap_representation[,1],
                   UMAP2 = umap_representation[,2])

  return(df)
}


