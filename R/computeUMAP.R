#' Uniform Manifold Approximate Projection (UMAP)
#'
#' @param x A targeted intensity dataframe (rows: pixels, cols: features)
#'
#' @return object of class umap, containing at least a component with an embedding and a component with configuration settings
#'
#' @import ggplot2
#' @importFrom umap umap
#' @export
#'s

computeUMAP <- function(x, plot = FALSE){

  umap_representation <- umap::umap(x)

  if(plot == TRUE){

    umap_coordinates <- umap_representation$layout

    df <- data.frame(UMAP1 = umap_coordinates[,1],
                     UMAP2 = umap_coordinates[,2])

    scatter_plot <- ggplot(df) +
      geom_point(aes(x = UMAP1, y = UMAP2), alpha = 0.9) +
      scale_color_viridis_c() +
      theme_minimal() +
      theme(legend.position = "right") +
      labs(title = "UMAP on MALDI-imaging data",
           x = "UMAP Dimension 1",
           y = "UMAP Dimension 2")

    print(scatter_plot)
  }

  return(umap_representation)
}


