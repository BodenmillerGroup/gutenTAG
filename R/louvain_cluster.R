#' Title Perform Louvain clustering
#'
#' @param x Output of getIntensityDF function.
#' @param metric The distance metric to be used (angular, L2)
#'
#' @return The cluster membership for each pixel.
#'
#' @importFrom N2R Knn
#' @importFrom igraph graph_from_adjacency_matrix
#' @importFrom igraph cluster_louvain
#' @importFrom igraph membership
#' @export
#'
#' @examples
#' louvain_cluster(final)
louvain_cluster <- function(x, k = 30, metric = "angular", plot = FALSE){

  intensity <- x$IntensityDF

  # 1. Build adjacency matrix
  KNN_graph_matrix <- N2R::Knn(as.matrix(intensity), k = k, verbose = FALSE, indexType = metric)

  KNN_graph_matrix <- KNN_graph_matrix + t(KNN_graph_matrix)

  # 2. Build graph from adjacency matrix
  graph <- igraph::graph_from_adjacency_matrix(KNN_graph_matrix, mode = 'undirected', weighted = TRUE)

  # 3. Perform louvain clustering
  clustering <- igraph::cluster_louvain(graph, resolution = 1)

  # get membership for each pixel
  membership <- as.character(igraph::membership(clustering))

  # bind cluster membership to pixel coordinates
  coords <- x$SpatialCoords
  cluster_coords <- cbind(coords, membership)

  if(plot == TRUE){
    # ggplot spatial distribution of clusters
    ggplot(cluster_coords, aes(x = x, y = y, color = membership)) +
      geom_point(size = 1) +
      scale_color_discrete(name="Cluster Membership") +
      theme_minimal() +
      labs(title = "Cluster Memberships on Spatial Coordinates",
           x = "X Coordinate",
           y = "Y Coordinate")
  }

  print(length(unique(membership)))
  return(membership)

}
