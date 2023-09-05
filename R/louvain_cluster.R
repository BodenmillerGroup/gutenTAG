#' Title Perform Louvain clustering
#'
#' @param x Output of getIntensityDF function.
#' @param k Number of nearest neighbours when building adjacency matrix.
#' @param coords If only providing an intensity dataframe for argument x, spatial coordinates must also be provided in the coords argument.
#' @param metric The distance metric to be used (angular, L2)
#' @param plot Plot the cluster membership over the spatial coordinates.
#' @param resolution Clustering parameter to determine size of communities. Small values produce larger communities, larger values produce smaller communities.
#'
#' @return The cluster membership for each pixel.
#'
#' @importFrom N2R Knn
#' @importFrom igraph graph_from_adjacency_matrix
#' @importFrom igraph cluster_louvain
#' @importFrom igraph membership
#' @importFrom Matrix t
#' @export
#'
#' @examples
#' louvain_cluster(final)
louvain_cluster <- function(x, coords, k = 300, metric = "angular", plot = FALSE, resolution = 0.5){

  if(is.null(dim(x))){
    intensity <- x$IntensityDF
  }else{
    intensity <- x
  }

  # 1. Build adjacency matrix
  KNN_graph_matrix <- N2R::Knn(as.matrix(intensity), k = k, verbose = FALSE, indexType = metric)

  KNN_graph_matrix <- as(KNN_graph_matrix, "dgCMatrix")

  # transpose matrix and make sparse
  t_KNN_graph_matrix <- Matrix::t(KNN_graph_matrix)
  #t_KNN_graph_matrix <- as(t_KNN_graph_matrix, "dgCMatrix")

  KNN_graph_matrix <- KNN_graph_matrix + t_KNN_graph_matrix

  # 2. Build graph from adjacency matrix
  graph <- igraph::graph_from_adjacency_matrix(KNN_graph_matrix, mode = 'undirected', weighted = TRUE)

  # 3. Perform louvain clustering)
  clustering <- igraph::cluster_louvain(graph, resolution = resolution)

  # get membership for each pixel
  membership <- as.character(igraph::membership(clustering))

  # bind cluster membership to pixel coordinates
  coords <- x$SpatialCoords
  cluster_coords <- cbind(coords, membership)

  if(plot == TRUE){
    # ggplot spatial distribution of clusters
    ggplot(cluster_coords, aes(x = x, y = y, color = membership)) +
      geom_point(size = 1) +
      #scale_color_discrete(name = "Cluster Membership") +
      scale_color_brewer(palette = "Set1", name = "Cluster Membership") +
      theme_minimal() +
      labs(title = "Cluster Memberships on Spatial Coordinates",
           x = "X Coordinate",
           y = "Y Coordinate")
  }

  return(membership)

}
