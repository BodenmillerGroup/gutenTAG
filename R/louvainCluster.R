#' Title Perform Louvain clustering
#'
#' @param x Output of assignMetapeaks function.
#' @param k Number of nearest neighbours when building adjacency matrix.
#' @param coords If only providing an intensity dataframe for argument x, spatial coordinates must also be provided in the coords argument.
#' @param metric The distance metric to be used (angular, L2)
#' @param resolution Clustering parameter to determine size of communities. Small values produce larger communities, larger values produce smaller communities.
#'
#' @return The cluster membership for each pixel.
#'
#' @importFrom N2R Knn
#' @importFrom igraph graph_from_adjacency_matrix
#' @importFrom igraph cluster_louvain
#' @importFrom igraph membership
#' @importFrom Matrix t
#' @export louvainCluster

louvainCluster <- function(x, coords, k = 300, metric = "angular", resolution = 0.5){

  # validity checks
  .valid.louvainCluster(x, coords, k, metric, resolution)

  # if x is a list, take out the intensity df, if it is the intensity df itself, leave it
  if(is.null(dim(x))){
    intensity <- x$IntensityDF
  }else{
    intensity <- x
  }

  # 1. Build adjacency matrix
  KNN_graph_matrix <- Knn(as.matrix(intensity), k = k, verbose = FALSE, indexType = metric)

  KNN_graph_matrix <- as(KNN_graph_matrix, "dgCMatrix")

  # transpose matrix and make sparse
  t_KNN_graph_matrix <- t(KNN_graph_matrix)
  #t_KNN_graph_matrix <- as(t_KNN_graph_matrix, "dgCMatrix")

  KNN_graph_matrix <- KNN_graph_matrix + t_KNN_graph_matrix

  # 2. Build graph from adjacency matrix
  graph <- graph_from_adjacency_matrix(KNN_graph_matrix, mode = 'undirected', weighted = TRUE)

  # 3. Perform louvain clustering)
  clustering <- cluster_louvain(graph, resolution = resolution)

  # get membership for each pixel
  membership <- as.character(membership(clustering))

  # bind cluster membership to pixel coordinates
  coords <- x$SpatialCoords
  cluster_coords <- cbind(coords, membership)

  #if(plot == TRUE){
  #  # ggplot spatial distribution of clusters
  #  ggplot(cluster_coords, aes(x = x, y = y, color = membership)) +
  #    geom_point(size = 1) +
  #    #scale_color_discrete(name = "Cluster Membership") +
  #    scale_color_brewer(palette = "Set1", name = "Cluster Membership") +
  #    theme_minimal() +
  #    labs(title = "Cluster Memberships on Spatial Coordinates",
  #         x = "X Coordinate",
  #         y = "Y Coordinate")
  #}

  return(membership)

}
