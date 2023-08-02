#' Title Perform Louvain clustering
#'
#' @param x Output of getIntensityDF() function.
#'
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
louvain_cluster <- function(x, metric = "angular"){

  intensity <- x$IntensityDF

  # 1. Build adjacency matrix
  print("step1")
  KNN_graph_matrix <- N2R::Knn(as.matrix(intensity), k = 30, verbose = FALSE, indexType = metric)
  print("step2")
  KNN_graph_matrix <- KNN_graph_matrix + t(KNN_graph_matrix)

  # 2. Build graph from adjacency matrix
  graph <- igraph::graph_from_adjacency_matrix(KNN_graph_matrix, mode = 'undirected', weighted = TRUE)

  # 3. Perform louvain clustering
  clustering <- igraph::cluster_louvain(graph)

  # get membership for each pixel
  membership <- as.character(igraph::membership(clustering))

  return(membership)

}
