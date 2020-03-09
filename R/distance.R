##' @rdname make_distance
##' @aliases make_relationship
##'
##' @title Generate Distance Matrix
##'
##' @description Compute the distance matrix of using shortest paths of a (directed) \code{\link[igraph]{igraph}} structure, normalising by the diameter of the network, preserving node/column/row names (and direction).
##'
##' @param mat precomputed adjacency or commonlink matrix.
##' @param graph An \code{\link[igraph]{igraph}} object. May be directed or weighted.
##' @param directed logical. Whether directed information is passed to the distance matrix.
##' @param absolute logical. Whether distances are scaled as the absolute difference from the diameter (maximum possible). Defaults to TRUE. The alternative is to calculate a relative difference from the diameter for a geometric decay in distance.
##' @keywords graph network igraph adjacency
##' @import igraph
##' @examples 
##' 
##' library("igraph")
##' graph_test_edges <- rbind(c("A", "B"), c("B", "C"), c("B", "D"))
##' graph_test <- graph.edgelist(graph_test_edges, directed = TRUE)
##' adjacency_matrix <- make_adjmatrix_graph(graph_test)
##' distance_matrix <- make_distance_adjmat(adjacency_matrix)
##' 
##' @return A numeric matrix of values in the range [0, 1] where lower values are closer
##' 
##' @export
make_distance_graph <- function(graph, directed = TRUE, absolute = FALSE){
  if(directed == FALSE) graph <- as.undirected(graph)
  diam <- diameter(graph)
  if (absolute){
    mat <- (diam-shortest.paths(graph))/diam
  } else {
    mat <- 1^-diam/(diam*shortest.paths(graph))
    diag(mat) <- 1
  }
  rownames(mat) <- colnames(mat) <- names(V(graph))
  return(mat)
}

##' @rdname make_distance
##' @importFrom igraph graph_from_adjacency_matrix
##' @export
make_distance_adjmat <- function(mat, directed = TRUE, absolute = FALSE){
  diag(mat) <- 0
  graph <- graph_from_adjacency_matrix(mat, weighted = NULL, mode = "undirected")
  diam <- diameter(graph)
  if (absolute){
    mat <- (diam-shortest.paths(graph))/diam
  } else {
    mat <- 1^-diam/(diam*shortest.paths(graph))
    diag(mat) <- 1
  }
  rownames(mat) <- colnames(mat) <- names(V(graph))
  return(mat)
}
