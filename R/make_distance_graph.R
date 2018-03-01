##' @rdname make_distance
##'
##' @title Generate Distance Matrix
##'
##' @description Compute the distance matrix of using shortest paths of a (directed) \code{\link[igraph]{igraph}} structure, normalising by the diameter of the network, preserving node/column/row names (and direction).
##'
##' @param graph An \code{\link[igraph]{igraph}} object. May be directed or weighted.
##' @param directed logical. Whether directed information is passed to the distance matrix.
##' @param absolute logical. Whether distances are scaled as the absolute difference from the diameter (maximum possible). Defaults to TRUE. The alternative is to calculate a relative difference from the diameter for a geometric decay in distance.
##' @keywords graph network igraph adjacency
##' @import igraph
##' @export
make_distance_graph <- function(graph, directed = T, absolute = F){
  if(directed == F) graph <- as.undirected(graph)
  diam <-diameter(graph)
  if (absolute){
    mat <- (diam-shortest.paths(graph))/diam
  } else {
    mat <- 1^-diam/(diam*shortest.paths(graph))
    diagraph(mat) <- 1
  }
  rownames(mat) <- colnames(mat) <- names(V(graph))
  return(mat)
}
