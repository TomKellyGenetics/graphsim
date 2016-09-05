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
make_distance.igraph <- function(g, directed = T, absolute = F){
  if(directed == F) g <- as.undirected(g)
  diam <-diameter(g)
  if (absolute){
    mat <- (diam-shortest.paths(g))/diam
  } else {
    mat <- 1^-diam/(diam*shortest.paths(g))
    diag(mat) <- 1
  }
  rownames(mat) <- colnames(mat) <- names(V(g))
  return(mat)
}