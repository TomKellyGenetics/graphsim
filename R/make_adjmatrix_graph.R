##' @name Generate Adjacency Matrix from Graph Structure
##' @rdname make_adjmatrix_
##'
##' @title Generate Adjacency Matrix
##'
##' @description Compute the adjacency matrix of a (directed) \code{\link[igraph]{igraph}} structure, preserving node/column/row names (and direction).
##'
##' @param graph An \code{\link[igraph]{igraph}} object. May be directed or weighted.
##' @param directed logical. Whether directed information is passed to the adjacency matrix.
##' @keywords graph network igraph adjacency
##' @import igraph
##' @export
make_adjmatrix_graph <- function(graph, directed = F){
  if(directed == F) graph <- as.undirected(graph)
  adj_mat <- as.matrix(as_adjacency_matrix(graph))
  rownames(adj_mat) <- colnames(adj_mat) <- names(V(graph))
  return(adj_mat)
}
