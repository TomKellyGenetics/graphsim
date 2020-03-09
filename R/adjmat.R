##' @name make_adjmatrix
##' @aliases make_adjmatrix_graph
##' @rdname make_adjmatrix_graph
##'
##' @title Generate Adjacency Matrix
##' 
##' @description Compute the adjacency matrix of a (directed) \code{\link[igraph]{igraph}} structure, preserving node/column/row names (and direction).
##'
##' @param graph An \code{\link[igraph]{igraph}} object. May be directed or weighted.
##' @param directed logical. Whether directed information is passed to the adjacency matrix.
##' @keywords graph network igraph adjacency
##' @import igraph
##' @examples 
##' 
##' library("igraph")
##' graph_test_edges <- rbind(c("A", "B"), c("B", "C"), c("B", "D"))
##' graph_test <- graph.edgelist(graph_test_edges, directed = TRUE)
##' adjacency_matrix <- make_adjmatrix_graph(graph_test)
##' 
##' @return An adjacency matrix compatible with generating an expression matrix
##' 
##' @export make_adjmatrix_graph
make_adjmatrix_graph <- function(graph, directed = FALSE){
  if(directed == FALSE) graph <- as.undirected(graph)
  adj_mat <- as.matrix(as_adjacency_matrix(graph))
  rownames(adj_mat) <- colnames(adj_mat) <- names(V(graph))
  return(adj_mat)
}