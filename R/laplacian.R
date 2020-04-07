##' @name make_laplacian
##' @aliases make_laplacian_graph
##' @rdname make_laplacian_graph
##'
##' @title Generate Laplacian Matrix
##' 
##' @description Compute the Laplacian matrix of a (directed) \code{\link[igraph]{igraph}} structure, preserving node/column/row names (and direction).
##'
##' @param graph An \code{\link[igraph]{igraph}} object. May be directed or weighted.
##' @param directed logical. Whether directed information is passed to the Laplacian matrix.
##' @keywords graph network igraph Laplacian
##' @importFrom igraph laplacian_matrix graph.edgelist as.undirected
##' @import igraph
##' @importFrom Matrix Matrix
##' @importClassesFrom Matrix dgCMatrix
##' @examples 
##' 
##' library("igraph")
##' graph_test_edges <- rbind(c("A", "B"), c("B", "C"), c("B", "D"))
##' graph_test <- graph.edgelist(graph_test_edges, directed = TRUE)
##' laplacian_matrix <- make_laplacian_graph(graph_test)
##' 
##' @return An Laplacian matrix compatible with generating an expression matrix
##' 
##' @export make_laplacian_graph
make_laplacian_graph <- function(graph, directed = FALSE){
  if(directed == FALSE) graph <- as.undirected(graph)
  laplacian <- as.matrix(laplacian_matrix(graph))
  rownames(laplacian) <- colnames(laplacian) <- names(V(graph))
  return(laplacian)
}
