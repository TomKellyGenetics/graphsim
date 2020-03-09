##' @name make_commonlink
##' @rdname make_commonlink
##'
##' @title Generate Common Link Matrix
##'
##' @description Compute the common link matrix of a (directed) \code{\link[igraph]{igraph}} structure, preserving node / column / row names (and direction).
##'
##' @param adj_mat precomputed adjacency matrix.
##' @param graph An \code{\link[igraph]{igraph}} object. May be directed or weighted.
##' @param directed logical. Whether directed information is passed to the adjacency matrix.
##' @keywords graph network igraph neighbourhood
##' @import igraph
##' @examples 
##' 
##' library("igraph")
##' graph_test_edges <- rbind(c("A", "B"), c("B", "C"), c("B", "D"))
##' graph_test <- graph.edgelist(graph_test_edges, directed = TRUE)
##' adjacency_matrix <- make_adjmatrix_graph(graph_test)
##' common_link_matrix <- make_commonlink_adjmat(adjacency_matrix)
##' 
##' @return An integer matrix of number of links shared between nodes
##' 
##' @export
make_commonlink_adjmat <- function(adj_mat){
  comm_mat <- matrix(NA, nrow(adj_mat), ncol(adj_mat))
  for(ii in 1:nrow(adj_mat)){
    for(jj in 1:ncol(adj_mat)){
      comm_mat[ii, jj] <- sum(adj_mat[ii,] & adj_mat[,jj])
    }
  }
  rownames(comm_mat) <- rownames(adj_mat)
  colnames(comm_mat) <- colnames(adj_mat)
  return(comm_mat)
}

##' @rdname make_commonlink
##' @export
make_commonlink_graph <- function(graph, directed = FALSE){
  adj_mat <- make_adjmatrix_graph(graph, directed = directed)
  comm_mat <- make_commonlink_adjmat(adj_mat)
}
