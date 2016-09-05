##' @name Generate Common Link Matrix from Graph Structure
##' @rdname make_commonlink
##'
##' @title Generate Common Link Matrix
##'
##' @description Compute the common link matrix of a (directed) \code{\link[igraph]{igraph}} structure, preserving node/column/row names (and direction).
##'
##' @param adj_mat precomputed adjacency matrix.
##' @param graph An \code{\link[igraph]{igraph}} object. May be directed or weighted.
##' @param directed logical. Whether directed information is passed to the adjacency matrix.
##' @keywords graph network igraph neighbohood
##' @import igraph
make_commonlink.matrix <- function(adj_mat){
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

make_commonlink_graph <- function(g, directed = F){
  adj_mat <- make_adjmatrix_graph(g, directed = directed)
  comm_mat <- make_commonlink_adjmat(adj_mat)
}
