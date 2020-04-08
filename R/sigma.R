##' @name make_sigma
##' @rdname make_sigma
##'
##' @title Generate Sigma Matrix
##'
##' @description Compute the Sigma matrix from an \code{\link[igraph]{igraph}} structure 
##' or pre-computed matrix. These are compatible with \code{\link[mvtnorm]{rmvnorm}} and
##' \code{\link[graphsim]{generate_expression}}.
##'
##' @param mat precomputed adjacency, commonlink, or scaled distance matrix.
##' @param graph An \code{\link[igraph]{igraph}} object. May be directed or weighted.
##' @param cor numeric. Simulated maximum correlation/covariance of two adjacent nodes. Default to 0.8.
##' @param directed logical. Whether directed information is passed to the distance matrix.
##' @param comm logical whether a common link matrix is used to compute sigma. Defaults to FALSE (adjacency matrix).
##' @param laplacian logical whether a Laplacian matrix is used to compute sigma. Defaults to FALSE (adjacency matrix).
##' @param absolute logical. Whether distances are scaled as the absolute difference from
##' the diameter (maximum possible). Defaults to TRUE. The alternative is to calculate a
##' relative difference from the diameter for a geometric decay in distance.
##' @keywords graph network igraph mvtnorm
##' @importFrom igraph as_adjacency_matrix
##' @import igraph
##' @examples 
##' 
##' library("igraph")
##' graph_test_edges <- rbind(c("A", "B"), c("B", "C"), c("B", "D"))
##' graph_test <- graph.edgelist(graph_test_edges, directed = TRUE)
##' adjacency_matrix <- make_adjmatrix_graph(graph_test)
##' sigma_matrix <- make_sigma_mat_adjmat(adjacency_matrix, cor = 0.8)
##' 
##' @return a numeric covariance matrix of values in the range [-1, 1]
##' @export
make_sigma_mat_adjmat <- function(mat, cor = 0.8){
  sig <- ifelse(mat > 0, cor, 0)
  diag(sig) <- 1
  rownames(sig) <- rownames(mat)
  colnames(sig) <- colnames(mat)
  return(sig)
}

##' @rdname make_sigma
##' @export
make_sigma_mat_laplacian <- function(mat, cor = 0.8){
  sig <- ifelse(mat < 0, mat*cor/min(mat), 0)
  diag(sig) <- 1
  rownames(sig) <- rownames(mat)
  colnames(sig) <- colnames(mat)
  return(sig)
}

##' @rdname make_sigma
##' @export
make_sigma_mat_graph <- function(graph, cor = 0.8, comm = FALSE, laplacian = FALSE, directed = FALSE){
  if(comm && laplacian){
    warning("Error: only one of commonlink or laplacian can be used")
    stop()
  }
  if(!comm && !laplacian) mat <- make_adjmatrix_graph(graph, directed = directed)
  if(comm) mat <- make_commonlink_adjmat(mat)
  if(laplacian){
    mat <- make_laplacian_graph(graph, directed = directed)
    mat <- abs(mat)
    diag(mat)[diag(mat) == 0] <- 1
    mat <- apply(mat, 1, function(x) x/max(x))
    mat <- apply(mat, 2, function(x) x/max(x))
  } else{
    diag(mat) <- 1
  }
  sig <- ifelse(mat>0, cor*mat/max(mat), 0)
  diag(sig) <- 1
  rownames(sig) <- rownames(mat)
  colnames(sig) <- colnames(mat)
  return(sig)
}

##' @rdname make_sigma
##' @export
make_sigma_mat_dist_adjmat <- function(mat, cor = 0.8, absolute = FALSE){
  if(!(all(diag(mat) == 1))) stop("distance matrix must have diagonal of zero")
  if(!(max(mat[mat != 1]) > 0) || !(max(mat[mat!=1]) <= 1)) stop("distance matrix expected, not adjacency matrix")
  sig <- mat/max(mat[mat != 1]) * cor
  sig <- ifelse(sig > 0, sig, 0)
  diag(sig) <- 1
  rownames(sig) <- rownames(mat)
  colnames(sig) <- colnames(mat)
  return(sig)
}


##' @rdname make_sigma
##' @export
make_sigma_mat_dist_graph <- function(graph, cor = 0.8, absolute = FALSE){
  mat <- make_distance_graph(graph, absolute = absolute)
  sig <- make_sigma_mat_dist_adjmat(mat, cor)
  return(sig)
}
