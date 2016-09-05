##' @rdname make_sigma
##'
##' @title Generate Sigma Matrix
##'
##' @description Compute the Sigma matrix from an \code{\link[igraph]{igraph}} structure or pre-computed matrix. These are compatible with \code{\link[mvtnorm]{rmvnorm}} and \code{\link[graphsim]{generate_expression}}.
##'
##' @param mat precomputed adjacency, commonlink, or scaled distance matrix.
##' @param graph An \code{\link[igraph]{igraph}} object. May be directed or weighted.
##' @param cor numeric. Simulated maximum correlation/covariance of two adjacent nodes. Default to 0.8.
##' @param directed logical. Whether directed information is passed to the distance matrix.
##' @param comm logical whether a common link matrix is used to compute sigma. Defaults to FALSE (adjacency matrix).
##' @param absolute logical. Whether distances are scaled as the absolute difference from the diameter (maximum possible). Defaults to TRUE. The alternative is to calculate a relative difference from the diameter for a geometric decay in distance.
##' @keywords graph network igraph mvtnorm
##' @import igraph mvtnorm
##' @export
make_sigma_mat.matrix <- function(mat, cor){
  sig <- ifelse(mat>0, cor, 0)
  diag(sig)<-1
  rownames(sig) <- rownames(mat)
  colnames(sig) <- colnames(mat)
  return(sig)
}

##' @rdname make_sigma
##' @export
make_sigma_mat.igraph <- function(graph, cor, comm = F, directed = F){
  mat <- make_adjmatrix.igraph(graph, directed = directed)
  if(comm) mat <- make_commonlink.matrix(mat)
  diag(mat) <- 1
  sig <- ifelse(mat>0, cor*mat/max(mat), 0)
  diag(sig)<-1
  rownames(sig) <- rownames(mat)
  colnames(sig) <- colnames(mat)
  return(sig)
}

##' @rdname make_sigma
##' @export
make_sigma_mat_dist.matrix <- function(mat, cor, absolute = F){
  sig <- mat/max(mat[mat!=1]) * cor
  sig <- ifelse(sig>0, sig, 0)
  diag(sig)<-1
  rownames(sig) <- rownames(mat)
  colnames(sig) <- colnames(mat)
  return(sig)
}

##' @rdname make_sigma
##' @export
make_sigma_mat_dist.igraph <- function(graph, cor, absolute = F){
  mat <- make_distance.igraph(graph, absolute = absolute)
  sig <- make_sigma_mat_dist.matrix(mat, cor)
  return(sig)
}
