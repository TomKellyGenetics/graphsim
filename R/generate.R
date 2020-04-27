##' @name generate_expression
##' @rdname generate_expression
##'
##' @title Generate Simulated Expression
##'
##' @description Compute simulated continuous expression data from a graph network structure. Requires an \code{\link[igraph]{igraph}} pathway structure and a matrix of states (1 for activating and -1 for inhibiting) for link signed correlations, from a vector of edge states to a signed adjacency matrix for use in \code{\link[graphsim]{generate_expression}}. Uses graph structure to pass a sigma covariance matrix from \code{\link[graphsim]{make_sigma_mat_dist_graph}} or \code{\link[graphsim]{make_sigma_mat_graph}} on to \code{\link[mvtnorm]{rmvnorm}}.
##'
##' @param n number of observations (simulated samples).
##' @param mat precomputed adjacency, laplacian, commonlink, or scaled distance matrix.
##' @param graph An \code{\link[igraph]{igraph}} object. May must be directed if states are used.
##' @param state numeric vector. Vector of length E(graph). Sign used to calculate state matrix, may be an integer state or inferred directly from expected correlations for each edge. May be applied a scalar across all edges or as a vector for each edge respectively. May also be entered as text for "activating" or "inhibiting" or as integers for activating (0,1) or inhibiting (-1,2). Compatible with inputs for \code{\link[graphsim]{plot_directed}}. Also takes a pre-computed state matrix from \code{\link[graphsim]{make_state_matrix}} if applied to the same graph multiple times.
##' @param cor numeric. Simulated maximum correlation/covariance of two adjacent nodes. Default to 0.8.
##' @param mean mean value of each simulated gene. Defaults to 0. May be entered as a scalar applying to all genes or a vector with a separate value for each.
##' @param dist logical. Whether a graph distance (\code{\link[graphsim]{make_sigma_mat_dist_graph}}) or derived matrix (\code{\link[graphsim]{make_sigma_mat_graph}}) is used to compute the sigma matrix.
##' @param comm,absolute,laplacian logical. Parameters for Sigma matrix generation. Passed on to \code{\link[graphsim]{make_sigma_mat_dist_graph}} or \code{\link[graphsim]{make_sigma_mat_graph}}.
##' @keywords graph network igraph mvtnorm simulation
##' @importFrom  mvtnorm rmvnorm
##' @importFrom igraph as_adjacency_matrix graph.edgelist
##' @import igraph
##' @importFrom igraph is_igraph
##' @importFrom Matrix nearPD
##' @importFrom matrixcalc is.symmetric.matrix is.positive.definite
##' @importFrom gplots heatmap.2
##' @examples
##' 
##' library("igraph")
##' graph_test_edges <- rbind(c("A", "B"), c("B", "C"), c("B", "D"))
##' graph_test <- graph.edgelist(graph_test_edges, directed = TRUE)
##' n <- 100
##' generate_expression(n, graph_test, cor = 0.8)
##' 
##' adjacency_matrix <- make_adjmatrix_graph(graph_test)
##' generate_expression_mat(n, adjacency_matrix, cor = 0.8)
##' 
##' @return numeric matrix of simulated data (log-normalised counts)
##' 
##' @export
generate_expression <- function(n, graph, state = NULL, cor = 0.8, mean = 0, comm = FALSE, dist = FALSE, absolute = FALSE, laplacian = FALSE){
  if(!is.integer(n)){
    if(is.numeric(n)){
      if(floor(n) == n){
        n <- as.integer(n)
      } else{
        n <- floor(n)
        print(paste("rounding to sample size", n))
      }
    } else{
      stop("sample size n must be an integer of length 1")
    }
  }
  if(length(n) > 1) stop("sample size n must be an integer of length 1")
  if(!is_igraph(graph)) stop("graph must be an igraph class")
  if(is.vector(state) || length(state) == 1) state <- make_state_matrix(graph, state)
  if(!(is.vector(mean)) || length(mean) == 1 ) mean <- rep(mean,length(V(graph)))
  if(dist){
    sig <- make_sigma_mat_dist_graph(graph, cor, absolute = absolute)
  } else {
    sig <- make_sigma_mat_graph(graph, cor, comm = comm, laplacian = laplacian)
  }
  ## migrate state to calling sigma ##
  if(!(is.null(state))) sig <- state * sig
  if(is.symmetric.matrix(sig) == FALSE) {
    warning("sigma matrix was not positive definite, nearest approximation used.")
    sig <- as.matrix(nearPD(sig, corr=T, keepDiag = TRUE)$mat) #postive definite correction
  }
  if(is.positive.definite(sig) == FALSE) {
    warning("sigma matrix was not positive definite, nearest approximation used.")
    sig <- as.matrix(nearPD(sig, corr=T, keepDiag = TRUE)$mat) #postive definite correction
  }
  expr_mat <- t(rmvnorm(n,mean=mean, sigma=sig))
  rownames(expr_mat) <- names(V(graph))
  colnames(expr_mat) <- paste0("sample_", 1:n)
  return(expr_mat)
}


##' @rdname generate_expression
##' @export
generate_expression_mat <- function(n, mat, state = NULL, cor = 0.8, mean = 0, comm = FALSE, dist = FALSE, absolute = FALSE, laplacian = FALSE){
  if(!is.integer(n)){
    if(is.numeric(n)){
      if(floor(n) == n){
        n <- as.integer(n)
      } else{
        n <- floor(n)
        print(paste("rounding to sample size", n))
      }
    } else{
      stop("sample size n must be an integer of length 1")
    }
  }
  if(length(n) > 1) stop("sample size n must be an integer of length 1")
  if(!is.matrix(mat)) stop("graph must be an igraph class")
  if(is.vector(state) || length(state) == 1) state <- make_state_matrix(graph, state)
  if(!(is.vector(mean)) || length(mean) == 1 ) mean <- rep(mean,ncol(mat))
  if(dist){
    distmat <- make_distance_adjmat(mat)
    sig <- make_sigma_mat_dist_adjmat(distmat, cor, absolute = absolute)
  } else {
    if(comm && laplacian){
      warning("Error: only one of commonlink or laplacian can be used")
      stop()
    }
    if(!comm && !laplacian) sig <- make_sigma_mat_adjmat(mat, cor)
    if(comm){
      mat <- make_commonlink_adjmat(mat)
      sig <- make_sigma_mat_comm(mat, cor)
    }
    if(laplacian){
      mat <- make_laplacian_adjmat(mat)
      sig <- make_sigma_mat_laplacian(mat, cor)
    }
  }
  if(!(is.null(state))) sig <- state * sig
  if(is.symmetric.matrix(sig) == FALSE) {
    warning("sigma matrix was not positive definite, nearest approximation used.")
    sig <- as.matrix(nearPD(sig, corr=T, keepDiag = TRUE)$mat) #postive definite correction
  }
  if(is.positive.definite(sig) == FALSE) {
    warning("sigma matrix was not positive definite, nearest approximation used.")
    sig <- as.matrix(nearPD(sig, corr=T, keepDiag = TRUE)$mat) #postive definite correction
  }
  expr_mat <- t(rmvnorm(n,mean=mean, sigma=sig))
  rownames(expr_mat) <- names(colnames(mat))
  colnames(expr_mat) <- paste0("sample_", 1:n)
  return(expr_mat)
}
