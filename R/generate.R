##' @name generate_expression
##' @rdname generate_expression
##'
##' @title Generate Simulated Expression
##'
##' @description Compute simulated continuous expression data from a graph network structure. Requires an \code{\link[igraph]{igraph}} pathway structure and a matrix of states (1 for activating and -1 for inhibiting) for link signed correlations, from a vector of edge states to a signed adjacency matrix for use in \code{\link[graphsim]{generate_expression}}. Uses graph structure to pass a sigma covariance matrix from \code{\link[graphsim]{make_sigma_mat_dist_graph}} or \code{\link[graphsim]{make_sigma_mat_graph}} on to \code{\link[mvtnorm]{rmvnorm}}.
##'
##' @param n number of observations (simulated samples).
##' @param graph An \code{\link[igraph]{igraph}} object. May must be directed if states are used.
##' @param state numeric vector. Vector of length E(graph). Sign used to calculate state matrix, may be an integer state or inferred directly from expected correlations for each edge. May be applied a scalar across all edges or as a vector for each edge respectively. May also be entered as text for "activating" or "inhibiting" or as integers for activating (0,1) or inhibiting (-1,2). Compatible with inputs for \code{\link[graphsim]{plot_directed}}. Also takes a pre-computed state matrix from \code{\link[graphsim]{make_state_matrix}} if applied to the same graph multiple times.
##' @param cor numeric. Simulated maximum correlation/covariance of two adjacent nodes. Default to 0.8.
##' @param mean mean value of each simulated gene. Defaults to 0. May be entered as a scalar applying to all genes or a vector with a separate value for each.
##' @param dist logical. Whether a graph distance (\code{\link[graphsim]{make_sigma_mat_dist_graph}}) or derived matrix (\code{\link[graphsim]{make_sigma_mat_graph}}) is used to compute the sigma matrix.
##' @param comm,absolute logical. Parameters for Sigma matrix generation. Passed on to \code{\link[graphsim]{make_sigma_mat_dist_graph}} or \code{\link[graphsim]{make_sigma_mat_graph}}.
##' @keywords graph network igraph mvtnorm simulation
##' @import igraph mvtnorm
##' @importFrom igraph is_igraph
##' @importFrom Matrix nearPD
##' @importFrom matrixcalc is.symmetric.matrix is.positive.definite
##' @examples
##' 
##' library("igraph")
##' graph_test_edges <- rbind(c("A", "B"), c("B", "C"), c("B", "D"))
##' graph_test <- graph.edgelist(graph_test_edges, directed = TRUE)
##' adjacency_matrix <- make_adjmatrix_graph(graph_test)
##' n <- 100
##' generate_expression(n, graph_test, cor = 0.8)
##' 
##' @return numeric matrix of simulated data (log-normalised counts)
##' 
##' @export
generate_expression <- function(n, graph, state = NULL, cor = 0.8, mean = 0, comm = FALSE, dist = FALSE, absolute = FALSE){
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
    sig <- make_sigma_mat_graph(graph, cor, comm = comm)
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
  rownames(expr_mat) <-names(V(graph))
  colnames(expr_mat) <-paste0("sample_", 1:n)
  return(expr_mat)
}
