#' The graphsim package
#' 
#' graphsim is a package to simulate normalised expression data from networks 
#' for biological pathways using \sQuote{\code{igraph}} objects and multivariate
#' normal distributions. 
#' 
#' @rdname graphsim-package
#' @name graphsim-package
#' @aliases graphsim-package graphsim
#' @docType package
#'
#' @section Introduction:
#' This package enables the generation of simulated gene expression datasets 
#' containing pathway relationships from a known underlying network.
#' These simulated datasets can be used to evaluate various bioinformatics 
#' methodologies, including statistical and network inference procedures.
#' 
#' These are computed by 1) resolving inhibitory states to derive a consistent
#' matrix of positive and negative edges, 2) inferring relationships between
#' nodes from paths in the graph, 3) weighting these in a Sigma (\eqn{\Sigma}) 
#' covariance matrix and 4) using this to sample a multivariate normal 
#' distribution.
#' 
NULL