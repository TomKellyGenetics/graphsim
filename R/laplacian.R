##' @name make_laplacian
##' @aliases make_laplacian_graph
##' @rdname make_laplacian
##'
##' @title Generate Laplacian Matrix
##' 
##' @description Compute the Laplacian matrix of a (directed) \code{\link[igraph:aaa-igraph-package]{igraph}}
##' structure, preserving node/column/row names (and direction).
##'
##' @param mat precomputed adjacency matrix.
##' @param graph An \code{\link[igraph:aaa-igraph-package]{igraph}} object. May be directed or weighted.
##' @param directed logical. Whether directed information is passed to the Laplacian matrix.
##' @keywords graph network igraph Laplacian
##' @importFrom igraph laplacian_matrix graph.edgelist as.undirected graph_from_adjacency_matrix
##' @import igraph
##' @importFrom Matrix Matrix
##' @importClassesFrom Matrix dgCMatrix
##' @examples 
##' 
##' # construct a synthetic graph module
##' library("igraph") 
##' graph_test_edges <- rbind(c("A", "B"), c("B", "C"), c("B", "D"))
##' graph_test <- graph.edgelist(graph_test_edges, directed = TRUE)
##' # compute Laplacian matrix for toy example
##' laplacian_matrix <- make_laplacian_graph(graph_test)
##' laplacian_matrix
##' 
##' # compute Laplacian matrix from adjacency matrix
##' adjacency_matrix <- make_adjmatrix_graph(graph_test)
##' laplacian_matrix <- make_laplacian_adjmat(adjacency_matrix)
##' laplacian_matrix
##' 
##' # construct a synthetic graph network
##' graph_structure_edges <- rbind(c("A", "C"), c("B", "C"), c("C", "D"), c("D", "E"),
##'                                c("D", "F"), c("F", "G"), c("F", "I"), c("H", "I"))
##' graph_structure <- graph.edgelist(graph_structure_edges, directed = TRUE)
##' # compute Laplacian matrix for toy network
##' graph_structure_laplacian_matrix <- make_laplacian_graph(graph_structure)
##' graph_structure_laplacian_matrix
##'  
##' # import graph from package for reactome pathway
##' # TGF-\eqn{\Beta} receptor signaling activates SMADs (R-HSA-2173789)
##' TGFBeta_Smad_graph <- identity(TGFBeta_Smad_graph)
##' 
##' # compute Laplacian matrix for TGF-\eqn{\Beta} receptor signaling activates SMADs
##' TGFBeta_Smad_laplacian_matrix <- make_laplacian_graph(TGFBeta_Smad_graph)
##' dim(TGFBeta_Smad_laplacian_matrix)
##' TGFBeta_Smad_laplacian_matrix[1:12, 1:12]
##' # visualise matrix
##' library("gplots")
##' heatmap.2(TGFBeta_Smad_laplacian_matrix, scale = "none", trace = "none",
##'           col = colorpanel(50, "blue", "white", "red"))
##' 
##' @return An Laplacian matrix compatible with generating an expression matrix
##' 
##' @export
make_laplacian_adjmat <- function(mat, directed = FALSE){
  graph <- graph_from_adjacency_matrix(mat, weighted = TRUE, mode = ifelse(directed, "directed", "undirected"))
  laplacian <- as.matrix(laplacian_matrix(graph))
  rownames(laplacian) <- colnames(laplacian) <- names(V(graph))
  return(laplacian)
}

##' @rdname make_laplacian
##' @export
make_laplacian_graph <- function(graph, directed = FALSE){
  if(directed == FALSE) graph <- as.undirected(graph)
  laplacian <- as.matrix(laplacian_matrix(graph))
  rownames(laplacian) <- colnames(laplacian) <- names(V(graph))
  return(laplacian)
}
