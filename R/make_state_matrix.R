##' @name Make State Matrix
##' @rdname make_state_matrix
##'
##' @title Make State Matrix
##'
##' @description Functions to compute the matrix of states (1 for activating and -1 for inhibiting) for link signed correlations, from a vector of edge states to a signed adjacency matrix for use in \code{\link[graphsim]{generate_expression}}.
##'
##' @param graph An \code{\link[igraph]{igraph}} object. May be directed or weighted as long as a shortest path can be computed.
##' @param state numeric vector. Vector of length E(graph). Sign used to calculate state matrix, may be an integer state or inferrfed directly from expect correlations for each edge.
##' @keywords graph network igraph mvtnorm simulation
##' @import igraph
##' @export
make_state_matrix <- function(graph, state){
  state <- sign(state) # coerce to vector or 1 and -1 if not already
  edges <- as.matrix(get.edgelist(graph)[grep(-1, state),])
  if(length(grep(-1, state))==1) edges <- t(edges)
  state_mat <- matrix(1, length(V(graph)), length(V(graph)))
  rows <- apply(edges, 1, function(x) grep(unlist(as.list(x)[1]), names(V(graph))))
  cols <- apply(edges, 1, function(x) grep(unlist(as.list(x)[2]), names(V(graph))))
  for(ii in 1:length(rows)){
    state_mat[rows[ii], cols[ii]] <- state_mat[rows[ii], cols[ii]] * -1
    sub_edge <- get.edgelist(graph)[setdiff(c(1:length(E(graph))),intersect(grep(names(V(graph))[rows[ii]], as.matrix(get.edgelist(graph))[,1]),  grep(names(V(graph))[cols[ii]], as.matrix(get.edgelist(graph))[,2]))),]
    sub_graph <- graph.edgelist(sub_edge)
    clust <- clusters(sub_graph)
    down_ind <- grep(clust$membership[match(names(V(graph))[cols[ii]], names(clust$membership))], clust$membership) #find downstream genes
    if(all(is.na(down_ind))==F){ #if downstream genes
      down_genes <- names(clust$membership[down_ind])
      ind <- match(down_genes, names(V(graph)))
      state_mat[ind, ind] <- state_mat[ind, ind] * -1
    }
  }
  return(state_mat)
}
