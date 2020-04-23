##' @name make_state_matrix
##' @rdname make_state_matrix
##'
##' @title Make State Matrix
##' 
##' @description Functions to compute the matrix of states (1 for activating and -1 for inhibiting) 
##' for link signed correlations, from a vector of edge states to a signed adjacency matrix for use
##' in \code{\link[graphsim]{generate_expression}}.
##' @param graph An \code{\link[igraph]{igraph}} object. May be directed or weighted as long as
##' a shortest path can be computed.
##' @param state numeric vector. Vector of length E(graph). Sign used to calculate state matrix,
##' may be an integer state or inferred directly from expected correlations for each edge. May be
##' applied a scalar across all edges or as a vector for each edge respectively. May also be
##' entered as text for "activating" or "inhibiting" or as integers for activating (0,1) or
##' inhibiting (-1,2). Compatible with inputs for \code{\link[graphsim]{plot_directed}}. Vector 
##' input is supported either directly calling the function with a value for each edge in 
##'  \code{E(graph)} or as an edge "attribute" in the igraph object (using 
##'  \code{E(g)$state <- states}).
##' @keywords graph network igraph mvtnorm simulation
##' @import igraph
##' @examples 
##' 
##' library("igraph")
##' graph_test_edges <- rbind(c("A", "B"), c("B", "C"), c("B", "D"))
##' graph_test <- graph.edgelist(graph_test_edges, directed = TRUE)
##' state_matrix <- make_state_matrix(graph_test)
##' 
##' @return An integer matrix indicating the resolved state
##' (activating or inhibiting for each edge or path between nodes)
##' @export
make_state_matrix <- function(graph, state = NULL){
  if(!is.null(get.edge.attribute(graph, "state"))){
    state <- get.edge.attribute(graph, "state")
  } else {
    # add default state if not specified
    if(is.null(state)){
      state <- "activating"
    }
  }
  if(length(state) == 1) state <- rep(state, length(E(graph)))
  state[state == 2] <- -1
  state[state == 1] <- 1
  if(is.character(state)){
    state <- as.list(state)
    state[grep("activating", state)] <- 1
    state[grep("activation", state)] <- 1
    state[grep("activate", state)] <- 1
    state[grep("active", state)] <- 1
    state[grep("positive", state)] <- 1
    state[grep("inhibiting", state)] <- -1
    state[grep("inhibition", state)] <- -1
    state[grep("inhibitory", state)] <- -1
    state[grep("inhibit", state)] <- -1
    state[grep("negative", state)] <- -1
    if(is.character(state)){
      warning("Please give state as a scalar or vector of length(E(graph)): input must be 'activating', 'inhibiting' or an integer")
    }
    state <- unlist(as.numeric(state))
  }
  if(!all(state %in% -1:2)){
    state <- sign(state) # coerce to vector or 1 and -1 if not already
    warning("State inferred from non-integer weighted edges: Please give numeric states as integers: 0 or 1 for activating, -1 or 2 for inhibiting")
  }
  # define shortest paths to all possible nodes
  #check if connected or use connected subgraphs
  if(is.connected(graph)){
    #define paths from 1st node
    paths <- shortest_paths(graph, V(graph)$name[1])$vpath
  } else {
    subgraphs <- decompose(graph)
    nodes <- sapply(subgraphs, function(subgraph) V(subgraph)$name[1])
    paths <- list()
    jj <- 1
    for(ii in 1:length(nodes)){
      subpaths <- shortest_paths(subgraphs[[ii]], V(subgraphs[[ii]])$name[1])$vpath
      paths[[jj:(jj+length(subpaths))]] <- subpaths
      jj <- jj + length(subpaths)
    }
    paths <- lapply(decompose(graph), function(subgraph){
      shortest_paths(subgraph, V(subgraph)$name[1])$vpath}
    )
  }
  
  edges <- as.matrix(get.edgelist(graph)[grep(-1, state),])
  if(length(grep(-1, state))==1) edges <- t(edges)
  state_mat <- matrix(1, length(V(graph)), length(V(graph)))
  if(length(edges) > 0){
    rows <- unlist(apply(edges, 1, function(x) grep(unlist(as.list(x)[1], use.names = FALSE), names(V(graph)))))
    cols <- unlist(apply(edges, 1, function(x) grep(unlist(as.list(x)[2], use.names = FALSE), names(V(graph)))))
    for(ii in 1:length(rows)){
      state_mat[rows[ii], cols[ii]] <- state_mat[rows[ii], cols[ii]] * -1
      sub_edge <- get.edgelist(graph)[setdiff(c(1:length(E(graph))),
                                              intersect(grep(names(V(graph))[rows[ii]],as.matrix(get.edgelist(graph))[,1]), 
                                                        grep(names(V(graph))[cols[ii]], as.matrix(get.edgelist(graph))[,2]))),]
      sub_graph <- graph.edgelist(sub_edge)
      clust <- clusters(sub_graph)
      down_ind <- grep(clust$membership[match(names(V(graph))[cols[ii]], names(clust$membership))], clust$membership) #find downstream genes
      if(all(is.na(down_ind))==F){ #if downstream genes
        down_genes <- names(clust$membership[down_ind])
        ind <- match(down_genes, names(V(graph)))
        state_mat[ind, ind] <- state_mat[ind, ind] * -1
      }
    }
  }
  #ensure symmetric matrix
  state_mat <- state_mat * t(state_mat)
  rownames(state_mat) <- colnames(state_mat) <- names(V(graph))
  return(state_mat)
}
