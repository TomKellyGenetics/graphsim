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
##' @importFrom igraph graph.edgelist get.edge.attribute get.edgelist
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
  graph <- as.undirected(graph, mode = "collapse")
  E(graph)$state <-  state
  # remove duplicate edges (and corresponding states)
  graph <- as.directed(graph, mode = "arbitrary")
  state <- get.edge.attribute(graph, "state")
  if(length(state) == 1) state <- rep(state, length(E(graph)))
  state[state == 2] <- -1
  state[state == 1 | state == 0] <- 1
  if(is.numeric(state)) state[!(state %in% -1:2)] <- sign(state[!(state %in% -1:2)])
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
  #define starting states
  state_mat <- matrix(1, length(V(graph)), length(V(graph)))
  rownames(state_mat) <- colnames(state_mat) <- names(V(graph))
  # check for inhibiting edges
  if(all(state == 1)){
    #return without resolving for activations
    return(state_mat)
  } else if (length(V(graph)) == 1){
    state_mat <- 1
    return(state_mat)
  } else if (length(V(graph)) == 2){
    if(length(state) == 1){
      state_mat[row(state_mat) != col(state_mat)] <- state # one edge
      return(state_mat)
    } else{
      warning("only one edge state expected")
    }
  } else {
    #resolve inhibitions downstream
    # define shortest paths to all possible nodes
    #check if connected or use connected subgraphs
    compute_paths <- function(graph){
      if(is.connected(graph)){
        #define paths from 1st node
        paths <- shortest_paths(as.undirected(graph), V(graph)$name[1])$vpath
      } else {
        subgraphs <- decompose(graph)
        nodes <- sapply(subgraphs, function(subgraph) V(subgraph)$name[1])
        paths <- as.list(rep(NA, length(V(graph))))
        jj <- 0
        for(ii in 1:length(nodes)){
          subpaths <- shortest_paths(as.undirected(subgraphs[[ii]]), V(subgraphs[[ii]])$name[1])$vpath
          paths[jj+1:length(subpaths)] <- subpaths
          jj <- jj + length(subpaths)
        }
      }
      paths
    }
    paths <- compute_paths(graph)
    #check for cycles
    is.cyclic <- function(paths){
      cylic_paths <- sapply(paths, function(path){
        if(length(path) <= 1){
          FALSE
        } else if (path[1] == path[length(path)]){
          TRUE
        } else {
          FALSE
        }
      })
      any(cylic_paths)
    }
    if(is.cyclic(paths)){
      warning("Graph contains a cycle, computing minimal spanning tree. This may result in unresolved inhibitions.")
      tree <- minimum.spanning.tree(graph)
      paths <- compute_paths(tree)
      if(is.cyclic(paths)){
        warning("Graph contains cycles and cannot compute a minimal spanning tree. This will result in unresolved inhibitions.")
        stop()
      }
    }
    #sort paths into longest to shortest
    paths <- paths[order(sapply(paths, length), decreasing = TRUE)]
    #remove subpaths
    for(ii in 2:length(paths)){
      remove <- FALSE
      for(jj in 1:(ii-1)){
        # test if all nodes in path in a larger path above
        if(all(names(paths[[ii]]) %in% names(paths[[jj]]))){
          remove <- TRUE
        }
      }
      if(remove){
        paths[[ii]] <- NA
      }
    }
    na.omit.list <- function(y) { return(y[!sapply(y, function(x) all(is.na(x)))]) }
    paths <- na.omit.list(paths)
    # check inhibition state through graph paths
    edges <- as.matrix(get.edgelist(graph))
    for(ii in 1:length(paths)){
      if(length(paths[[ii]]) > 1){
        state_path <- c(1, rep(NA, length(paths[[ii]])-1))
        for(jj in 2:length(paths[[ii]])){
          #find edges
          kk <- which(
          edges[,1] == names(paths[[ii]])[jj-1] & edges[,2] == names(paths[[ii]])[jj] |
          edges[,2] == names(paths[[ii]])[jj-1] & edges[,1] == names(paths[[ii]])[jj]
          )
          #check for multiple edges
          if(length(kk) > 1){
            #check for same edges
            if(all(apply(edges[kk,], 2, function(x) all(x == x[1])))){
              #check for same state
              if(all(state[kk] == state[kk][1])){
                state[kk] <- state[kk][1]
                kk <- kk[1]
              } else {
                Mode <- function(x) unique(x)[which.max(tabulate(match(x, unique(x))))][1]
                state[kk] <- Mode(state[kk])
                kk <- kk[1]
              }
            }
          }
          #state for edges
          if(length(state_path[jj]) == length(state[kk])){
            state_path[jj] <- state[kk]
          } else {
            print(ii)
            print(jj)
            stop("State not compatible with graph paths")
          }
        }
        #cumulative product changes state along a path (*-1)
        state_path <- cumprod(state_path)
        #skip if activating path
        if(any(state_path != 1)){
          #cross-product produces a matrix form
          state_cross <- state_path %*% t(state_path)
          #add to state_matrix (all adjusted to start of path)
          state_mat[names(paths[[ii]]), names(paths[[ii]])] <- state_cross
        }
      }
    }
    #check for negative and positive links not in paths (against the 1st node)
    neg_nodes <- colnames(state_mat)[which(state_mat[1,] == -1)]
    pos_nodes <- colnames(state_mat)[2:ncol(state_mat)][which(state_mat[1,2:ncol(state_mat)] == 1)] # can assume at least 3 nodes due to checks above
    #check for if not computed in shorted paths already
    ##paths <- compute_paths(tree) # not minimal spanning tree
    for(aa in neg_nodes){
      for(bb in pos_nodes){
        aa_in_path <- sapply(paths, function(path) aa %in% names(path))
        bb_in_path <- sapply(paths, function(path) bb %in% names(path))
        #check if in ANY same path
        if(!any(aa_in_path & bb_in_path)){
          #if not invert
          state_mat[aa, bb] <- state_mat[bb, aa] <- -1
        }
      }
    }
    # check if same sign (to do: allow disconnected)
    for(aa in pos_nodes){
      for(bb in pos_nodes){
        aa_in_path <- sapply(paths, function(path) aa %in% names(path))
        bb_in_path <- sapply(paths, function(path) bb %in% names(path))
        #check if in ANY same path
        if(!any(aa_in_path & bb_in_path)){
          #if not invert
          state_mat[aa, bb] <- state_mat[bb, aa] <- 1
        }
      }
    }
    #ensure symmetric matrix
    state_mat[t(state_mat) == -1] <- -1
    return(state_mat)
  }
}
