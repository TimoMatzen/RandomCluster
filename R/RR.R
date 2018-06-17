#' Randomness Recycler 
#' 
#' Randomness recycler method for sampling from the Random Cluster model. Introduced by
#' Fill and Huber (2000).
#' 
#' @param n.nodes The number of nodes in the network
#' @param p The wiring probability
#' @param q The clustering coefficient of the Random Cluster model
#' 
#' @return A adjacency matrix of a graph sampled from the Random Cluster model.
#' 
#' @export
#' 


RandomnessRecycler <- function(n.nodes, p, q) {
  
  i <- 0
  
  # First create the empty set
  x.mat <- matrix(0, n.nodes, n.nodes)
  
  # Than create Et
  Et <- matrix(NA,nrow = 1, ncol = 2)
  
  # Create E
  E <- which(upper.tri(x.mat) == TRUE, arr.ind=T)


  # Start while loop
  #!all(apply(E, 1, paste, collapse = "") %in% 
   #      apply(Et, 1, paste, collapse = ""))
  while (!all(apply(E, 1, paste, collapse = "") %in% 
                    apply(Et, 1, paste, collapse = ""))) {
    i <- i + 1
    # Sample an edge
   
    e <- sample(which(!apply(E, 1, paste, collapse = "") %in% 
                        apply(Et, 1, paste, collapse = "")),1)
    
    # Make graph object
    graph <- graph_from_adjacency_matrix(x.mat, mode = "undirected")
    
    # Sample state
    state <- sample(c(0,1),1, prob = c((1-p),p))
    # Set state
    x.mat[E[e,1], E[e,2]] <- x.mat[E[e,2], E[e,1]] <- state
  
    # Get in which component the nodes are
    comp <- components(graph)$membership
    
    if (state == 0) {
      # Accept state
      Et <- rbind(Et, E[e, ])
    }
    
    if (state == 1 ) {
      if (comp[E[e,1]] == comp[E[e,2]]) {
        # If connected, accept the proposed edge
        Et <- rbind(Et, E[e, ])
      } else if (((1/q) >= runif(1))) {
        # Unconnected with prob 1/q Accept the edge
        Et <- rbind(Et, E[e, ])
      } else {
        # Reject; but we have information about the  clusters
        
        # Sample w
        w <- sample(c(E[e,1], E[e,2]),1 )
        # Get all vertices in same component W
        comp.w <- which(comp == comp[w], arr.ind = T)
        
        # Delete from Et all edges in or adjacent to W
        # Colour them zero in x
        # Set edge to zero
        x.mat[Et[apply(Et, 1, function(r) any(r %in% as.vector(comp.w))),, drop = FALSE]] <- 0
        x.mat <- SymMat(x.mat)
        #x.mat[Et[apply(Et, 1, function(r) any(r %in% comp.w)),, drop = FALSE][,c(2,1)]] <- 0
        # Set edge to zero
        #x.mat[E[e,1], E[e,2]] <- 0
        #x.mat[E[e,2], E[e,1]] <- 0
        
        # Delete the edges leading in and out of the component
        #Et <- Et[!apply(Et, 1, function(r) any(r %in% w)),, drop = FALSE]
        
        # Create a spanning tree of the nodes that are deleted
        
        # Make a edgelist
        #x.edge <- as_edgelist(graph)
        #a <- Et[apply(Et, 1, function(r) any(r %in% w)),, drop = FALSE]
        #a <-  Et[apply(Et, 1, function(r) any(r %in% comp.w)),, drop = FALSE]
        #a <- Et[apply(Et, 1, function(r) any(r %in% comp.w)),, drop = FALSE]
        #a <- x.edge[apply(x.edge, 1, function(r) any(r %in% comp.w)),, drop = FALSE]
        # Add edge (v,w) to the spanning tree
        #a <- rbind(a, E[e,])
        # Create graph
        #a <- graph_from_edgelist(a, directed = FALSE)
        # use mst to create uinique spanning tree of the edges
        #a <- mst(a)
        # Change back into edge list
        #a <- as_edgelist(a)
        
        # Graph Et
        #graphEt <- graph_from_edgelist(Et[-1, ], directed = FALSE)
        
        
        # Delete all edges adjacent or in C
        Et <- Et[!apply(Et, 1, function(r) any(r %in% w)),, drop = FALSE]
        Et <- Et[!apply(Et, 1, function(r) any(r %in% as.vector(comp.w))),, drop = FALSE]
        
        # Sample edges from the spanning tree
        #rho <- p/q
        #t <- sample(c(1,0), size = nrow(a), replace = TRUE, 
        #         prob = c(rho/(1-p+rho),(1-p)/(1-p+rho)))
        
        # Set edge states
        #x.mat[a] <- t
        #x.mat[a[,c(2,1)]] <- t
        
        # Add the spanning tree back into Et
        #Et <- rbind(Et, a)
        # Set edge to zero
        x.mat[E[e,1], E[e,2]] <- 0
        x.mat[E[e,2], E[e,1]] <- 0
        
      }
    }
    
    if(!isSymmetric(x.mat)) {
      print("not symmetric at")
      print(i)
    }
  }
  
  return(x.mat)
}


  