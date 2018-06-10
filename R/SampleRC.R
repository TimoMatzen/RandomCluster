# TODO: Add possibility to also sample with Randomness Recycler

#' Sample Random Cluster
#' 
#' Function to sample from the Random Cluster model. Currently this function works with 
#' the coupling from the past algorithm that is used in Grimmit (2006).
#' 
#' @param p The probability of forming an edge
#' @param q The cluster coefficient from the Random Cluster model
#' @param n.edges The number of edges in the model
#' @param n.nodes The number of nodes in the network
#' @param n.iter The maximum number of iterations to run CFTP
#' @param verbose Whether to print the progress.
#' 
#' @return Returns a matrix with samples from the random cluster model.
#' 
#' @export
#' 
SampleRC <- function(p, q, n.edges, n.nodes, n.iter, verbose = TRUE) {
  
  # Starting at a higher number of runs when network is large.
  if(n.nodes < 5 ) {
    M <- 2^(0:10)
  } else {
    M <- 2^(5:10)
  }
  
  # Matrix to save samples
  output <- matrix(1,ncol = n.edges, nrow = n.iter)
  for (iter in 1:n.iter) {
    # Set empty set
    s1 <- rep(0, n.edges)
    # Set complete set
    s2 <- rep(1, n.edges)
    # Create vector of edge samples to save the edge sample
    e.samps <- vector(mode = "numeric", length = 0)
    # Create a uniform samples vector
    U.samps <- vector(mode = "numeric", length = 0)
    for (n in 1:length(M)) {
      # Go back further in time
      m <- M[n]
      
      # Pick a random sample of the edges for each state
      if (n == 1) {
        # If at it 1 take only 1 sample
        e.samp <- sample(1:n.edges, M[n], replace = TRUE)
        # Pick a random uniform sample
        U.samp <- runif(M[n])
        
      } else {
        # If it != 1 take samples iter[p] - iter[p-1] 
        e.samp <- sample(1:n.edges, (M[n] - M[n-1]), replace = TRUE)
        # Pick a random uniform sample
        U.samp <- runif((M[n] - M[n-1]))
      }
      
      # Append e.samps with the new sample
      e.samps <- append(e.samps, e.samp)
      # Append U.samps with the new sample
      U.samps <- append(U.samps, U.samp)
      
      for (i in m:1) {
        # Make transitions
        s1 <- transition(n.nodes, s1, U.samps[i], e.samps[i], p, q)
        s2 <- transition(n.nodes, s2, U.samps[i], e.samps[i], p, q)
      }
      if (!all(s1 == s2)) {
        # Reset the vectors to the empty/full set
        s1 <- rep(0, n.edges)
        s2 <- rep(1, n.edges)
      } else {
        break
      }
      
    }
    # Save output
    output[iter,] <- s1
    # Track progress
    if (verbose == TRUE) {
      print(iter)
    }
  }
  return(output)
}