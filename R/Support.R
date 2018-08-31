# Calculate n.edges from n.nodes
#
# x: n.nodes
# 
# @keyword internal
#' @export
ed <- function(x) {
  # Takes as input the number of nodes
  # Outputs the number of edges
  y <- x*((x-1)/ 2)
  return(y)
}

# Transition function (CFTP)
#
# n.nodes: number of nodes
# s: current state
# U: random uniform number on (0,1)
# r: random edge to switch state
# p: probability for an edge
# q: cluster coefficient
#
# @keyword internal
transition <- function(n.nodes, s, U, r ,p, q) {
  if (s[r] == 1) {
    s2 <- s
    # Change random edge to 0
    s2[r] <- 0
    # Apply transition function
    out <- RC(n.nodes, s2, p, q)/(RC(n.nodes, s, p, q) + RC(n.nodes, s2, p, q))
    
    # Check whether to accept change
    s[r] <- 1 * !(U <= out)
    
  } else {
    s2 <- s
    # Change random edge to 1
    s2[r] <- 1
    # Apply transition function
    out <- RC(n.nodes, s, p, q)/(RC(n.nodes, s, p, q) + RC(n.nodes, s2, p, q))
    # Check whether to accept change
    s[r] <- 1 * (U > out)
  }
  
  # Return the adjusted state
  return(s)
}

# Curie Weiss partition function
# 
# Can be used as partition function for Random Cluster model with q is 2 as it you can
# proof that these are the same.
# 
# @param n.nodes The number of nodes in the network
# @param p The probability for an edge to be present
# 
# @return Returns the normalizing constant for the Random Cluster model with q = 2.
# 
# @keyword internal
#' @export
ZC <- function(n.nodes, p) {
  # Eerst moet sigma berekend worden, afgeleid van Maartens document
  sig <- (log(1-p))/-2
  # Z eerst gelijk aan 0
  Z <- 0
  # Vervolgens loopen over alle nodes.
  for (i in 0:n.nodes){
    Z <- Z + (choose(n.nodes, i) * exp((.5*sig) * (2*i - n.nodes)^2))
  }
  # Om de partitie functie van curie-weiss gelijk te krijgen aan die van het RC nog
  # vermenigvuldigen met een term die afhangt van sigma en het aantal nodes.
  Z <- Z * exp((-.5*sig*n.nodes^2))
  return(Z)
}
# Product function numerator Random Cluster function
#
# x: n.edges
# p: edge probability
#
# @keyword internal
b <- function(x,p){
  p^(x)*(1-p)^(1-x)
}

# Function for copying upper tri matrix to lower tri
#
# m: Matrix
#
# return: Symmetric matrix
#' @export
SymMat <- function(m) {
  m[lower.tri(m)] <- t(m)[lower.tri(m)]
  m
}

# Function to just calculate the RC function
#
# n.nodes: Number of nodes
# s: State of system
# p: Probability for forming an edge
# q: Clustering measure from the Random Cluster model
#
# @keyword internal
#' @export 

RC <- function(n.nodes, s, p, q) {
  # Create empty adjacency matrix
  mat <- matrix(0,n.nodes, n.nodes)
  # Fill upper tri with the state
  mat[upper.tri(mat)] <- s
  # Make matrix symmetric
  mat <- SymMat(mat)
  # Create the edge list
  net <- graph_from_adjacency_matrix(mat, mode = "undirected")
  # Get number of components
  k <- components(net)$no
  # Calculate RC measure
  out <- p^sum(s) * (1-p)^(length(s)-sum(s)) * q^k
  return(out)
}