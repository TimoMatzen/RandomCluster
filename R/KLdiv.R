#' Kullback Leibler divergence 
#' Calculates the KL divergence between the random cluster model and the erdos-renyi model
#
#' @param x Matrix with sampled configurations
#' @param n.nodes The number of nodes in the model
#' @param p.true The p vthat was sampled with in the erdos renyi model
#' @param p.mod The p value of the model
#' 
#' @return The kl divergence between the Erdos-Renyi model and the Random Cluster model
#' 
#' @export



KLdivergence <- function (x, n.nodes, p.true, p.mod) {
  
  # First calculate the Random Cluster partition function
 
  n.edges <- ed(n.nodes)
  # Calculate the KL divergence if p.mod!=p.true
  if (p.true != p.mod) {
  # Calculte the KL divergence
    div <- 0
    for (i in 1:nrow(x)) {
      Per <- RCProb(x[i,], p.true, q = 1, n.nodes)
      Prc <- RCProb(x[i, ], p.mod, q = 2, n.nodes, Curie = TRUE)
      div <- div + ( Per * log(Per/Prc))
     
    }
  } else {
    # Get partition with curie-weiss approximation
    Z <- ZC(n.nodes, p.true)
    
    # Set div to 0
    div <- 0
      for (i in 1:nrow(x)) {
      
     
      # Calculate components
      mat <- matrix(NA,n.nodes, n.nodes)
      mat[upper.tri(mat)] <- x[i,]
      mat <- SymMat(mat)
      net <- graph_from_adjacency_matrix(mat, mode = "undirected")
      k <- components(net)$no
      
      # Calculate ER prob
      Per <- (RCProb(x[i, ], p.true, q = 1, n.nodes)) # Normalize the probability value
      
      # Calculate KL divergence with rewritten equation when p.true == p.mod
      div <- div +  (k  * Per)
      
      
      }
    # Multiply div by log(z) - log(2)
    div <- log(Z) - (log(2) * div)
  }
  
  # Return the kl divergence between the models
  return(div)
}