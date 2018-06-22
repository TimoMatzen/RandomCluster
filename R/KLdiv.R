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
  # With q = 2 the Curie Weiss partition function can be used
  #Zrc <- ZC(n.nodes, p.mod)
  #mat <- matrix(0, n.nodes, n.nodes)
  # Loop over every configuration
  #div <- 0
  #for (i in 1:nrow(x)) {
   # Create adjacency matrix
   #mat[upper.tri(mat)] <- x[i,]
   #mat <- SymMat(mat)
    
   # Make graph object from adjacency matrix
   #graph <- graph_from_adjacency_matrix(mat, mode = "undirected")
    
  # Calculate the number of components in the graph
   #k <- components(graph)$no
   # Calculate Erdos-Renyi prob; q = 1
   #Pr <- RCProb(x[i,], p.true, 1,  n.nodes)
   #Prc <- prod(apply(x[i,,drop = FALSE],2,b, p = p.mod), 2^k)
   #Prc <- RCProb(x[i, ], p.mod, q = 2, n.nodes, Curie = TRUE)
    # Calculate KL divergence
   #div <- div +  ((log(Pr) + log(Zrc) - 
                    # log(prod(apply(x[i,,drop = FALSE],2,b, p = p.mod))) + log(2^k))  * Pr)
   #div <- div + (log(Pr) + log(Zrc) - log())
   #div <- div + ( (log(Pr) + log(Zrc) - log(prod(apply(x[i,,drop = FALSE],2,b, p = p.mod))) + log(2^k)) * Pr)
     
   #}
  
  # Multiply by log(Zrc)-log(2) to get KLdivergence
   #div <- (log(Zrc)) - log(2) * div
  
  # Calculte the KL divergence
  div <- 0
  for (i in 1:nrow(x)) {
    Per <- RCProb(x[i,], p.true, q = 1, n.nodes)
    Prc <- RCProb(x[i, ], p.mod, q = 2, n.nodes, Curie = TRUE)
    div <- div + ( Per * log(Per/Prc))
   
  }
  
  # Return the kl divergence between the models
  return(div)
}