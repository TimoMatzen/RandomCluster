#' Gibbs Sampler for the Random Cluster model
#' 
#' Function to sample configurations from the Random cluster model using a Gibbs sampler.
#' 
#' @param t The number of iterations to run the algorithm
#' @param p The wiring probability
#' @param q The clustering coefficient of the Random Cluster model
#' @param n.nodes The number of nodes in the configuration
#' @param burnin Burnin period for the markov chain
#' 
#' @return The sampled configurations; minus the burnin
#' 
#' @export
#' 
GibbsRC <- function(t,p,q,n.nodes, burnin = 50) {
  n.edges <- ed(n.nodes)
  #Initialise the matrix to save all the samples
  mat1 <- matrix(NA, nrow = t, ncol = n.edges)
  
  # Initialiseer de eerste sample, starten bij alles op 0
  mat <- matrix(NA,n.nodes, n.nodes)
  mat[upper.tri(mat)] <- seq(1:n.edges)
  edge.list<- which(mat > 0, arr.ind = T )
  mat1[1,] <- rep(0,n.edges)
  
  for(i in 1:(t-1)){
    # Sample an edge of the configuration
    j <- sample(seq(1,n.edges,by = 1),1, prob = rep((1/n.edges),n.edges))
    
    # Sla de vector op
    ome <- mat1[i,]
    if(ome[j]==0){
      ome1 <- ome
      ome1[j] <- 1
    }else{
      ome1 <- ome
      ome1[j] <- 0
    }
    # Estimate network if the edge is 0
    net0 <- graph_from_data_frame(d = matrix(edge.list[ome != 0,],ncol = 2),
                                  vertices = seq(1:n.nodes), directed = F)
    k0 <- components(net0)$no
    
    # Estimate the network if the edge is 1
    net1 <- graph_from_data_frame(d = matrix(edge.list[ome1 != 0,],ncol = 2),vertices = seq(1:n.nodes),
                                  directed = F)
    k1 <- components(net1)$no
    
    # Estimate probabilities
    p0 <- (prod(apply(matrix(ome,nrow = 1),2,b,p = p))*q^k0)
    p1 <- (prod(apply(matrix(ome1,nrow = 1),2,b,p = p))*q^k1)
    
    # Sample a random uniform number
    U <- runif(1)
    
    # Calculate change statistic
    if(ome[j] == 0){
      P <- p0/(p0+p1)
      if(U > P){
        ome[j] <- 1
      }
    }else{
      P <- p1/(p0+p1)
      if(U <= P){
        ome[j] <- 0
      }
    }
    mat1[i+1,] <- ome 
  }
  
  # Return the samples minus the burnin
  return(mat1[-(1:burnin),])
}
