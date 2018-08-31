#' Random Cluster Probability
#' 
#' Calculates the probability of a edge configuration given parameters p and q.
#' 
#' @param conf A binary vector with the edge configuration
#' @param p Probability of forming an edge between two nodes
#' @param q Clustering coefficient RC model
#' @param n.nodes The number of nodes in the network
#' @param Curie Whether to use the Curie-Weiss partition function, is only possible when q = 2.
#' 
#' 
#' @return Returns the configuration probability.
#' 
#' @export

RCProb <- function(conf, p,q,n.nodes, Curie = F){
  if(Curie == T & q != 2){
    stop('Can not use Curie-Weiss partition function when q != 2')
  }
  # Calculate maximum number of edges.
  n.edges <- ed(n.nodes)
  # Check whether configuration and network are compatible
  if(n.edges != length(conf)){
    stop('Edges are not as expected from the number of nodes')
  }
  # Create edgelist
  #mat <- matrix(NA, n.nodes, n.nodes)
  #mat[upper.tri(mat)] <- conf
  #mat <- SymMat(mat)
  #diag(mat) <- 0
  #edge.list<- which(mat > 0, arr.ind = T )
  # Creeren edgelist
  mat <- matrix(NA,n.nodes, n.nodes)
  mat[upper.tri(mat)] <- seq(1:n.edges)
  edge.list<- which(mat > 0, arr.ind = T )
  
  # Create Omega
  if(q != 1 & Curie == F){
    l <- rep(list(0:1), n.edges)
    omega <- as.matrix(expand.grid(l))
    
  }
  # Initialiseer Zrc op 0
  Zrc <- 0
  # Wanneer q is 1, gelijkstellen aan erdos renyi en p teruggeven
  if(q == 1){
    p <- prod(apply(matrix(conf,nrow = 1),2,b,p = p))
    
    # Log likelihood
    #p <- (sum(conf)*log(p)+((length(conf)) - sum(conf))*log(1-p))
    return(p)
  }
  if(Curie == T){ # Use the Curie Weiss partition function
    Zrc <- ZC(n.nodes,p)
  }
  if(q != 1 & Curie == F){
    p1 <- rep(NA, nrow(omega))
    for(i in 1:nrow(omega)){
      
      net <- graph_from_data_frame(d = matrix(edge.list[omega[i,]!=0,], ncol = 2),
                                   vertices = seq(1:n.nodes), directed = F)
      k <- components(net)$no
      p1[i] <- (prod(apply(omega[i,, drop = FALSE],2,b,p = p))*q^k)
      Zrc <- Zrc + p1[i]
    }
  }
  
  # Current configuration
  conf.mat <- matrix(0, n.nodes, n.nodes)
  conf.mat[upper.tri(conf.mat)] <- conf
  conf.mat <- SymMat(conf.mat)
  
  net <- graph_from_adjacency_matrix(conf.mat, mode = "undirected")
  
  #net <- graph_from_data_frame(d = matrix(edge.list[conf != 0,], ncol = 2),
                               #vertices = seq(1:n.nodes), directed = F)
  k <- components(net)$no
  
  # Normalize
  prob <- (prod(apply(matrix(conf,nrow = 1),2,b,p = p))*q^k)/(Zrc)
   #prob <- (sum(conf)*log(p)+((length(conf)) - sum(conf))*log(1-p) + (k * log(q)))-log(Zrc)

  
  #prob <- (prod(conf, p)*q^k)/Zrc
  #prob <- (prod(apply(matrix(conf,nrow = 1),2,prod,p = p))*q^k)/(Zrc)
  #p1 <- p1/Zrc
  #dat <- cbind(omega, p1)
  return(prob)
}


