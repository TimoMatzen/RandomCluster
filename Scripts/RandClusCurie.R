# Laden packages
library(igraph)


# Dit is de partitie functie voor q == 2 vanuit de afleiding van Maarten.
ZC <- function(n, p) {
  # Eerst moet sigma berekend worden, afgeleid van Maartens document
  sig <- (log(1-p))/-2
  # Z eerst gelijk aan 0
  Z <- 0
  # Vervolgens loopen over alle nodes.
  for (i in 0:n){
    Z <- Z + (choose(n,i)*exp((.5*sig)*(2*i-n)^2))
  }
  # Om de partitie functie van curie-weiss gelijk te krijgen aan die van het RC nog
  # vermenigvuldigen met een term die afhangt van sigma en het aantal nodes.
  Z <- Z*exp((-.5*sig*n^2))
  return(Z)
}

# Functie voor het eerste product
b <- function(x,p){
  p^(x)*(1-p)^(1-x)
}

# Function for copying upper tri matrix to lower tri
SymMat <- function(m) {
  m[lower.tri(m)] <- t(m)[lower.tri(m)]
  m
}

# Function to just calculate the RC function
RC <- function(n, s, p, q) {
  # n: Number of nodes
  # s: State of system
  # p: Probability
  # q: Clustering measure
  
  # Create empty adjacency matrix
  mat <- matrix(0,n,n)
  # Fill upper tri with the state
  mat[upper.tri(mat)] <- s
  # Make matrix symmetric
  mat <- SymMat(mat)
  # Create the edge list
  net <- graph_from_adjacency_matrix(mat, mode = "undirected")
  # Get number of components
  k <- components(net)$no
  # Calculate RC measure
  out <- p^sum(s) * (1-p)^(1-sum(s)) * q^k
  return(out)
  
}

# Random cluster functie
RCProb <- function(p,q,n,conf, Curie = F){
  if(Curie == T & q != 2){
    stop('Can not use Curie-Weiss partition function when q != 2')
  }
  # Alle mogelijke edges.
  n.edges <- n*((n-1)/2)
  if(n.edges != length(conf)){
    stop('Edges are not as expected from the number of nodes')
  }
  # Creeren edgelist
  mat <- matrix(NA,n,n)
  mat[upper.tri(mat)] <- seq(1:n.edges) 
  edge.list<- which(mat > 0, arr.ind = T )
  
  # Omega creeren
  if(q != 1 & Curie == F){
    l <- rep(list(0:1), n.edges); omega <- as.matrix(expand.grid(l))
  }
  # Initialiseer Zrc op 0
  Zrc <- 0
  # Wanneer q is 1, gelijkstellen aan erdos renyi en p teruggeven
  if(q == 1){
    p <- prod(apply(matrix(conf,nrow = 1),2,b,p = p))
    return(p)
  }
  if(Curie == T){
    Zrc <- ZC(n,p)
  }
  if(q != 1 & Curie == F){
    p1 <- rep(NA, nrow(omega))
    for(i in 1:nrow(omega)){
      
      net <- graph_from_data_frame(d = matrix(edge.list[omega[i,]!=0,],ncol = 2),vertices = seq(1:n), directed = F)
      k <- components(net)$no
      p1[i] <- prod(apply(matrix(omega[i,], nrow = 1),2,b,p=p))*q^k
      Zrc <- Zrc + p1[i]
    }
  }
  
  net <- graph_from_data_frame(d = matrix(edge.list[conf != 0,], ncol = 2),vertices = seq(1:n), directed = F)
  k <- components(net)$no
  prob <- (prod(apply(matrix(conf,nrow = 1),2,b,p = p))*q^k)/(Zrc)
  #p1 <- p1/Zrc
  #dat <- cbind(omega, p1)
  return(prob)
}
