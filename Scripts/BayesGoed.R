#########################################################################
##########  Code for Bayesian model comparison RC and ER  ###############
##########             Code by Timo Matzen                ###############
##########    All the functions used are part of the      ###############
## RandomCluster package: https://github.com/TimoMatzen/RandomCluster ###
#########################################################################

# Load packages
library(RandomCluster)
library(igraph)

# Set parameters
p.grid <- seq(0.1,0.9,.05)
bayes <- rep(0,length(p.grid))
n.iter <- 60
n.nodes <- 15

# Define function to integrate
integrand <- function(x,p, q) {
  if (q == 2) {
    Pr <- dbeta(p,1,1) * RCProb(x,p, q = q, n.nodes = n.nodes, Curie = TRUE)
    #return(log(Pr))    
    #print(Pr * 10^18)
    return(Pr * (1.55*10^15))
    
  } else {
    Pr <- dbeta(p,1,1) * RCProb(x,p, q = q, n.nodes = n.nodes)
    
    return(Pr * (1.55*10^15) )
    
  }
}

# starting value power
power <- 13

# Simulate from random cluster model and calculate bayes factor
for (i in 1:length(p.grid)) {

  # Simulate data
  data <- matrix(NA, nrow = n.iter, ncol = ed(n.nodes) )
  
  
  # Sample Random Cluster Data
  for (j in 1:n.iter) {
    print(j)
    mat <- RandomnessRecycler(n.nodes, p.grid[i], q = 2)
    data[j, ] <- mat[upper.tri(mat)]
    
  }
  # Rest H1 and H0 to 0
  H1 <- H0 <- 0
  # Raise constant in integrand untill H1 and  H0 != 0
  while (H1 == 0 | H0 == 0) {
    
    
    integrand <- function(x,p, q) {
      if (q == 2) {
        Pr <- dbeta(p,1,1) * RCProb(x,p, q = q, n.nodes = n.nodes, Curie = TRUE)
        #return(log(Pr))    
        #print(Pr * 10^18)
        return(Pr * (1.55*10^power))
        
      } else {
        Pr <- dbeta(p,1,1) * RCProb(x,p, q = q, n.nodes = n.nodes)
        #return(log(Pr))
        #
        #print(Pr * 10^16)
        return(Pr * (1.55*10^power) )
        
      }
    }
    
    
    # Calculate hypotheses
    H1 <- integrate(Vectorize(function(x) 
      prod(apply(data, 1, integrand, p = x, q = 2))),0,1)$value
    
    H0 <- integrate(Vectorize(function(x) 
      prod(apply(data, 1, integrand, p = x, q = 1))) ,0,1)$value
    
    power <- power + 1
  }
  
  # Calculate bayes factor
  a <- H1/H0
  print(a)
  # Save the bayes factor
  bayes[i] <- a
  
  # Increase power after each iteration
 
}

# Save data
save(bayes, file = "bayes15nod.Rdata")

plot(y = log(bayes), x = p.grid, xaxt = "n",type = "l",
     ylim = c(round(-1 * max(log(bayes))), round(max(log(bayes)))), yaxt = "n",
     main = paste0("Bayes factor for ", n.iter, " samples from the RC model with ", n.nodes, " nodes. h0 = ER; h1 = RC"),
     xlab = "P", bty = "n")
axis(2, at = seq(round(-1 * max(log(bayes)), 0), round(max(log(bayes)),0), by= 2), las = 1)
axis(1, at = seq(0, p.grid[length(p.grid)],.05 ), las = 2)
abline(0,0)


p.grid <- seq(0.1,0.9,.05)
bayes <- rep(0,length(p.grid))
n.iter <- 60
n.nodes <- seq(5,10,1)
par(mfrow = c(3,2))
for (n in 1:length(n.nodes)) {
  # Integral function
  integrand <- function(x,p, q) {
    if (q == 2) {
      Pr <- dbeta(p,1,1) * RCProb(x,p, q = q, n.nodes = n.nodes[n], Curie = TRUE)
      #return(log(Pr))   
      #print(Pr * 10^9)
      return(Pr * 10^(n.nodes[n] -1))
    } else {
      Pr <- dbeta(p,1,1) * RCProb(x, p, q = q, n.nodes = n.nodes[n])
      #return(log(Pr))
      #
      #print(Pr * 10^9)
      return(Pr * 10^(n.nodes[n]-1 ))
    }
    
  }
    
  # Simulate from random cluster model and calculate bayes factor
  for (i in 1:length(p.grid)) {
    # Simulate data
    data <- matrix(NA, nrow = n.iter, ncol = ed(n.nodes[n]) )
    
    
    # Sample Random Cluster Data
    # Sample erdos renyi data
    for (j in 1:n.iter) {
      print(j)
      #mat <- matrix(as_adjacency_matrix(erdos.renyi.game(n.nodes, p.grid[j])), 
          #          nrow = n.nodes
      #data[j,] <- sample(c(0,1), ed(n.nodes[n]), prob = c((1-p.grid[i]), p.grid[i]), replace = T)
  
      mat <- RandomnessRecycler(n.nodes[n], p.grid[i], q = 2)
      data[j, ] <- mat[upper.tri(mat)]
    }
    
    # Calculate hypotheses
    H1 <- integrate(Vectorize(function(x) 
      prod(apply(data, 1, integrand, p = x, q = 2))),0,1)$value
    
    H0 <- integrate(Vectorize(function(x) 
      prod(apply(data, 1, integrand, p = x, q = 1))) ,0,1)$value
    
    # Calculate bayes factor
    a <- H1/H0
    print(a)
    # Save the bayes factor
    bayes[i] <- as.numeric(a)
    
    if(is.nan(a)) {
      break
    }
  }
  plot(y = log(bayes), x = p.grid, type = "l", xaxt = "n",
       ylim = c(round(-1 * max(log(bayes))), round(max(log(bayes)))), yaxt = "n",
       main = paste0("Bayes factor for samples from the RC model with ", n.nodes[n], " nodes. h0 = ER; h1 = RC"),
       xlab = "P")
  axis(2, at = seq(round(-1 * max(log(bayes)), 0), round(max(log(bayes)),0), by= 2), las = 1)
  axis(1, at = seq(0, p.grid[length(p.grid)],.05 ), las = 2)
  abline(0,0)
  #plot(y = log(bayes), x = p.grid, type = "l", xaxt = "n",
  #     ylim = c(round( min(log(bayes))), round(-1 * min(log(bayes)))), yaxt = "n",
  #     main = paste0("Bayes factor for ", n.iter, " samples from the ER model with ", n.nodes[n], " nodes. h0 = ER; h1 = RC"),
  #     xlab = "P")
  #axis(2, at = seq(round( min(log(bayes)), 0), round(-1 * min(log(bayes)),0), by= 2), las = 1)
  #axis(1, at = seq(0, p.grid[length(p.grid)],.05 ), las = 2)
  #abline(0,0)
}





H1 <- integrate(Vectorize(function(x) 
  prod(apply(data, 1, integrand, p = x, q = 2))),0,1)$value

as.numeric(a)
10^400

prod(apply(data, 1, integrand, p = .2, q = 1) )

RCProb(data, p = .2, q = 1, n.nodes = 10) #* dbeta(.2,1,1)

integraal <- function(conf,n.nodes, q = 1) {
  n <- length(conf)
  k <- sum(conf)
  
  # Calculate components
  mat <- matrix(0,n.nodes, n.nodes)
  mat[upper.tri(mat)] <- conf
  mat <- SymMat(mat)
  graph <- graph_from_adjacency_matrix(mat, mode = "undirected")
  comps <- components(graph)$no
  
  Z <- RC(n.nodes, conf, .3, q=2)
  
  # Calculate integral 
  integral <- (((factorial(k)*factorial(n-k))/factorial(n+1)) * q^comps)
  return(integral)
  
}

integrate(Vectorize(function(x)RCProb(data, x, 2, 5, Curie = TRUE )), 0, 1)$value
H0 <- integraal(data, 5)
H1 <- integraal(data, 5,  q=2)/

H1/H0

integrate(function(x) dbinom( sum(data), length(data), x), 0, 1)$value
