#########################################################################
##########  Code for Bayesian model comparison RC and ER  ###############
##########             Code by Timo Matzen                ###############
##########    All the functions used are part of the      ###############
## RandomCluster package: https://github.com/TimoMatzen/RandomCluster ###
#########################################################################

# Load packages
library(RandomCluster)
library(igraph)

load("bayes5nod.Rdata")

# Define number of nodes, adjust for different network sizes
n.nodes <- 5

# Define grid for different sample sizes
grid <- matrix(c(seq(1,110,10), seq(1, 330, 30), seq(1, 660, 60)), ncol = 11, 
               nrow = 3, byrow = TRUE )

# Array to save bayes factors, third dimensions is the different sample size
bayes <- array(NA, dim = c(10, dim(dat5nod)[3], 3))

# Define integrand
integrand <- function(x,p, q) {
  if (q == 2) {
    Pr <- dbeta(p,1,1) * apply(x, 1, RCProb, p = p, q = q , 
                               n.nodes = n.nodes, Curie = TRUE)
    
    Pr <- sum(log(Pr))
    
    return(Pr)
    
  } else {
    Pr <- dbeta(p,1,1) * apply(x, 1, RCProb, p = p, q = q , 
                               n.nodes = n.nodes)
    Pr <- sum(log(Pr))
    return(Pr)
    
  }
}



# Set probs, for getting macimum for preventing underflow
probs <- seq(0.1, 0.9, .1)


# Loop over different sample sizes
for (k in 1:3) {
  # Print progress
  print(paste("Sample at ", k))
  # Loop over pvalues
  for (i in 1:dim(dat5nod)[3]) {
    
    
    # Loop over samples
    for (g in 1:(length(grid[k,])-1)) {
      
      # Select data to calculate BF
      X <- dat15nod[grid[k, g]:(grid[ k, g+1]-1),,i]
      
      # Caculate max sum log RC
      maxRC <- function(p, x) {
        a <- sapply(p, integrand, x=x, q=2 )
        return(max(a))
      }
      
      bRC <- maxRC(p = probs, x = X)
      
      # Define integral RC
      integrandRC <- function(p,x) {
        
        a <- sapply(p, integrand, x=x, q= 2 )
        return(exp(a-bRC))
        
      }
      
      # Define integral ER
      integrandER <- function(p,x) {
        
        a <- sapply(p, integrand, x=x, q=1 )
        return(exp(a-bRC))
        
      }
      
      
      # Calculate hypotheses
      H1 <- integrate(function(x) integrandRC(x, X),0, 1)$value
      
      H0 <- integrate(function(x) integrandER(x, X),0, 1)$value
      
      # Calculate bayes factor
      a <- H1/H0
      
      # Save the bayes factor
      bayes[g, i, k] <- a
     }
    
  
  }
}

bayes15nod <- bayes
  
save(bayes15nod, file = "bayes15nod.Rdata")
