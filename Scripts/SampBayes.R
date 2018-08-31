#########################################################################
######## Code for sample Bayesian model comparison RC and ER ############
##########             Code by Timo Matzen                ###############
##########    All the functions used are part of the      ###############
## RandomCluster package: https://github.com/TimoMatzen/RandomCluster ###
#########################################################################

# Import necessary packages
library(RandomCluster)
library(igraph)

###################################################
################Sample for 5 nodes#################
###################################################

# P grid for sampling
p.grid <- seq(0, 0.9,.05)

# Sample ten times largest data size, for variance
n.iter <- 600
# Degine number of nodes, can be adjusted for larger networks
n.nodes <- 5

# Create data matrix to store values
dat5nod <- array(0,dim=c(n.iter, ed(n.nodes),length(p.grid)))

# Simulate from random cluster model 
for (i in 1:length(p.grid)) {
  
  # Sample Random Cluster Data
  for (j in 1:n.iter) {
    print(j)
    mat <- RandomnessRecycler(n.nodes, p.grid[i], q = 2)
    # Save data 
    dat5nod[j, ,i] <- mat[upper.tri(mat)]
    
  }
  # Print iter
  print(paste0("Sampling with 5 nodes at iteration: ", i, "from ", length(p.grid)))
 
}

# Save data
save(dat5nod, file = "dat5nod.Rdata")

###################################################
################Sample for 10 nodes#################
###################################################

# Change number of nodes
n.nodes <- 10

# Create data matrix to store values
dat10nod <- array(0,dim=c(n.iter, ed(n.nodes),length(p.grid)))


# Simulate from random cluster model 
for (i in 1:length(p.grid)) {
  
  # Sample Random Cluster Data
  for (j in 1:n.iter) {
    print(j)
    mat <- RandomnessRecycler(n.nodes, p.grid[i], q = 2)
    dat10nod[j, ,i] <- mat[upper.tri(mat)]
    
  }
  # Print iter
  print(paste0("Sampling with 10 nodes at iteration: ", i, "from ", length(p.grid)))
  
}

# Save data
save(dat10nod, file = "dat10nod.Rdata")

###################################################
################Sample for 15 nodes#################
###################################################

# Change number of nodes
n.nodes <- 15

# Create data matrix to store values
dat15nod <- array(0,dim=c(n.iter, ed(n.nodes),length(p.grid)))


# Simulate from random cluster model 
for (i in 1:length(p.grid)) {
  
  # Sample Random Cluster Data
  for (j in 1:n.iter) {
    print(j)
    mat <- RandomnessRecycler(n.nodes, p.grid[i], q = 2)
    dat15nod[j,,i ] <- mat[upper.tri(mat)]
    
  }
  # Print iter
  print(paste0("Sampling with 15 nodes at iteration: ", i, "from ", length(p.grid)))
  
}

# Save data
save(dat15nod, file = "dat15nod.Rdata")



