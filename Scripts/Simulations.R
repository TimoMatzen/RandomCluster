### Investigate the RC model degree distribution ##
library(igraph)
library(qgraph)


### Random Cluster model


### Simulate the degree distribution

# Set number of iterations
n.iter <- 500
# Set parameters
p <- seq(.1,.9,.1)
# Ising model 
q <- 2
# Set number of nodes
n.nodes <- 15
# Vector to save mean degrees
degree.RC <-matrix(NA, nrow = n.iter, ncol = length(p)) 
for (j in 1:length(p)) {
  print(paste0("Sampling with p: ", p[j]))
  for (i in 1:n.iter){
    # Take a network sample 
    RC.samp1 <- RandomnessRecycler( n.nodes,p[j], 2 )
    
    # Save the degree of node 1
    degree.RC[i,j] <- rowSums(RC.samp1)[1]
    
    print(i)
  }
}

## Erdos Renyi random graph
# Erdos Renyi model 

degree.ER <- matrix(NA, nrow = n.iter, ncol = length(p))

for (j in 1:length(p)) {
  print(paste0("Sampling with p: ", p[j]))
  for (i in 1:n.iter) {
    ER.samp1 <- RandomnessRecycler( n.nodes,p[j], 1 )
    degree.ER[i, j] <-  rowSums(ER.samp1)[1]
    
    print(i)
  }
}


par(mfrow = c(3,3))
for (i in 1:length(p)) {
  boxplot(degree.RC[,i], degree.ER[,i], main = paste0("RC vs ER wirh p = ",p[i],  " with 1000 iterations and ",n.nodes, " nodes."),
          col = c(5,6),  names = c('RC', 'ER'), yaxt = "n")
  axis(2, seq(0,20,1))
}

par(mfrow = c(1,1))
qqplot(degree.ER[,1], degree.RC[, 1])

