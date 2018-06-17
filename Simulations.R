### Investigate the RC model degree distribution ##
library(igraph)
library(qgraph)


### Random Cluster model
## Simulate mean degrees with changing p

# Set parameters
p <- seq(.1,.7,.01)
# Ising model 
q <- 2
# Set number of nodes
n.nodes <- 10
# Vector to save mean degrees
degree.RC <- rep(0,length(p))
for (i in 1:length(p)){
  par (mfrow = c(1,2))
  # Take a network sample 
  RC.samp1 <- SampleRC(p, q, n.nodes, 1, mat = TRUE )
  qgraph(RC.samp1)
  # Save the mean degrees
  degree.RC[i] <- mean(rowSums(RC.samp1))
  

  # Plot the structure of the network
  #RC <- qgraph(RC.samp1)
}
plot(x = p, y = degree.RC, ylab = "Mean Degree", xlab= "p",
     main = "Random Cluster Degree with q = 2", type = 'l')
# Plot the degree distribution
#hist(degree.RC)

## Erdos Renyi random graph
# Erdos Renyi model 

degree.ER <- rep(0, length(p))
for (i in 1:length(p)) {
  # Take a network sample erdos renyi, q = 1
  ER.samp1 <- SampleRC(p[i], 1, n.nodes, 1, mat = TRUE )

  # Get the degree of each node
  degree.ER[i] <- mean(rowSums(ER.samp1))
}

plot(x = p, y = degree.ER, ylab = "Mean Degree", xlab= "p",
     main = "Erdos-Renyi Degree", type = 'l')

# Plot the structure of the network
qgraph(ER.samp1, layout = RC$layout)


# Plot the degree distribution
#hist(degree.ER)


### Simulate the degree distribution

sqrt(2)/(1 + sqrt(2))
# Set number of iterations
n.iter <- 1000
# Set parameters
p <- .4
# Ising model 
q <- 2
# Set number of nodes
n.nodes <- 10
# Vector to save mean degrees
degree.RC <- rep(0,n.iter)
for (i in 1:n.iter){
  # Take a network sample 
  RC.samp1 <- SampleRC(p, q, n.nodes, 1, mat = TRUE )
  
  # Save the degree of node 1
  degree.RC[i] <- rowSums(RC.samp1)[1]
}


## Erdos Renyi random graph
# Erdos Renyi model 

degree.ER <- rep(0, n.iter)

for (i in 1:n.iter) {
  graph <- erdos.renyi.game(n.nodes,p)
  degree.ER[i] <- rowSums(matrix(get.adjacency(graph), n.nodes, n.nodes))[1]

}


par(mfrow = c(1,1))
boxplot(degree.RC, degree.ER, main = "RC(p=.4, q=2) vs ER(p=.4) with 1000 iterations and 10 nodes.",
        col = c(5,6),  names = c('RC', 'ER'))
q/
