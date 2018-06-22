### Investigate the RC model degree distribution ##
library(igraph)
library(qgraph)


### Random Cluster model
## Simulate mean degrees with changing p

# Set parameters
p <- seq(.1,.9, .05)
# Ising model 
q <- 2
# Set number of nodes
n.nodes <- 10
# Vector to save mean degrees
degree.RC <- rep(NA,length(p))
for (i in 1:length(p)){
  print(paste0("Now sampling with ", p[i]))
  deg <- rep(NA, 1000)
  for (j in 1:1000){
    
    par (mfrow = c(1,2))
    # Take a network sample 
    RC.samp1 <- RandomnessRecycler(n.nodes, p[i],  2)
    # Save the mean degrees
    deg[j] <- mean(rowSums(RC.samp1))
  
    # Plot the structure of the network
    #RC <- qgraph(RC.samp1)
    print(j)
  }
  degree.RC[i] <- mean(deg)
  
}


## Erdos Renyi random graph
# Erdos Renyi model 

degree.ER <- rep(NA, length(p))
for (i in 1:length(p)) {
  print(paste0("Now sampling with ", p[i]))
  deg <- rep(NA, 1000)
  for (j in 1:1000) {
    # Take a network sample erdos renyi, q = 1
    #ER.samp1 <- RandomnessRecycler(n.nodes, p[i], 1)
    ER.samp1 <- erdos.renyi.game(n.nodes, p[i])
    ER.samp1 <- matrix(get.adjacency(ER.samp1), n.nodes, n.nodes)
    
    # Get the degree of each node
    deg[j] <- mean(rowSums(ER.samp1))
    print(j)
  }
  degree.ER[i] <- mean(deg)
}

par(mfrow = c(1,1))
plot(x = p, y = degree.RC, ylab = "Mean Degree", xlab= "p",
     main = "Mean degrees Random Cluster vs Erdos Renyi", type = 'l',
     xaxt = "n", yaxt = "n")
axis(1, seq(0.1,0.9,0.05), las = 2)
axis(2, seq(min(degree.RC), max(degree.RC), 1), las = 2)

lines(x = p, y = degree.ER, ylab = "Mean Degree", xlab= "p",
     main = "Erdos-Renyi Degree", type = 'l', col = 'red')
legend("bottomright", col = c("black", "red"),legend = c("RC (q=2)", "ER (q=1)"), 
       pch = 1)

### Simulate the degree distribution

# Set number of iterations
n.iter <- 1000
# Set parameters
p <- seq(.1,.5,.1)
# Ising model 
q <- 2
# Set number of nodes
n.nodes <- 10
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
    graph <- erdos.renyi.game(n.nodes, p[j])
    degree.ER[i, j] <- rowSums(matrix(get.adjacency(graph), n.nodes, n.nodes))[1]
  
    print(i)
  }
}


par(mfrow = c(1,1))
boxplot(degree.RC[,5], degree.ER[,5], main = "RC(p=.5, q=2) vs ER(p=.5) with 1000 iterations and 10 nodes.",
        col = c(5,6),  names = c('RC', 'ER'), yaxt = "n")
axis(2, seq(0,20,1))



