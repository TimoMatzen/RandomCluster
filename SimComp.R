### Components Simulations ###
library(igraph)
library(RandomCluster)
# Different values of p 
p <-0.5
# Number of nodes
n.nodes <- 10
# Number of iterations
n.iter <- 1000
# Vector to save components RC
RC.k <- rep(0, n.iter)
# Vector to save components ER
ER.k <- rep(0, n.iter)
for (i in 1:n.iter) {
  # Simulate RC; q = 2
  RC <- SampleRC(p, 2, n.nodes, 1, mat = T, verbose = FALSE)
  # Create edgelist RC
  net.RC <- graph_from_adjacency_matrix(RC, mode = "undirected")
  # Get number of components
  RC.k[i] <- components(net.RC)$no
  # Simulate ER
  net.ER <-  erdos.renyi.game(n.nodes, p)
  # Get number of components
  ER.k[i] <- components(net.ER)$no
  # Print the progress
  print(i)
}

par(mfrow = c(1,2))
plot(table(RC.k), main = "Comp RC model for p=0.5; q=2")
plot(table(ER.k), main = "Comp ER model for p=0.5; q=1")
