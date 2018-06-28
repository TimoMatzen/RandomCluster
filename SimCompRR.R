##### Sample and look at components
library(igraph)
n.iter <- 1000
p <- seq(0.1,0.5,0.1)
n.nodes <- 10

# Vector for saving the number of components
comps.RC <- matrix(rep(NA, n.iter) * length(p), nrow = n.iter, ncol = length(p))
comps.ER <- matrix(rep(NA,n.iter) * length(p), nrow = n.iter, ncol = length(p))

for (j in 1:length(p)) {
  print(paste0("Starting sampling with p: ", p[j]))
  for (i in 1:n.iter) {
    mat.RC <- RandomnessRecycler(n.nodes, p[j], 2)
    graph <- graph_from_adjacency_matrix(mat.RC, mode = "undirected")
    # Get the number of components
    comps.RC[i, j] <- components(graph)$no
    # Sample from erdos renyu
    erdos <- erdos.renyi.game(n.nodes, p[j])
    # Get components
    comps.ER[i, j] <- components(erdos)$no
    print(i)
  
  }
  
}

# Add column names to the matrices
colnames(comps.RC) <- colnames(comps.ER) <- p



# Plot the components for each p value
for (j in 1:length(p)) {
  par(mfrow = c(1, 2))
  plot(table(comps.ER[, j]), 
       main = paste0( n.nodes, " nodes, ", n.iter, " iterations, p: ",
                     p[j], ", q = 1. (ER)"),
       xlab = "Number of components", ylab = "Frequency")
  plot(table(comps.RC[, j]), 
       main = paste0( n.nodes, ", nodes ", n.iter, " iterations, p: ", 
                     p[j], ", q = 2. (RC)"),
       xlab = "Number of components", ylab = "Frequency")
  
}

