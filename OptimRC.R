## Using Optim to fit a RC model on Erdos Renyi simulated data
library(igraph)
library(RandomCluster)

# Define number of nodes in model
n.nodes = 10

# First define the function that you want to minimalize
# Wrapper for RCprob
RCfn <- function(data, par) {
  # Remember to adjust n.nodes argument within this wrapper when using bigger networks
  sum(-log(apply(data,1, RCProb, p = par, q = 2, n.nodes = n.nodes, Curie = T)))
}

# Wrapper voor RCprob Erdos Renyi
ERfn <- function(data, par) {
  # Remember to adjust n.nodes argument within this wrapper when using bigger networks
  sum(-log(apply(data, 1, RCProb, p = par, q = 1, n.nodes = n.nodes)))
}

# Simulate data with erdos renyi
n.iter <- 100
p.grid <- seq(0.1, 0.9, by = 0.05)
# Dataframe with optimized parameters
par <- rep(NA, length(p.grid))
for(j in 1:length(p.grid)) {
  data <- matrix(NA, nrow = n.iter, ncol = ed(n.nodes) )
  for (i in 1:n.iter) {
    #mat <- matrix(as_adjacency_matrix(erdos.renyi.game(n.nodes, p.grid[j])), 
                   # nrow = n.nodes)
    mat <- RandomnessRecycler(n.nodes, p.grid[j], q = 2)
    data[i, ] <- matrix(mat[upper.tri(mat)])
  }
  # Run optim
  result <- optim(par = 0, fn = ERfn, data = data, lower = .01, upper = .9, method = "Brent")

  # Save result
  par[j] <- result$par
  print(j)
}
 

# Make a plot of the optimalized values vs ER sample values

qqplot( par, p.grid,  ylim = range(p.grid,par), xlim = range(p.grid,par),
       xaxt="n", yaxt = "n", 
       main = paste0("Optim fit for RC model (q = 2) on ER sampled (q = 1) data. (", n.nodes, " nodes)"),
       xlab = "Optimized p-values (optim) for RC model",
       ylab = "ER sample p values (1000 iters)"
       )
axis(2, p.grid, las = 2)
axis(1, p.grid, las = 2)
abline(0,1)
