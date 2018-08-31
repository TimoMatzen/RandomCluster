#########################################################################
########   Code for MLE Simulation for the RC and ER model  #############
##########             Code by Timo Matzen                ###############
##########    All the functions used are part of the      ###############
## RandomCluster package: https://github.com/TimoMatzen/RandomCluster ###
#########################################################################

# Load packages
library(igraph)
library(RandomCluster)

# Define number of nodes in model
n.nodes <- 15

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
p.grid <- seq(0.1, 0.9, by = 0.01)

# Dataframe with optimized parameters
par <- rep(NA, length(p.grid))
for(j in 1:length(p.grid)) {
  data <- matrix(NA, nrow = n.iter, ncol = ed(n.nodes) )
  for (i in 1:n.iter) {
    #mat <- matrix(as_adjacency_matrix(erdos.renyi.game(n.nodes, p.grid[j])), 
                   # nrow = n.nodes)
    mat <- RandomnessRecycler(n.nodes, p.grid[j], q = 1)
    data[i, ] <- matrix(mat[upper.tri(mat)])
  }
  # Run optim
  result <- optim(par = 0, fn = RCfn, data = data,
                  lower = .01, upper = .9, method = "Brent")

  # Save result
  par[j] <- result$par
  print(j)
}

# Save as dataframe
MLE15nod <- data.frame("True"= p.grid, "MLE" = par)

# Set plotting parameters
par(cex.main = 1.3, mar = c(4.5, 6, 4, 7) + 0.1, mgp = c(3, 1, 0), cex.lab = 1.3, 
    font.lab = 2, cex.axis = 1.3, las = 1)
# PLot MLE estimates
qqplot(MLE5nod$MLE, p.grid,  ylim = range(p.grid,MLE5nod$MLE), 
       xlim = range(p.grid,MLE5nod$MLE),
       xaxt="n", yaxt = "n", 
       main = paste("MLE fit RC model on ER data (n = 5)", sep = " "),
       xlab = "",
       ylab = "", bty = "n", lwd = 2, cex.main = 2, cex.lab = 1.3, cex.axis = 1.3
       )
axis(2, at = seq(0.1,.9,.05), labels = seq(0.1,.9,.05))
axis(1, at = seq(0.1,.9,.05), labels = seq(0.1,.9,.05), las = 2)
ablineclip(a = 0,b = 1, lty = 2)
# Add text and legend
mtext("True p-value (ER)", side = 2, line = 3.5, las = 0, cex = 1.7)
mtext("MLE fit p-value (RC)", side = 1, 
      line = 3.5, las = 1, cex = 1.3)
legend("topleft", legend =  "True p == MLE p", lty = 2, lwd = 3, col = "grey",
       bty = "n",  cex = 1.4)
