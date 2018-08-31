## Bayes Factor calculation voor the RC model (Q=2) and the ER model (Q=1)
library(RandomCluster)
library(igraph)

n.iter <- 1000
n.nodes <- 10
data <- matrix(NA, nrow = n.iter, ncol = ed(n.nodes) )

# Sample erdos renyi data
for (i in 1:n.iter) {
  p <- .2
  mat <- matrix(as_adjacency_matrix(erdos.renyi.game(n.nodes, p)), 
                nrow = n.nodes)
  data[i,] <- mat[upper.tri(mat)]
}

# Sample Random Cluster Data
for (i in 1:n.iter) {
  print(i)
  p <- runif(1,0,1)
  mat <- RandomnessRecycler(n.nodes, p, q = 2)
  data[i, ] <- mat[upper.tri(mat)]
  
}

function(data,p)

integrate(Vectorize(function(p) RCProb(data[,], p = p, q = 1, n.nodes = n.nodes, Curie = FALSE)), 0, 1)$value



RCfn(.3)
# 

integrate()

# First set the bayes factor to one a.k.a no evidence for either model.
BF <- 1

BFs <- rep(NA, n.iter)

# Loop over the simulated data
for (i in 1:n.iter) {
  
  # Wrapper for RCprob with q = 2
  RCfn <- function(par) {
    RCProb(data[i,], p = par, q = 2, n.nodes = n.nodes, Curie = TRUE)
  }
  
  # Wrapper for RCprob with q = 1
  ERfn <- function(par) {
    RCProb(data[i,], p = par, q = 1, n.nodes = n.nodes, Curie = FALSE)
  }
  
  # Integrate over RCfn, from theta 0 to 1
  RCmod <- integrate(Vectorize(RCfn), 0, 1)$value
  ERmod <- integrate(Vectorize(ERfn), 0, 1)$value
  
  # Calculate the update for the bayes factor
  update <- RCmod/ERmod
  print(update)
  
  # Update the bayes factor
  BF <- BF * update
  
  BFs[i] <- BF
  
}


# Plot RC sampled data
plot(y = log(BFs), x = seq(1,n.iter), type = "l", xaxt = "n",
     ylim = c(round(-1 * max(log(BFs))), round(max(log(BFs)))), yaxt = "n",
     main = paste0("Bayes factor for ", n.iter, " samples from the RC model with ", n.nodes, " nodes. h0 = ER; h1 = RC"),
     xlab = "Sample n")
axis(2, at = seq(round(-1 * max(log(BFs)), 0), round(max(log(BFs)),0), by= 2), las = 1)
axis(1, at = seq(0, n.iter, 20), las = 2)
abline(0,0)

# PLot Er sampled
plot(y = log(BFs), x = seq(1,n.iter), type = "l", xaxt = "n",
     ylim = c(round(min(log(BFs)),2), round(-1 * min(log(BFs)),2)), yaxt = "n",
     main = paste0("Bayes factor for ", n.iter, " samples from the ER model with ", n.nodes, " nodes. h0 = ER; h1 = RC"),
     xlab = "Sample n")
axis(2, at = seq(round(min(log(BFs)), 0),round(-1 * min(log(BFs)),0), by =5), las = 1)
axis(1, at = seq(0, n.iter, 20), las = 2)
abline(0,0)
legend("topright", legend = "Under is evidence h0: ER; above is evidence h1: RC", pch='l')

## Poging om op een andere manier te integreren

