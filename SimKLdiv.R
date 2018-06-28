### Calculating KL divergences ##
library(RandomCluster)
library(igraph)

# Set true parameters
n.nodes <- 5
p.true <- .5
n.iter <- 30

# Get omega 

mat <- matrix(NA,n.nodes, n.nodes)
mat[upper.tri(mat)] <- seq(1:n.edges)
edge.list<- which(mat > 0, arr.ind = T )
l <- rep(list(0:1), n.edges)
omega <- as.matrix(expand.grid(l))

# Calculate the kl divergence, try p = .2
# Create a grid of p values
p.mod <- seq(0.01,.9, .01)

# Create matrix to save divergences between the two models
kl.div <- rep(NA,length(p.mod))

# Go over all the p values
for (i in 1:length(p.mod)){
  kl.div[i] <- KLdivergence(omega, n.nodes, p.true, p.mod[i])
}
par(mfrow = c(1,1))
plot(x = p.mod, y= kl.div, xaxt = "n", 
     main = paste0("KL divergence between ER & RC with p = ", p.true),
     xlab = "P value RC", ylab = "KL divergence")
axis(1, at = seq(0, .9, .05), las = 2)

abline(v = p.mod[which.min(kl.div)] )

# Keep p the same for the kldiv
p.mod <- p.true <- seq(0.01,.9, .01)

for (i in 1:length(p.mod)){
  kl.div[i] <- KLdivergence(data, n.nodes, p.true[i], p.mod[i])
}
par(mfrow = c(1,1))
plot(x = p.mod, y= kl.div, xaxt = "n", main = "KL-divergence between RC and ER with 4 nodes; equal p ",
     xlab = "P value", ylab = "KL divergence")
axis(1, at = seq(0, .9, .05), las = 2)
legend("topright", legend = paste0("Maximum KL divergence at p: ", p.mod[which.max(kl.div)]  ))

abline(v = p.mod[which.max(kl.div)] )








### Bigger network with grid of p-values
n.iter <- 1000

# Create a grid of p values
p.mod <- p.true <- seq(0,.99, .02)


# define the number of nodes
n.nodes <- 20

kl.div <- rep(NA,length(p.mod))

for(j in 1:length(p.mod)) {
  data <- matrix(NA, nrow = n.iter, ncol = ed(n.nodes) )
  
  # Take 1000 samples each time
  for (i in 1:n.iter) {
    mat <- matrix(as_adjacency_matrix(erdos.renyi.game(n.nodes, p.mod[j])), 
                  nrow = n.nodes)
    data[i,] <- mat[upper.tri(mat)]
  }
  
  # Calculate KL divergence
  kl.div[j] <- KLdivergence(data, n.nodes, p.true[j], p.mod[j])
  print(j)
 
}

par(mfrow = c(1,1))
plot(x = p.mod, y= kl.div, xaxt = "n",yaxt = "n", main = paste0("KL-divergence between RC and ER with ", n.nodes, " nodes; equal p ", n.iter, " samples ER."),
     xlab = "P value", ylab = "KL divergence")
axis(1, at = seq(0, 1, .05), las = 2)
axis(2, at= seq(0, max(kl.div),1), las = 2)
legend("topright", legend = paste0("Maximum KL divergence at p: ", p.mod[which.max(kl.div)]  ))

abline(v = p.mod[which.max(kl.div)] )
