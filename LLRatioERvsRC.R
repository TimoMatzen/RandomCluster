### Likelihood Ratio ###

# Sample from both models and calculate the likelihood

n.iter <- 100
p <- seq(0.1,0.5,0.005)
n.nodes <- 13

# Create matrices to output the data
probs.RC <- rep(NA, length(p))
probs.ER <- rep(NA, length(p))

# Configurations
conf.RC <- matrix(NA, ncol = length(p), nrow =  edges(13))

# Sample from each model and calculate probabilities
for (j in 1:length(p)) {
  print(paste0("Starting sampling with p: ", p[j]))
  
  mat.RC <- RandomnessRecycler(n.nodes = n.nodes, p[j], 2)
  # Get the configuration
  conf.RC[ ,j] <- mat.RC[upper.tri(mat.RC)]
    
  # Sample from erdos renyi
  #erdos <- erdos.renyi.game(n.nodes, p[j])
  #mat.ER <- matrix(get.adjacency(erdos), n.nodes, n.nodes)
  #conf.ER <- mat.ER[upper.tri(mat.ER)]
    
  # Get RC measure q = 2
  probs.RC[j] <- RCProb(conf.RC[,j], p[j], q = 2, n.nodes, Curie = TRUE)
  probs.ER[ j] <- RCProb(conf.RC[,j], p[j], n.nodes, q = 1)
  
  print(i)
    
  
}


llr.RC <- probs.RC/probs.ER
par(mfrow = c(1,1))
plot(p,llr.RC , ylab = "Likelihoor ratio (Q=2 / Q=1)", 
     main = "Samples from the random cluster model (q=2) with 13 nodes.", xaxt = "n",
     yaxt = "n")
axis(1, at = seq(.1, .5, by = .025), las=2)
axis(2, at = seq(0, 25, by = 1), las=2)

abline(1,0)

which.max(llr.RC)


## Samples from the erdos renyi model

p <- seq(0.1,0.5,0.005)
n.nodes <- 13

# Create matrices to output the data
probs.RC <- rep(NA, length(p))
probs.ER <- rep(NA, length(p))

# Add labels
colnames(probs.RC) <- colnames(probs.ER) <- p

# Sample from each model and calculate probabilities
for (j in 1:length(p)) {
  print(paste0("Starting sampling with p: ", p[j]))
  
  # Sample from erdos renyi
  erdos <- erdos.renyi.game(n.nodes, p[j])
  mat.ER <- matrix(get.adjacency(erdos), n.nodes, n.nodes)
  conf.ER <- mat.ER[upper.tri(mat.ER)]
  
  # Get RC measure q = 2
  probs.RC[j] <- RCProb(conf.ER, p[j], q = 2, n.nodes, Curie = TRUE)
  probs.ER[ j] <- RCProb(conf.ER, p[j], n.nodes, q = 1)
  
  print(i)
  
  
}


llr.ER <- probs.ER/probs.RC
par(mfrow = c(1,1))
plot(p,llr.ER , ylab = "Likelihood ratio (Q=1/Q=2)", 
     main = "Samples from the Erdos-Renyi model with 13 nodes.",
     xaxt = "n", yaxt = "n")
axis(1, at = seq(.1, .9, by = .025), las=2)
axis(2, at = seq(0, 15, by = 1), las=2)

abline(1,0)

