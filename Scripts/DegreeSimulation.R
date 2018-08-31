#########################################################################
######## Code for degree distribution  comparison RC and ER #############
##########             Code by Timo Matzen                ###############
##########    All the functions used are part of the      ###############
## RandomCluster package: https://github.com/TimoMatzen/RandomCluster ###
#########################################################################

# Load packages
library(RandomCluster)
library(igraph)


# Set sampling parameters
n.nodes <- 30
n.iter <- 100

# Critical value chosen so that the graph will be connected
# Can be changed to disconnected value from thesis
p <- .05

# Matrices to save the degrees
degreesRC <- matrix(NA, nrow = n.iter, ncol = n.nodes)
degreesER <- matrix(NA, nrow = n.iter, ncol = n.nodes)

## Vectors to keep track of unconnected networks
# RC not connected
ncRC <- 0
# ER not connected
ncER <- 0

# Iter not connected
ncItRC <- vector()
ncItER <- vector()

# Start sampling 
for (i in 1:n.iter) {
  
  # Keep track of iteration
  print(paste0("Sampling at iteration ", i))
  
  ## Sample RC
  RC.mat <-RandomnessRecycler(n.nodes, p, q = 2)
  ## Sample ER
  ER.mat <- RandomnessRecycler(n.nodes, p, q = 1)
  
  # Create graph object
  RC <- graph_from_adjacency_matrix(RC.mat, mode = "undirected")
  ER <- graph_from_adjacency_matrix(ER.mat, mode = "undirected")
  
  # Get the degrees
  degreesRC[i, ] <- degree(RC)
  degreesER[i, ] <- degree(ER)
  
  if (components(RC)$no != 1) {
    # Count how many configurations are not connected
    ncRC <- ncRC +1
    # Save the iteration
    ncItRC <- append(ncItRC, i)
    
  }
  
  if (components(ER)$no != 1) {
    # Count how many configurations are not connected
    ncER <- ncER +1
    # Save the iteration
    ncItER <- append(ncItER, i)
    
  }
  

}

# Table the degree distributions
degreesRC <- table(degreesRC)
degreesER <- table(degreesER)



# binomial sample
bin <- table(rbinom(n.nodes * n.iter, n.nodes-1, p))

## Plot the data 
# Set plotting parameters
par(cex.main = 1.3, mar = c(4.5, 6, 4, 7) + 0.1, mgp = c(3, 1, 0), cex.lab = 1.3, 
    font.lab = 2, cex.axis = 1.3, las = 1)
# Plot first lines
plot(names(degreesRC), as.numeric(degreesRC), col = "red", lwd = 2, 
     bty = "n", main = "Degree Distribution Disconnected", xaxt = "none", 
     xlim = c(0, 10), ylim = c(0,1400),type = "l", cex.main = 2, xlab = "", ylab = "")
axis(side = 1, at = seq(0, 10, 1), las = 2)
# Add lines
lines(names(bin), as.numeric(bin), col = "black", lwd = 2, type = "b")
lines(names(degreesER), as.numeric(degreesER), col = "blue", 
      lwd = 2, type = "l")
#lines(names(pois), as.numeric(pois), lwd = 2, type = "b", col)
# Add legend and texts
legend("topright", legend = c("Random Cluster model", "Erdös-Renyi model", 
                              "Binomial Distribution"),
       col=c("red", "blue", "black"), lty = c(1, 1, 1), lwd = 3, 
       cex = 1.2, bty = "n")
mtext("Degree", side = 1, 
      line = 2.5, las = 1, cex = 1.3)
mtext("Frequencies", side = 2, 
      line = 3.5, las = 0, cex = 1.7)

#########################################################################
####################Sample Disconnected##################################
#########################################################################



# Set sampling parameters
n.nodes <- 30
n.iter <- 100

# Critical value chosen so that the graph will be connected
p <- .14

# Matrices to save the degrees
degreesRC <- matrix(NA, nrow = n.iter, ncol = n.nodes)
degreesER <- matrix(NA, nrow = n.iter, ncol = n.nodes)

# RC not connected
ncRC <- 0
# ER not connected
ncER <- 0

# Iter not connected
ncItRC <- vector()
ncItER <- vector()

# Start sampling 
for (i in 1:n.iter) {
  
  # Keep track of iteration
  print(paste0("Sampling at iteration ", i))
  
  ## Sample RC
  RC.mat <-RandomnessRecycler(n.nodes, p, q = 2)
  ## Sample ER
  ER.mat <- RandomnessRecycler(n.nodes, p, q = 1)
  
  # Create graph object
  RC <- graph_from_adjacency_matrix(RC.mat, mode = "undirected")
  ER <- graph_from_adjacency_matrix(ER.mat, mode = "undirected")
  
  # Get the degrees
  degreesRC[i, ] <- degree(RC)
  degreesER[i, ] <- degree(ER)
  
  if (components(RC)$no != 1) {
    # Count how many configurations are not connected
    ncRC <- ncRC +1
    # Save the iteration
    ncItRC <- append(ncItRC, i)
    
  }
  
  if (components(ER)$no != 1) {
    # Count how many configurations are not connected
    ncER <- ncER +1
    # Save the iteration
    ncItER <- append(ncItER, i)
    
  }
  
  
}

# Table the degree distributions
degreesRC <- table(degreesRC)
degreesER <- table(degreesER)

# Save the data
save(degreesRC, file = "RCdeg30dis.Rdata")
save(degreesER, file = "ERdeg30dis.Rdata")

# Plot the data
par(mfrow = c(1,1))

# binomial sample
bin <- table(rbinom(n.nodes * n.iter, n.nodes-1, .14))

# Plot the data 
par(cex.main = 1.3, mar = c(4.5, 6, 4, 7) + 0.1, mgp = c(3, 1, 0), cex.lab = 1.3, 
    font.lab = 2, cex.axis = 1.3, las = 1)
plot(names(degreesRC), as.numeric(degreesRC), col = "red", lwd = 2, 
     bty = "n", main = "Degree Distribution Disconnected", xaxt = "none", 
     xlim = c(0, 14), type = "l", cex.main = 2, xlab = "", ylab = "")
axis(side = 1, at = seq(0, 14, 1), las = 2)
lines(names(bin), as.numeric(bin), col = "black", lwd = 2, type = "b")
lines(names(degreesER), as.numeric(degreesER), col = "blue", 
      lwd = 2, type = "l")
legend("topright", legend = c("Random Cluster model", "Erdös-Renyi model", "Binomial Distribution"),
       col=c("red", "blue", "black"), lty = c(1, 1, 1), lwd = 3, cex = 1.2, bty = "n")
mtext("Degree", side = 1, 
      line = 2.5, las = 1, cex = 1.3)
mtext("Frequencies", side = 2, 
      line = 3.5, las = 0, cex = 1.7)
