#########################################################################
###### Code for Fragmentation Simulation for the RC and ER model  #######
##########             Code by Timo Matzen                ###############
##########    All the functions used are part of the      ###############
## RandomCluster package: https://github.com/TimoMatzen/RandomCluster ###
#########################################################################

# Load packages
library(igraph)
library(RandomCluster)


######################################
###### Sampling With 4 nodes #########
######################################

# Code can be adjusted to sample with different network sizes

# Different values of p, max div, p_c, super critical
# Needs to be adjusted for n.nodes != 4
p.grid <- c(.32, (2/4), (2/4) + .2)

# Number of nodes
n.nodes <- 4

# Number of iterations
n.iter <- 1000

# Matrix to save components RC
RC.k.four <- matrix(NA, n.iter, length(p.grid))

# Vector to save components ER
ER.k.four <-  matrix(NA, n.iter, length(p.grid))

# Start the simulation
for(p in 1:length(p.grid)) {
  # Print progress
  print(paste0("Sampling at p ", p.grid[p]))
  for (i in 1:n.iter) {
    
    # Simulate RC; q = 2
    RC <- RandomnessRecycler(n.nodes,p.grid[p], 2)
    # Create graph object RC
    graph.RC <- graph_from_adjacency_matrix(RC, mode = "undirected")
    
    # Get number of components
    RC.k.four[i, p] <- components(graph.RC)$no
    
    # Simulate ER
    ER <- RandomnessRecycler(n.nodes, p.grid[p], 1)
    # Create graph object ER
    graph.ER <- graph_from_adjacency_matrix(ER, mode = "undirected")
    # Get number of components
    ER.k.four[i, p] <- components(graph.ER)$no
    
    
  }
}

# Add colnames
colnames(RC.k.four) <- c("maxKL", "Pcrit", "Psupercrit")
colnames(ER.k.four) <- c("maxKL", "Pcrit", "Psupercrit")


# Create plot

# Adjust p to plot for different values of p.grid
p <- 1

# Set the plotting parameters
par(cex.main = 1.3, mar = c(4.5, 6, 4, 7) + 0.1, mgp = c(3, 1, 0), cex.lab = 1.3, 
    font.lab = 2, cex.axis = 1.3, las = 1)
# Create plot
plot(names(table(RC.k.four[,p])),as.numeric(table(RC.k.four[,p])), type = "b", frame = FALSE, pch = 19, 
     col = "red", xlab = "", ylab = "", main = expression(paste("Fragmentation for max KL-div,  p =  0.32"), sep = ""), lty = 1,
     lwd = 2, cex.main = 2, yaxt = "n", xaxt = "n", cex.lab = 1.3, cex.axis = 1.3, ylim = c(0, 1000))
axis(2, at = seq(0, 1000, 100), labels = seq(0, 1000, 100), las = 2)
axis(1, at = seq(1,4,1), labels = seq(1,4,1) )
# Add a second line
lines(names(table(ER.k.four[,p])), as.numeric(table(ER.k.four[,p])), pch = 18, col = "blue", type = "b",
      lty = 2, lwd = 2)
# Add a legend to the plot
legend("top", legend=c("Random Cluster model", "Erdös-Renyi model"),
       col=c("red", "blue"), lty = 1:2, cex=1.2, bty = "n", lwd = 2)
mtext("Frequencies", side = 2, line = 3.2, las = 0, cex = 1.7)
mtext("Components", side = 1, 
      line = 2.5, las = 1, cex = 1.3)


######################################
###### Sampling With 5 nodes #########
######################################

# Different values of p 
p.grid <- c(0.27, (2/5), (2/5)+.2)
# Number of nodes
n.nodes <- 5
# Number of iterations
n.iter <- 1000
# Matrix to save components RC
RC.k.five <- matrix(NA, n.iter, length(p.grid))
# Vector to save components ER
ER.k.five <-  matrix(NA, n.iter, length(p.grid))

# Start the simulation
for(p in 1:length(p.grid)) {
  # Print progress
  print(paste0("Sampling at p ", p.grid[p]))
  for (i in 1:n.iter) {
    print(i)
    # Simulate RC; q = 2
    RC <- RandomnessRecycler(n.nodes,p.grid[p], 2)
    # Create graph object RC
    graph.RC <- graph_from_adjacency_matrix(RC, mode = "undirected")
    # Get number of components
    RC.k.five[i, p] <- components(graph.RC)$no
    
    # Simulate ER
    ER <- RandomnessRecycler(n.nodes, p.grid[p], 1)
    # Create graph object ER
    graph.ER <- graph_from_adjacency_matrix(ER, mode = "undirected")
    # Get number of components
    ER.k.five[i, p] <- components(graph.ER)$no
    
    
  }
}

# Add colnames
colnames(RC.k.five) <- p.grid
colnames(ER.k.five) <- p.grid


# Save the data 
save(RC.k.five, file = "CompFiveRC.Rdata")
save(ER.k.five, file = "CompFiveER.Rdata")


p <- 3

par(cex.main = 1.3, mar = c(4.5, 6, 4, 7) + 0.1, mgp = c(3, 1, 0), cex.lab = 1.3, 
    font.lab = 2, cex.axis = 1.3, las = 1)
# Create a first line
plot(names(table(RC.k.five[,p])),as.numeric(table(RC.k.five[,p])), type = "b", frame = FALSE, pch = 19, 
     col = "red", xlab = "", ylab = "", main = expression(paste("Fragmentation for ", p > p[c] ,", p =  0.6"), sep = ""), lty = 1,
     lwd = 2, cex.main = 2, yaxt = "n", xaxt = "n", cex.lab = 1.3, cex.axis = 1.3, ylim = c(0, 1000))
axis(2, at = seq(0, 1000, 100), labels = seq(0, 1000, 100), las = 2)
axis(1, at = seq(1,5,1), labels = seq(1,5,1) )
# Add a second line
lines(names(table(ER.k.five[,p])), as.numeric(table(ER.k.five[,p])), pch = 18, col = "blue", type = "b",
      lty = 2, lwd = 2)
# Add a legend to the plot
legend("top", legend=c("Random Cluster model", "Erdös-Renyi model"),
       col=c("red", "blue"), lty = 1:2, cex=1.2, bty = "n", lwd = 2)
mtext("Frequencies", side = 2, line = 3.2, las = 0, cex = 1.7)
mtext("Components", side = 1, 
      line = 2.5, las = 1, cex = 1.3)


############################################
##########Sampling with 6 nodes#############
############################################

# Different values of p 
p.grid <- c(.23, (2/6), (2/6)+.2)
# Number of nodes
n.nodes <- 6
# Number of iterations
n.iter <- 1000
# Matrix to save components RC
RC.k.six <- matrix(NA, n.iter, length(p.grid))
# Vector to save components ER
ER.k.six <-  matrix(NA, n.iter, length(p.grid))

# Start the simulation
for(p in 1:length(p.grid)) {
  # Print progress
  print(paste0("Sampling at p ", p.grid[p]))
  for (i in 1:n.iter) {
    print(i)
    # Simulate RC; q = 2
    RC <- RandomnessRecycler(n.nodes,p.grid[p], 2)
    # Create graph object RC
    graph.RC <- graph_from_adjacency_matrix(RC, mode = "undirected")
    # Get number of components
    RC.k.six[i, p] <- components(graph.RC)$no
    
    # Simulate ER
    ER <- RandomnessRecycler(n.nodes, p.grid[p], 1)
    # Create graph object ER
    graph.ER <- graph_from_adjacency_matrix(ER, mode = "undirected")
    # Get number of components
    ER.k.six[i, p] <- components(graph.ER)$no
    
    
  }
}

# Add colnames
colnames(RC.k.six) <- p.grid
colnames(ER.k.six) <- p.grid


# Save the data 
save(RC.k.six, file = "CompSixRC.Rdata")
save(ER.k.six, file = "CompSixER.Rdata")

#par(mfrow = c(1,2))
#plot(table(RC.k.five[,8]), main = "Comp RC model for p=0.5; q=2")
#plot(table(ER.k.five[,8]), main = "Comp ER model for p=0.5; q=1")

p <- 3

par(cex.main = 1.3, mar = c(4.5, 6, 4, 7) + 0.1, mgp = c(3, 1, 0), cex.lab = 1.3, 
    font.lab = 2, cex.axis = 1.3, las = 1)
# Create a first line
plot(names(table(RC.k.six[,p])),as.numeric(table(RC.k.six[,p])), type = "b", frame = FALSE, pch = 19, 
     col = "red", xlab = "", ylab = "", main = expression(paste("Fragmentation for ", p > p[c] ,", p =  0.53"), sep = ""), lty = 1,
     lwd = 2, cex.main = 2, yaxt = "n", xaxt = "n", cex.lab = 1.3, cex.axis = 1.3, ylim = c(0, 1000))
axis(2, at = seq(0, 1000, 100), labels = seq(0, 1000, 100), las = 2)
axis(1, at = seq(1,6,1), labels = seq(1,6,1) )
# Add a second line
lines(names(table(ER.k.six[,p])), as.numeric(table(ER.k.six[,p])), pch = 18, col = "blue", type = "b",
      lty = 2, lwd = 2)
# Add a legend to the plot
legend("top", legend=c("Random Cluster model", "Erdös-Renyi model"),
       col=c("red", "blue"), lty = 1:2, cex=1.2, bty = "n", lwd = 2)
mtext("Frequencies", side = 2, line = 3.2, las = 0, cex = 1.7)
mtext("Components", side = 1, 
      line = 2.5, las = 1, cex = 1.3)


