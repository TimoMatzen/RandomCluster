#########################################################################
## Code for Calculating the KL divergences between the RC and ER model ##
##########             Code by Timo Matzen                ###############
##########    All the functions used are part of the      ###############
## RandomCluster package: https://github.com/TimoMatzen/RandomCluster ###
#########################################################################

library(RandomCluster)
library(igraph)

# Set number of nodes to sample
n.nodes <- 6

# Make matrix with all possible configurations
n.edges <- ed(n.nodes)
mat <- matrix(NA,n.nodes, n.nodes)
mat[upper.tri(mat)] <- seq(1:n.edges)
edge.list<- which(mat > 0, arr.ind = T )
l <- rep(list(0:1), n.edges)
omega <- as.matrix(expand.grid(l))


# Create matrix to save divergences between the two models
kl.div <- rep(NA,length(p.mod))

# Keep p the same for the kldiv
p.mod <- p.true <- seq(0.01,.9, .01)

# Calculate the KLdivergence for all values of the p.grid
for (i in 1:length(p.mod)){
  kl.div[i] <- KLdivergence(omega, n.nodes, p.true[i], p.mod[i])
}

# Set plotting parameters
par(cex.main = 1.3, mar = c(4.5, 6, 4, 7) + 0.1, mgp = c(3, 1, 0), cex.lab = 1.3, 
    font.lab = 2, cex.axis = 1.3, las = 1)
# Plot the KL-divergence
plot(x = p.mod, y= kl.div, xaxt = "n", main = paste0("KL divergence ER/RC with ", n.nodes," nodes") ,
     cex.main = 2, xlab = "", ylab = "", bty = "n")
# Add polygon around Pc
polygon(c(.23, seq(.23, .43,.01), .43), c(0, kl.div[23:43],0), col='grey') 
axis(1, at = seq(0, .9, .05), las = 2)
# Add labels and text
mtext("P value", side = 1, 
      line = 3.5, las = 1, cex = 1.3)
mtext("KL divergence", side = 2, 
      line = 3.5, las = 0, cex = 1.7)
# Add abline for maximum divergence and critical value Pc
segments(x0 = p.mod[which.max(kl.div)], y0 = 0, y1 = kl.div[27], lwd = 2 )
segments(x0 = .33, y0 = 0, y1= kl.div[33], col = "red" , lwd = 2)
legend("topright", legend = c(paste0("Max divergence at p: ", p.mod[which.max(kl.div)]  ), paste0("Critical p-value at ", round(q/n.nodes, 2)) ),
       col = c("black", "red"), 
       lty = 1, bty = "n", lwd = 3, cex = 1.2)






