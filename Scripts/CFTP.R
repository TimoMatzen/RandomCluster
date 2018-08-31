#### Coupling From the Past (CFTP) #####
library(qgraph)
library("animation")
library("magick")
# Function for calculating the number of edges
edges <- function(x) {
  # Takes as input the number of nodes
  # Outputs the number of edges
  y <- x*((x-1)/ 2)
  return(y)
}

# Define the transition function
transition <- function(n, s, U, r ,p, q) {
  # Calculate the transition probability for the system
  # n: number of nodes
  # s: current state
  # U: random uniform number on (0,1)
  # r: random edge to switch state
  # p: probability for an edge
  # q: cluster coefficient
  if (s[r] == 1) {
    s2 <- s
    # Change random edge to 0
    s2[r] <- 0
    # Apply transition function
    out <- RC(n, s2, p, q)/(RC(n, s, p, q) + RC(n, s2, p, q))
    # Check whether to accept change
    s[r] <- 1 * !(U <= out)
    
  } else {
    s2 <- s
    # Change random edge to 1
    s2[r] <- 1
    # Apply transition function
    out <- RC(n, s, p, q)/(RC(n, s, p, q) + RC(n, s2, p, q))
    # Check whether to accept change
    s[r] <- 1 * (U > out)
  }
  
  # Return the adjusted state
  return(s)
}

# Define the size of the graph (number of nodes), 
# first with a small network
n.nodes <- 4

# Calculate n.edges
n.edges <- edges(n.nodes)

# Set p and q to a value
p <- .6
# q to 2; ising representation
q <- 2

# Create a sequence to loop over

n.iter <- 1000

# Function for sampling from the random cluster model
SampleRC <- function(n.edges, n.nodes, p, q,  n.iter, verbose = TRUE, adjacency = FALSE) {
  
  if(n.nodes < 5 ) {
    M <- 2^(0:10)
  } else {
    M <- 2^(5:10)
  }
  # Matrix to save samples
  output <- matrix(1,ncol = n.edges, nrow = n.iter)
  for (iter in 1:n.iter) {
    # Set states
    s1 <- rep(0, n.edges)
    s2 <- rep(1, n.edges)
    # Create vector of edge samples to save the edge sample
    e.samps <- vector(mode = "numeric", length = 0)
    # Create a uniform samples vector
    U.samps <- vector(mode = "numeric", length = 0)
    for (n in 1:length(M)) {
      # Go back further in time
      m <- M[n]
      
      # Pick a random sample of the edges for each state
      if (n == 1) {
        # If at it 1 take only 1 sample
        e.samp <- sample(1:n.edges, M[n], replace = TRUE)
        # Pick a random uniform sample
        U.samp <- runif(M[n])
       
      } else {
        # If it != 1 take samples iter[p] - iter[p-1] 
        e.samp <- sample(1:n.edges, (M[n] - M[n-1]), replace = TRUE)
        # Pick a random uniform sample
        U.samp <- runif((M[n] - M[n-1]))
      }
      
      # Append e.samps with the new sample
      e.samps <- append(e.samps, e.samp)
      # Append U.samps with the new sample
      U.samps <- append(U.samps, U.samp)
      
      for (i in m:1) {
        # Make transitions
        s1 <- transition(n.nodes, s1, U.samps[i], e.samps[i], p, q)
        s2 <- transition(n.nodes, s2, U.samps[i], e.samps[i], p, q)
      }
      if (!all(s1 == s2)) {
        # Reset the vectors to the empty/full set
        s1 <- rep(0, n.edges)
        s2 <- rep(1, n.edges)
      } else {
        break
      }
    
    }
    # Save output
    output[iter,] <- s1
    # Track progress
    if (verbose == TRUE) {
      print(iter)
    }
  }
  if (adjacency == FALSE) {
    return(output)
  } else {
    # Create adjacency matrix
    adjacency <- matrix(NA,n.nodes, n.nodes)
    # Set diagonal to 0
    diag(adjacency) <- 0
    # Save setup in matrix
    adjacency[upper.tri(adjacency)] <- s1
    # Make matrix symmetrical
    adjacency <- SymMat(adjacency)
    return(adjacency)
  }
  
}

output <- SampleRC(n.edges , n.nodes, 1000 )


# Get probabilities from CFTP sample
frequency <- as.data.frame(table(apply(output, 1, paste, collapse = "")))
prop <- frequency[,2]/1000
dat <- matrix(0, nrow = nrow(frequency), ncol = n.edges)
for (i in 1:nrow(frequency)) {
  dat[i,] <- as.numeric(strsplit(as.character(frequency[i,1]), "")[[1]])
}
propRC <- apply(dat, 1, RCProb,p = p, q = q, n.nodes = n.nodes, Curie = T )

# Calculate difference with actual model
diff <- prop - propRC
# Plot the difference
plot(x = 1:nrow(dat), y = diff, ylim = c(-.1, .1))

# Sample erdos reny
p <- seq(0.01,.3, .01)
q <- 2
# Define the size of the graph (number of nodes), 
# first with a small network
n.nodes <- 20
# Calculate n.edges
n.edges <- edges(n.nodes)

saveGIF({ 
  for (i in 1:length(p)) {
    par(mfrow = c(1,2))
   # Sample from randomcluster model with different p each tim
    ad <- SampleRC(n.edges, n.nodes, p[i], q, 1, adjacency = TRUE )
    qgraph(ad, title = paste0("Random-Cluster with p = ", p[i], " and q = ", q))
    
    # Sample from randomcluster model with different p each tim
    erdos <- erdos.renyi.game(n.nodes,p[i])
    erdos.mat <- as.matrix(get.adjacency(erdos))
    qgraph(erdos.mat, title = paste0("Erdos-Renyi with p = ", p[i]))
  }
}, movie.name = 'RCvsRG.gif',
ani.width = 800, ani.height = 600)

saveGIF({ 
  for (i in 1:length(p)) {
    # Sample from randomcluster model with different p each tim
    erdos <- erdos.renyi.game(n.nodes,p[i])
    erdos.mat <- as.matrix(get.adjacency(erdos))
    qgraph(erdos.mat, title = paste0("p = ", p[i]))
  }
}, movie.name = 'Erdos.gif',
ani.width = 600, ani.height = 600)



