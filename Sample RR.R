## Simulation with the Randomness Recycler
n.iter <- 5000
n.nodes <- 5
p <- .8
n.edges <- edges(n.nodes)
samp <- matrix(0,nrow = n.iter, ncol = n.edges)
for (i in 1:n.iter) {
  out <- RandomnessRecycler(n.nodes, p, 2)
  samp[i, ] <- out[upper.tri(out)]
  print(i)
}
for (i in 1:n.iter) {
  graph <- erdos.renyi.game(n.nodes,p)
  mat <- matrix(get.adjacency(graph), n.nodes, n.nodes)
  samp[i,] <- mat[upper.tri(mat)]
  print(i)
}


samp2 <- apply(samp, 1, paste, collapse = "")

dat <- data.frame(setup = names(table(samp2)/ n.iter), 
                  prob = as.vector(table(samp2)/ n.iter) )


comp <- matrix(0, nrow = nrow(dat), ncol = 2)
for(i in 1:nrow(dat)) {
    comp[i,1] <- RCProb(as.numeric(strsplit(as.character(dat[i,1]), "")[[1]]), p, 2, n.nodes, Curie = T)
}




# Add sampled probs
comp[,2] <- dat$prob
diff <-  comp[,2] - comp[,1]
plot(x = 1:nrow(dat), y = diff, ylim = c(-.01, .01))
abline(0,0)
sum(diff)

which.max(diff)


