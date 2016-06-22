################################################################################
# Script   : clustering.R
# Descrip. : Research on clustering algorithms
################################################################################
# Author   : (c) Miquel Torrens, 2016.06.19
# Modified :     -
################################################################################

# Dependencies
require(igraph)
require(mvtnorm)

# Load sample matrices
source('~/Desktop/community/matrices.R')
source('~/Desktop/community/functions.R')

# Matrices
set.seed(666)
XI <- mvtnorm::rmvnorm(mean = rep(0, 9), n = 1e5, sigma = WI)
X0 <- mvtnorm::rmvnorm(mean = rep(0, 9), n = 1e5, sigma = W0)
X1 <- mvtnorm::rmvnorm(mean = rep(0, 9), n = 1e5, sigma = W1)
X2 <- mvtnorm::rmvnorm(mean = rep(0, 9), n = 1e5, sigma = W2)
X3 <- mvtnorm::rmvnorm(mean = rep(0, 9), n = 1e5, sigma = W3)
X4 <- mvtnorm::rmvnorm(mean = rep(0, 21), n = 1e5, sigma = W4)
X5 <- mvtnorm::rmvnorm(mean = rep(0, 9), n = 1e5, sigma = W5)
X6 <- mvtnorm::rmvnorm(mean = rep(0, 18), n = 1e5, sigma = W6)
X7 <- mvtnorm::rmvnorm(mean = rep(0, 18), n = 1e5, sigma = W7)

# Graph objects
w <- TRUE
m <- 'undirected'
GI <- graph_from_adjacency_matrix(abs(cor(XI)), weighted = w, mode = m)
G0 <- graph_from_adjacency_matrix(abs(cor(X0)), weighted = w, mode = m)
G1 <- graph_from_adjacency_matrix(abs(cor(X1)), weighted = w, mode = m)
G2 <- graph_from_adjacency_matrix(abs(cor(X2)), weighted = w, mode = m)
G3 <- graph_from_adjacency_matrix(abs(cor(X3)), weighted = w, mode = m)
G4 <- graph_from_adjacency_matrix(abs(cor(X4)), weighted = w, mode = m)
G5 <- graph_from_adjacency_matrix(abs(cor(X5)), weighted = w, mode = m)
G6 <- graph_from_adjacency_matrix(abs(cor(X6)), weighted = w, mode = m)
G7 <- graph_from_adjacency_matrix(abs(cor(X7)), weighted = w, mode = m)

# Try several methods of community detection
repeat {
  aux <- try(comm.detection(XI, GI, 8, short = TRUE), silent = TRUE)
  if (class(aux) != 'try-error') { break }
}
repeat {
  aux <- try(comm.detection(X0, G0, 3, short = TRUE), silent = TRUE)
  if (class(aux) != 'try-error') { break }
}
comm.detection(X1, G1, 3, short = TRUE)
comm.detection(X2, G2, 3, 5, short = TRUE)
comm.detection(X3, G3, 3, short = TRUE)
comm.detection(X4, G4, 4, 7, short = TRUE)
comm.detection(X5, G5, 3, 7, short = TRUE)
comm.detection(X6, G6, 3, short = TRUE)
comm.detection(X7, G7, 3, 7, short = TRUE)

# Second try (rounded numbers have an effect)
GI <- graph_from_adjacency_matrix(round(abs(cov(XI)), 2), weighted = w, mode = m)
G0 <- graph_from_adjacency_matrix(round(abs(cov(X0)), 2), weighted = w, mode = m)
G1 <- graph_from_adjacency_matrix(round(abs(cov(X1)), 2), weighted = w, mode = m)
G2 <- graph_from_adjacency_matrix(round(abs(cov(X2)), 2), weighted = w, mode = m)
G3 <- graph_from_adjacency_matrix(round(abs(cov(X3)), 2), weighted = w, mode = m)
G4 <- graph_from_adjacency_matrix(round(abs(cov(X4)), 2), weighted = w, mode = m)
G5 <- graph_from_adjacency_matrix(round(abs(cov(X5)), 2), weighted = w, mode = m)
G6 <- graph_from_adjacency_matrix(round(abs(cov(X6)), 2), weighted = w, mode = m)
G7 <- graph_from_adjacency_matrix(round(abs(cov(X7)), 2), weighted = w, mode = m)

repeat {
  aux <- try(comm.detection(XI, GI, 8), silent = TRUE)
  if (class(aux) != 'try-error') { break }
}
repeat {
  aux <- try(comm.detection(X0, G0, 3), silent = TRUE)
  if (class(aux) != 'try-error') { break }
}
comm.detection(X1, G1, 3, short = TRUE)
comm.detection(X2, G2, 3, short = TRUE)
comm.detection(X3, G3, 3, short = TRUE)
comm.detection(X4, G4, 4, 7, short = TRUE)
comm.detection(X5, G5, 3, 7, short = TRUE)
comm.detection(X6, G6, 3, short = TRUE)
comm.detection(X7, G7, 3, 7, short = TRUE)














comm.detection(X4, G4, 4, 7)
comm.detection(X4, G4, 4, 7, short = TRUE)

comm.detection(X4, G4, 4, nneighbors = 5)






