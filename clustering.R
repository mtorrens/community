################################################################################
# Script   : clustering.R
# Descrip. : Research on clustering algorithms
################################################################################
# Author   : (c) Miquel Torrens, 2016.06.19
# Modified :     -
################################################################################
# source('~/Desktop/community/clustering.R')
################################################################################

################################################################################
# Parameters
rounded <- TRUE

# Dependencies
require(igraph)
require(mvtnorm)

# Load sample matrices
source('~/Desktop/community/functions.R')
source('~/Desktop/community/matrices.R')
gfam <- graph_from_adjacency_matrix

# Matrices
set.seed(666)
X0 <- mvtnorm::rmvnorm(mean = rep(0, 9), n = 1e5, sigma = W0)
X1 <- mvtnorm::rmvnorm(mean = rep(0, 9), n = 1e5, sigma = W1)
X2 <- mvtnorm::rmvnorm(mean = rep(0, 9), n = 1e5, sigma = W2)
X3 <- mvtnorm::rmvnorm(mean = rep(0, 9), n = 1e5, sigma = W3)
X4 <- mvtnorm::rmvnorm(mean = rep(0, 9), n = 1e5, sigma = W4)
X5 <- mvtnorm::rmvnorm(mean = rep(0, 9), n = 1e5, sigma = W5)
X6 <- mvtnorm::rmvnorm(mean = rep(0, 18), n = 1e5, sigma = W6)
X7 <- mvtnorm::rmvnorm(mean = rep(0, 18), n = 1e5, sigma = W7)
X8 <- mvtnorm::rmvnorm(mean = rep(0, 21), n = 1e5, sigma = W8)

# Different theoretical memberships
mem1 <- c(1, 1, 1, 2, 2, 2, 3, 3, 3)
mem2 <- c(1, 1, 1, 2, 2, 2, 3, 4, 5)
mem3 <- c(1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 6, 7, 7, 7)
mem4 <- c(1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 3, 3, 3, 2, 2, 2, 1, 1, 1)
mem5 <- c(1, 1, 1, 2:7)
mem6 <- c(rep(1, 6), rep(2, 6), rep(3, 6))
################################################################################

################################################################################
# Community detection
# Choose if the sample correlations should be rounded
if (rounded == TRUE) {
  # Graph objects
  w <- TRUE
  m <- 'undirected'
  G0 <- gfam(abs(cor(X0)), weighted = w, mode = m)
  G1 <- gfam(abs(cor(X1)), weighted = w, mode = m)
  G2 <- gfam(abs(cor(X2)), weighted = w, mode = m)
  G3 <- gfam(abs(cor(X3)), weighted = w, mode = m)
  G4 <- gfam(abs(cor(X4)), weighted = w, mode = m)
  G5 <- gfam(abs(cor(X5)), weighted = w, mode = m)
  G6 <- gfam(abs(cor(X6)), weighted = w, mode = m)
  G7 <- gfam(abs(cor(X7)), weighted = w, mode = m)
  G8 <- gfam(abs(cor(X8)), weighted = w, mode = m)

  # Try several methods of community detection
  repeat {
    aux <- try(comm.detection(X0, G0, 8, short = TRUE, truth = 1:9),
               silent = TRUE)
    if (class(aux) != 'try-error') { break }
  }
  repeat {
    aux <- try(comm.detection(X1, G1, 3, short = TRUE, truth = rep(1, 9)),
               silent = TRUE)
    if (class(aux) != 'try-error') { break }
  }
  comm.detection(X2, G2, 3, short = TRUE, truth = mem1)
  comm.detection(X3, G3, 3, short = TRUE, truth = mem1)
  comm.detection(X4, G4, 3, 5, short = TRUE, truth = mem1, optb = mem2)
  comm.detection(X5, G5, 3, 7, short = TRUE, truth = mem1, optb = mem5)
  comm.detection(X6, G6, 3, short = TRUE, truth = mem6)
  comm.detection(X7, G7, 3, short = TRUE, truth = mem6)
  comm.detection(X8, G8, 4, 7, short = TRUE, truth = mem3, optb = mem4)
} else {
  # Second try (rounded numbers have an effect)
  G0 <- gfam(round(abs(cor(X0)), 2), weighted = w, mode = m)
  G1 <- gfam(round(abs(cor(X1)), 2), weighted = w, mode = m)
  G2 <- gfam(round(abs(cor(X2)), 2), weighted = w, mode = m)
  G3 <- gfam(round(abs(cor(X3)), 2), weighted = w, mode = m)
  G4 <- gfam(round(abs(cor(X4)), 2), weighted = w, mode = m)
  G5 <- gfam(round(abs(cor(X5)), 2), weighted = w, mode = m)
  G6 <- gfam(round(abs(cor(X6)), 2), weighted = w, mode = m)
  G7 <- gfam(round(abs(cor(X7)), 2), weighted = w, mode = m)
  G8 <- gfam(round(abs(cor(X8)), 2), weighted = w, mode = m)

  comm.detection(X2, G2, 3, short = TRUE, truth = mem1)
  comm.detection(X3, G3, 3, short = TRUE, truth = mem1)
  comm.detection(X4, G4, 3, 5, short = TRUE, truth = mem1, optb = mem2)
  comm.detection(X5, G5, 3, 7, short = TRUE, truth = mem1, optb = mem5)
  comm.detection(X6, G6, 3, short = TRUE, truth = mem6)
  comm.detection(X7, G7, 3, short = TRUE, truth = mem6)
  comm.detection(X8, G8, 4, 7, short = TRUE, truth = mem3, optb = mem4)
}
################################################################################
# END OF SCRIPT
