################################################################################
# Script   : clustering.R
# Descrip. : Research on clustering algorithms
################################################################################
# Author   : (c) Miquel Torrens, 2016.06.19
# Modified :     -
################################################################################

################################################################################
# Dependencies
require(igraph)
require(mvtnorm)

# Load sample matrices
source('~/Desktop/matrices.R')
source('~/Desktop/functions.R')
################################################################################

# Matrices
set.seed(666)
X0 <- mvtnorm::rmvnorm(mean = rep(0, 9), n = 1e5, sigma = W0)
X1 <- mvtnorm::rmvnorm(mean = rep(0, 9), n = 1e5, sigma = W1)
X4 <- mvtnorm::rmvnorm(mean = rep(0, 21), n = 1e5, sigma = W4)

# Try several methods of community detection
comm.detection(X4, G4, 4, 7)

################################################################################
comm.detection <- function(X, G, k1, k2 = NULL) {
################################################################################
  # Algorithm names
  algs <- c(paste('SC Unnormalised (k = ', k2,')', sep = ''),
            paste('SC Shi-Malik (k = ', k2,')', sep = ''),
            paste('SC Ng-Weiss-Joran (k = ', k2,')', sep = ''),
            paste('SC Unnormalised (k = ', k1,')', sep = ''),
            paste('SC Shi-Malik (k = ', k1,')', sep = ''),
            paste('SC Ng-Weiss-Joran (k = ', k1,')', sep = ''),
            'Girvan-Newman', 'Fast Greedy', 'Infomap', 'Label propagation',
            'Modularity maximisation', 'Louvain', 'Spinglass', 'Walktrap',
            'Optimal')

  # Spectral clustering
  if (! is.null(k2)) {
    res01 <- sclust(X, k = k2, method = 'unnormalised')
    res02 <- sclust(X, k = k2, method = 'shi')
    res03 <- sclust(X, k = k2, method = 'ng')
  }
  res04 <- sclust(X, k = k1, method = 'unnormalised')
  res05 <- sclust(X, k = k1, method = 'shi')
  res06 <- sclust(X, k = k1, method = 'ng')

  # Algorithms from igraph
  clu07 <- cluster_edge_betweenness(G)  # Girvan-Newman
  clu08 <- cluster_fast_greedy(G)
  clu09 <- cluster_infomap(G)
  clu10 <- cluster_label_prop(G)
  clu11 <- cluster_leading_eigen(G)  # Modularity maximisation
  clu12 <- cluster_louvain(G)
  clu13 <- cluster_spinglass(G)
  clu14 <- cluster_walktrap(G)
  clu15 <- cluster_optimal(G)
  for (i in sprintf('%02.0f', 7:15)) {
    assign(paste('res', i, sep = ''),
           membership(get(paste('clu', i, sep = ''))))
    assign(paste('com', i, sep = ''),
           length(unique(get(paste('res', i, sep = '')))))
  }

  # Print results
  idxs <- unique(c(ifelse(rep(! is.null(k2), 3), 1:3, rep(4, 3)), 4:15))
  for (i in sprintf('%02.0f', idxs)) {
    id <- as.numeric(i)
    mod <- modularity(G4, membership = get(paste('res', i, sep = '')))
    assign(paste('mod', i, sep = ''), mod)
    cat('Method ', i, ': ', round(mod, 4), ' [', algs[id], ']\n', sep = '')
    cat(sort.clusters(get(paste('res', i, sep = ''))), '\n')
  }
}





res01 <- sclust(X4, k = 7, method = 'unnormalised')
res02 <- sclust(X4, k = 7, method = 'shi')
res03 <- sclust(X4, k = 7, method = 'ng')
res04 <- sclust(X4, k = 4, method = 'unnormalised')
res05 <- sclust(X4, k = 4, method = 'shi')
res06 <- sclust(X4, k = 4, method = 'ng')

clu07 <- cluster_edge_betweenness(G4)  # Girvan-Newman
clu08 <- cluster_fast_greedy(G4)  # Needs undirected graph
clu09 <- cluster_infomap(G4)
clu10 <- cluster_label_prop(G4)
clu11 <- cluster_leading_eigen(G4)
clu12 <- cluster_louvain(G4)
clu13 <- cluster_spinglass(G4)
clu14 <- cluster_walktrap(G4)
clu15 <- cluster_optimal(G4)
for (i in sprintf('%02.0f', 7:15)) {
  assign(paste('res', i, sep = ''),
         membership(get(paste('clu', i, sep = ''))))
  assign(paste('com', i, sep = ''),
         length(unique(get(paste('res', i, sep = '')))))
}


compare(res13, res14)


for (i in sprintf('%02.0f', 1:15)) {
  mod <- modularity(G4, membership = get(paste('res', i, sep = '')))
  #mod <- modularity(G4, membership = get(paste('res', i, sep = '')), weights = E(G4)$weights)
  assign(paste('mod', i, sep = ''), mod)
  cat('Method ', i, ': ', round(mod, 4), '\n', sep = '')
  cat(sort.clusters(get(paste('res', i, sep = ''))), '\n')
}

# modularity(G4, membership = res01)
# modularity(G4, membership = res02)
# modularity(G4, membership = res03)
# modularity(G4, membership = res04)
# modularity(G4, membership = res05)
# modularity(G4, membership = res06)
# modularity(G4, membership = res07)

