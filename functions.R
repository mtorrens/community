################################################################################
# Sort the indexes of the clusters resulting from k-means
sort.clusters <- function(clusts, desc = FALSE) {
################################################################################
  noms <- unique(clusts)
  if (desc == TRUE) {
    idx <- length(noms):1
  } else {
    idx <- 1:length(noms)
  }
  new.clust <- rep(NA, length(clusts))
  for (n in noms) {
    new.clust[which(clusts == n)] <- idx[which(noms == n)]
  }
  return(new.clust)
}

################################################################################
# Spectral clustering for predictors
sclust <- function(X, k, method = 'unnormalised', similarity.method = 'covm',
                   center = TRUE, scale = TRUE, empty.diag = TRUE,
                   nneighbors = NULL, epsilon = NULL, sort.result = TRUE,
                   seed = 666, ...) {
################################################################################
  # Center the matrix
  X <- scale(X, center = center, scale = scale)

  # Compute a similarity matrix
  S <- abs(t(X) %*% X)

  # Compute the adjacency graph according to some similarity criterion
  if (similarity.method == 'covm') {
    W <- matrix(rep(1, nrow(S) ** 2), nrow = nrow(S))
  } else if (similarity.method == 'knn') {
    if (is.null(nneighbors)) {
      nneighbors <- ceiling(log(ncol(X)))
    }
    nns <- apply(S, 2, function(x) { order(x)[1:nneighbors] })
    W <- matrix(rep(0, ncol(S) ** 2), nrow = ncol(S))
    for (j in 1:ncol(W)) {
      W[as.numeric(nns[, j]), j] <- 1
    }
    W <- sign(W + t(W))  # Make the graph undirected (not MUTUAL)
  } else if (similarity.method == 'e-neigh') {
    if (is.null(epsilon)) {
      #epsilon <- quantile(S, 0.8)
      n <- nrow(X)
      d <- ncol(X)
      epsilon <- log(n) ** d / n
    }
    W <- matrix(rep(0, ncol(S) ** 2), nrow = ncol(S))
    for (j in 1:ncol(W)) {
      W[which(S[, j] > epsilon), j] <- 1
    }
  }

  # Final matrix for the algorithm
  W <- W * S

  # Elements of the model
  ones <- matrix(rep(1, ncol(W)), ncol = 1)
  d <- as.numeric(W %*% ones)
  D <- diag(d)

  # Compute the Laplacian matrix
  if (method == 'unnormalised') {
    L <- D - W
  } else if (method == 'shi') {
    L <- diag(d ** (-1)) %*% (D - W)
  } else {
    L <- diag(d ** (-1 / 2)) %*% W %*% diag(d ** (-1 / 2))
  }

  # Eigen decomposition of the Laplacian
  ev <- eigen(L)
  evecs <- ev$vectors  # Last eigenvector always (multiple of) 1 vector
  evals <- ev$values  # Last eigenvalue always = 0

  # Compute kmeans on last k eigenvectors
  idx <- (length(evals) - k + 1):length(evals)  # Luxburg
  if (method == 'ng') {
    # Build Y with normalized rows of X
    Z <- evecs[, idx]
    Y <- Z
    for (j in 1:ncol(Z)) {
      vlength <- sum(Z[j, ] ** 2) ** (1 / 2)
      Y[j, ] <- Z[j, ] / vlength
    }
  } else {
    Y <- evecs[, idx]
  }

  # Perform k-means
  set.seed(seed)
  clust <- stats::kmeans(Y, centers = k)

  # End
  if (sort.result == TRUE) {
    return(sort.clusters(clust[['cluster']]))
  } else {
    return(clust[['cluster']])
  }
}

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

