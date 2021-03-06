################################################################################
# Script   : functions.R
# Descrip. : Community detection functions
################################################################################
# Author   : (c) Miquel Torrens, 2016.06.19
# Modified :     Miquel Torrens, 2016.06.25
################################################################################
# source('~/Desktop/community/functions.R')
################################################################################

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
      if (epsilon > max(S) | epsilon < min(S)) {
        epsilon <- quantile(S, 1 / k)
      }
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
    idx <- 1:k
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
  Y[is.nan(Y)] <- 0

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
zhang.newman <- function(A, k, max.iter = 1e3, no.empty = TRUE, verbose = TRUE,
                         seed = 666) {
################################################################################
  # Construction of B matrix
  d <- A %*% rep(1, nrow(A))
  m <- sum(A) / 2
  B <- matrix(rep(NA, nrow(A) ** 2), nrow = nrow(A), ncol = nrow(A))
  for (i in 1:nrow(B)) {
    for (j in 1:ncol(B)) {
      B[i, j] <- A[i, j] - (d[i] * d[j]) / (2 * m)
    }
  }

  # Eigendecomposition
  U <- eigen(B)$vectors
  l <- eigen(B)$values

  # Rank p approximation
  U <- U[, l > 0]
  l <- l[l > 0]

  # R vectors
  R <- matrix(rep(NA, nrow(U) * ncol(U)), nrow = nrow(U), ncol = ncol(U))
  for (i in 1:nrow(R)) {
    for (el in 1:ncol(R)) {
      R[i, el] <- sqrt(l[el]) * U[i, el]
    }
  }
  R <- t(R)

  # Random initialisation
  split <- nrow(A) %/% k
  remain <- nrow(A) %% k
  set.seed(seed)
  pool <- c(as.numeric(sapply(1:k, rep, split)), sample(1:k, remain))  
  group <- sample(pool, length(pool), replace = FALSE)

  # Loop
  iter <- 1
  repeat {
    new.group <- rep(0, nrow(A))
    for (n in 1:ncol(R)) {
      mine <- group[n]
      scores <- c()
      for (w in 1:k) {
        aux <- cbind(R[, which(group == w)], rep(0, nrow(R)))
        if (w == mine) {
          scores <- c(scores, t(rowSums(aux) - R[, n]) %*% R[, n])
        } else {
          scores <- c(scores, t(rowSums(aux)) %*% R[, n])
        }
      }
      new.group[n] <- which.max(scores)
    }

    # Do not leave any cluster behind
    while (no.empty == TRUE & length(unique(new.group)) != k) {
      alone <- which(! 1:k %in% unique(new.group))
      for (a in alone) {
        new.group[which(group == a)[1]] <- a
      }
    }

    # See if results have changed
    if (identical(group, new.group)) {
      group <- new.group
      break
    } else {
      group <- new.group
      if (iter >= max.iter) {
        if (verbose == TRUE) {
          cat('Warning: The algorithm did not converge.\n')
        }
        break
      }
      iter <- iter + 1
    }
  }

  if (verbose == TRUE) {
    cat('Number of iterations:', iter, '\n')
  }
  return(sort.clusters(group))
}

################################################################################
comm.detection <- function(X, G, k1, k2 = NULL, short = FALSE, truth = NULL,
                           optb = NULL, ...) {
################################################################################
  # Algorithm names
  algs <- c(paste('SC Unnormalised (k = ', k2,')', sep = ''),
            paste('SC Shi-Malik (k = ', k2,')', sep = ''),
            paste('SC Ng-Weiss-Jordan (k = ', k2,')', sep = ''),
            paste('SC Unnormalised (k = ', k1,')', sep = ''),
            paste('SC Shi-Malik (k = ', k1,')', sep = ''),
            paste('SC Ng-Weiss-Jordan (k = ', k1,')', sep = ''),
            'Girvan-Newman', 'Fast Greedy', 'Infomap', 'Label propagation',
            'Modularity maximisation', 'Louvain', 'Spinglass', 'Walktrap',
            'Optimal', paste('SC Unnormalised KNN (k = ', k1,')', sep = ''),
            paste('SC Unnormalised Eps-N (k = ', k1,')', sep = ''),
            paste('SC Shi-Malik KNN (k = ', k1,')', sep = ''),
            paste('SC Shi-Malik Eps-N (k = ', k1,')', sep = ''),
            paste('SC Ng-Weiss-Jordan KNN (k = ', k1,')', sep = ''),
            paste('SC Ng-Weiss-Jordan Eps-N (k = ', k1, ')', sep = ''),
            paste('Zhang-Newman (k = ', k1, ')', sep = ''),
            paste('Zhang-Newman (Forced k = ', k1, ')', sep = ''),
            paste('Zhang-Newman (k = ', k2, ')', sep = ''),
            paste('Zhang-Newman (Forced k = ', k2, ')', sep = ''))

  # Spectral clustering
  if (! is.null(k2)) {
    res01 <- sclust(X, k = k2, method = 'unnormalised')
    res02 <- sclust(X, k = k2, method = 'shi')
    res03 <- sclust(X, k = k2, method = 'ng')
  }
  res04 <- sclust(X, k = k1, method = 'unnormalised')
  res05 <- sclust(X, k = k1, method = 'shi')
  res06 <- sclust(X, k = k1, method = 'ng')

  if (short == FALSE) {
    res16 <- sclust(X, k1, method = 'unnormalised', similarity.method = 'knn', ...)
    res17 <- sclust(X, k1, method = 'unnormalised', similarity.method = 'e-neigh', ...)
    res18 <- sclust(X, k1, method = 'shi', similarity.method = 'knn', ...)
    res19 <- sclust(X, k1, method = 'shi', similarity.method = 'e-neigh', ...)
    res20 <- sclust(X, k1, method = 'ng', similarity.method = 'knn', ...)
    res21 <- sclust(X, k1, method = 'ng', similarity.method = 'e-neigh', ...)
  }

  # Zhang-Newman algorithm
  S <- round(cor(X), 2)
  res22 <- zhang.newman(A = S, k = k1, verbose = FALSE, no.empty = FALSE)
  res23 <- zhang.newman(A = S, k = k1, verbose = FALSE, no.empty = TRUE)
  if (! is.null(k2)) {
    res24 <- zhang.newman(A = S, k = k2, verbose = FALSE, no.empty = FALSE)
    res25 <- zhang.newman(A = S, k = k2, verbose = FALSE, no.empty = TRUE)
  }

  # Algorithms from igraph
  idxs <- c(7, 11)
  clu07 <- cluster_edge_betweenness(G)  # Girvan-Newman
  clu11 <- try(cluster_leading_eigen(G))  # Modularity maximisation
  if (short == FALSE) {
    idxs <- 7:15
    clu08 <- cluster_fast_greedy(G)
    clu09 <- cluster_infomap(G)
    clu10 <- cluster_label_prop(G)
    clu12 <- cluster_louvain(G)
    clu13 <- cluster_spinglass(G)
    clu14 <- cluster_walktrap(G)
    clu15 <- cluster_optimal(G)
  }
  for (i in sprintf('%02.0f', idxs)) {
    assign(paste('res', i, sep = ''),
           membership(get(paste('clu', i, sep = ''))))
    assign(paste('com', i, sep = ''),
           length(unique(get(paste('res', i, sep = '')))))
  }

  # Print results
  idxs <- unique(c(ifelse(rep(! is.null(k2), 3), 1:3, rep(4, 3)), 4:23))
  if (! is.null(k2)) { idxs <- c(idxs, 24:25) }
  if (short) {
    idxs <- setdiff(idxs, c(8:10, 12:21))
  }
  for (i in sprintf('%02.0f', idxs)) {
    id <- as.numeric(i)
    mod <- modularity(G, membership = get(paste('res', i, sep = '')))
    assign(paste('mod', i, sep = ''), mod)
    cat('Method ', i, ': ', round(mod, 4), ' [', algs[id], ']\n', sep = '')
    cat(sort.clusters(get(paste('res', i, sep = ''))), '\n')
  }
  if (! is.null(truth)) {
    mod <- modularity(G, membership = truth)
    cat('TRUTH:', round(mod, 4), '\n')
    cat(sort.clusters(truth), '\n')
  }
  if (! is.null(optb)) {
    mod <- modularity(G, membership = optb)
    cat('TRUTH (OPT. B):', round(mod, 4), '\n')
    cat(sort.clusters(optb), '\n')
  }
}
# END OF SCRIPT
