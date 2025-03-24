#' Generator of simulation setting 1
#'
#' @param n sample size
#' @param G number of genes
#' @param C group size (e.g., number of celltypes)
#' @param seed an integer to specify the random seed
#'
#' @return
#' \describe{
#'  \item{X}{data matrix of dimension n by (GC)}
#'  \item{group_label}{a (GC)-dimensional vector, group membership}
#'  \item{signal_indices}{a vector, provides the indices of signals in PC}
#'  \item{pc1}{a (GC)-dimensional vector, true PC1 by design}
#' }
#' @export
simulator1 <- function(n = 500, G = 2000, C = 3, seed = 123){

  set.seed(seed)

  group_label <- rep(1: G, each = C)

  pc1 <- rep(0, G * C)

  group_indices <- c(1: floor(C * G * 0.05))

  signal_indices <- sample(group_indices, floor(0.8 * length(group_indices)), replace = F)

  pc1[signal_indices] <- rbeta(length(signal_indices), 2, 2)

  pc1 <- pc1 / norm(pc1, type = "2")

  rho <- 5

  z <- rnorm(n)

  X <- rho * z %*% t(pc1) + matrix(rnorm(n * G * C), nrow = n, ncol = G * C)

  return(list(X = X, group_label = group_label, signal_indices = signal_indices, pc1 = pc1))

}


#' Generator of simulation setting 2
#'
#' @param n sample size
#' @param G number of genes
#' @param cond sprsity condition (1: group sparsity, 2: individual sparsity, 3: group and individual sparsity)
#' @param seed an integer to specify the random seed
#'
#' @return
#' \describe{
#'  \item{X}{data matrix of dimension n by (10G)}
#'  \item{group_label}{a (10G)-dimensional vector, group membership}
#'  \item{signal_indices}{a vector, provides the indices of signals in PC}
#'  \item{pc1}{a (10G)-dimensional vector, true PC1 by design}
#' }
#' @export
simulator2 <- function(n = 500, G = 2000, cond = 3, seed = 123){

  set.seed(seed)

  C <- 10

  group_label <- rep(1: G, each = C)

  pc1 <- rep(0, G * C)

  if (cond == 1){
    signal_indices <- sample(1: G * C, floor(0.05 * G * C))
  } else if (cond == 2){
    signal_indices <- c(1: floor(0.05 * G * C))
  } else if (cond == 3){
    group_indices <- c(1: floor(0.05 * G * C))
    signal_indices <- sample(group_indices, 0.8 * length(group_indices))
  }

  pc1[signal_indices] <- rbeta(length(signal_indices), 2, 2)

  pc1 <- pc1 / norm(pc1, type = "2")

  rho <- 5

  z <- rnorm(n)

  X <- rho * z %*% t(pc1) + matrix(rnorm(n * G * C), nrow = n, ncol = G * C)

  return(list(X = X, group_label = group_label, signal_indices = signal_indices, pc1 = pc1))
}

#' Generator of simulation setting 3
#'
#' @param n sample size
#' @param G number of genes
#' @param seed an integer to specify the random seed
#'
#' @return
#' \describe{
#'  \item{X}{data matrix of dimension n by (10G)}
#'  \item{group_label}{a (10G)-dimensional vector, group membership}
#'  \item{signal_indices1}{a vector, provides the indices of signals in PC1}
#'  \item{signal_indices2}{a vector, provides the indices of signals in PC2}
#'  \item{signal_indices3}{a vector, provides the indices of signals in PC3}
#'  \item{pc1}{a (10G)-dimensional vector, true PC1 by design}
#'  \item{pc2}{a (10G)-dimensional vector, true PC2 by design}
#'  \item{pc3}{a (10G)-dimensional vector, true PC3 by design}
#' }
#' @export
simulator3 <- function(n = 500, G = 2000, seed = 123){

  set.seed(seed)

  C = 10

  group_label <- rep(1: G, each = C)

  pc1 <- rep(0, G * C)
  pc2 <- rep(0, G * C)
  pc3 <- rep(0, G * C)

  group_indices1 <- c(1: floor(C * G * 0.05))
  group_indices2 <- c((floor(C * G * 0.05) + 1): floor(C * G * 0.1))
  group_indices3 <- c((floor(C * G * 0.15) + 1): floor(C * G * 0.2))

  signal_indices1 <- sample(group_indices1, floor(0.8 * length(group_indices1)))
  signal_indices2 <- sample(group_indices2, floor(0.8 * length(group_indices2)))
  signal_indices3 <- sample(group_indices3, floor(0.8 * length(group_indices3)))

  pc1[signal_indices1] <- rbeta(length(signal_indices1), 2, 2)
  pc2[signal_indices2] <- rbeta(length(signal_indices2), 2, 2)
  pc3[signal_indices3] <- rbeta(length(signal_indices3), 2, 2)

  pc1 <- pc1 / norm(pc1, type = "2")
  pc2 <- pc2 / norm(pc2, type = "2")
  pc3 <- pc3 / norm(pc3, type = "2")

  rho1 <- 20
  rho2 <- 5
  rho3 <- 3

  z <- matrix(rnorm(n * 3), nrow = n)

  X <- z %*% diag(c(rho1, rho2, rho3)) %*% t(cbind(pc1, pc2, pc3))

  return(list(X = X, group_label = group_label, signal_indices1 = signal_indices1,
              signal_indices2 = signal_indices2, signal_indices3 = signal_indices3,
              pc1 = pc1, pc2 = pc2, pc3 = pc3))

}
