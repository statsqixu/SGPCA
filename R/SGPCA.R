#' Sparse Group Principal Component Analysis for one single PC estimation
#' @importFrom stats rbeta rnorm
#' @importFrom plotly %>%
#'
#'
#' @param X a data matrix
#' @param group_label a vector for group membership labels
#' @param tau_reg a float to specify the individual thresholding level
#' @param eta_reg a float to specify the group thresholding level
#' @param max_iter a positive integer to specify the maximum iterations to operate
#' @param tol a float to specify the stopping threshold
#'
#' @return A list of with the following components:
#' \describe{
#'  \item{v}{a P-dimensional vector, estimated the first PC}
#'  \item{u}{a N-dimensional vector, scores of the first PC}
#'  \item{d}{a float, variance explained by the first PC}
#' }
#' @export
SGPCA_1pc <- function(X, group_label, tau_reg = 1, eta_reg = 1, max_iter = 100, tol = 1e-5){

  group_label = as.integer(group_label) # convert group labels to integers, from 1 to G

  X = scale(X, center = TRUE, scale = FALSE) # center the data matrix

  N <- nrow(X)
  P <- ncol(X)
  G <- length(unique(group_label))

  v <- rep(0, P) # PC
  u <- rep(0, N) # score of PC projection
  d <- 0 # variance captured by PC

  ## initialization using SPC from PMA package
  spc.tuning <- PMA::SPC.cv(X, nfolds = 3, sumabsvs = seq(1, sqrt(P), length.out = 5), trace = FALSE)
  spc.result <- PMA::SPC(X, K = 1, sumabsv = spc.tuning$bestsumabsv, trace = FALSE)
  init_v <- spc.result$v

  v <- SGPCA_cpp(X, group_label, tau_reg, eta_reg, max_iter, tol, init_v)
  u <- X %*% v / norm(X %*% v, type = "2")
  d <- norm(X %*% v, type = "2")

  output <- list(v = v, u = u, d = d)

  return(output)
}

#' Sparse Group Principal Component Analysis for multiple PCs estimation
#'
#' @param X a data matrix
#' @param group_label a vector for group membership labels
#' @param J number of principal components to be estimated
#' @param tau_reg a float to specify the individual thresholding level
#' @param eta_reg a float to specify the group thresholding level
#' @param max_iter a positive integer to specify the maximum iterations to operate
#' @param tol a float to specify the stopping threshold
#'
#' @return a list of the following components
#' \describe{
#'  \item{V}{a P * J matrix, estimated J PCs}
#'  \item{U}{a N * J matrix, scores on the estimated J PCs}
#'  \item{D}{a J-dimensional vector, variance explained by J PCs}
#' }
#' @export
SGPCA <- function(X, group_label, J = 1, tau_reg = 1, eta_reg = 1,
                  max_iter = 100, tol = 1e-5){

  ## tau_reg, eta_reg specify individual and group thresholding levels
  ## the length of tau_reg and eta_reg should be 1 (uniform for all PCs)
  ## or J (each PC has its own specification)

  if (J > 1){

    if (length(tau_reg) == 1){tau_reg = rep(tau_reg, J)}
    else if (length(tau_reg) != J){stop("The length of tau_reg should be 1 or J")}

    if (length(eta_reg) == 1){eta_reg = rep(eta_reg, J)}
    else if (length(eta_reg) != J){stop("The length of eta_reg should be 1 or J")}
  }

  else{

    if (length(tau_reg) > 1){stop("The length of tau_reg should be 1")}
    if (length(eta_reg) > 1){stop("The length of eta_reg should be 1")}
  }

  group_label = as.integer(group_label) # convert group labels to integers, from 1 to G

  X = scale(X, center = TRUE, scale = FALSE) # center the data matrix

  N <- nrow(X)
  P <- ncol(X)
  G <- length(unique(group_label))

  V <- matrix(0, P, J) # PC
  U <- matrix(0, N, J) # score of PC projection
  D <- rep(0, J) # square root ofvariance captured by PC

  for (j in 1: J){

    pcj_fit <- SGPCA_1pc(X, group_label, tau_reg[j], eta_reg[j], max_iter, tol)
    V[, j] <- pcj_fit$v
    U[, j] <- pcj_fit$u
    D[j] <- pcj_fit$d

    X <- X - D[j] * U[, j]%*% t(V[, j])
  }

  return(list(V = V, U = U, D = D))
}

#' Tuning parameter for SGPCA using resampling-based method and alignment criterion
#'
#' @param X a data matrix
#' @param group_label a vector for group membership labels
#' @param J number of principal components to be estimated
#' @param B number of resampling datasets
#' @param rho resampling ratio
#' @param tau_range search range for tau
#' @param eta_range search range for eta
#' @param max_iter a positive integer to specify the maximum iterations to operate
#' @param tol a float to specify the stopping threshold
#' @param verbose TRUE/FALSE, output algorithm info (TRUE) or not (FALSE)
#' @param mode string, "auto" or "manual". "auto" mode return the tuning parameter pairs that maximize the alignment score. "manual" mode output a figure "alignment v.s. mean support size", user can specify the tuning parameters interactively.
#'
#' @return a list of the following components
#' \describe{
#'  \item{V}{a P * J matrix, estimated J PCs}
#'  \item{U}{a N * J matrix, scores on the estimated J PCs}
#'  \item{D}{a J-dimensional vector, variance explained by J PCs}
#' }
#' @export
SGPCA.rs <- function(X, group_label, J = 1, B = 20, rho = 0.5,
                    tau_range = 10 ^ seq(-2, 2, length.out = 20),
                    eta_range = 10 ^ seq(-2, 2, length.out = 20),
                    max_iter = 100, tol = 1e-5, verbose = TRUE,
                    mode = "auto"){

  group_label = as.integer(group_label) # convert group labels to integers, from 1 to G

  X = scale(X, center = TRUE, scale = FALSE) # center the data matrix

  N <- nrow(X)
  P <- ncol(X)
  G <- length(unique(group_label))

  V <- matrix(0, P, J) # PC
  U <- matrix(0, N, J) # score of PC projection
  D <- rep(0, J) # square root of variance captured by PC

  for (j in 1: J){

    if (verbose){cat("Processing PC", j, "...\n")}

    spc.tuning <- PMA::SPC.cv(X, nfolds = 3, sumabsvs = seq(1, sqrt(P), length.out = 5), trace = FALSE)
    spc.result <- PMA::SPC(X, K = 1, sumabsv = spc.tuning$bestsumabsv, trace = FALSE)
    init_v <- spc.result$v

    if (verbose){cat("Tuning parameters for PC", j, "...\n")}

    set.seed(1234)
    tune_result <- SGPCA_rs_cpp(X, group_label, B, rho, tau_range, eta_range, max_iter, tol, mode, init_v)

    tune_result$param_grid <- expand.grid(tau = tau_range, eta = eta_range)

    tune_df <- data.frame(tau = tune_result$param_grid$tau,
                          eta = tune_result$param_grid$eta,
                          mean_support_size = tune_result$mean_support_size,
                          alignment = tune_result$alignment)

    if (mode == "manual"){
      options(viewer = rstudioapi::viewer)
      plot <- plotly::plot_ly(tune_df) %>%
              plotly::add_trace(x = ~mean_support_size,
                          y = ~alignment,
                          type = "scatter",
                          mode = "markers",
                          name = "Alignment",
                          text = ~paste("tau = ", tau, "\n",
                                        "eta = ", eta, "\n",
                                        "mean_support_size = ", mean_support_size, "\n",
                                        "alignment = ", alignment),
                          hoverinfo = "text") %>%
              plotly::layout(title = "Alignment v.s. mean support size",
                       xaxis = list(title = "Mean Support Size"),
                       yaxis = list(title = "Alignment"))
      print(plot)

      readline(prompt = "Press Enter to continue...")
      repeat{
        tau <- as.numeric(readline(prompt = "Enter the tau value:"))
        if (tau < 0){cat("Invalid tau value. Please enter a positive number.\n")}
        else{break}
      }
      repeat{
        eta <- as.numeric(readline(prompt = "Enter the eta value:"))
        if (eta < 0){cat("Invalid eta value. Please enter a positive number.\n")}
        else{break}
      }
    } else if (mode == "auto"){
      best_idx <- which.max(tune_result$alignment)
      tau <- tune_result$param_grid$tau[best_idx]
      eta <- tune_result$param_grid$eta[best_idx]
      if (verbose){
        cat("Best parameters found:\n")
        cat("tau = ", tau, "\n")
        cat("eta = ", eta, "\n")
        cat("Alignment = ", tune_df$alignment[best_idx], "\n")
        cat("Mean Support Size = ", tune_df$mean_support_size[best_idx], "\n")
      }
    }

    if (verbose){cat("Estimating PC", j, "with selected parameters...\n")}

    pcj_fit <- SGPCA_1pc(X, group_label, tau * sqrt(rho), eta * sqrt(rho), max_iter, tol)
    V[, j] <- pcj_fit$v
    U[, j] <- pcj_fit$u
    D[j] <- pcj_fit$d

    X <- X - D[j] * U[, j] %*% t(V[, j])

    if (verbose){cat("Estimating PC", j, "completed!\n")}
  }

  return(list(V = V, U = U, D = D))
}
