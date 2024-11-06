simulateData <- function(n,  # number of subjects (scalar)
                         SS, # number of functional time points (scalar)
                         TT  # number of longitudinal time points (scalar)
){
  ###################################################################################################
  ## Description: Function for generating densely observed longitudinal functional data as described 
  ##              in Supplementary Materials Appendix C. In particular, the simulated data contain  
  ##              no outliers (Case 1) and are generated from the standard product FPCA model detailed 
  ##              in Section 2.1.
  ## Args:        (see above)
  ## Returns:     List, length 2
  ##               y: simulated functional data list of length n
  ##                   each list contains a vector of Y_i(s, t) (list, (SS x TT) x 1 x n)
  ##               truth_funcs: truth function list of length 8
  ##                   mu: true mean function (vector, (SS x TT) x 1)
  ##                   psi: list of 2
  ##                     1st true longitudinal eigenfunction (vector, SS x 1)
  ##                     2nd true longitudinal eigenfunction (vector, SS x 1)
  ##                   phi: list of 2
  ##                     1st true functional eigenfunction (vector, TT x 1)
  ##                     2nd true functional eigenfunction (vector, TT x 1)
  ##                   tau: true longitudinal eigenvalues (vector, J x 1)
  ##                   vartheta: true functional eigenvalues (vector, K x 1)
  ##                   K_s: true marginal longitudinal covariance (matrix, SS x SS)
  ##                   K_s: true marginal functional covariance (matrix, TT x TT)
  ##                   K: true covariance kernel (matrix, (SSxTT) x (SSxTT))
  ## simulateData Outline:
  ##              1. Define model components
  ##              2. Generate longitudinal functional data
  #############################################################################
  
  #############################################################################
  # 1. Define model components
  #############################################################################
  # Install missing packages
  list.of.packages <- c("Rcpp", "pracma", "LFBayes")
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if (length(new.packages)) install.packages(new.packages) 
  
  # Load libraries
  library(Rcpp); library(pracma); library(LFBayes)
  
  # Longitudinal and functional time points
  t <- seq(from = 0, to = 1, length.out = TT) # vector, TT x 1
  s <- seq(from = 0, to = 1, length.out = SS) # vector, SS x 1
  
  # Longitudinal and functional eigen-components
  J <- 2; K <- 2
  
  # Mean function to generate components for the mean function 
   ## Input: longitudinal time points s (vector, TT x 1)
   ##        functional time points s (vector, SS x 1)  
   ## Output: a vector of mean function with dimension (TT x SS) x 1
  mu <- function(s, t){
    mu <- matrix(0, nrow = length(t), ncol = length(s))
    for(i in 1:length(t)){
      for(j in 1:length(s)){
        mu[i,j] <- 20 * sqrt(1 - (s[j] - 0.5)^2 - (t[i] - 0.5)^2)
      }
    }
    c(mu)
  }
  
  # Longitudinal eigenfunctions, psi_j(s) (vector, SS x 1)
  psi_1 <- -2^{1/2} * cos(2 * pi * s)
  psi_2 <- 2^{1/2} * sin(2 * pi * s)
  
  # Functional eigenfunctions, phi_k(t) (vector, TT x 1)
  phi_1 <- -2^{1/2} * cos(3 * pi * t)
  phi_2 <- -2^{1/2} * sin(3 * pi * t)

  # Marginal longitudinal eigenvalues (vector, J x 1)
  tau <- 2*(exp(1.2 * (3 - c(1, 2)))/sum(exp(1.2 * c(1, 2))))
  FVE_tau <- tau/sum(tau)
  
  # Marginal functional eigenvalues (vector, K x 1)
  vartheta <- 2*(exp(1.6 * (3 - c(1, 2)))/sum(exp(1.6 * c(1, 2))))
  FVE_vartheta <- vartheta/sum(vartheta)
  
  # Covariance matrix of subject-specific product scores (matrix, J x K)
  V <- diag(c(tau %*% t(vartheta)))
  
  # Subject-specific product scores (vector, (J+K) x 1)
  chi <- mvrnorm(n = 1, mu = rep(0, 4), Sigma = V)
  
  # Calculate the true underlying longitudinal/functional eigenfunctions
  # Two-dimensional eigenfunctions
  gamma_1 <- kronecker(psi_1, phi_1)
  gamma_2 <- kronecker(psi_2, phi_1)
  gamma_3 <- kronecker(psi_1, phi_2)
  gamma_4 <- kronecker(psi_2, phi_2)
  gamma <- cbind(gamma_1, gamma_2, gamma_3, gamma_4)
  
  cov_kernel <- gamma %*% V %*% t(gamma) 
  K_s <- LFBayes::get_marginal_long(cov_kernel, SS, TT)
  tau_truth <- eigen(K_s)$values[1:2] * (s[2] - s[1])
  psi_truth <- eigen(K_s)$vectors[,1:2]
  K_t <- LFBayes::get_marginal_func(t(cov_kernel), SS, TT)/(t[2] - t[1])
  vartheta_truth <- eigen(K_t)$values[1:2] * (t[2] - t[1])
  phi_truth <- eigen(K_t)$vectors[,1:2]
  
  # True underlying function, marginal covariances and covariance kernel
  truth_funcs <- list(mu = mu(s, t),
                      psi = list(psi_truth[,1], -1 * psi_truth[,2]),
                      phi = list(phi_truth[,1], phi_truth[,2]),
                      tau = tau_truth,
                      vartheta = vartheta_truth,
                      K_s = K_s, 
                      K_t = K_t,
                      K = cov_kernel)
  
  # Product FPCA components (vector, (SS x TT)  x 1)
  # Function to generate product FPCA components 
    ## Input: longitudinal time points s (vector, TT x 1)
    ##        functional time points s (vector, SS x 1)  
    ##        covariance matrix V (matrix, J x K)
    ## Output: components for the product FPCA part (vector, (SS x TT)  x 1)
  signal <- function(s, t, V){
    chi <- mvrnorm(n = 1, mu = rep(0, 4), Sigma = V)
    
    signal <- (chi[1] * kronecker(psi_1, phi_1)) + (chi[2] * kronecker(psi_2, phi_1)) +
      (chi[3] * kronecker(psi_1, phi_2)) + (chi[4] * kronecker(psi_2, phi_2))
    
    c(signal)
  }

  #############################################################################
  # 2. Generate functional data
  #############################################################################
  # Generate the outcome c(Y_i(s, t)) in list format for all subjects
  # Each list for one subject containing (SS x TT)  x 1 outcome components
  y <- list()
  for(i in 1:n){
    y[[i]] <- mu(s, t) + signal(s, t, V) + rnorm(n = SS*TT, mean = 0, sd = sqrt(5))
  }
  
  return(list(y = y, truth_funcs = truth_funcs))
}
