MCMC <- function(Y,           # List of simulated functional data, each list contains a vector of Y_i(s, t) (list, (SS x TT) x 1 x n )
                 n,           # number of subjects (scalar)
                 SS,          # number of functional time points (scalar)
                 TT,          # number of longitudinal time points (scalar)
                 p1,          # Number of longitudinal basis functions (integer)
                 p2,          # Number of functional basis functions (integer)
                 q1,          # Number of longitudinal latent functions (integer)
                 q2,          # Number of functional latent functions (integer)
                 nchain,      # Number of MCMC chains (integer)
                 iter,        # Total number of iterations in each MCMC chain (integer)
                 burnin,      # Number of iterations used for burn-in (integer, 0-iter)
                 thin         # Number of thinning (integer)
){
  #############################################################################
  ## Description: Function for performing posterior estimation using a Gibbs sampler 
  ##              described in "Central posterior envelopes for Bayesian longitudinal 
  ##              functional principal component analysis" by Boland et al. (2024). 
  ##              The sampling model, prior distributions, and posterior distributions
  ##              employed and calculated are described in Section 2.1 and Supplementary 
  ##              Materials Appendix A of the manuscript. Posterior sample object of size C 
  ##              using the B-LFPCA model are returned from this function. 
  ## Args:        see above
  ## Returns:     List of posterior sample object returned from `run_mcmc` function in the `LFBayes` package
  ## MCMC Outline:
  ##              1. Define global variables 
  ##              2. Employ MCMC Chain using the `LFBayes` package 
  ##              3. Collect posterior model components
  #############################################################################
  
  # Install missing packages
  list.of.packages <- c("MASS", "Matrix", "data.table", "magrittr", "splines",
                        "RcppProgress", "RcppArmadillo", "LFBayes")
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if (length(new.packages)) install.packages(new.packages) 
  
  # Load packages
  library(MASS)
  library(Matrix)
  library(data.table)
  library(magrittr)
  library(splines)
  library(RcppProgress)
  library(RcppArmadillo)
  library(LFBayes)
  
  #############################################################################
  # 1. Define global variables 
  #############################################################################
  # Number of longitudinal and functional time points
  t <- seq(from = 0, to = 1, length.out = TT)
  s <- seq(from = 0, to = 1, length.out = SS)
  
  # Number of eigencomponents in longitudinal and functional dimensions
  J <- 2; K <- 2
  
  # Number of subjects
  n <- 30
  
  # Number of basis functions
  p1 <- p1
  p2 <- p2
  
  # Basis functions for longitudinal and functional dimensions
  tt <- list(); tt[[1]] <- 1:(TT*SS); tt <- rep(tt, n)
  Bs <- bs(s, df = p1, knots = c(1/8, 1/4, 3/8, 1/2, 5/8, 3/4, 7/8), intercept = TRUE)
  Bt <- bs(t, df = p2, knots = c(1/8, 1/4, 3/8, 1/2, 5/8, 3/4, 7/8), intercept = TRUE)
  
  # Define the missing values in y in list format
  missing <- list()
  for(ii in 1:n){
    missing[[ii]] <- numeric(0)
  }
  
  # Define the design matrix 
  X <- cbind(rep(1, n))
  
  #############################################################################
  # 2. MCMC procedure using the `LFBayes` package
  #############################################################################
  ### MCMC control parameters
  iter <- iter # Number of iterations
  burnin <- burnin # Burnin iterations
  thin <- thin # Thinning for each chain
  nchain <- nchain # Number of chains
  q1 <- q1 # Number of latent factors for longitudinal dimension
  q2 <- q2 # Number of latent factors for functional dimension
  
  ### Input longitudinal functional response
  y <- Y
  
  ### Running MCMC
  posterior_samples <- LFBayes::run_mcmc(y, missing, X, Bs, Bt,
                                q1, q2, iter, thin, burnin, nchain)
  return(posterior_samples)
}