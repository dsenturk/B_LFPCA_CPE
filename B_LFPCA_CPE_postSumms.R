postSumms <- function(posterior_samples, # output object from MCMC() of FDpostSumms_BLFPCA_MCMC.R
                      neig,              # number of eigencomponents (integer) 
                      n,                 # number of subjects (scalar)
                      SS,                # number of functional time points (scalar)
                      TT,                # number of longitudinal time points (scalar)
                      p1,                # number of longitudinal basis functions (integer)
                      p2,                # number of functional basis functions (integer)
                      q1,                # number of longitudinal latent functions (integer)
                      q2,                # number of functional latent functions (integer)
                      nchain,            # number of MCMC chains (integer)
                      iter,              # total number of iterations in each MCMC chain (integer)
                      burnin,            # number of iterations used for burn-in (integer, 0-iter)
                      thin,              # number of thinning (integer)
                      alpha              # alpha-level of (1 - alpha)100% credible intervals and central posterior envelopes (scalar, 0-1)
                      ){
  #############################################################################
  ## Description: Function for calculating the traditional and proposed functional 
  ##              depth posterior summaries for the model components described in 
  ##              "Central Posterior Envelopes for Bayesian Functional Principal 
  ##              Component Analysis" by Boland et al. (2022). Calculations included 
  ##              are the traditional posterior summaries for the mean function, 
  ##              covariance function, eigenvalues, and eigenfunctions detailed in 
  ##              Section 2.1 and the proposed posterior summaries for the mean function
  ##              and eigenfunctions detailed in Section 3. Specifically, traditional 
  ##              point estimates are calculated for all the model components, and 
  ##              both traditional credible intervals and proposed central posterior 
  ##              envelopes (CPEs) are calculated for the mean function and eigenfunctions. 
  ##              For quick reference, Table 1 in the manuscript contains a summary 
  ##              of the nomenclature and notation of all the posterior summaries.
  ## Args:        (see above)
  ## Returns:     list()
  ##              mu: (list, length 2)
  ##                point_ests: (list, length 2)
  ##                  mu.hat: estimated mean, \hat{\mu}(s,t) (vector, (SSxTT) x 1)
  ##                  mu.med.hat: MVD median mean, \tilde{m}{\mu(s,t)} (vector, (SSxTT) x 1)
  ##                credible_ints_cpes: (list, length 2)
  ##                    simultaneous_parametric: P^s_{1 - \alpha}{\mu(s,t)} (list, length 2)
  ##                      lower: lower bound (vector, (SSxTT) x 1)
  ##                      upper: upper bound (vector, (SSxTT) x 1)
  ##                    MVD: D^{\star}_{1 - \alpha}{\mu(s,t)} (list, length 2)
  ##                      lower: lower bound (vector, (SSxTT) x 1)
  ##                      upper: upper bound (vector, (SSxTT) x 1)
  ##              K_s: (list, length 1)
  ##                point_ests: (list, length 4)
  ##                  K_s.tilde: estimated marginal longitudinal covariance, \tilde{K}_s(s, s') (matrix, SS x SS)
  ##                  K_s.bar: estimated marginal longitudinal covariance, \bar{K}_s(s, s') (matrix, SS x SS)
  ##                  K_s.med.tilde: median of MVD marginal longitudinal covariance, \tilde{m}{K_s(s, s')} (matrix, SS x SS)
  ##                  K_s.med.bar: median of kernel marginal longitudinal covariance, \bar{m}{K_s(s, s')} (matrix, SS x SS)
  ##              K_t: (list, length 1)
  ##                point_ests: (list, length 4)
  ##                  K_t.tilde: estimated marginal functional covariance, \tilde{K}_t(t, t') (matrix, TT x TT)
  ##                  K_t.bar: estimated marginal functional covariance, \bar{K}_t(t, t') (matrix, TT x TT)
  ##                  K_t.med.tilde: median of MVD marginal functional covariance, \tilde{m}{K_t(t, t')} (matrix, TT x TT)
  ##                  K_t.med.bar: median of kernel marginal functional covariance, \bar{m}{K_t(t, t')} (matrix, TT x TT)
  ##              K: (list, length 1)
  ##                point_ests: (list, length 2)
  ##                  K.bar: estimated covariance kernel, \hat{K}{(s, t),(s', t')} (matrix, (SSxTT) x (SSxTT))
  ##                  K.med.bar: MVD median covariance kernel, \hat{m}[K{(s, t),(s', t')}] (matrix, (SSxTT) x (SSxTT))
  ##              tau_j: (list, length J)
  ##                tau_1: (list, length 1)
  ##                ...
  ##                tau_J: (list, length 1)
  ##                 point_ests: (list, length 3)
  ##                  tau_j.hat: estimated longitudinal eigenvalue, \hat{\tau}_2 (scalar)
  ##                  tau_j.tilde: estimated longitudinal eigenvalue, \tilde{\tau}_2 (scalar)
  ##                  tau_j.bar: estimated longitudinal eigenvalue, \bar{\tau}_2 (scalar)
  ##              vartheta_k: (list, length 2)
  ##                vartheta_1: (list, length 1)
  ##                ...
  ##                vartheta_K: (list, length 1)
  ##                 point_ests: (list, length 3)
  ##                  vartheta_k.hat: estimated functional eigenvalue, \hat{vartheta}_2 (scalar)
  ##                  vartheta_k.tilde: estimated functional eigenvalue, \tilde{vartheta}_2 (scalar)
  ##                  vartheta_k.bar: estimated functional eigenvalue, \bar{vartheta}_2 (scalar)
  ##              psi_j: (list, length J)
  ##                psi_1: (list, length 2)
  ##                ...
  ##                psi_J: (list, length 2)
  ##                  point_ests: (list, length 6)
  ##                    psi_j.hat: estimated longitudinal eigenfunction, \hat{\psi}_k(s) (vector, SS x 1)
  ##                    psi_j.tilde: estimated longitudinal eigenfunction via covariance, \tilde{\psi}_k(s) (vector, SS x 1)
  ##                    psi_j.bar: estimated longitudinal eigenfunction via covariance kernel, \bar{\psi}_k(s) (vector, SS x 1)
  ##                    psi_j.med.hat: MBD median longitudinal eigenfunction, \hat{m}{\psi_k(s)} (vector, SS x 1)
  ##                    psi_j.med.tilde: MVD median longitudinal eigenfunction, \tilde{m}{\psi_k(s)} (vector, SS x 1)
  ##                    psi_j.med.bar: kernel MVD median longitudinal eigenfunction, \bar{m}{\psi_k(s)} (vector, SS x 1)
  ##                  credible_ints_cpes: (list, length 4)
  ##                    simultaneous_parametric: P^s_{1 - \alpha}{\psi_j(s)} (list, length 2)
  ##                      lower: lower bound (vector, SS x 1)
  ##                      upper: upper bound (vector, SS x 1)
  ##                    MBD: D_{1 - \alpha}{\psi_j(s)} (list, length 2)
  ##                      lower: lower bound (vector, SS x 1)
  ##                      upper: upper bound (vector, SS x 1)
  ##                    MVD.marg: D^{\star}_{1 - \alpha}{\psi_j(s)} (list, length 2)
  ##                      lower: lower bound (vector, SS x 1)
  ##                      upper: upper bound (vector, SS x 1)
  ##                    MVD.kernel: D^{\dager}_{1 - \alpha}{\psi_j(s)} (list, length 2)
  ##                      lower: lower bound (vector, SS x 1)
  ##                      upper: upper bound (vector, SS x 1)
  ##              phi_k: (list, length K)
  ##                phi_1: (list, length 2)
  ##                ...
  ##                phi_K: (list, length 2)
  ##                  point_ests: (list, length 6)
  ##                    phi_k.hat: estimated functional eigenfunction, \hat{\phi}_k(t) (vector, TT x 1)
  ##                    phi_k.tilde: estimated functional eigenfunction via covariance, \tilde{\phi}_k(t) (vector, TT x 1)
  ##                    phi_k.bar: estimated functional eigenfunction via covariance kernel, \bar{\phi}_k(t) (vector, TT x 1)
  ##                    phi_k.med.hat: MBD median functional eigenfunction, \hat{m}{\phi_k(t)} (vector, TT x 1)
  ##                    phi_k.med.tilde: MVD median functional eigenfunction, \tilde{m}{\phi_k(t)} (vector, TT x 1)
  ##                    phi_k.med.bar: kernel MVD median functional eigenfunction, \bar{m}{\phi_k(t)} (vector, TT x 1)
  ##                  credible_ints_cpes: (list, length 4)
  ##                    simultaneous_parametric: P^s_{1 - \alpha}{\phi_j(t)} (list, length 2)
  ##                      lower: lower bound (vector, TT x 1)
  ##                      upper: upper bound (vector, TT x 1)
  ##                    MBD: D_{1 - \alpha}{\phi_j(t)} (list, length 2)
  ##                      lower: lower bound (vector, TT x 1)
  ##                      upper: upper bound (vector, TT x 1)
  ##                    MVD.marg: D^{\star}_{1 - \alpha}{\phi_j(t)} (list, length 2)
  ##                      lower: lower bound (vector, TT x 1)
  ##                      upper: upper bound (vector, TT x 1)
  ##                    MVD.kernel: D^{\dager}_{1 - \alpha}{\phi_j(t)} (list, length 2)
  ##                      lower: lower bound (vector, TT x 1)
  ##                      upper: upper bound (vector, TT x 1)
  ##                ranks: (list, length 6)
  ##                  mu_MVD: matrix(3 x C): id   # number of iteration
  ##                                         mbd  # mbd values
  ##                                         rank # ranks according to mbd values
  ##                  MVD.marginal_cov_func: matrix(3 x C): id   # number of iteration
  ##                                                        mvd  # mvd values for marginal functional covariance
  ##                                                        rank # ranks according to mvd values
  ##                  MVD.marginal_cov_long: matrix(3 x C): id   # number of iteration
  ##                                                        mvd  # mvd values for marginal longitudinal covariance
  ##                                                        rank # ranks according to mvd values
  ##                  MVD.cov_kernel: matrix(3 x C): id   # number of iteration
  ##                                                 mvd  # mvd values for covariance kernel
  ##                                                 rank # ranks according to mvd values
  ##                  MBD.psi_j: (list, length J)
  ##                    In each list: matrix(3 x C): id   # number of iteration
  ##                                                 mbd  # mbd values for longitudunal eigenfunction
  ##                                                 rank # ranks according to mbd values
  ##                  MBD.phi_k: (list, length K)
  ##                    In each list: matrix(3 x C): id   # number of iteration
  ##                                                 mbd  # mbd values for functional eigenfunction
  ##                                                 rank # ranks according to mbd values
  ##                mu_sample: matrix of posterior mean function samples (matrix, (SS x TT) x C)
  ##                tau_long: matrix of posterior longitudinal eigenvalue samples (matrix, J x C)
  ##                vartheta_func: matrix of posterior functional eigenvalue samples (matrix, K x C)
  ##                Psi_long: array of posterior longitudinal eigenfunction samples (array, SS x J x C)
  ##                Phi_func: array of posterior functional eigenfunction samples (array, TT x K x C)
  ## postSumms Outline:
  ##              1. Traditional posterior summaries
  ##              2. Traditional credible intervals
  ##              3. Proposed functional depth based point estimates
  ##              4. Proposed functional depth based credible intervals
  #############################################################################
  
  #############################################################################
  # 1. Preliminaries
  #############################################################################
  
  # Install missing packages
  list.of.packages <- c("splines", "Rcpp", "acid", "pracma", "expm", "LFBayes")
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if (length(new.packages)) install.packages(new.packages) 
  
  # Load libraries
  library(splines)
  library(Rcpp)
  library(acid)
  library(pracma)
  library(expm)
  library(LFBayes)
  
  # Function to calculate number of combinations C_n^k
  # Input: n: # number of total subjects (integer)
  #        k: # number of chosen subjects (integer)
  # Output: combinat: # number of combinations C_n^k (integer)
  combinat <- function(n, p) 
  { # tota
    if (n < p) {combinat = 0}
    else {combinat = exp(lfactorial(n) - (lfactorial(p) + lfactorial(n-p)))}
  }
  
  # Function to calculatre MBD values
  # Input: Transposed functional posterior estimates, matrix (T x C)
  # Output: data.frame with columns c("id", "mbd", "mbd.rank")
  #        id: row index, 1,...,C
  #        mbd: MBD value, MBD_{C, 2}{X^{(c)}(t)}
  #        mbd.rank: rank of the MBD value, MBD_{C, 2}{X^{[c]}(t)} for function X^{(c)}(t)
  MBD <- function(data  # Transposed functional posterior estimates, matrix (T x C)
  ){
    p = dim(data)[1]  
    n = dim(data)[2]  
    rmat = apply(data, 1, rank)
    down = rmat - 1
    up = n - rmat
    mbd <- (rowSums(up * down) / p + n - 1) / combinat(n, 2)  # Calculate MBD values 
    mbd.rank <- rank(mbd)  # Rank the MBD values from smallest to largest (1 -> C)
    
    # Return data.frame with columns c("id", "mbd", "mbd.rank")
    #        id: row index, 1,...,C
    #        mbd: MBD value, MBD_{C, 2}{X^{(c)}(t)}
    #        mbd.rank: rank of the MBD value, MBD_{C, 2}{X^{[c]}(t)} for function X^{(c)}(t)
    return(data.frame(id = 1:n, mbd = mbd, rank = mbd.rank))
  }
  
  # Function to calculate MVD of a sample of covariance functions
  # Input: functional posterior covariance estimates (list of length C)
  # Output: data.frame with columns c("id", "mvd", "mvd.rank")
  #        id: row index, 1,...,C
  #        mvd: MBD value, MVD_{C, 2}{X^{(c)}(s,t)}
  #        mvd.rank: rank of the MVD value, MVD_{C, 2}{X^{[c]}(s,t)} for function X^{(c)}(s,t)
  MVD <- function(list  # Functional posterior covariance estimates, (list, length C)
  ){
    data <- do.call(cbind, lapply(1:length(list), function(m){  # Vectorize the covariance surfaces
      as.vector(list[[m]])
    }))
    
    p = dim(data)[1]  
    n = dim(data)[2]  
    rmat = apply(data, 1, rank)
    down = rmat - 1
    up = n - rmat
    mvd <- (rowSums(up * down) / p + n - 1) / combinat(n, 2)  # Calculate MVD values 
    mvd.rank <- rank(mvd)  # Rank the MVD values from smallest to largest (1 -> M)
    
    # Return data.frame with columns c("id", "mvd", "mvd.rank")
    #        id: row index, 1,...,C
    #        mvd: MBD value, MVD_{C, 2}{X^{(c)}(s,t)}
    #        mvd.rank: rank of the MVD value, MVD_{C, 2}{X^{[c]}(s,t)} for function X^{(c)}(s,t)
    return(data.frame(id = 1:n, mvd = mvd, rank = mvd.rank))
  }
  
  # Function to calculate simultaneous parametric credible intervals
  # Input: data,  # Sample of posterior estimates (matrix, C x T)
  #        alpha  # Alpha-level of credible interval (scalar, 0-1)
  # Output: list of length 2
  #            upper: upper bounds (vector, T x 1) 
  #            lower: lower bounds (vector, T x 1) 
  simultaneous_parametric <- function(data,  # Sample of posterior estimates, (matrix, C x T)
                                      alpha  # Alpha-level of credible interval, (scalar, 0-1)
  ){
    point.mean <- colMeans(data)  # Calculate pointwise mean, (vector, T x 1)
    point.sd <- apply(data, 2, sd)  # Calculate pointwise standard deviation, (vector, T x 1)
    
    # Calculate sample quantile of the absolute standardized deviation, c_alpha
    c_alpha.m <- unlist(lapply(1:nrow(data), function(m){
      max(abs((data[m,] - point.mean))/point.sd)
    }))
    c_alpha <- quantile(c_alpha.m, probs = 1 - alpha)
    
    # Return list() containing upper (vector, T x 1) and lower (vector, T x 1) bounds
    return(list(lower = point.mean - c_alpha * point.sd,
                upper = point.mean + c_alpha * point.sd))
  }
  
  # Function to calculate MVD or MBD credible interval
  # Input: rank,  # Data.frame returned from the MBD() or MVD() functions
  #        data,  # Sample of posterior estimates (matrix, C x T)
  #        alpha  # Alpha-level of credible interval (scalar, 0-1)
  # Output: list of length 2
  #            upper: upper bounds (vector, T x 1) 
  #            lower: lower bounds (vector, T x 1) 
  depth_credible_interval <- function(rank,  # Data.frame returned from the MBD() or MVD() functions
                                      data,  # Posterior estimates (matrix, C x T)
                                      alpha  # Alpha-level of credible intervals
  ){
    
    M <- nrow(data)  # Number of posterior estimates, M
    
    # Obtain the row indices for the posterior estimates with the largest MBD/MVD values
    row.indices <- rank$id[rank$rank %in% ((alpha * M) + 1):M]
    
    # Subset of posterior estimates with largest MBD or MVD to form credible interval 
    data.subset <- data[row.indices, ]
    
    # Return list of lower (T x 1) and upper (T x 1) bounds
    return(list(lower = apply(data.subset, 2, min),
                upper = apply(data.subset, 2, max)))
  }
  
  # Function to calculate integrals for latent functions
  # Input: latent: latent functions (vector of length of time points)
  #        times: times (vector of length of time points)
  # Output: latent_int: integral results for latent functions (vector of latent function dimensions)
  integrated_latent <- function(latent, times){
    latent_dim <- ncol(latent)
    latent_int <- numeric(latent_dim)
    for(i in 1:latent_dim){
      latent_int[i] <- pracma::trapz(times, latent[, i] * latent[, i])
    }
    return(latent_int)
  }
  
  # Function to calculate integrals for B-spline basis functions
  # Input: spline: B-spline basis functions (vector of length of time points)
  #        times: times (vector of length of time points)
  # Output: latent_int: integral results for B-spline basis functions (matrix of column(B-spline basis function) x column(B-spline basis function))
  integrated <- function(spline, times){
    spline_dim <- ncol(spline)
    spline_int <- matrix(nrow = spline_dim, ncol = spline_dim)
    for(i in 1:spline_dim){
      for(j in 1:spline_dim){
        spline_int[i, j] <- pracma::trapz(times, spline[, i] * spline[, j])
      }
    }
    return(spline_int)
  }
  
  # Function to extract covariance 
  extract_covlatent2 <- function(latent, S.sub, H.sub, spline, spline_int_sqrt, spline_int_sqrt_inv,
                                 spline_int, latent_trapz){
    spline_dim <- ncol(spline)
    time_size <- nrow(spline)
    HD <- diag(c(H.sub %*% latent_trapz))
    cov_latent <- spline %*% (diag(c(S.sub %*% diag(spline_int))) + latent %*% HD %*% t(latent)) %*% t(spline)
    # cov_latent <- spline %*% spline_int_sqrt_inv %*% cov_latent %*% t(spline %*% spline_int_sqrt_inv)
    return(cov_latent)
  }
  
  # Function to extract eigenfunctions
  extract_eigenfn <- function(latent, S.sub, H.sub, spline, spline_int_sqrt, spline_int_sqrt_inv,
                              spline_int, latent_trapz, numeig){
    spline_dim <- ncol(spline)
    time_size <- nrow(spline)
    HD <- diag(c(H.sub %*% latent_trapz))
    cov_latent <- spline_int_sqrt %*% (diag(c(S.sub %*% diag(spline_int))) + latent %*% HD %*% t(latent)) %*% spline_int_sqrt
    eig.cov_latent <- eigen(cov_latent)
    eigenfn <- eig.cov_latent$vectors; eigval <- eig.cov_latent$values
    eigen_spline <- matrix(nrow = time_size, ncol = 2)
    for(i in 1:2){
      eigen_spline[, i] <- spline %*% spline_int_sqrt_inv %*% eigenfn[, i]
    }
    for(i in 1:2){
      eigen_spline[, i] <-  eigen_spline[, i]/norm(eigen_spline[, i], type = "2")
    }
    return(list(eigenfn = eigen_spline, eigenval = eigval))
  }
  
  # Align the eigenfunctions
  # Input: eigenfn: eigenfunctions that need to be aligned
  # Output: eigenfn: eigenfunctions that finish alignment
  eigenfn_align <- function(eigenfn){
    eigenfn_mean <- rep(0, nrow(eigenfn))
    for(i in 1:ncol(eigenfn)){
      if(norm(eigenfn[, i] - -1 * eigenfn_mean, type = "2") < 
         norm(eigenfn[, i] - eigenfn_mean, type = "2")){
        eigenfn[, i] <- -1 * eigenfn[, i]
      }
      eigenfn_mean <- (eigenfn[, i] - eigenfn_mean) / i + eigenfn_mean
    }
    return(eigenfn)
  }
  
  # Function to calculate IMSE
  # Input: est: estimated values (can be scalar or vector)
  #        truth: underlying values (can be scalar of vector)
  IMSE <- function(est, truth){
    sum((est - truth)^2)/sum(truth^2)
  }
  
  # Function to calculate IMSE when functions have been aligned
  # Input: est: estimated functions (vector)
  #        truth: underlying true functions (vector)
  IMSE.flip_eigen <- function(est, truth){
    if(as.numeric(IMSE(est, truth)) > as.numeric(IMSE(-1 * est, truth))){
      -1 * est
    } else {
      est
    }
  }
  
  # Function to record the index of functions that have been aligned
  # Input: est: estimated functions (vector)
  #        truth: underlying true functions (vector)
  IMSE.flip_eigen_index <- function(est, truth){
    if(as.numeric(IMSE(est, truth)) > as.numeric(IMSE(-1 * est, truth))){
      FALSE
    } else {
      TRUE
    }
  }
  
  # Define the eigencomponents that needs to be retained
  K <- neig; J <- neig

  #############################################################################
  # 2. Obtain posterior estimates needed for calculations
  #############################################################################
  # Create a null list to store all the summary results
  post_summs <- NULL
  post_summs <- list(
      mu = list(point_ests = list(mu.hat = NA, 
                                  mu.med.hat = NA), 
                credible_ints = list(simultaneous_parametric = NA,
                                     MVD = NA)),
      K_s = list(point_ests = list(K_s.tilde = NA,
                                   K_s.bar = NA,
                                   K_s.med.tilde = NA,
                                   K_s.med.bar = NA)),
      K_t = list(point_ests = list(K_t.tilde = NA,
                                   K_t.bar = NA,
                                   K_t.med.tilde = NA,
                                   K_t.med.bar = NA)),
      K = list(point_ests = list(K.bar = NA,
                                 K.med.bar = NA)),
      tau_j = list(), vartheta_k = list(),
      psi_j = list(), phi_k = list())
    for(j in 1:J){
      post_summs$tau_j[[j]] <- list(point_ests = list(tau_j.hat = NA,
                                                      tau_j.tilde = NA,
                                                      tau_j.bar = NA))
      
      post_summs$psi_j[[j]] <- list(point_ests = list(psi_j.hat = NA,
                                                      psi_j.tilde = NA,
                                                      psi_j.bar = NA, 
                                                      psi_j.med.hat = NA,
                                                      psi_j.med.tilde = NA,
                                                      psi_j.med.bar = NA),
                                    credible_ints = list(simultaneous_parametric = NA,
                                                         MBD = NA, 
                                                         MVD.marg = NA,
                                                         MVD.kernel = NA))
    }
    for(k in 1:K){
      post_summs$vartheta_k[[k]] <- list(point_ests = list(vartheta_k.hat = NA,
                                                           vartheta_k.tilde = NA,
                                                           vartheta_k.bar = NA))
      
      post_summs$phi_k[[k]] <- list(point_ests = list(phi_k.hat = NA,
                                                      phi_k.tilde = NA,
                                                      phi_k.bar = NA,
                                                      phi_k.med.hat = NA,
                                                      phi_k.med.tilde = NA,
                                                      phi_k.med.bar = NA),
                                    credible_ints = list(simultaneous_parametric = NA,
                                                         MBD = NA, 
                                                         MVD.marg = NA,
                                                         MVD.kernel = NA))
    }
    names(post_summs$tau_j) <- paste("tau_", 1:J, sep = "")
    names(post_summs$vartheta_k) <- paste("vartheta_", 1:K, sep = "")
    names(post_summs$psi_j) <- paste("psi_", 1:J, sep = "")
    names(post_summs$phi_k) <- paste("phi_", 1:K, sep = "")
    
    # Define time points
    t <- seq(from = 0, to = 1, length.out = TT)
    s <- seq(from = 0, to = 1, length.out = SS)
    
    # Define splines
    Bs1 <- bs(s, df = p1, knots = c(1/8, 1/4, 3/8, 1/2, 5/8, 3/4, 7/8), intercept = TRUE)
    Bt1 <- bs(t, df = p2, knots = c(1/8, 1/4, 3/8, 1/2, 5/8, 3/4, 7/8), intercept = TRUE)
    splineT <- Bt1
    splineS <- Bs1
    
    gridsize_t <- t[2] - t[1]; gridsize_s <- s[2] - s[1]
    
    m_index <- seq(burnin + 1, iter, by = thin)
    
    # Extract Lambda results from posterior samples
    Lambda <- lapply(X = 1:nchain, function(X){
      Lambda.sub <- posterior_samples$Lambda[[X]]
      Lambda.sub <- Lambda.sub[ , , m_index]
      Lambda.sub
    })
    
    # Extract Gamma results from posterior samples
    Gamma <- lapply(X = 1:nchain, function(X){
      Gamma.sub <- posterior_samples$Gamma[[X]]
      Gamma.sub <- Gamma.sub[ , , m_index]
      Gamma.sub
    })
    
    # Extract Beta results from posterior samples
    Beta <- lapply(X = 1:nchain, function(X){
      Beta.sub <- posterior_samples$Beta[[X]]
      Beta.sub <- Beta.sub[ , , m_index]
      Beta.sub
    })
    
    # Extract Sigma results from posterior samples
    Sigma <- lapply(X = 1:nchain, function(X){
      Sigma.sub <- posterior_samples$Sigma[[X]]
      Sigma.sub <- Sigma.sub[ , , m_index]
      Sigma.sub
    })
    
    # Extract H results from posterior samples
    H <- lapply(X = 1:nchain, function(X){
      H.sub <- posterior_samples$H[[X]]
      H.sub <- H.sub[ , , m_index]
      H.sub
    })
    
    # Save memory
    rm(posterior_samples); gc()
    
    #############################################################################
    # 3. Calculate posterior summaries of mean functions \mu(s, t)
    #############################################################################
    # Caluate mean functions
    mu <- do.call(cbind, lapply(X = 1:1500, FUN = function(X){
      kronecker(Bs1 %*% Gamma[[1]][, , X], Bt1 %*% Lambda[[1]][, , X]) %*% (Beta[[1]][, X])
    }))
    for(Y in 2:nchain){
      mu <- cbind(mu, do.call(cbind, lapply(X = 1:1500, FUN = function(X){
        kronecker(Bs1 %*% Gamma[[Y]][, , X], Bt1 %*% Lambda[[Y]][, , X]) %*% (Beta[[Y]][, X])
      })))
    }
    
    # Summaries related with mean functions
    post_summs$mu$point_ests$mu.hat <- rowMeans(mu)

    post_summs$mu$credible_ints$simultaneous_parametric <- simultaneous_parametric(data = t(mu), alpha = alpha)
    
    mu.MVD <- MBD(mu)
    mu.MVD <- mu.MVD[order(-mu.MVD$rank),]
    post_summs$ranks$mu_MVD <- mu.MVD
    
    post_summs$mu$point_ests$mu.med.hat <- (mu[, mu.MVD$id[1]])
    post_summs$mu$credible_ints$MVD <- depth_credible_interval(mu.MVD, t(mu), alpha)
    
    #############################################################################
    # 4. Calculate Covariance Kernels 
    #############################################################################
    splineS_int <- integrated(splineS, s)
    splineT_int <- integrated(splineT, t)
    splineS_int_sqrt <- pracma::sqrtm(splineS_int); splineT_int_sqrt <- pracma::sqrtm(splineT_int)
    splineS_int_sqrt_inv <- splineS_int_sqrt$Binv; splineT_int_sqrt_inv <- splineT_int_sqrt$Binv
    splineS_int_sqrt <- splineS_int_sqrt$B; splineT_int_sqrt <- splineT_int_sqrt$B
    
    # Calculate marginal functional covariance
    marginal_cov_func <- lapply(X = 1:1500, function(X){
      Lambda.sub <- Lambda[[1]][, , X]
      Gamma.sub <- Gamma[[1]][, , X]
      Psi <- splineS %*% Gamma.sub
      Phi <- splineT %*% Lambda.sub
      Psi_trapz <- integrated_latent(Psi, s)
      Phi_trapz <- integrated_latent(Phi, t)
      H.sub <- H[[1]][, , X]
      H.sub <- matrix(1/diag(H.sub), nrow = q1, ncol = q2)
      S.sub <- 1/(Sigma[[1]][, , X])
      
      cov_func <- extract_covlatent2(
        Lambda.sub, S.sub, H.sub, splineT, splineT_int_sqrt,
        splineT_int_sqrt_inv, splineS_int, Psi_trapz)
      
      return(cov_func)
    })
    for(Y in 2:nchain){
      marginal_cov_func <- append(marginal_cov_func, lapply(X = 1:1500, function(X){
        Lambda.sub <- Lambda[[Y]][, , X]
        Gamma.sub <- Gamma[[Y]][, , X]
        Psi <- splineS %*% Gamma.sub
        Phi <- splineT %*% Lambda.sub
        Psi_trapz <- integrated_latent(Psi, s)
        Phi_trapz <- integrated_latent(Phi, t)
        H.sub <- H[[Y]][, , X]
        H.sub <- matrix(1/diag(H.sub), nrow = q1, ncol = q2)
        S.sub <- 1/(Sigma[[Y]][, , X])
        
        cov_func <- extract_covlatent2(
          Lambda.sub, S.sub, H.sub, splineT, splineT_int_sqrt,
          splineT_int_sqrt_inv, splineS_int, Psi_trapz)
        
        return(cov_func)}))
    }
    
    # Summaries related with marginal functional covariance
    post_summs$K_t$point_ests$K_t.tilde <-  Reduce("+", marginal_cov_func)/length(marginal_cov_func)
    MVD.marginal_cov_func <- MVD(marginal_cov_func); gc()
    MVD.marginal_cov_func <- MVD.marginal_cov_func[order(-MVD.marginal_cov_func$rank),]
    post_summs$ranks$MVD.marginal_cov_func <- MVD.marginal_cov_func
    
    post_summs$K_t$point_ests$K_t.med.tilde <- marginal_cov_func[[MVD.marginal_cov_func$id[[1]]]]
    
    # Extract functional eigenfunctions and eigenvalues from mariginal functional covariance
    ## Extract functional eigenfunctions
    eigen_func <- lapply(X = 1:6000, function(X){
      K <- marginal_cov_func[[X]]
      eigen.K <- eigen(K)
      vartheta_m <- eigen.K$values[1:2] * gridsize_t
      phi_m <- eigen.K$vectors[, c(1:2)]
      return(list(vartheta_m = vartheta_m, phi_m = phi_m))
    })
    
    Phi_func <- array(dim = c(TT, K, 6000))
    for(Y in 1:2){
      phi_data <- do.call(cbind, lapply(1:6000, function(X){
        eigen_func[[X]]$phi_m[, Y]
      }))
      phi_data <- eigenfn_align(phi_data)
      Phi_func[ , Y, ] <- phi_data
    }
    # Extract functional eigenvalues
    vartheta_func <- do.call(cbind, lapply(X = 1:6000, function(X){
      eigen_func[[X]]$vartheta_m
    }))
    
    # Save memory
    rm(eigen_func); gc()
    
    # Calculate marginal longitudinal covariance
    marginal_cov_long <- lapply(X = 1:1500, function(X){
      Lambda.sub <- Lambda[[1]][, , X]
      Gamma.sub <- Gamma[[1]][, , X]
      Psi <- splineS %*% Gamma.sub
      Phi <- splineT %*% Lambda.sub
      Psi_trapz <- integrated_latent(Psi, s)
      Phi_trapz <- integrated_latent(Phi, t)
      H.sub <- H[[1]][, , X]
      H.sub <- matrix(1/diag(H.sub), nrow = q1, ncol = q2)
      S.sub <- 1/(Sigma[[1]][, , X])
      
      cov_long <- extract_covlatent2(
        Gamma.sub, t(S.sub), t(H.sub), splineS, splineS_int_sqrt,
        splineS_int_sqrt_inv, splineT_int, Phi_trapz)
      return(cov_long)
    })
    for(Y in 2:nchain){
      marginal_cov_long <- append(marginal_cov_long, lapply(X = 1:1500, function(X){
        Lambda.sub <- Lambda[[Y]][, , X]
        Gamma.sub <- Gamma[[Y]][, , X]
        Psi <- splineS %*% Gamma.sub
        Phi <- splineT %*% Lambda.sub
        Psi_trapz <- integrated_latent(Psi, s)
        Phi_trapz <- integrated_latent(Phi, t)
        H.sub <- H[[Y]][, , X]
        H.sub <- matrix(1/diag(H.sub), nrow = q1, ncol = q2)
        S.sub <- 1/(Sigma[[Y]][, , X])
        
        cov_long <- extract_covlatent2(
          Gamma.sub, t(S.sub), t(H.sub), splineS, splineS_int_sqrt,
          splineS_int_sqrt_inv, splineT_int, Phi_trapz)
        return(cov_long)}))
    }
    
    # Posteriors summaries related with marginal longitudinal covariance
    post_summs$K_s$point_ests$K_s.tilde <- Reduce("+", marginal_cov_long)/length(marginal_cov_long)
    MVD.marginal_cov_long <- MVD(marginal_cov_long); gc()
    MVD.marginal_cov_long <- MVD.marginal_cov_long[order(-MVD.marginal_cov_long$rank),]
    post_summs$ranks$MVD.marginal_cov_long <- MVD.marginal_cov_long
    
    post_summs$K_s$point_ests$K_s.med.tilde <- marginal_cov_long[[MVD.marginal_cov_long$id[[1]]]]
    
    # Extract marginal longitudinal eigenfunctions
    eigen_long <- lapply(X = 1:6000, function(X){
      K <- marginal_cov_long[[X]]
      eigen.K <- eigen(K)
      tau_m <- eigen.K$values[1:2] * gridsize_s
      psi_m <- eigen.K$vectors[, c(1:2)]
      return(list(tau_m = tau_m, psi_m = psi_m))
    })
    
    Psi_long <- array(dim = c(SS, J, 6000))
    for(Y in 1:2){
      psi_data <- do.call(cbind, lapply(1:6000, function(X){
        eigen_long[[X]]$psi_m[, Y]
      }))
      psi_data <- eigenfn_align(psi_data)
      Psi_long[ , Y, ] <- psi_data
    }
    # Extract marginal longitudinal eigenvalues
    tau_long <- do.call(cbind, lapply(X = 1:6000, function(X){
      eigen_long[[X]]$tau_m
    }))
    
    # Save memory
    rm(eigen_long); gc()
    
    # Calculate covariance kernel
    tmp3 <- kronecker(Bs1, Bt1)
    cov_kernel <- lapply(X = 1:1500, function(X){
      tmp1 <- kronecker(Gamma[[1]][ , , X], Lambda[[1]][, , X])
      tmp2 <- tmp1 %*% tcrossprod(diag(diag(1/(H[[1]][, , X]))), tmp1)
      Omega <- tmp2 + diag(as.vector(1/(Sigma[[1]][, , X])))
      K <- tmp3 %*% tcrossprod(Omega, tmp3)
      return(K)
    })
    for(Y in 2:nchain){
      cov_kernel <- append(cov_kernel, lapply(X = 1:1500, function(X){
        tmp1 <- kronecker(Gamma[[Y]][ , , X], Lambda[[Y]][, , X])
        tmp2 <- tmp1 %*% tcrossprod(diag(diag(1/(H[[Y]][, , X]))), tmp1)
        Omega <- tmp2 + diag(as.vector(1/(Sigma[[Y]][, , X])))
        K <- tmp3 %*% tcrossprod(Omega, tmp3)
        return(K)
      }))
    }
    gc()
    # Poseriors summaries related with covariance kernel
    start.time <- Sys.time()
    post_summs$K$point_ests$K.bar <- Reduce("+", cov_kernel)/length(cov_kernel)
    post_summs$K_t$point_ests$K_t.bar <- LFBayes::get_marginal_func(post_summs$K$point_ests$K.bar, SS, TT) 
    post_summs$K_s$point_ests$K_s.bar <- LFBayes::get_marginal_long(post_summs$K$point_ests$K.bar, SS, TT)
    
    MVD.cov_kernel <- MVD(cov_kernel); gc()
    end.time <- Sys.time()
    end.time - start.time
    # Note: 1.5 hours for this step on Windows server (simulation)
    MVD.cov_kernel <- MVD.cov_kernel[order(-MVD.cov_kernel$rank),]
    
    post_summs$K$point_ests$K.med.bar <- cov_kernel[[MVD.cov_kernel$id[[1]]]]
    post_summs$K_t$point_ests$K_t.med.bar <- marginal_cov_func[[MVD.cov_kernel$id[[1]]]]
    post_summs$K_s$point_ests$K_s.med.bar <- marginal_cov_long[[MVD.cov_kernel$id[[1]]]]
    post_summs$ranks$MVD.cov_kernel <- MVD.cov_kernel
    
    # Save memory
    rm(marginal_cov_func, marginal_cov_long, cov_kernel); gc()

    #############################################################################
    # 5. Eigencomponents summaries 
    #############################################################################
    psi_1 <- -2^{1/2} * cos(2 * pi * s); psi_1 <- psi_1/sqrt(sum(psi_1^2))
    psi_2 <- 2^{1/2} * sin(2 * pi * s); psi_2 <- psi_2/sqrt(sum(psi_2^2))
    psi.truth <- list(psi_1, psi_2) 
    phi_1 <- -2^{1/2} * cos(3 * pi * s); phi_1 <- phi_1/sqrt(sum(phi_1^2))
    phi_2 <- -2^{1/2} * sin(3 * pi * s); phi_2 <- phi_2/sqrt(sum(phi_2^2))
    phi.truth <- list(phi_1, phi_2)
    
    eigen.K_s.tilde <- eigen(post_summs$K_s$point_ests$K_s.tilde)
    eigen.K_t.tilde <- eigen(post_summs$K_t$point_ests$K_t.tilde)
    # eigen.K_s.bar <- eigen(splineS_int_sqrt %*% post_summs$K_s$point_ests$K_s.bar %*% splineS_int_sqrt)
    eigen.K_s.bar <- eigen(post_summs$K_s$point_ests$K_s.bar)
    eigen.K_s.bar$values <- (eigen.K_s.bar$values)
    eigen.K_t.bar <- eigen(post_summs$K_t$point_ests$K_t.bar)
    eigen.K_t.bar$values <- eigen.K_t.bar$values
    
    for(j in 1:J){
      ## a. Point estimates
      data <- Re(Psi_long[, j, ])
      if(IMSE.flip_eigen_index(-1 * rowMeans(data), psi.truth[[j]])){
        data <- -1 * data
        Psi_long[, j, ] <- data
      } 
      MBD.psi_j <- MBD(data)
      MBD.psi_j <- MBD.psi_j[order(-MBD.psi_j$rank),]
      post_summs$ranks$MBD.psi_j[[j]] <- MBD.psi_j
      post_summs$psi_j[[j]]$point_ests$psi_j.hat <- IMSE.flip_eigen(rowMeans(data), psi.truth[[j]])
      post_summs$psi_j[[j]]$point_ests$psi_j.med.hat <- IMSE.flip_eigen(data[ , MBD.psi_j$id[[1]]], psi.truth[[j]])
      post_summs$psi_j[[j]]$point_ests$psi_j.med.tilde <- IMSE.flip_eigen(data[ , MVD.marginal_cov_long$id[[1]]], psi.truth[[j]])
      post_summs$psi_j[[j]]$point_ests$psi_j.med.bar <- IMSE.flip_eigen(data[ , MVD.cov_kernel$id[[1]]], psi.truth[[j]])
      
      eigen.K_s.tilde$vectors[, j] <- IMSE.flip_eigen(eigen.K_s.tilde$vectors[, j], psi.truth[[j]])
      post_summs$psi_j[[j]]$point_ests$psi_j.tilde <- eigen.K_s.tilde$vectors[, j]
      eigen.K_s.bar$vectors[, j] <- IMSE.flip_eigen(eigen.K_s.bar$vectors[, j], psi.truth[[j]])
      post_summs$psi_j[[j]]$point_ests$psi_j.bar <- eigen.K_s.bar$vectors[, j]
      
      post_summs$tau_j[[j]]$point_ests$tau_j.hat <- mean(tau_long[j, ])
      post_summs$tau_j[[j]]$point_ests$tau_j.tilde <- eigen.K_s.tilde$values[j]
      post_summs$tau_j[[j]]$point_ests$tau_j.bar <- eigen.K_s.bar$values[j]
      
      ## b. Credible intervals
      data <- t(data)
      post_summs$psi_j[[j]]$credible_ints$simultaneous_parametric <- simultaneous_parametric(data,
                                                                                             alpha = alpha)
      post_summs$psi_j[[j]]$credible_ints$MBD <- depth_credible_interval(rank = MBD.psi_j, 
                                                                         data = data,
                                                                         alpha = alpha)
      post_summs$psi_j[[j]]$credible_ints$MVD.marg <- depth_credible_interval(rank = MVD.marginal_cov_long, 
                                                                              data = data,
                                                                              alpha = alpha)
      post_summs$psi_j[[j]]$credible_ints$MVD.kernel <- depth_credible_interval(rank = MVD.cov_kernel, 
                                                                                data = data,
                                                                                alpha = alpha)
      
    }
    
    for(k in 1:K){
      ## a. Point estimates
      data <- Re(Phi_func[, k ,])
      if(IMSE.flip_eigen_index(-1 * rowMeans(data), phi.truth[[k]])){
        data <- -1 * data
        Phi_func[, k ,] <- data
      } 
      MBD.phi_k <- MBD(data)
      MBD.phi_k <- MBD.phi_k[order(-MBD.phi_k$rank),]
      post_summs$ranks$MBD.phi_k[[k]] <- MBD.phi_k
      post_summs$phi_k[[k]]$point_ests$phi_k.hat <- IMSE.flip_eigen(rowMeans(data), phi.truth[[k]])
      post_summs$phi_k[[k]]$point_ests$phi_k.med.hat <- IMSE.flip_eigen(data[ , MBD.phi_k$id[[1]]], phi.truth[[k]])
      post_summs$phi_k[[k]]$point_ests$phi_k.med.tilde <- IMSE.flip_eigen(data[ , MVD.marginal_cov_func$id[[1]]], phi.truth[[k]])
      post_summs$phi_k[[k]]$point_ests$phi_k.med.bar <- IMSE.flip_eigen(data[ , MVD.cov_kernel$id[[1]]], phi.truth[[k]])
      
      eigen.K_t.tilde$vectors[, k] <- IMSE.flip_eigen(eigen.K_t.tilde$vectors[, k], phi.truth[[k]])
      post_summs$phi_k[[k]]$point_ests$phi_k.tilde <- eigen.K_t.tilde$vectors[, k]
      eigen.K_t.bar$vectors[, k] <- IMSE.flip_eigen(eigen.K_t.bar$vectors[, k], phi.truth[[k]])
      post_summs$phi_k[[k]]$point_ests$phi_k.bar <- eigen.K_t.bar$vectors[, k]
      
      post_summs$vartheta_k[[k]]$point_ests$vartheta_k.hat <- mean(vartheta_func[k, ])
      post_summs$vartheta_k[[k]]$point_ests$vartheta_k.tilde <- eigen.K_t.tilde$values[k]
      post_summs$vartheta_k[[k]]$point_ests$vartheta_k.bar <- eigen.K_t.bar$values[k]
      
      ## b. Credible intervals
      post_summs$phi_k[[k]]$credible_ints$simultaneous_parametric <- simultaneous_parametric(t(data),
                                                                                             alpha = alpha)
      
      post_summs$phi_k[[k]]$credible_ints$MBD <- depth_credible_interval(rank = MBD.phi_k, 
                                                                         data = t(data),
                                                                         alpha = alpha)
      post_summs$phi_k[[k]]$credible_ints$MVD.marg <- depth_credible_interval(rank = MVD.marginal_cov_func, 
                                                                              data = t(data),
                                                                              alpha = alpha)
      post_summs$phi_k[[k]]$credible_ints$MVD.kernel <- depth_credible_interval(rank = MVD.cov_kernel, 
                                                                                data = t(data),
                                                                                alpha = alpha)
      
    }
    
    post_summs <- append(post_summs, list(mu_sample = mu, 
                                          tau_long = tau_long, 
                                          vartheta_func = vartheta_func, 
                                          Psi_long = Psi_long, 
                                          Phi_func = Phi_func))
  }
    
  
  