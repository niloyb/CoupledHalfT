# Functions for a blocked Gibbs sampler with half-t priors
# require(Rcpp)
# require(inline)
# require(RcppEigen)

## Rcpp functions ##
# crossprodCpp <- '
# using Eigen::Map;
# using Eigen::MatrixXd;
# using Eigen::Lower;
# const Map<MatrixXd> A(as<Map<MatrixXd> >(AA));
# const int m(A.rows()), n(A.cols());
# MatrixXd AtA(MatrixXd(n, n).setZero().selfadjointView<Lower>().rankUpdate(A.adjoint()));
# return wrap(AtA);
# '
# #' fcprd
# #' @description Calculates t(AA)*AA using RcppEigen
# #' @param AA matrix
# #' @return t(AA)*AA
# fcprd <- cxxfunction(signature(AA = "matrix"), crossprodCpp, "RcppEigen")

## Matrix calculation helper functions ##
#' X_eta_tX
#' @description Calculates X * Diag(1/eta) * t(X_transpose)
#' @param eta vector of length p
#' @param X matrix of length n by p
#' @param X_transpose Pre-calculated transpose of X
#' @return X Diag(1/eta) t(X_transpose)
#' @export
X_eta_tX <- function(eta, X, X_transpose){
  # NOTE: (1) For MacOS with veclib BLAS, crossprod is fast via multit-hreading
  # return(crossprod(X_transpose*c(1/eta)^0.5))
  return(fcprd(X_transpose*c(1/eta)^0.5))
}
#' M_matrix
#' @description Calculates I + X * Diag(1/eta) * t(X_transpose)
#' @param xi positive scalar
#' @param eta vector of length p
#' @param X_eta_tX_matrix X * Diag(1/eta) * t(X_transpose), matrix n by n
#' @param n positive integer
#' @return I + X Diag(1/eta) t(X_transpose)
#' @export
M_matrix <- function(xi, eta, X_eta_tX_matrix, n){
  if(length(eta)==0) return(diag(n))
  return(diag(n) + (xi^-1)*X_eta_tX_matrix)
}

## xi update given eta ##
# Unnormalized posterior pdf of log(xi)
#' log_ratio
#' @description Unnormalized posterior pdf of log(xi)
#' @param xi positive scalar
#' @param eta length p vector
#' @param X_eta_tX_matrix X * Diag(1/eta) * t(X_transpose), matrix n by n
#' @param y length n vector
#' @param a0 positive scalar
#' @param b0 positive scalar
#' @param xi_interval support of prior distribution of xi
#' @return Unnormalized posterior pdf of log(xi)
#' @export
log_ratio <- function(xi, eta, X_eta_tX_matrix, y, a0, b0, xi_interval)
{
  n <- length(y)
  M <- M_matrix(xi,eta,X_eta_tX_matrix,n)
  chol_M <- chol(M)
  log_det_M <- 2*sum(log(diag(chol_M)))
  M_inverse <- chol2inv(chol_M)
  ssr <- b0 + t(y)%*%((M_inverse)%*%y)
  log_likelihood <- -0.5*log_det_M -0.5*(n+a0)*log(ssr)
  
  if((xi_interval[1]<=xi)&(xi<=xi_interval[2])){
    log_prob <- -log(sqrt(xi)*(1+xi))
  } else {
    log_prob <- log(0)
  }
  
  return(list('log_likelihood'=log_likelihood+log_prob,
              'ssr' = ssr, 'M_matrix_inverse' = M_inverse))
}
# Unnormalized posterior pdf of log(xi) for approximate MCMC
#' log_ratio_approx
#' @description Unnormalized posterior pdf of log(xi) for approximate MCMC
#' @param xi positive scalar
#' @param eta length p vector
#' @param X n by p matrix
#' @param X_transpose Pre-calculated transpose of X
#' @param y length n vector
#' @param a0 positive scalar
#' @param b0 positive scalar
#' @param xi_interval support of prior distribution of xi
#' @param active_set length p vector of booleans
#' @return Unnormalized posterior pdf of log(xi)
#' @export
log_ratio_approx <- function(xi, eta, X, X_transpose, y, a0, b0, active_set,
                             xi_interval)
{
  n <- length(y)
  
  if((xi_interval[1]<=xi)&(xi<=xi_interval[2])){
    log_prob <- -log(sqrt(xi)*(1+xi))
  } else {
    log_prob <- log(0)
  }
  
  if (sum(active_set)==0)
  {
    M_inverse <- diag(n)
    ssr <- b0 + sum(y^2)
    log_likelihood <- -0.5*(n+a0)*log(ssr)
  } else{
    eta <- eta[active_set]
    X <- X[,active_set,drop=F]
    X_transpose <- X_transpose[active_set,,drop=F]
    if(sum(active_set)==1){
      woodbury_matrix_part <- xi*eta + cpp_prod(X_transpose, X)
    } else{
      woodbury_matrix_part <- xi*diag(eta) + cpp_prod(X_transpose, X)
    }
    woodbury_matrix_part_inverse <- chol2inv(chol(woodbury_matrix_part))
    M_inverse <- diag(n) -
      cpp_prod(cpp_prod(X,woodbury_matrix_part_inverse),X_transpose)
    log_det_M <- sum( log( xi^(-1)*(svd((eta^(-0.5))*X_transpose)$d)^2 + 1 ) )
    ssr <- b0 + t(y)%*%((M_inverse)%*%y)
    log_likelihood <- -0.5*log_det_M -0.5*(n+a0)*log(ssr)
  }
  return(list('log_likelihood'=log_likelihood+log_prob, 'ssr' = ssr,
              'M_matrix_inverse' = M_inverse))
}
# Metropolis-Hastings update of xi given eta
#' xi_update
#' @description Metropolis-Hastings update of xi
#' @param current_xi positive scalar
#' @param eta length p vector
#' @param X n by p matrix
#' @param X_transpose Pre-calculated transpose of X
#' @param y length n vector
#' @param a0 positive scalar
#' @param b0 positive scalar
#' @param xi_interval support of prior distribution of xi
#' @param std_MH standard deviation of log-normal MH proposal
#' @param approximate_algo_delta approximate MCMC error (non-negative scalar)
#' @return Metropolis-Hastings update of xi given eta
#' @export
xi_update <- function(current_xi, eta, X, X_transpose, y, a0, b0, std_MH,
                      approximate_algo_delta=0, fixed=FALSE, xi_interval)
{
  n <- length(y)
  if (fixed==TRUE){
    min_xi <- current_xi
    active_set <- ((min_xi*eta)^(-1) > approximate_algo_delta)
    
    # Matrix calculations
    if (sum(active_set)>n) # Woodbury inversion when sum(active_set)>n
    {
      X_eta_tX_matrix <- X_eta_tX(eta[active_set], X[ ,active_set,drop=F],
                                  X_transpose[active_set, ,drop=F])
      log_ratio_current_and_ssr <- log_ratio(current_xi, eta[active_set],
                                             X_eta_tX_matrix, y, a0, b0, xi_interval)
    } else { # Woodbury inversion when sum(active_set)<n
      log_ratio_current_and_ssr <-
        log_ratio_approx(current_xi, eta, X, X_transpose, y, a0, b0, active_set, xi_interval)
    }
  } else {
    proposed_xi <- exp(rnorm(1, log(current_xi), std_MH))
    min_xi <- min(current_xi, proposed_xi)
    active_set <- ((min_xi*eta)^(-1) > approximate_algo_delta)
    
    # Matrix calculations
    if (sum(active_set)>n) # Woodbury inversion when sum(active_set)>n
    {
      X_eta_tX_matrix <- X_eta_tX(eta[active_set],X[ ,active_set,drop=F],
                                  X_transpose[active_set, ,drop=F])
      log_ratio_current_and_ssr <- log_ratio(current_xi, eta[active_set],
                                             X_eta_tX_matrix, y, a0, b0, xi_interval)
      log_ratio_proposed_and_ssr <- log_ratio(proposed_xi, eta[active_set],
                                              X_eta_tX_matrix, y, a0, b0, xi_interval)
    } else { # Woodbury inversion when sum(active_set)<n
      log_ratio_current_and_ssr <-
        log_ratio_approx(current_xi, eta, X, X_transpose, y, a0, b0, active_set, xi_interval)
      log_ratio_proposed_and_ssr <-
        log_ratio_approx(proposed_xi, eta, X, X_transpose, y, a0, b0, active_set, xi_interval)
    }
    
    # MH accept- reject step
    log_accept_prob <-
      (log_ratio_proposed_and_ssr$log_likelihood -
         log_ratio_current_and_ssr$log_likelihood) +
      (log(proposed_xi)-log(current_xi))
    
    if (log(runif(1))<log_accept_prob)
    {
      current_xi <- proposed_xi
      log_ratio_current_and_ssr <- log_ratio_proposed_and_ssr
    }
  }
  return(list('xi'=current_xi, 'ssr' = log_ratio_current_and_ssr$ssr,
              'M_matrix' = log_ratio_current_and_ssr$M_matrix_inverse,
              'active_set' = active_set))
}

## sigma^2 update step given eta and xi
#' sigma2_update
#' @description sigma^2 update step given eta and xi
#' @param xi positive scalar
#' @param eta length p vector
#' @param n positive integer
#' @param ssr positive scalar
#' @param a0 positive scalar
#' @param b0 positive scalar
#' @param sigma2_fixed_value NULL for sigma2 varying, positive scalar when fixed
#' @return sigma^2 update step given eta and xi
#' @export
sigma2_update <- function(xi, eta, n, ssr, a0, b0, sigma2_fixed_value=NULL)
{
  if(!is.null(sigma2_fixed_value)) return(sigma2_fixed_value)
  return(1/(rgamma(1, shape = (n+a0)/2, rate = (ssr)/2)))
}

## beta update given xi, eta, sigma2
#' beta_update
#' @description beta given xi, eta, sigma2 update using algo of Bhattacharya
#' @param xi positive scalar
#' @param sigma2 positive scalar
#' @param eta length p vector
#' @param X n by p matrix
#' @param X_transpose Pre-calculated transpose of X
#' @param y length n vector
#' @param M_matrix_inverse n by n matrix
#' @param active_set length p vector of booleans
#' @return beta update given xi, eta, sigma2
#' @export
beta_update <- function(xi, sigma2, eta, X, X_transpose, y,
                        M_matrix_inverse, active_set)
{
  p <- length(eta)
  n <- length(y)
  u = rnorm(p, 0, 1)
  u = (sqrt(xi*eta)^-1)*u
  v = X%*%u + c(rnorm(n,0,1))
  v_star <- M_matrix_inverse%*%(y/sqrt(sigma2) - v)
  U <- (xi^-1)* ( ((eta[active_set])^(-1))*(X_transpose[active_set, ,drop=F]) )
  u[active_set] <- u[active_set] + U%*%v_star
  beta <- sqrt(sigma2)*(u)
  return(beta)
}


## Full blocked Gibbs samplers ##
#' half_t_kernel
#' @description Full blocked Gibbs kernel for half-t priors
#' @param X n by p matrix
#' @param X_transpose Pre-calculated transpose of X
#' @param y length n vector
#' @param a0 positive scalar
#' @param b0 positive scalar
#' @param xi_interval support of prior distribution of xi
#' @param std_MH standard deviation of log-normal MH proposal
#' @param xi_current current xi value (positive scalar)
#' @param sigma2_current current sigma2 value (positive scalar)
#' @param beta_current current beta value (vector of length p)
#' @param eta_current current eta value (vector of length p)
#' @param approximate_algo_delta approximate MCMC error (non-negative scalar)
#' @param nrepeats_eta number of slice sampling steps
#' @param verbose boolean for printing/ not printing run time
#' @param xi_fixed boolean for fixing / not fixing xi
#' @param sigma2_fixed boolean for fixing / not fixing sigma2
#' @param t_dist_df degree of freedom v for Half-t(v). Take v >=1.
#' @return (beta, eta, sigma2, xi) sampled from the half-t prior blocked Gibbs kernel
#' @export
half_t_kernel <-
  function(X, X_transpose, y, a0=1, b0=1, std_MH=0.8,
           xi_current, sigma2_current,
           beta_current, eta_current, approximate_algo_delta,
           nrepeats_eta = 1, verbose = FALSE,
           xi_fixed=FALSE, sigma2_fixed=FALSE, t_dist_df, 
           xi_interval)
  {
    if (verbose) ptm <- proc.time()
    n <- dim(X)[1]
    p <- dim(X)[2]
    # xi update
    xi_new <- xi_update(xi_current, eta_current, X, X_transpose, y, a0, b0, std_MH,
                        approximate_algo_delta, fixed=xi_fixed, xi_interval=xi_interval)
    # sigma2 update
    sigma2_fixed_value <- NULL
    if(sigma2_fixed==TRUE) sigma2_fixed_value <- sigma2_current
    sigma2_new <- sigma2_update(xi_new$xi, eta_current, n, xi_new$ssr, a0, b0,
                                sigma2_fixed_value)
    # beta update
    beta_new <- beta_update(xi_new$xi, sigma2_new, eta_current, X, X_transpose, y,
                            xi_new$M_matrix, xi_new$active_set)
    # eta update
    eta_new <- eta_update_half_t(xi_new$xi, sigma2_new, beta_new, eta_current,
                                 t_dist_df, nrepeats_eta)
    
    if (verbose) print(proc.time()[3]-ptm[3])
    output <- list( 'beta_samples'=beta_new, 'eta_samples'=eta_new,
                    'sigma2_samples'=sigma2_new, 'xi_samples'=xi_new$xi)
    return(output)
  }

#' half_t_mcmc
#' @description Blocked Gibbs MCMC for Bayesian shrinkage with half-t priors
#' @param chain_length length of Markov chain
#' @param burnin Markov chain burn in
#' @param X n by p matrix
#' @param X_transpose Pre-calculated transpose of X
#' @param y length n vector
#' @param a0 positive scalar
#' @param b0 positive scalar
#' @param xi_interval support of prior distribution of xi
#' @param rinit Initial distribution
#' @param xi_current current xi value (positive scalar)
#' @param sigma2_current current sigma2 value (positive scalar)
#' @param beta_current current beta value (vector of length p)
#' @param eta_current current eta value (vector of length p)
#' @param approximate_algo_delta approximate MCMC error (non-negative scalar)
#' @param nrepeats_eta number of slice sampling steps
#' @param verbose boolean for printing/ not printing run time
#' @param xi_fixed boolean for fixing / not fixing xi
#' @param sigma2_fixed boolean for fixing / not fixing sigma2
#' @param t_dist_df degree of freedom v for Half-t(v). Take v >=1.
#' @return (beta, eta, sigma2, xi) sampled from the half-t prior blocked Gibbs kernel
#'@export
half_t_mcmc <-
  function(chain_length, burnin, X, X_transpose, y, a0=1, b0=1, std_MH=0.8,
           rinit=NULL, approximate_algo_delta=0, nrepeats_eta = 1,
           verbose = FALSE, xi_fixed=FALSE, sigma2_fixed=FALSE, t_dist_df,
           xi_interval=c(0,Inf))
  {
    n <- dim(X)[1]
    p <- dim(X)[2]
    if(is.null(rinit)){
      # Initializing from the prior
      rinit <- function(){
        # xi <- (1/rt(1, df=1))^2
        xi <- tan(runif(1)*(atan(xi_interval[2]^0.5)-atan(xi_interval[1]^0.5))+
                    atan(xi_interval[1]^0.5))^2
        sigma2 <- 1/rgamma(1, shape = a0/2, rate = b0/2)
        eta <- (1/rt(p, df=t_dist_df))^2
        beta <- rnorm(p)*sqrt(sigma2/(xi*eta))
        return(list(xi = xi, sigma2 = sigma2, beta = beta, eta = eta))
      }
    }
    # Initializing first chain
    chain <- rinit()
    xi_samples <- rep(NA, (chain_length-burnin))
    sigma2_samples <- rep(NA, (chain_length-burnin))
    beta_samples <- matrix(NA, nrow=(chain_length-burnin), ncol=p)
    eta_samples <- matrix(NA, nrow=(chain_length-burnin), ncol=p)
    xi_samples[1] <- xi_current <- chain$xi
    sigma2_samples[1] <- sigma2_current <- chain$sigma2
    beta_samples[1,] <- beta_current <- chain$beta
    eta_samples[1,] <- eta_current <- chain$eta
    
    i <- 1
    while (i <= chain_length)
    {
      output <-
        half_t_kernel(X, X_transpose, y, a0=a0, b0=b0, std_MH=std_MH,
                      xi_current, sigma2_current, beta_current, eta_current,
                      approximate_algo_delta, nrepeats_eta = nrepeats_eta,
                      verbose = verbose, xi_fixed=xi_fixed,
                      sigma2_fixed=sigma2_fixed, t_dist_df, xi_interval=xi_interval)
      xi_current <- output$xi_samples
      sigma2_current <- output$sigma2_samples
      beta_current <- output$beta_samples
      eta_current <- output$eta_samples
      if(i > burnin)
      {
        xi_samples[(i-burnin)] <- xi_current
        sigma2_samples[(i-burnin)] <- sigma2_current
        beta_samples[(i-burnin),] <- beta_current
        eta_samples[(i-burnin),] <- eta_current
      }
      if (verbose) print(i)
      i <- i + 1
    }
    return(list( 'beta_samples'=beta_samples, 'eta_samples'=eta_samples,
                 'sigma2_samples'=sigma2_samples, 'xi_samples'=xi_samples))
  }





