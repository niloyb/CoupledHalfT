# Functions for the coupled blocked Gibbs sampler with half-t priors
require(lubridate)

## Couplings for xi update given eta ##
# Coupled Metropolis-Hastings update of xi given eta
#' xi_update
#' @description Coupled Metropolis-Hastings update of xi given eta
#' @param current_xi_1,current_xi_2 current xi values (positive scalars)
#' @param eta_1,eta_2 current eta values (vector length p)
#' @param X n by p matrix
#' @param X_transpose Pre-calculated transpose of X
#' @param y length n vector
#' @param a0 positive scalar
#' @param b0 positive scalar
#' @param xi_interval support of prior distribution of xi
#' @param std_MH standard deviation of log-normal MH proposal
#' @param approximate_algo_delta approximate MCMC error (non-negative scalar)
#' @param epsilon_xi threshold (when above/below use common random numbers/ maximal coupling of proposals)
#' @return Coupled Metropolis-Hastings update of xi given eta
#' @export
crn_max_xi_coupling <- function(current_xi_1, eta_1, current_xi_2, eta_2,
                                X, X_transpose, y, a0, b0, std_MH,
                                approximate_algo_delta=0, epsilon_xi=0, fixed=FALSE,
                                xi_interval)
{
  if(fixed==TRUE){ # When xi is fixed in the Gibbs sampler
    min_xi_1 <- current_xi_1
    min_xi_2 <- current_xi_2
    active_set_1 <- ((min_xi_1*eta_1)^(-1) > approximate_algo_delta)
    active_set_2 <- ((min_xi_2*eta_2)^(-1) > approximate_algo_delta)

    if (sum(active_set_1)>n)
    {
      X_eta_tX_matrix_1 <-
        X_eta_tX(eta_1[active_set_1],X[, active_set_1, drop=F],
                 X_transpose[active_set_1, , drop=F])
      log_ratio_current_ssr_matrixinv_1 <-
        log_ratio(current_xi_1, eta_1[active_set_1], X_eta_tX_matrix_1, y, a0, b0, xi_interval)
    } else {
      log_ratio_current_ssr_matrixinv_1 <-
        log_ratio_approx(current_xi_1, eta_1, X, X_transpose, y, a0, b0, active_set_1, xi_interval)
    }

    if (sum(active_set_2)>n)
    {
      X_eta_tX_matrix_2 <-
        X_eta_tX(eta_2[active_set_2],X[, active_set_2, , drop=F], X_transpose[active_set_2, , drop=F])
      log_ratio_current_ssr_matrixinv_2 <-
        log_ratio(current_xi_2, eta_2[active_set_2], X_eta_tX_matrix_2, y, a0, b0, xi_interval)
    } else {
      log_ratio_current_ssr_matrixinv_2 <-
        log_ratio_approx(current_xi_2, eta_2, X, X_transpose, y, a0, b0, active_set_2, xi_interval)
    }
  } else { # When xi is varying in the Gibbs sampler
    standard_normal <- rnorm(1, mean = 0, sd = 1)
    log_proposed_xi_1 <- standard_normal*sqrt(std_MH) + log(current_xi_1)

    relative_error_delta <- abs(log(current_xi_1)-log(current_xi_2))
    # print(relative_error_delta)
    if ((0 < relative_error_delta) & (relative_error_delta < epsilon_xi)) {
      # using max coupling to get the proposal in the MH-algo
      if( dnorm(log_proposed_xi_1, mean = log(current_xi_1), sd = std_MH, log = TRUE) +
          log(runif(1)) < dnorm(log_proposed_xi_1, mean = log(current_xi_2), sd = std_MH, log = TRUE) )
      {
        log_proposed_xi_2 <- log_proposed_xi_1
      } else {
        reject <- TRUE
        y_proposal <- NA
        attempts <- 0
        while(reject){
          attempts <- attempts + 1
          # print(attempts)
          y_proposal <- rnorm(1, mean = log(current_xi_2), sd = std_MH)
          reject <- ( dnorm(y_proposal, mean = log(current_xi_2), sd = std_MH, log = TRUE) +
                        log(runif(1)) < dnorm(y_proposal, mean = log(current_xi_1), sd = std_MH, log = TRUE) )
        }
        log_proposed_xi_2 <- y_proposal
      }
    } else {
      # using common random numbers to get the proposal in the MH-algo
      log_proposed_xi_2 <- standard_normal*sqrt(std_MH) + log(current_xi_2)
    }

    proposed_xi_1 <- exp(log_proposed_xi_1)
    proposed_xi_2 <- exp(log_proposed_xi_2)

    min_xi_1 <- min(current_xi_1, proposed_xi_1)
    min_xi_2 <- min(current_xi_2, proposed_xi_2)
    active_set_1 <- ((min_xi_1*eta_1)^(-1) > approximate_algo_delta)
    active_set_2 <- ((min_xi_2*eta_2)^(-1) > approximate_algo_delta)

    if (sum(active_set_1)>n)
    {
      X_eta_tX_matrix_1 <- X_eta_tX(eta_1[active_set_1],X[, active_set_1, drop=F], X_transpose[active_set_1, , drop=F])
      log_ratio_current_ssr_matrixinv_1 <- log_ratio(current_xi_1, eta_1[active_set_1], X_eta_tX_matrix_1, y, a0, b0, xi_interval)
      log_ratio_proposed_ssr_matrixinv_1 <- log_ratio(proposed_xi_1, eta_1[active_set_1], X_eta_tX_matrix_1, y, a0, b0, xi_interval)
    } else {
      log_ratio_current_ssr_matrixinv_1 <- log_ratio_approx(current_xi_1, eta_1, X, X_transpose, y, a0, b0, active_set_1, xi_interval)
      log_ratio_proposed_ssr_matrixinv_1 <- log_ratio_approx(proposed_xi_1, eta_1, X, X_transpose, y, a0, b0, active_set_1, xi_interval)
    }

    if (sum(active_set_2)>n)
    {
      X_eta_tX_matrix_2 <- X_eta_tX(eta_2[active_set_2],X[, active_set_2, , drop=F], X_transpose[active_set_2, , drop=F])
      log_ratio_current_ssr_matrixinv_2 <- log_ratio(current_xi_2, eta_2[active_set_2], X_eta_tX_matrix_2, y, a0, b0, xi_interval)
      log_ratio_proposed_ssr_matrixinv_2 <- log_ratio(proposed_xi_2, eta_2[active_set_2], X_eta_tX_matrix_2, y, a0, b0, xi_interval)
    } else {
      log_ratio_current_ssr_matrixinv_2 <- log_ratio_approx(current_xi_2, eta_2, X, X_transpose, y, a0, b0, active_set_2, xi_interval)
      log_ratio_proposed_ssr_matrixinv_2 <- log_ratio_approx(proposed_xi_2, eta_2, X, X_transpose, y, a0, b0, active_set_2, xi_interval)
    }

    log_u <- log(runif(1))

    log_accept_prob_1 <- (log_ratio_proposed_ssr_matrixinv_1$log_likelihood - log_ratio_current_ssr_matrixinv_1$log_likelihood) + (log(proposed_xi_1)-log(current_xi_1))
    if (log_u<log_accept_prob_1){
      current_xi_1 <- proposed_xi_1
      log_ratio_current_ssr_matrixinv_1 <- log_ratio_proposed_ssr_matrixinv_1
    }

    log_accept_prob_2 <- (log_ratio_proposed_ssr_matrixinv_2$log_likelihood - log_ratio_current_ssr_matrixinv_2$log_likelihood) + (log(proposed_xi_2)-log(current_xi_2))
    if (log_u<log_accept_prob_2){
      current_xi_2 <- proposed_xi_2
      log_ratio_current_ssr_matrixinv_2 <- log_ratio_proposed_ssr_matrixinv_2
    }
  }
  return(list('xi_values'=c(current_xi_1, current_xi_2), 'log_ratio_ssr_matrix_inv_1'=log_ratio_current_ssr_matrixinv_1,  'log_ratio_ssr_matrix_inv_2'=log_ratio_current_ssr_matrixinv_2, 'active_set_1'= active_set_1, 'active_set_2'= active_set_2))
}

## Couplings for sigma^2 update given eta ##
digamma <- function(x, alpha, beta){
  return(alpha * log(beta) - lgamma(alpha) - (alpha+1) * log(x) - beta / x)
}
rigamma <- function(n, alpha, beta){
  return(1/rgamma(n = n, shape = alpha, rate = beta))
}
rigamma_coupled <- function(alpha1, alpha2, beta1, beta2){
  x <- rigamma(1, alpha1, beta1)
  if (digamma(x, alpha1, beta1) + log(runif(1)) < digamma(x, alpha2, beta2)){
    return(c(x,x))
  } else {
    reject <- TRUE
    y <- NA
    while (reject){
      y <- rigamma(1, alpha2, beta2)
      reject <- (digamma(y, alpha2, beta2) + log(runif(1)) < digamma(y, alpha1, beta1))
    }
    return(c(x,y))
  }
}
### sigma^2 maximal coupling update step given eta and xi
sigma2_update_maximal_coupling <- function(xi_1, eta_1, xi_2, eta_2, n, ssr_1, ssr_2, a0, b0, sigma2_fixed_value=NULL)
{
  if(!is.null(sigma2_fixed_value)){
    sample <- c(sigma2_fixed_value, sigma2_fixed_value)
  } else {
    sample <- rigamma_coupled(((n+a0)/2), ((n+a0)/2), (ssr_1/2), (ssr_2/2))
  }
  return(sample)
}
### sigma^2 maximal coupling update step given eta and xi
sigma2_update_crn <- function(xi_1, eta_1, xi_2, eta_2, n, ssr_1, ssr_2, a0, b0,
                              sigma2_fixed_value=NULL)
{
  if(!is.null(sigma2_fixed_value)){
    sample <- c(1/sigma2_fixed_value, 1/sigma2_fixed_value)
  } else {
    # Common random number gamma draw
    crn_gamma <-  rgamma(1, shape = (n+a0)/2, rate = 1)
    sample <- crn_gamma / c((ssr_1/2), (ssr_2/2))
  }
  return(1/sample)
}
# Coupled update of sigma2 given eta
#' xi_update
#' @description Coupled update of sigma2 given eta
#' @param xi_1,xi_2 current xi values (positive scalars)
#' @param eta_1,eta_2 current eta values (vector length p)
#' @param n number of data points
#' @param ssr_1,ssr_2 postiive scalars
#' @param a0 positive scalar
#' @param b0 positive scalar
#' @param epsilon_sigma2 threshold (when above/below use common random numbers/ maximal coupling of proposals)
#' @param sigma2_fixed_value boolean (when TRUE/FALSE sigma2 fixed/varying)
#' @return Coupled update of sigma2 given eta
#' @export
crn_max_sigma2_coupling <- function(xi_1, eta_1, xi_2, eta_2, n, ssr_1, ssr_2, a0, b0,
                                    epsilon_sigma2=0, sigma2_fixed_value=NULL){
  relative_error_delta <- abs(ssr_1-ssr_2)
  if((0 < relative_error_delta) & (relative_error_delta < epsilon_sigma2)){
    output <- sigma2_update_maximal_coupling(xi_1, eta_1, xi_2, eta_2, n, ssr_1, ssr_2, a0, b0, sigma2_fixed_value)
  } else {
    output <- sigma2_update_crn(xi_1, eta_1, xi_2, eta_2, n, ssr_1, ssr_2, a0, b0, sigma2_fixed_value)
  }
  return(output)
}

## Common random numbers coupling of beta given xi, eta, sigma2 ##
#'@export
crn_joint_beta_update <-
  function(xi_1, sigma2_1, eta_1, xi_2, sigma2_2, eta_2,
           X, X_transpose, y, M_matrix_inverse_1, M_matrix_inverse_2,
           active_set_1, active_set_2)
{
  n <- nrow(X)
  p <- nrow(X_transpose)
  # Using same common random numbers for draws on two chains
  random_u <- rnorm(p, 0, 1)
  random_delta <- c(rnorm(n,0,1))
  u_1 = (sqrt(xi_1*eta_1)^-1)*random_u
  v_1 = X%*%u_1 + random_delta
  v_star_1 <- M_matrix_inverse_1%*%(y/sqrt(sigma2_1) - v_1)
  if(sum(active_set_1)>0){
    U_1 = (xi_1^-1)*((eta_1[active_set_1]^(-1))*(X_transpose[active_set_1, ,drop=F]))
    u_1[active_set_1] <- u_1[active_set_1] + U_1%*%v_star_1
  }
  beta_parameter_1 <- sqrt(sigma2_1)*(u_1)
  u_2 = (sqrt(xi_2*eta_2)^-1)*random_u
  v_2 = X%*%u_2 + random_delta
  v_star_2 <- M_matrix_inverse_2%*%(y/sqrt(sigma2_2) - v_2)
  if(sum(active_set_2)>0){
    U_2 = (xi_2^-1)*((eta_2[active_set_2]^(-1))*(X_transpose[active_set_2, ,drop=F]))
    u_2[active_set_2] <- u_2[active_set_2] + U_2%*%v_star_2
  }
  beta_parameter_2 <- sqrt(sigma2_2)*(u_2)
  return(cbind(beta_parameter_1, beta_parameter_2))
}

## Full coupled blocked Gibbs samplers ##
#' coupled_half_t_kernel
#' @description Coupled blocked Gibbs kernel for half-t priors
#' @param X n by p matrix
#' @param X_transpose Pre-calculated transpose of X
#' @param y length n vector
#' @param a0 positive scalar
#' @param b0 positive scalar
#' @param xi_interval support of prior distribution of xi
#' @param std_MH standard deviation of log-normal MH proposal
#' @param xi_1_current,xi_2_current current xi values (positive scalar)
#' @param sigma2_1_current,sigma2_2_current current sigma2 values (positive scalar)
#' @param beta_1_current,beta_2_current current beta values (vector of length p)
#' @param eta_1_current,eta_2_current current eta values (vector of length p)
#' @param approximate_algo_delta approximate MCMC error (non-negative scalar)
#' @param epsilon_eta eta common random numbers/ maximal coupling coupling threshold
#' @param epsilon_xi xi common random numbers/ maximal coupling coupling threshold
#' @param epsilon_sigma2 sigma2 common random numbers/ maximal coupling coupling threshold
#' @param nrepeats_eta number of slice sampling steps
#' @param verbose boolean for printing/ not printing run time
#' @param xi_fixed boolean for fixing / not fixing xi
#' @param sigma2_fixed boolean for fixing / not fixing sigma2
#' @param t_dist_df degree of freedom v>=1 for Half-t(v).
#' @return (beta1, beta2, eta1, eta2, sigma2_1, sigma2_2, xi1, xi2, metric_d)
#' @export
coupled_half_t_kernel <-
  function(X, X_transpose, y, a0=1, b0=1, std_MH=0.8,
           xi_1_current, xi_2_current, sigma2_1_current, sigma2_2_current,
           beta_1_current, beta_2_current, eta_1_current, eta_2_current,
           approximate_algo_delta=0, epsilon_eta = 0.5,
           epsilon_xi = Inf, epsilon_sigma2=Inf, nrepeats_eta=1,
           verbose = FALSE, xi_fixed=FALSE, sigma2_fixed=FALSE, t_dist_df,
           two_scale=TRUE, xi_interval)
  {
    n <- dim(X)[1]
    p <- dim(X)[2]

    if (verbose) ptm <- proc.time()
    
    if(two_scale){
      
      # (nrepeats_eta-1) steps of slice sampling
      if (is.infinite(nrepeats_eta)){ stop("Number of slice sampling must be finite") }
      if (nrepeats_eta>1){
        for (i in 1:(nrepeats_eta-1)) {
          eta_sample <-
            eta_update_half_t_crn_couple(xi_1_current, beta_1_current, eta_1_current, sigma2_1_current,
                                         xi_2_current, beta_2_current, eta_2_current, sigma2_2_current,
                                         t_dist_df)
          eta_1_current <- eta_sample[,1]
          eta_2_current <- eta_sample[,2]
        }
      }
      
      # Calculating the metric
      # When epsilon_eta >=1, relative_error_delta <= epsilon_eta always.
      if(epsilon_eta >= 1){
        relative_error_delta <- 1
      } else {
        relative_error_delta <-
          half_t_max_couple_prob(xi_1_current, beta_1_current, eta_1_current, sigma2_1_current,
                                 xi_2_current, beta_2_current, eta_2_current, sigma2_2_current,
                                 t_dist_df, iterations=1)
        # Checking overflow
        # if(typeof(relative_error_delta)=='list') return(relative_error_delta)
      }
      if (verbose) print(relative_error_delta)
      
      # Using max coupling of 1-step slice sampling when close, CRN of 1-step slice sampling when far away
      if (relative_error_delta <= epsilon_eta){ 
        eta_sample <-
          eta_update_half_t_max_couple(xi_1_current, beta_1_current, eta_1_current, sigma2_1_current,
                                       xi_2_current, beta_2_current, eta_2_current, sigma2_2_current,
                                       t_dist_df)
      } else {
        eta_sample <-
          eta_update_half_t_crn_couple(xi_1_current, beta_1_current, eta_1_current, sigma2_1_current,
                                       xi_2_current, beta_2_current, eta_2_current, sigma2_2_current,
                                       t_dist_df)
      }
    } else {
      # Switch-to-CRN coupling: now no need to calculate the metric
      relative_error_delta <- 1 
      eta_sample <-
        eta_update_half_t_max_couple_till_you_miss(xi_1_current, beta_1_current, eta_1_current, sigma2_1_current,
                                                   xi_2_current, beta_2_current, eta_2_current, sigma2_2_current,
                                                   t_dist_df)
    }
    
    eta_1_new <- eta_sample[,1]
    eta_2_new <- eta_sample[,2]
    if (verbose) print(proc.time()[3]-ptm[3])

    if (verbose) print(c('Eta components coupled:', sum(eta_1_new==eta_2_new)))
    
    if (verbose) ptm <- proc.time()
    xi_sample <-
      crn_max_xi_coupling(xi_1_current, eta_1_new, xi_2_current, eta_2_new,
                          X, X_transpose, y, a0, b0, std_MH,
                          approximate_algo_delta, epsilon_xi, fixed=xi_fixed, xi_interval)
    xi_1_new <- xi_sample$xi_values[1]
    xi_2_new <- xi_sample$xi_values[2]
    if (verbose) print(proc.time()[3]-ptm[3])

    if (verbose) ptm <- proc.time()
    sigma2_fixed_value <- NULL
    if(sigma2_fixed==TRUE){sigma2_fixed_value <- sigma2_1_current}
    sigma2_sample <-
      crn_max_sigma2_coupling(xi_1_new,eta_1_new,xi_2_new,eta_2_new,n,
                              (xi_sample$log_ratio_ssr_matrix_inv_1)$ssr,
                              (xi_sample$log_ratio_ssr_matrix_inv_2)$ssr,a0,b0,
                              epsilon_sigma2,sigma2_fixed_value)
    sigma2_1_new <- sigma2_sample[1]
    sigma2_2_new <- sigma2_sample[2]
    if (verbose) print(proc.time()[3]-ptm[3])

    if (verbose) ptm <- proc.time()
    M_inverse_1 <- (xi_sample$log_ratio_ssr_matrix_inv_1)$M_matrix_inverse
    M_inverse_2 <- (xi_sample$log_ratio_ssr_matrix_inv_2)$M_matrix_inverse
    active_set_1 <- xi_sample$active_set_1
    active_set_2 <- xi_sample$active_set_2
    beta_samples <- crn_joint_beta_update(xi_1_new, sigma2_1_new, eta_1_new,
                                          xi_2_new, sigma2_2_new, eta_2_new,
                                          X, X_transpose, y, M_inverse_1, M_inverse_2,
                                          active_set_1, active_set_2)
    beta_1_new <- beta_samples[,1]
    beta_2_new <- beta_samples[,2]
    if (verbose) print(proc.time()[3]-ptm[3])

    output <- list('beta_1_samples'=beta_1_new, 'beta_2_samples'=beta_2_new,
                   'eta_1_samples'=eta_1_new, 'eta_2_samples'=eta_2_new,
                   'sigma2_1_samples'=sigma2_1_new, 'sigma2_2_samples'=sigma2_2_new,
                   'xi_1_samples'=xi_1_new, 'xi_2_samples'=xi_2_new,
                   'metric_d'=relative_error_delta)
    return(output)
  }


#' coupled_half_t_mcmc
#' @description Coupled Blocked Gibbs MCMC for Bayesian shrinkage with half-t priors
#' @param mc_chain_size minimum length of Markov chain
#' @param X n by p matrix
#' @param X_transpose Pre-calculated transpose of X
#' @param y length n vector
#' @param a0 positive scalar
#' @param b0 positive scalar
#' @param xi_interval support of prior distribution of xi
#' @param rinit Initial distribution
#' @param approximate_algo_delta approximate MCMC error (non-negative scalar)
#' @param epsilon_eta eta common random numbers/ maximal coupling coupling threshold
#' @param epsilon_xi xi common random numbers/ maximal coupling coupling threshold
#' @param epsilon_sigma2 sigma2 common random numbers/ maximal coupling coupling threshold
#' @param nrepeats_eta number of slice sampling steps
#' @param verbose boolean for printing/ not printing run time
#' @param preallocate pre-allocated memory
#' @param max_iterations number of maximum iterations
#' @param totalduration maximum duration
#' @param xi_fixed boolean for fixing / not fixing xi
#' @param sigma2_fixed boolean for fixing / not fixing sigma2
#' @param lag Lag
#' @param t_dist_df degree of freedom v>=1 for Half-t(v)
#' @return (beta, eta, sigma2, xi) sampled from the half-t prior blocked Gibbs kernel
#' @export
coupled_half_t_mcmc <-
  function(mc_chain_size, X, X_transpose, y, a0=1, b0=1, std_MH=0.8, rinit=NULL,
           approximate_algo_delta=0, epsilon_eta=0.5, epsilon_xi=Inf, epsilon_sigma2=Inf,
           nrepeats_eta=1, verbose = FALSE, preallocate = 100, max_iterations = Inf,
           totalduration = Inf, xi_fixed=FALSE, sigma2_fixed=FALSE, lag=1, t_dist_df,
           two_scale=TRUE, xi_interval=c(0,Inf))
{
  starttime <- Sys.time() # record starting time
  n <- dim(X)[1]
  p <- dim(X)[2]
  #
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
  #
  ## Initializing chains
  m <- mc_chain_size
  xi_samples1 <-     matrix(nrow = m+preallocate+1, ncol = 1)
  sigma2_samples1 <- matrix(nrow = m+preallocate+1, ncol = 1)
  beta_samples1 <-   matrix(nrow = m+preallocate+1, ncol = p)
  eta_samples1 <-    matrix(nrow = m+preallocate+1, ncol = p)
  xi_samples2 <-     matrix(nrow = m+preallocate, ncol = 1)
  sigma2_samples2 <- matrix(nrow = m+preallocate, ncol = 1)
  beta_samples2 <-   matrix(nrow = m+preallocate, ncol = p)
  eta_samples2 <-    matrix(nrow = m+preallocate, ncol = p)
  metric_d <- matrix(nrow = m+preallocate, ncol = 1)
  #
  nrowsamples1 <- m+preallocate+1
  # drawing initial states
  chain1 <- rinit()
  chain2 <- rinit()
  xi_samples1[1,]  <-    xi_1_current <-     chain1$xi
  sigma2_samples1[1,] <- sigma2_1_current <- chain1$sigma2
  beta_samples1[1,]  <-  beta_1_current <-   chain1$beta
  eta_samples1[1,] <-    eta_1_current <-    chain1$eta
  xi_samples2[1,]  <-    xi_2_current <-     chain2$xi
  sigma2_samples2[1,] <- sigma2_2_current <- chain2$sigma2
  beta_samples2[1,]  <-  beta_2_current <-   chain2$beta
  eta_samples2[1,] <-    eta_2_current <-    chain2$eta
  current_nsamples1 <- 1
  #
  # Intializing the first chain to have the L^th marginal at start
  for (l in 1:lag)
  {
    first_chain_update <-
      half_t_kernel(X, X_transpose, y, a0, b0, std_MH,
                    xi_1_current, sigma2_1_current, beta_1_current, eta_1_current,
                    approximate_algo_delta=approximate_algo_delta,
                    nrepeats_eta = nrepeats_eta, verbose = verbose,
                    xi_fixed=xi_fixed, sigma2_fixed=sigma2_fixed, t_dist_df, xi_interval)
    xi_1_current <- first_chain_update$xi_samples
    sigma2_1_current <- first_chain_update$sigma2_samples
    beta_1_current <- first_chain_update$beta_samples
    eta_1_current <- first_chain_update$eta_samples
  }
  #
  current_nsamples1 <- current_nsamples1 + 1
  xi_samples1[current_nsamples1,]  <-    xi_1_current
  sigma2_samples1[current_nsamples1,] <- sigma2_1_current
  beta_samples1[current_nsamples1,]  <-  beta_1_current
  eta_samples1[current_nsamples1,] <-    eta_1_current

  # Setting up coupled chain
  iter <- lag
  meet <- FALSE
  finished <- FALSE
  meetingtime <- Inf
  while (!finished && iter < max_iterations){
    currentime <- Sys.time()
    elapsedtime <-
      as.numeric(lubridate::as.duration(lubridate::ymd_hms(currentime) -
                                          lubridate::ymd_hms(starttime)), "seconds")
    if (elapsedtime > totalduration){
      # time is up, interrupt function
      return(list(finished = FALSE, message = "interrupted because time is up"))
    }
    if (meet){
      # chain_state1 <- single_kernel(chain_state1, ...)
      first_chain_update <-
        half_t_kernel(X, X_transpose, y, a0,b0, std_MH,
                      xi_1_current, sigma2_1_current, beta_1_current, eta_1_current,
                      approximate_algo_delta=approximate_algo_delta,
                      nrepeats_eta = nrepeats_eta, verbose = verbose,
                      xi_fixed=xi_fixed, sigma2_fixed=sigma2_fixed, t_dist_df, xi_interval)
      # chain_state2 <- chain_state1
      xi_1_current <- first_chain_update$xi_samples
      sigma2_1_current <- first_chain_update$sigma2_samples
      beta_1_current <- first_chain_update$beta_samples
      eta_1_current <- first_chain_update$eta_samples
      xi_2_current <- xi_1_current
      sigma2_2_current <- sigma2_1_current
      beta_2_current <- beta_1_current
      eta_2_current <- eta_1_current
      metric_d_current <- 0
    } else {
      # res_coupled_kernel <- coupled_kernel(chain_state1, chain_state2, ...)
      output <-
        coupled_half_t_kernel(X, X_transpose, y, a0, b0, std_MH,
                              xi_1_current, xi_2_current, sigma2_1_current, sigma2_2_current,
                              beta_1_current, beta_2_current, eta_1_current, eta_2_current,
                              approximate_algo_delta = approximate_algo_delta,
                              epsilon_eta = epsilon_eta, epsilon_xi = epsilon_xi, epsilon_sigma2 = epsilon_sigma2,
                              nrepeats_eta = nrepeats_eta, verbose = verbose,
                              xi_fixed=xi_fixed, sigma2_fixed=sigma2_fixed, t_dist_df,
                              two_scale=two_scale, xi_interval)
      xi_1_current <- output$xi_1_samples
      sigma2_1_current <- output$sigma2_1_samples
      beta_1_current <- output$beta_1_samples
      eta_1_current <- output$eta_1_samples
      xi_2_current <- output$xi_2_samples
      sigma2_2_current <- output$sigma2_2_samples
      beta_2_current <- output$beta_2_samples
      eta_2_current <- output$eta_2_samples
      metric_d_current <- output$metric_d
    }
    if ((current_nsamples1+1) > nrowsamples1){
      # print('increase nrow')
      new_rows <- nrowsamples1-1
      nrowsamples1 <- nrowsamples1 + new_rows
      xi_samples1 <- rbind(xi_samples1, matrix(NA, nrow = new_rows, ncol = 1))
      sigma2_samples1 <- rbind(sigma2_samples1, matrix(NA, nrow = new_rows, ncol = 1))
      beta_samples1 <- rbind(beta_samples1, matrix(NA, nrow = new_rows, ncol = p))
      eta_samples1 <- rbind(eta_samples1, matrix(NA, nrow = new_rows, ncol = p))
      xi_samples2 <- rbind(xi_samples2, matrix(NA, nrow = new_rows, ncol = 1))
      sigma2_samples2 <- rbind(sigma2_samples2, matrix(NA, nrow = new_rows, ncol = 1))
      beta_samples2 <- rbind(beta_samples2, matrix(NA, nrow = new_rows, ncol = p))
      eta_samples2 <- rbind(eta_samples2, matrix(NA, nrow = new_rows, ncol = p))
      metric_d <- rbind(metric_d, matrix(NA, nrow = new_rows, ncol = 1))
    }
    xi_samples1[current_nsamples1+1,] <- xi_1_current
    sigma2_samples1[current_nsamples1+1,] <- sigma2_1_current
    beta_samples1[current_nsamples1+1,] <- beta_1_current
    eta_samples1[current_nsamples1+1,] <- eta_1_current
    xi_samples2[current_nsamples1,] <- xi_2_current
    sigma2_samples2[current_nsamples1,] <- sigma2_2_current
    beta_samples2[current_nsamples1,] <- beta_2_current
    eta_samples2[current_nsamples1,] <- eta_2_current
    metric_d[current_nsamples1,] <- metric_d_current

    current_nsamples1 <- current_nsamples1 + 1
    if (all(eta_1_current==eta_2_current) && (xi_1_current==xi_2_current)
        && (sigma2_1_current==sigma2_2_current) && !meet){
      # recording meeting time tau
      meet <- TRUE
      meetingtime <- (iter+1)
    }
    iter <- iter + 1
    # stop after max(m, tau) steps
    if (iter >= (lag + max(meetingtime, mc_chain_size))){
      finished <- TRUE
    }
    if(verbose) print(iter)
  }
  xi_samples1 <-    xi_samples1[1:current_nsamples1,,drop=F]
  xi_samples2 <-    xi_samples2[1:(current_nsamples1-1),,drop=F]
  sigma2_samples1 <- sigma2_samples1[1:current_nsamples1,,drop=F]
  sigma2_samples2 <- sigma2_samples2[1:(current_nsamples1-1),,drop=F]
  beta_samples1 <-  beta_samples1[1:current_nsamples1,,drop=F]
  beta_samples2 <-  beta_samples2[1:(current_nsamples1-1),,drop=F]
  eta_samples1 <-   eta_samples1[1:current_nsamples1,,drop=F]
  eta_samples2 <-   eta_samples2[1:(current_nsamples1-1),,drop=F]
  metric_d_samples <- metric_d[1:(current_nsamples1-1),,drop=F]
  final_output <- list('beta_samples1'=beta_samples1, 'beta_samples2'=beta_samples2,
                       'eta_samples1'=eta_samples1, 'eta_samples2'=eta_samples2,
                       'sigma2_samples1'=sigma2_samples1, 'sigma2_samples2'=sigma2_samples2,
                       'xi_samples1'=xi_samples1, 'xi_samples2'=xi_samples2,
                       'metric_d'=metric_d_samples, meetingtime = meetingtime,
                       finished = TRUE)
  return(final_output)
}



#' meetingtime_half_t
#' @description Function which just returns the meeting time of a coupled chain
#' @param X n by p matrix
#' @param X_transpose Pre-calculated transpose of X
#' @param y length n vector
#' @param a0 positive scalar
#' @param b0 positive scalar
#' @param xi_interval support of prior distribution of xi
#' @param rinit Initial distribution
#' @param approximate_algo_delta approximate MCMC error (non-negative scalar)
#' @param epsilon_eta eta common random numbers/ maximal coupling coupling threshold
#' @param epsilon_xi xi common random numbers/ maximal coupling coupling threshold
#' @param epsilon_sigma2 sigma2 common random numbers/ maximal coupling coupling threshold
#' @param nrepeats_eta number of slice sampling steps
#' @param verbose boolean for printing/ not printing run time
#' @param max_iterations number of maximum iterations
#' @param totalduration maximum duration
#' @param xi_fixed boolean for fixing / not fixing xi
#' @param sigma2_fixed boolean for fixing / not fixing sigma2
#' @param lag Lag
#' @param t_dist_df degree of freedom v>=1 for Half-t(v)
#' @return (beta, eta, sigma2, xi) sampled from the half-t prior blocked Gibbs kernel
#' @export
meetingtime_half_t <-
  function(X, X_transpose, y, a0=1, b0=1, std_MH=0.8, rinit=NULL, approximate_algo_delta=0,
           epsilon_eta=0.5, epsilon_xi=Inf, epsilon_sigma2=Inf, nrepeats_eta=1,
           verbose = FALSE, max_iterations = Inf, totalduration = Inf, lag=1,
           xi_fixed=FALSE, sigma2_fixed=FALSE, t_dist_df, two_scale=TRUE, xi_interval=c(0,Inf)){
  starttime <- Sys.time() # record starting time
  n <- dim(X)[1]
  p <- dim(X)[2]
  #
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
  #
  ## Initialising chains
  chain1 <- rinit()
  chain2 <- rinit()
  xi_1_current <-     chain1$xi
  sigma2_1_current <- chain1$sigma2
  beta_1_current <-   chain1$beta
  eta_1_current <-    chain1$eta
  xi_2_current <-     chain2$xi
  sigma2_2_current <- chain2$sigma2
  beta_2_current <-   chain2$beta
  eta_2_current <-    chain2$eta

  # Intializing the first chain to have the L^th marginal at start
  for (l in 1:lag)
  {
    first_chain_update <-
      half_t_kernel(X, X_transpose, y, a0, b0, std_MH,
                    xi_1_current, sigma2_1_current, beta_1_current, eta_1_current,
                    approximate_algo_delta=approximate_algo_delta,
                    nrepeats_eta = nrepeats_eta, verbose = verbose,
                    xi_fixed=xi_fixed, sigma2_fixed=sigma2_fixed, t_dist_df, xi_interval)
    xi_1_current <- first_chain_update$xi_samples
    sigma2_1_current <- first_chain_update$sigma2_samples
    beta_1_current <- first_chain_update$beta_samples
    eta_1_current <- first_chain_update$eta_samples
  }

  # Number of eta components coupled
  iter <- lag
  meetingtime <- Inf
  while (is.infinite(meetingtime) && iter < max_iterations){
    # Check if time is up already
    currentime <- Sys.time()
    elapsedtime <-
      as.numeric(lubridate::as.duration(lubridate::ymd_hms(currentime) -
                                          lubridate::ymd_hms(starttime)), "seconds")
    if (elapsedtime > totalduration){
      # time is up, interrupt function
      return(list(finished = FALSE, message = "interrupted because time is up"))
    }

    output <-
      coupled_half_t_kernel(X, X_transpose, y, a0, b0, std_MH,
                            xi_1_current, xi_2_current, sigma2_1_current, sigma2_2_current,
                            beta_1_current, beta_2_current, eta_1_current, eta_2_current,
                            approximate_algo_delta = approximate_algo_delta,
                            epsilon_eta = epsilon_eta, epsilon_xi = epsilon_xi, epsilon_sigma2 = epsilon_sigma2,
                            nrepeats_eta = nrepeats_eta, verbose = verbose,
                            xi_fixed=xi_fixed, sigma2_fixed=sigma2_fixed, t_dist_df,
                            two_scale=two_scale, xi_interval)
    # Checking overflow
    # if (!is.null(output$coupled_u)) return(output)
    xi_1_current <- output$xi_1_samples
    sigma2_1_current <- output$sigma2_1_samples
    beta_1_current <- output$beta_1_samples
    eta_1_current <- output$eta_1_samples
    xi_2_current <- output$xi_2_samples
    sigma2_2_current <- output$sigma2_2_samples
    beta_2_current <- output$beta_2_samples
    eta_2_current <- output$eta_2_samples

    if (all(eta_1_current==eta_2_current)){
      # recording meeting time
      meetingtime <- (iter+1)
    }
    iter <- iter + 1
  }
  return(list(finished = TRUE, meetingtime = meetingtime))
}

