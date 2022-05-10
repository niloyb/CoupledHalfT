###############################################################################
#### Functions for the eta updates for blocked Gibbs sampler with half-t priors
###############################################################################
#### Single chain MCMC functions
## eta updates given beta, xi, sigma2
# Incomplete Gamma Function
#' low_inc_gamma
#' @description Lower Incomplete Gamma Function.
#' @param rate positive scalar
#' @param upper_truncation positive scalar
#' @param log boolean to return log of lower incomplete gamma
#' @return log(pgamma(x,a)) or pgamma(x,a) where pgamma(x,a) := 1/Gamma(a) integral_0^x t^(a-1) exp(-t) dt.
#' @export
low_inc_gamma <- function(rate, upper_truncation, log = TRUE){
  return(pgamma(upper_truncation,rate, log.p = log))
}
# Inverse of Lower incomplete Gamma Function
#' low_inc_gamma_inv
#' @description Inverse of Lower incomplete Gamma Function.
#' @param rate positive scalar
#' @param y positive scalar
#' @param log boolean to return log of lower incomplete gamma
#' @return Returns y s.t. low_inc_gamma(a, x)=y
#' @details qgamma is the quantile function of the gamma distribution
#' @seealso \code{\link{low_inc_gamma}} and \code{\link[zipfR]{Rgamma.inv()}}
#' @export
low_inc_gamma_inv <- function(rate, y, log = TRUE){
  # If log=TRUE, the parameter y is taken to be on the natural logarithmic scale.
  if (log){
    if (any(y==0)){stop("Warnings: possible underflow")}
    if (any(y>0)){stop("Log second argument should be negative")}
  } else{
    if (any(y==1)){stop("Warnings: possible underflow")}
    if (any(y>1)){stop("Second argument should less than 1")}
  }
  return(qgamma(y, shape = rate, log.p = log))
}
# Confludent Hypergeomteric Function of the second kind function
#' conf_hyp
#' @description Inverse of Lower incomplete Gamma Function.
#' @param a positive scalar
#' @param b positive scalar
#' @param z positive scalar
#' @return Returns \eqn{U(a,b,z) := \frac{1}{\Gamma(a)} \int_0^\infty x^{a-1}(1+x)^{b-a-1} e^{-zx}dx}
#' @export
conf_hyp <- function(a, b, z){
  integrand_ <- function(x) x^(a-1) * (1+x)^(b-a-1) * exp(-z*x)
  integrate(f = integrand_, lower = 0, upper = +Inf)$val / gamma(a)
}
# Lower Incomplete Confludent Hypergeometric Function of the second kind function
#' low_inc_conf_hyp
#' @description Lower Incomplete Confluent Hypergeometric Function of the second kind
#' @param a positive scalar
#' @param b positive scalar
#' @param z positive scalar
#' @param t positive scalar
#' @return Returns \eqn{U(a,b,z,t) := \frac{1}{\Gamma(a)} \int_0^t x^{a-1}(1+x)^{b-a-1} e^{-zx}dx}
#' @export
low_inc_conf_hyp <- function(a, b, z, t, regularize=TRUE){
  integrand_ <- function(x) x^(a-1) * (1+x)^(b-a-1) * exp(-z*x)
  if (regularize){
    integrate(f = integrand_, lower = 0, upper = t)$val /
      integrate(f = integrand_, lower = 0, upper = +Inf)$val
  } else{
    integrate(f = integrand_, lower = 0, upper = t)$val / gamma(a)
  }
}
# Perfect sampling from univariate p(x) \prop x^(poly_exponent-1)*exp(-rate*x) on [0,trunc_upper]
#' r_trunc_poly_exp_crn
#' @description Perfect sampling from \eqn{p(x) \prop x^(poly_exponent-1)*exp(-rate*x) on [0,trunc_upper]}
#' @param poly_exponent positive scalar
#' @param rate positive vector (function is vectorized)
#' @param trunc_upper positive scalar
#' @param unif positive scalar
#' @return Returns vectorised samples from \eqn{p(x) \prop x^(poly_exponent-1)*exp(-rate*x) on [0,trunc_upper]}
#' @export
r_trunc_poly_exp_crn <- function(poly_exponent, rate, trunc_upper, unif){
  u <- unif*(low_inc_gamma(poly_exponent, rate*trunc_upper, log = FALSE))
  return(low_inc_gamma_inv(poly_exponent, u, log=FALSE)/rate)
}
#' @export
r_trunc_poly_exp <- function(poly_exponent, rate, trunc_upper){
  unif <- runif(length(rate))
  return(r_trunc_poly_exp_crn(poly_exponent, rate, trunc_upper, unif))
}
# Log Pdf of univariate p(x) \prop x^(poly_exponent-1)*exp(-rate*x) on [0,trunc_upper]
#' d_trunc_poly_exp
#' @description Log Pdf of \eqn{ p(x) \prop x^(poly_exponent-1)*exp(-rate*x) on [0,trunc_upper] }
#' @param poly_exponent positive scalar
#' @param rate positive vector
#' @param trunc_upper positive scalar
#' @param eta positive vector
#' @return Returns Log Pdf of \eqn{p(x) \prop x^(poly_exponent-1)*exp(-rate*x) on [0,trunc_upper]}
#' @export
d_trunc_poly_exp <- function(poly_exponent, rate, trunc_upper, eta){
  if(eta>trunc_upper) return(-Inf)
  log_d <- (poly_exponent-1)*log(eta)-rate*eta
  log_norm_constant <- -poly_exponent*log(rate)+lgamma(poly_exponent)+
    low_inc_gamma(poly_exponent, rate*trunc_upper, log = TRUE)
  return(log_d-log_norm_constant)
}
# Slice sampling from univariate p(x) \prop x^((v-1)/2)/(1+vx)^((v+1)/2)*exp(-mx) on [L,Inf]
#' eta_update_half_t
#' @description Slice sampler targeting \eqn{p(x) \prop x^((v-1)/2)/(1+vx)^((v+1)/2)*exp(-mx) on [0,Inf], }
#' where v is t_dist_df and m is rate.
#' @param xi positive scalar
#' @param sigma2 positive scalar
#' @param beta positive vector
#' @param eta positive vector
#' @param t_dist_df half-t degree of freedom
#' @param nrepeats number of slice sampling steps
#' @return Slice sampler targeting \eqn{p(x) \prop x^((v-1)/2)/(1+vx)^((v+1)/2)*exp(-mx) on [0,Inf], }
#' where v is t_dist_df and m is rate.
#' @export
eta_update_half_t <- function(xi, sigma2, beta, eta, t_dist_df, nrepeats=1)
{
  rate <- (beta^2)*(xi)/(2*sigma2)
  p <- length(eta)
  if (is.infinite(nrepeats)){
    stop('Perfect sampling not implemented for general t distribution shrinkage priors')
  } else {
    for (irepeat in 1:nrepeats){
      u <- runif(p)/(1 + t_dist_df*eta)^((1+t_dist_df)/2)
      eta <- r_trunc_poly_exp((1+t_dist_df)/2, rate,
                              (u^(-2/(1+t_dist_df))-1)/t_dist_df)
    }
  }
  return(eta)
}


###############################################################################
#### Coupled MCMC functions
## CRN Coupling of Eta Update
#' eta_update_half_t_crn_couple
#' @description Common random numbers coupling of Eta Update
#' @param xi_1,xi_2 xi values from the pair of chains
#' @param Beta_1,Beta_2 beta values (each vector of length p) from the pair of chains
#' @param eta_1,eta_2 eta values (each vector of length p) from the pair of chains
#' @param sigma2_1,sigma2_2 sigma values from the pair of chains
#' @param t_dist_df half-t degree of freedom
#' @return Returns (eta_1, eta_2) under common random numbers coupling
#' @export
eta_update_half_t_crn_couple <- function(xi_1, Beta_1, eta_1, sigma2_1,
                                         xi_2, Beta_2, eta_2, sigma2_2,
                                         t_dist_df){
  p <- length(eta_1)
  rate_1 <- (Beta_1^2)*(xi_1)/(2*sigma2_1)
  rate_2 <- (Beta_2^2)*(xi_2)/(2*sigma2_2)
  unif_crn_1 <- runif(p)
  u_1 <- unif_crn_1/(1 + t_dist_df*eta_1)^((1+t_dist_df)/2)
  u_2 <- unif_crn_1/(1 + t_dist_df*eta_2)^((1+t_dist_df)/2)
  unif_crn_2 <- runif(p)
  eta_1 <- r_trunc_poly_exp_crn((1+t_dist_df)/2, rate_1, (u_1^(-2/(1+t_dist_df))-1)/t_dist_df, unif_crn_2)
  eta_2 <- r_trunc_poly_exp_crn((1+t_dist_df)/2, rate_2, (u_2^(-2/(1+t_dist_df))-1)/t_dist_df, unif_crn_2)
  return(cbind(eta_1, eta_2))
}

## Maximal-coupling based Eta update
# Max couplings of univariates of form p(x) \prop x^(poly_exponent-1)*exp(-rate*x) on [0,trunc_upper]
trunc_poly_exp_max_couple <- function(poly_exponent, rate_1, rate_2,
                                      trunc_upper_1, trunc_upper_2){
  rp <- function() r_trunc_poly_exp(poly_exponent, rate_1, trunc_upper_1)
  rq <- function() r_trunc_poly_exp(poly_exponent, rate_2, trunc_upper_2)
  dp <- function(x){return(d_trunc_poly_exp(poly_exponent, rate_1, trunc_upper_1, x))}
  dq <- function(x){return(d_trunc_poly_exp(poly_exponent, rate_2, trunc_upper_2, x))}
  f <- get_max_coupling(rp, dp, rq, dq)
  return(f())
}
# Maximal-coupling based Eta update
eta_update_half_t_max_couple <- function(xi_1, Beta_1, eta_1, sigma2_1,
                                         xi_2, Beta_2, eta_2, sigma2_2,
                                         t_dist_df){
  p <- length(eta_1)
  rate_1 <- (Beta_1^2)*(xi_1)/(2*sigma2_1)
  rate_2 <- (Beta_2^2)*(xi_2)/(2*sigma2_2)
  etas_ <- matrix(0, nrow = p, ncol = 2)
  crn_unif <- runif(p)
  coupled_u <- cbind(crn_unif/(1 + t_dist_df*eta_1)^((1+t_dist_df)/2),
                     crn_unif/(1 + t_dist_df*eta_2)^((1+t_dist_df)/2))
  truncs_ <- ((coupled_u)^(-2/(1+t_dist_df))-1)/t_dist_df
  for (j in 1:p){
    etas_[j,] <-
      trunc_poly_exp_max_couple((1+t_dist_df)/2, rate_1[j], rate_2[j],
                                truncs_[j,1], truncs_[j,2])
  }
  eta_1 <- etas_[,1]
  eta_2 <- etas_[,2]
  return(cbind(eta_1, eta_2))
}

## Distance metric
# Total variation between two univariates p(x) \prop x^(r-1)*exp(-mx) on [0,T]
# with r poly_exponent; m rate_1, rate_2; T trunc_upper_1, trunc_upper_2.
# Note: function is vectorized.
trunc_poly_exp_tv <- function(poly_exponent, rate_1, rate_2,
                              trunc_upper_1, trunc_upper_2){
  # Probabilities of coupling each component
  coupling_probs <- rep(NA, length(rate_1))
  
  rate_index1 <- rate_1==rate_2
  rate_index2 <- rate_1<rate_2
  rate_index3 <- rate_1>rate_2
  
  K <- (-poly_exponent*log(rate_1)+low_inc_gamma(poly_exponent, rate_1*trunc_upper_1, log = TRUE))-
    (-poly_exponent*log(rate_2)+low_inc_gamma(poly_exponent, rate_2*trunc_upper_2, log = TRUE))
  K[rate_2!=rate_1] <- (K/(rate_2-rate_1))[rate_2!=rate_1]
  trunc_upper_min <- exp(((log(trunc_upper_1) + log(trunc_upper_2)) -
                            abs(log(trunc_upper_1)-log(trunc_upper_2)))/2)
  trunc_upper_max <- exp(((log(trunc_upper_1) + log(trunc_upper_2)) +
                            abs(log(trunc_upper_1)-log(trunc_upper_2)))/2)
  
  K_index1 <- K<=0
  K_index2 <- K>=trunc_upper_min
  K_index3 <- (0<K)&(K<trunc_upper_min)
  
  # Calculating coupling probabilities case by case for K_index and rate_index
  coupling_probs[rate_index1] <-
    exp(low_inc_gamma(poly_exponent,(rate_1*trunc_upper_min)[rate_2==rate_1], log = TRUE)-
          low_inc_gamma(poly_exponent,(rate_1*trunc_upper_max)[rate_2==rate_1], log = TRUE))
  coupling_probs[K_index1 & rate_index2] <-
    exp(low_inc_gamma(poly_exponent,(rate_2*trunc_upper_min)[K_index1 & rate_index2], log = TRUE)-
          low_inc_gamma(poly_exponent,(rate_2*trunc_upper_2)[K_index1 & rate_index2], log = TRUE))
  coupling_probs[K_index2 & rate_index2] <-
    exp(low_inc_gamma(poly_exponent,(rate_1*trunc_upper_min)[K_index2 & rate_index2], log = TRUE)-
          low_inc_gamma(poly_exponent,(rate_1*trunc_upper_1)[K_index2 & rate_index2], log = TRUE))
  coupling_probs[K_index3 & rate_index2] <-
    exp(low_inc_gamma(poly_exponent,(rate_1*K)[K_index3 & rate_index2], log = TRUE)-
          low_inc_gamma(poly_exponent,(rate_1*trunc_upper_1)[K_index3 & rate_index2], log = TRUE)) +
    exp(log(low_inc_gamma(poly_exponent,(rate_2*trunc_upper_min)[K_index3 & rate_index2], log = FALSE)-
       low_inc_gamma(poly_exponent,(rate_2*K)[K_index3 & rate_index2], log = FALSE))-
         low_inc_gamma(poly_exponent,(rate_2*trunc_upper_2)[K_index3 & rate_index2], log = TRUE))
  coupling_probs[K_index1 & rate_index3] <-
    exp(low_inc_gamma(poly_exponent,(rate_1*trunc_upper_min)[K_index1 & rate_index3], log = TRUE)-
          low_inc_gamma(poly_exponent,(rate_1*trunc_upper_1)[K_index1 & rate_index3], log = TRUE))
  coupling_probs[K_index2 & rate_index3] <-
    exp(low_inc_gamma(poly_exponent,(rate_2*trunc_upper_min)[K_index2 & rate_index3], log = TRUE)-
          low_inc_gamma(poly_exponent,(rate_2*trunc_upper_2)[K_index2 & rate_index3], log = TRUE))
  coupling_probs[K_index3 & rate_index3] <-
    exp(low_inc_gamma(poly_exponent,(rate_2*K)[K_index3 & rate_index3], log = TRUE)-
          low_inc_gamma(poly_exponent,(rate_2*trunc_upper_2)[K_index3 & rate_index3], log = TRUE)) +
    exp(log(low_inc_gamma(poly_exponent,(rate_1*trunc_upper_min)[K_index3 & rate_index3], log = FALSE)-
              low_inc_gamma(poly_exponent,(rate_1*K)[K_index3 & rate_index3], log = FALSE))-
          low_inc_gamma(poly_exponent,(rate_1*trunc_upper_1)[K_index3 & rate_index3], log = TRUE))
  
  # Correcting rounding errors
  coupling_probs[coupling_probs<0] <- 0
  coupling_probs[coupling_probs>1] <- 1
  
  if(is.null(nrow(coupling_probs))){return(1-exp(sum(log(coupling_probs))))}
  return(1-exp(rowSums(log(coupling_probs))))
}
# TV UB from componentwise coupling
half_t_max_couple_prob <- function(xi_1, Beta_1, eta_1, sigma2_1,
                                   xi_2, Beta_2, eta_2, sigma2_2,
                                   t_dist_df, iterations=1){
  p <- length(eta_1)
  rate_1 <- (Beta_1^2)*(xi_1)/(2*sigma2_1)
  rate_2 <- (Beta_2^2)*(xi_2)/(2*sigma2_2)
  
  tv_ub <- rep(0, iterations)
  for (i in 1:iterations){
    crn_unif <- runif(p)
    coupled_u <- cbind(crn_unif/(1 + t_dist_df*eta_1)^((1+t_dist_df)/2),
                       crn_unif/(1 + t_dist_df*eta_2)^((1+t_dist_df)/2))
    truncs_ <- ((coupled_u)^(-2/(1+t_dist_df))-1)/t_dist_df
    tv_ub[i] <-
      trunc_poly_exp_tv((1+t_dist_df)/2, rate_1, rate_2, truncs_[,1], truncs_[,2])
    
    if(is.na(tv_ub[i])) # Checking underflow
    {
      print(c("Warning: possible overflow", i))
      return(list(xi_1=xi_1, Beta_1=Beta_1, eta_1=eta_1, sigma2_1=sigma2_1,
                  xi_2=xi_2, Beta_2=Beta_2, eta_2=eta_2, sigma2_2=sigma2_2,
                  p=p, eta_lower_bound=eta_lower_bound, t_dist_df=t_dist_df,
                  iterations=iterations, coupled_u=coupled_u))
    }
  }
  return(mean(tv_ub))
}

##### New coupling scheme #####
## CRN Coupling of Eta Update
#' eta_update_half_t_max_couple_till_you_miss
#' @description Max Coupling componentwise until no meeting. Then CRN componentwise.
#' @param xi_1,xi_2 xi values from the pair of chains
#' @param Beta_1,Beta_2 beta values (each vector of length p) from the pair of chains
#' @param eta_1,eta_2 eta values (each vector of length p) from the pair of chains
#' @param sigma2_1,sigma2_2 sigma values from the pair of chains
#' @param t_dist_df half-t degree of freedom
#' @return Returns (eta_1, eta_2) under common random numbers coupling
#' @export
eta_update_half_t_max_couple_till_you_miss <- 
  function(xi_1, Beta_1, eta_1, sigma2_1, xi_2, Beta_2, eta_2, sigma2_2, t_dist_df){
    p <- length(eta_1)
    eta_sample <- cbind(rep(NA,p), rep(NA,p))
    ordered_components <- sample(c(1:p)) # c(1:p)
    max_couple <- TRUE
    for(i in 1:p){
      eta_sample[ordered_components[i],] <- 
        eta_update_half_t_max_couple(xi_1, Beta_1[ordered_components[i]], eta_1[ordered_components[i]], sigma2_1,
                                     xi_2, Beta_2[ordered_components[i]], eta_2[ordered_components[i]], sigma2_2, 
                                     t_dist_df)
      if(eta_sample[ordered_components[i],1]!=eta_sample[ordered_components[i],2]){break}
    }
    if(i<p){
      eta_sample[ordered_components[(i+1):p],] <- 
        eta_update_half_t_crn_couple(xi_1, Beta_1[ordered_components[(i+1):p]], eta_1[ordered_components[(i+1):p]], sigma2_1,
                                     xi_2, Beta_2[ordered_components[(i+1):p]], eta_2[ordered_components[(i+1):p]], sigma2_2, 
                                     t_dist_df)
    }
    return(eta_sample)
  }





