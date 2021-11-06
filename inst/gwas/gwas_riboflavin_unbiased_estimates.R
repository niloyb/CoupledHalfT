# Libraries
rm(list = ls())
set.seed(1)

library(doParallel)
library(doRNG)
registerDoParallel(cores = detectCores())

library(CoupledHalfT)
library(ggplot2)
library(latex2exp)
library(reshape2)
library(dplyr)
library(ggridges)

################################################################################
sample_unbiasedestimator <- 
  function(single_kernel, coupled_kernel, rinit, h = function(x) x, 
           k = 0, m = 1, lag = 1, max_iterations = Inf){
  if (k > m){
    print("error: k has to be less than m")
    return(NULL)
  }
  if (lag > m){
    print("error: lag has to be less than m")
    return(NULL)
  }
  starttime <- Sys.time()
  state1 <- rinit(); state2 <- rinit()
  # mcmcestimator computes the sum of h(X_t) for t=k,...,m
  mcmcestimator <- h(state1$chain_state)
  dimh <- length(mcmcestimator)
  if (k > 0){
    mcmcestimator <- rep(0, dimh)
  }
  # correction computes the sum of min(m-k+1, ceiling((t - k)/lag)) * (h(X_{t}) - h(Y_{t-lag})) for t=k+lag,..., tau - 1
  correction <- rep(0, dimh)
  time <- 0
  for (t in 1:lag){
    time <- time + 1
    state1 <- single_kernel(state1)
    if (time >= k){
      mcmcestimator <- mcmcestimator + h(state1$chain_state)
    }
  }
  if (time >= k + lag){
    correction <- correction + min(m-k+1, ceiling((time - k)/lag)) * (h(state1$chain_state) - h(state2$chain_state))
  }
  meetingtime <- Inf
  # time here is equal to lag; at this point we have X_lag,Y_0 and we are going to generate successively X_{t},Y_{t-lag} where time t is >= lag+1
  while ((time < max(meetingtime, m)) && (time < max_iterations)){
    time <- time + 1 # time is lag+1,lag+2,...
    if (is.finite(meetingtime)){
      state1 <- single_kernel(state1)
      state2 <- state1
      if (k <= time && time <= m){
        mcmcestimator <- mcmcestimator + h(state1$chain_state)
      }
    } else {
      res_coupled_kernel <- coupled_kernel(state1, state2)
      state1 <- res_coupled_kernel$state1
      state2 <- res_coupled_kernel$state2
      if (res_coupled_kernel$identical){
        meetingtime <- time
      }
      if (k <= time && time <= m){
        mcmcestimator <- mcmcestimator + h(state1$chain_state)
      }
      if (time >= k + lag){
        correction <- correction + min(m-k+1, ceiling((time - k)/lag)) * (h(state1$chain_state) - h(state2$chain_state))
      }
    }
  }
  uestimator <- mcmcestimator + correction
  cost <- lag + 2*(meetingtime - lag) + max(0, time - meetingtime)
  currentime <- Sys.time()
  elapsedtime <- as.numeric(lubridate::as.duration(lubridate::ymd_hms(currentime) - lubridate::ymd_hms(starttime)), "seconds")
  return(list(mcmcestimator = mcmcestimator / (m - k + 1), correction = correction / (m - k + 1), uestimator = uestimator / (m - k + 1),
              meetingtime = meetingtime, iteration = time, elapsedtime = elapsedtime, cost = cost))
  }


################################################################################
# Loading riboflavin dataset
data(riboflavin)

y <- as.vector(riboflavin$y)
X <- as.matrix(riboflavin$x)
colnames(X) <- NULL
rownames(X) <- NULL
# X <- matrix(scale(X), n, p)
X_transpose <- t(X)
n <- length(y)
p <-dim(X)[2]

# ################################################################################
# # Generate Synthetic data
# n <- 100
# p <- 1000
# s <- 10
# true_beta <- matrix(0,p,1)
# true_beta[1:s] = 2^(-(seq(s)/4-9/4))
# X <- matrix(rnorm(n*p), nrow = n, ncol = p)
# X_transpose <- t(X)
# #Error terms
# error_std <- 0.5
# error_terms = error_std*rnorm(n, mean = 0, sd = 1)
# y = X%*%true_beta + error_terms

t_dist_df <- 2
xi_interval <- c(0,Inf)
verbose=TRUE

################################################################################
single_kernel <- function(x){
  xi_current <- (x$chain_state)$xi_samples
  sigma2_current <- (x$chain_state)$sigma2_samples
  beta_current <- (x$chain_state)$beta_samples
  eta_current <- (x$chain_state)$eta_samples
  
  output <- half_t_kernel(X, X_transpose, y, a0=1, b0=1, std_MH=0.8,
                          xi_current, sigma2_current,
                          beta_current, eta_current, 
                          approximate_algo_delta=0,
                          t_dist_df=t_dist_df, xi_interval=xi_interval, verbose=verbose)
  return(list('chain_state'=output))
}


coupled_kernel <- function(x1, x2){
  xi_1_current <- (x1$chain_state)$xi_samples
  sigma2_1_current <- (x1$chain_state)$sigma2_samples
  beta_1_current <- (x1$chain_state)$beta_samples
  eta_1_current <- (x1$chain_state)$eta_samples
  
  xi_2_current <- (x2$chain_state)$xi_samples
  sigma2_2_current <- (x2$chain_state)$sigma2_samples
  beta_2_current <- (x2$chain_state)$beta_samples
  eta_2_current <- (x2$chain_state)$eta_samples
  
  output <- coupled_half_t_kernel(X, X_transpose, y, a0=1, b0=1, std_MH=0.8,
                                  xi_1_current, xi_2_current, sigma2_1_current, sigma2_2_current,
                                  beta_1_current, beta_2_current, eta_1_current, eta_2_current,
                                  t_dist_df=t_dist_df, xi_interval=xi_interval, verbose=verbose)
  
  state1 <- list('chain_state'=
                   list('beta_samples'=output$beta_1_samples,
                        'eta_samples'=output$eta_1_samples,
                        'sigma2_samples'=output$sigma2_1_samples,
                        'xi_samples'=output$xi_1_samples))
  state2 <- list('chain_state'=
                   list('beta_samples'=output$beta_2_samples,
                        'eta_samples'=output$eta_2_samples,
                        'sigma2_samples'=output$sigma2_2_samples,
                        'xi_samples'=output$xi_2_samples))
  
  if((max(abs(output$beta_1_samples-output$beta_2_samples))==0)&
     (max(abs(output$eta_1_samples-output$eta_2_samples))==0)&
     (max(abs(output$sigma2_1_samples-output$sigma2_2_samples))==0)&
     (max(abs(output$xi_1_samples-output$xi_2_samples))==0)){
    identical <- TRUE
  } else {
    identical <- FALSE
  }
  return(list('state1'=state1, 'state2'=state2, 'identical'=identical))
}

# Initializing from the prior
rinit <- function(a0=1, b0=1){
  # xi <- (1/rt(1, df=1))^2
  xi <- tan(runif(1)*(atan(xi_interval[2]^0.5)-atan(xi_interval[1]^0.5))+
              atan(xi_interval[1]^0.5))^2
  sigma2 <- 1/rgamma(1, shape = a0/2, rate = b0/2)
  eta <- (1/rt(p, df=t_dist_df))^2
  beta <- rnorm(p)*sqrt(sigma2/(xi*eta))
  return(list('chain_state'=
                list('xi_samples' = xi, 'sigma2_samples' = sigma2, 
                     'beta_samples' = beta, 'eta_samples' = eta)))
}


################################################################################
# Unbiased estimation
iterations <- 10
k <- 1000
m <- 2000
lag <- 750

# iterations <- 100
# k <- 250
# m <- 2000
# lag <- 200

unbiased_beta <- foreach(i = c(1:iterations), .combine = '+')%dopar%{
  output <- sample_unbiasedestimator(single_kernel, coupled_kernel, rinit, 
                                     h = function(x) x$beta_samples, 
                                     k = k, m = m, lag = lag)
  return(output$uestimator/iterations)
}

large_components <- order(abs(unbiased_beta), decreasing = TRUE)[1:5]

# MCMC estimate
mcmc_beta <- foreach(i = c(1:iterations), .combine = rbind)%dopar%{
  output <- half_t_mcmc(chain_length=m, burnin=k, X, X_transpose, y, t_dist_df=t_dist_df)
  return(output$beta_samples[,large_components, drop=FALSE])
}

# unbiased_riboflavin_data <- list('mcmc_beta'=mcmc_beta,'unbiased_beta'=unbiased_beta)
# save(unbiased_riboflavin_data, file='../BackUpFiles/CoupledHalfT_old/examples/big_data_examples/unbiased_riboflavin.RData')

# load('../BackUpFiles/CoupledHalfT_old/examples/big_data_examples/unbiased_riboflavin.RData')
# unbiased_beta <- unbiased_riboflavin_data$unbiased_beta
# mcmc_beta <- unbiased_riboflavin_data$mcmc_beta

# Comparison
# ymax <- max(mcmc_beta[,large_components],unbiased_beta[large_components])
# ymin <- min(mcmc_beta[,large_components],unbiased_beta[large_components])
# # matplot(mcmc_beta[,large_components], type='l', ylim = c(min(unbiased_beta[large_components]), max(unbiased_beta[large_components])))
# matplot(mcmc_beta[,large_components], type='l', ylim=c(ymin, ymax))
# abline(h=unbiased_beta[large_components])

# Traceplots with unbiased mean estimate
component <- 1
traceplot <- data.frame('iteration'=c((k+1):m), 'value'=mcmc_beta[,component])
plot <- ggplot(traceplot, aes(x = value)) + 
  geom_density(fill = "grey", colour = "grey", alpha = 0.5) + theme_classic(base_size = 14) +
  geom_vline(xintercept=unbiased_beta[large_components[component]], linetype="solid", color = "black", size=1) +
  xlab(TeX(paste0("$\\beta_{", large_components[component], "}$ posterior")))
plot
# ggsave(filename = "/Users/niloybiswas/Dropbox/horseshoe_coupling/Drafts/images/dataset_examples/unbiased_riboflavin1.pdf",
#        plot = plot, width = 3, height = 3) # base_size = 18

component <- 2
traceplot <- data.frame('iteration'=c((k+1):m), 'value'=mcmc_beta[,component])
plot <- ggplot(traceplot, aes(x = value)) + 
  geom_density(fill = "grey", colour = "grey", alpha = 0.5) + theme_classic(base_size = 14) +
  geom_vline(xintercept=unbiased_beta[large_components[component]], linetype="solid", color = "black", size=1) +
  xlab(TeX(paste0("$\\beta_{", large_components[component], "}$ posterior")))
plot
# ggsave(filename = "/Users/niloybiswas/Dropbox/horseshoe_coupling/Drafts/images/dataset_examples/unbiased_riboflavin2.pdf",
#        plot = plot, width = 3, height = 3) # base_size = 18

component <- 3
traceplot <- data.frame('iteration'=c((k+1):m), 'value'=mcmc_beta[,component])
plot <- ggplot(traceplot, aes(x = value)) + 
  geom_density(fill = "grey", colour = "grey", alpha = 0.5) + theme_classic(base_size = 14) +
  geom_vline(xintercept=unbiased_beta[large_components[component]], linetype="solid", color = "black", size=1) +
  xlab(TeX(paste0("$\\beta_{", large_components[component], "}$ posterior")))
plot
# ggsave(filename = "/Users/niloybiswas/Dropbox/horseshoe_coupling/Drafts/images/dataset_examples/unbiased_riboflavin3.pdf",
#        plot = plot, width = 3, height = 3) # base_size = 18











