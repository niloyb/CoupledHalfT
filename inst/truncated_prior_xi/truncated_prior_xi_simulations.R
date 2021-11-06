# Truncated prior for xi simulations
# Libraries
rm(list = ls())
set.seed(1)
library(CoupledHalfT)
library(dplyr)
library(ggplot2)
library(latex2exp)
library(gridExtra)
library(ggridges)
library(doParallel)
registerDoParallel(cores = detectCores())


## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
## One and two scale coupling simulations ## ## ## ## ## ## ## ## ## ## ## ##

# Simulation setup
iterations <- 100

# Horseshoe prior
t_dist_df <- 1

# Maximal coupling at every step
epsilon_eta <- Inf

# Fixed n varying p simulations
n <- 100
# meetingtimes_df <- data.frame()
meetingtimes_df <- 
  foreach(epsilon_eta = c(Inf,0.5), .combine = rbind) %:% 
  foreach(p = seq(50,150,20), .combine = rbind) %:% 
  foreach(i = c(1:iterations), .combine = rbind) %dopar% {
    ## Generate data
    s <- 10
    true_beta <- matrix(0,p,1)
    true_beta[1:s] = 2^(-(seq(s)/4-9/4))
    X <- matrix(rnorm(n*p), nrow = n, ncol = p)
    X_transpose <- t(X)
    #Error terms
    error_std <- 0.5
    error_terms = error_std*rnorm(n, mean = 0, sd = 1)
    y = X%*%true_beta + error_terms
    
    xi_interval <- c(1,p^2)
    
    mc_chain_size <- 1
    burnin <- 0
    max_iterations <- 1e+05
    meetingtime <- 
      meetingtime_half_t(X, X_transpose, y, a0 = 1, b0 = 1, std_MH = 0.8, 
                         epsilon_eta = epsilon_eta, 
                         max_iterations = max_iterations, t_dist_df=t_dist_df,
                         xi_interval=xi_interval)
    
    print(c(n,p,i,meetingtime$meetingtime)) # Will not print with %dopar%
    data.frame(n = n, p = p, meetingtime = meetingtime, epsilon_eta = epsilon_eta, 
               t_dist_df = t_dist_df, sparsity = s, error_std = error_std)
  }

# Save data
# save(meetingtimes_df, file = "examples/truncated_prior_xi/one_two_scale_coupling_meetings.RData")


## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
## Meeting time simulations ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

## ## ## ## ## 
## n<=p plot ## 
## ## ## ## ##

# Two scaling coupling simulation setup
iterations <- 100

epsilon_eta <- 0.5

# Varying p, nu simulations
n <- 100
tdist_meetingtimes_df <-
  foreach(t_dist_df = seq(2,1,-0.2), .combine = rbind) %:% 
  foreach(p = seq(200,1000,200), .combine = rbind) %:% 
  foreach(i = 1:iterations, .combine = rbind) %dopar% {
    ## Generate data
    s <- 20
    true_beta <- matrix(0,p,1)
    true_beta[1:s] = 2^(-(seq(s)/4-9/4))
    X <- matrix(rnorm(n*p), nrow = n, ncol = p)
    X_transpose <- t(X)
    #Error terms
    error_std <- 2
    error_terms = error_std*rnorm(n, mean = 0, sd = 1)
    y = X%*%true_beta + error_terms
    
    xi_interval <- c(1,p^2)
    max_iterations <- Inf
    meetingtime <- 
      meetingtime_half_t(X, X_transpose, y, a0 = 1, b0 = 1, std_MH = 0.8, 
                         epsilon_eta = epsilon_eta, 
                         max_iterations = max_iterations, t_dist_df=t_dist_df, 
                         xi_interval=xi_interval)
    print(c(t_dist_df, p, i, meetingtime$meetingtime))  # Will not print with %dopar%
    return(data.frame(n = n, p =p, sparsity = s, error_std = error_std, 
                      meetingtime = meetingtime, t_dist_df = t_dist_df))
  }

# save(tdist_meetingtimes_df, file = "examples/truncated_prior_xi/degree_of_freedom_meeting_times.RData") # meetingtimes_df


## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
## Statistical Performance Simulations ## ## ## ## ## ## ## ## ## ## ## ## ##
# Statistical Performance Functions
# MSE of Beta by Estimation of Time-Averages
BetaMSE <- function(true_beta, mc_beta_samples){
  return(mean((true_beta-colMeans(mc_beta_samples))^2))
}
# MSE of Beta by Estimation of Time-Averages
XBetaMSE <- function(true_beta, mc_beta_samples, X){
  return(mean((X%*%(true_beta-colMeans(mc_beta_samples)))^2))
}
ComponentwiseBetaMSE <- function(true_beta, mc_beta_samples){
  return((true_beta-colMeans(mc_beta_samples))^2)
}
# Covergae of credible intervals
BetaCoverageSymmetric <- function(true_beta, mc_beta_samples, qtl){
  componentwise_quantiles <- t(apply(mc_beta_samples, 2 , quantile , probs = c((1-qtl)/2,(1+qtl)/2)))
  coverage <- (componentwise_quantiles[,2]-true_beta)*(true_beta-componentwise_quantiles[,1])>0
  return(coverage)
}
# Covergae of credible intervals
BetaCoverageHD <- function(true_beta, mc_beta_samples, qtl){
  p <- length(true_beta)
  s <- max(which(true_beta!=0))
  coverage <- rep(NA,p)
  # Covering bimodalities w allowSplit
  for(i in 1:s){
    hdint <- HDInterval::hdi(density(mc_beta_samples[,i]), credMass = qtl, allowSplit=TRUE)
    coverage[i] <- any((hdint[,2]-true_beta[i])*(true_beta[i]-hdint[,1]) > 0)
  }
  componentwise_quantiles <- t(HDInterval::hdi(mc_beta_samples[,c((s+1):p)], credMass = qtl))
  coverage[(s+1):p] <- (componentwise_quantiles[,2]-true_beta[c((s+1):p)])*(true_beta[c((s+1):p)]-componentwise_quantiles[,1])>0
  return(coverage)
}

# Run simulations

## Generate data
n <- 100
p <- 500
s <- 20
true_beta <- matrix(0,p,1)
true_beta[1:s] = 2^(-(seq(s)/4-9/4))

xi_interval <- c(1,p^2)

# X <- matrix(rnorm(n*p, mean = 0, sd = 1), nrow = n, ncol = p)
# #Error terms
# error_std <- 2
# error_terms = error_std*rnorm(n, mean = 0, sd = 1)
# y = X%*%true_beta + error_terms
# X_transpose <- t(X)

stat_perform_df <- data.frame()
stat_perform_componentwise_df <- data.frame()
traceplots_df <- data.frame()
iterations <- 100
for (i in 1:iterations){
  # Generate data
  X <- matrix(rnorm(n*p, mean = 0, sd = 1), nrow = n, ncol = p)
  #Error terms
  error_std <- 2
  error_terms = error_std*rnorm(n, mean = 0, sd = 1)
  y = X%*%true_beta + error_terms
  X_transpose <- t(X)
  
  # Lasso Estimates
  cv.lasso <- glmnet::cv.glmnet(X, y, alpha = 1, family = 'gaussian', intercept = FALSE)
  lasso.estimates <- coef(cv.lasso, cv.lasso$lambda.min)[2:(p+1)]
  lasso.estimates <- data.frame(matrix(lasso.estimates, nrow = 1))
  colnames(lasso.estimates) <- paste("X", 1:p, sep = '')
  
  chain_df <- lasso.estimates[,c(1:(s+10))] %>%
    tidyr::gather(component, value, X1:X30) %>%
    dplyr::mutate(n = n, p =p, t_dist_df = 'lasso')
  traceplots_df <- rbind(traceplots_df, chain_df)
  stat_perform_df <-
    rbind(stat_perform_df,
          data.frame(n = n, p =p, t_dist_df = 'lasso',
                     BetaMSE=BetaMSE(true_beta, lasso.estimates),
                     XBetaMSE=XBetaMSE(true_beta, lasso.estimates, X),
                     BetaHDCoverage0.95=NA))
  stat_perform_componentwise_df <-
    rbind(stat_perform_componentwise_df,
          data.frame(n=n, p=p, t_dist_df = 'lasso', component = c(1:p),
                     BetaHDCoverage0.95 = NA,
                     ComponentwiseBetaMSE =
                       ComponentwiseBetaMSE(true_beta, lasso.estimates)))
  for (t_dist_df in c(seq(1,2,0.2),4,6,8)){
    if (t_dist_df==1){burnin <- 600} else {burnin <- 300}
    mc_chain_size <- burnin + 1000
    
    chain <- half_t_mcmc(mc_chain_size, burnin=0, X, X_transpose,y, 
                         a0=1, b0=1, std_MH=0.8, t_dist_df=t_dist_df, 
                         xi_interval=xi_interval)
    mc_beta_samples <- chain$beta_samples[c(burnin:mc_chain_size),]
    
    chain_df <- data.frame(mc_beta_samples[,c(1:(s+10))]) %>%
      tidyr::gather(component, value, X1:X30) %>%
      dplyr::mutate(n = n, p =p, t_dist_df = t_dist_df)
    traceplots_df <- rbind(traceplots_df, chain_df)
    
    stat_perform_df <-
      rbind(stat_perform_df,
            data.frame(n = n, p =p, t_dist_df = t_dist_df,
                       BetaMSE=BetaMSE(true_beta, mc_beta_samples),
                       XBetaMSE=XBetaMSE(true_beta, mc_beta_samples, X),
                       BetaHDCoverage0.95=mean(BetaCoverageHD(true_beta, mc_beta_samples, 0.95))))
    stat_perform_componentwise_df <-
      rbind(stat_perform_componentwise_df,
            data.frame(n=n, p=p, t_dist_df = t_dist_df, component = c(1:p),
                       BetaHDCoverage0.95 =
                         BetaCoverageHD(true_beta, mc_beta_samples, 0.95),
                       ComponentwiseBetaMSE =
                         ComponentwiseBetaMSE(true_beta, mc_beta_samples)))
    print(c(p, t_dist_df, i))
  }
}

traceplots_df <- traceplots_df %>% 
  dplyr::filter(component %in% c('X1','X9','X13','X25'))

# save(stat_perform_df, file = "../BackUpFiles/CoupledHalfT_old/examples/truncated_prior_xi/stat_performance_t_dist_prior.RData")
# save(stat_perform_componentwise_df, file = "../BackUpFiles/CoupledHalfT_old/examples/truncated_prior_xi/stat_perform_componentwise_df.RData")
# save(traceplots_df, file = "../BackUpFiles/CoupledHalfT_old/examples/truncated_prior_xi/traceplots_t_dist_prior.RData")


