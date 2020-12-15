# Plots of meeting times for two-scale coupling
# Note: depending on configuration %do% may be faster than %dopar%

# Libraries
rm(list = ls())
set.seed(1)
library(CoupledHalfT)
library(doParallel)
registerDoParallel(cores = detectCores()-2)

## ## ## ## ## 
## n<=p plot ## 
## ## ## ## ##

# Simulation setup
iterations <- 10
meetingtimes_df <- data.frame()

# Horseshoe prior
t_dist_df <- 1

# Two scaling coupling
epsilon_eta <- 0.5

# Fixed n varying p simulations
n <- 100
meetingtimes_df_p <- 
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
    
    mc_chain_size <- 1
    burnin <- 0
    max_iterations <- Inf
    meetingtime <- 
      meetingtime_half_t(X, X_transpose, y, a0 = 1, b0 = 1, std_MH = 0.8, 
                         epsilon_eta = epsilon_eta, 
                         max_iterations = max_iterations, t_dist_df=t_dist_df)
    print(c(n,p,s,error_std,i, meetingtime$meetingtime))  # Will not print with %dopar%
    gc()
    data.frame(n = n, p =p, meetingtime = meetingtime, 
               t_dist_df = t_dist_df, sparsity = s, error_std = error_std)
  }
meetingtimes_df <- rbind(meetingtimes_df, meetingtimes_df_p)

# Varying n simulations
p <- 150
meetingtimes_df_n <- 
  foreach(n=seq(120,200,20), .combine = rbind) %:% 
  foreach(i=c(1:iterations), .combine = rbind) %dopar% {
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
    
    mc_chain_size <- 1
    burnin <- 0
    max_iterations <- 1e+05
    meetingtime <- 
      meetingtime_half_t(X, X_transpose, y, a0 = 1, b0 = 1, std_MH = 0.8, 
                         epsilon_eta = epsilon_eta, 
                         max_iterations = max_iterations, t_dist_df=t_dist_df)
    print(c(n,p,s,error_std,i, meetingtime$meetingtime)) # Will not print with %dopar%
    gc()
    data.frame(n = n, p =p, meetingtime = meetingtime, 
               t_dist_df = t_dist_df, sparsity = s, error_std = error_std)
  }
meetingtimes_df <- rbind(meetingtimes_df, meetingtimes_df_n)

# Varying sparsity s simulations
n <- 100
p <- 150
meetingtimes_df_sparsity <- 
  foreach(s=seq(0,8,2), .combine = rbind) %:% 
    foreach(i=c(1:iterations), .combine = rbind) %dopar% {
      ## Generate data
      true_beta <- matrix(0,p,1)
      if (s>0){
        true_beta[1:s] = 2^(-(seq(s)/4-9/4))
      }
      X <- matrix(rnorm(n*p), nrow = n, ncol = p)
      X_transpose <- t(X)
      #Error terms
      error_std <- 0.5
      error_terms = error_std*rnorm(n, mean = 0, sd = 1)
      y = X%*%true_beta + error_terms
      
      mc_chain_size <- 1
      burnin <- 0
      max_iterations <- 1e+05
      meetingtime <- 
        meetingtime_half_t(X, X_transpose, y, a0 = 1, b0 = 1, std_MH = 0.8, 
                           epsilon_eta = epsilon_eta, 
                           max_iterations = max_iterations, t_dist_df=t_dist_df)
      print(c(n,p,s,error_std,i, meetingtime$meetingtime)) # Will not print with %dopar%
      gc()
      data.frame(n = n, p =p, meetingtime = meetingtime, 
                        t_dist_df = t_dist_df, sparsity = s, error_std = error_std)
  }
meetingtimes_df <- rbind(meetingtimes_df, meetingtimes_df_sparsity)

# Varying error_std simulations
n <- 100
p <- 150
meetingtimes_df_error_std <- 
  foreach(error_std=seq(0,0.4,0.1), .combine = rbind) %:%
    foreach(i=c(1:iterations), .combine = rbind) %dopar% {
      ## Generate data
      s <- 10
      true_beta <- matrix(0,p,1)
      true_beta[1:s] = 2^(-(seq(s)/4-9/4))
      X <- matrix(rnorm(n*p), nrow = n, ncol = p)
      X_transpose <- t(X)
      #Error terms
      error_terms = error_std*rnorm(n, mean = 0, sd = 1)
      y = X%*%true_beta + error_terms
      
      mc_chain_size <- 1
      burnin <- 0
      max_iterations <- 1e+05
      meetingtime <- 
        meetingtime_half_t(X, X_transpose, y, a0 = 1, b0 = 1, std_MH = 0.8, 
                           epsilon_eta = epsilon_eta, 
                           max_iterations = max_iterations, t_dist_df=t_dist_df)
      print(c(n,p,s,error_std,i, meetingtime$meetingtime)) # Will not print with %dopar%
      gc()
      data.frame(n = n, p =p, meetingtime = meetingtime, 
                        t_dist_df = t_dist_df, sparsity = s, error_std = error_std)
  }
meetingtimes_df <- rbind(meetingtimes_df, meetingtimes_df_error_std)

# Save data
# save(meetingtimes_df, file = "examples/two_scale_coupling/two_scale_coupling_plot.RData")

save(meetingtimes_df, file = "two_scale_coupling_plot.RData")
