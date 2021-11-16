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
# Number of chains
nchains <- 1000
chain_length <- 5e2
burnin <- 0

################################################################################
# Generate Synthetic data
n <- 20
p <- 2
s <- 10
true_beta <- matrix(0,p,1)
# true_beta[1:s] = 2^(-(seq(s)/4-9/4))
X <- matrix(rnorm(n*p), nrow = n, ncol = p)
X_transpose <- t(X)
#Error terms
error_std <- 0.5
error_terms = error_std*rnorm(n, mean = 0, sd = 1)
y = X%*%true_beta + error_terms

save(X, file='../BackUpFiles/CoupledHalfT_old/examples/big_data_examples/stein_X.RData')
save(y, file='../BackUpFiles/CoupledHalfT_old/examples/big_data_examples/stein_y.RData')

################################################################################
# Prior parameters
t_dist_df <- 2
xi_interval <- c(0,Inf)
verbose=FALSE

half_t_chains <- foreach(i=c(1:nchains))%dopar%{
  half_t_chain <- half_t_mcmc(chain_length,burnin,X,X_transpose,y, t_dist_df=t_dist_df)
}

save(half_t_chains, file='../BackUpFiles/CoupledHalfT_old/examples/big_data_examples/stein_chain_data.RData')






################################################################################
# For "exact" results, we will simulate many MCMC chains and 
# approximate the densities using those samples.
initial_value <- 10
sd_proposal <- 0.5

chain_length <- 300
number_of_chains <- 1000

all_chains <- foreach(i = 1:number_of_chains, .combine = rbind) %dorng% {
  X <- matrix(0, nrow = 1, ncol = chain_length)
  X[1] <- initial_value
  for (t in 2:chain_length) {
    proposal <- rnorm(1, mean = X[t - 1], sd = sd_proposal)
    if (log(runif(1)) < -0.5 * (proposal ^ 2 - X[t - 1] ^ 2)) {
      X[t] <- proposal
    } else {
      X[t] <- X[t - 1]
    }
  }
  X
}
all_chains

# save(all_chains, file='../BackUpFiles/CoupledHalfT_old/examples/big_data_examples/stein_rwmh_chain_data.RData')

