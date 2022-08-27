## R script for parallel computing on a cluster using Slurm

# Set working directory
setwd("/n/home11/nbiswas/coupledHalfT/")

# Set library path
.libPaths('/n/home11/nbiswas/r_packages_installed/')

# Libraries
rm(list = ls())
seed <- 1
set.seed(seed)
library(foreach)
library(parallel)
library(doParallel)
registerDoParallel(detectCores())

# Loading functions
function_scripts <- list.files("functions", full.names = TRUE)
sapply(function_scripts, source)

# Load GWAS data
load(file = "/n/home11/nbiswas/datasets/maize/design_matrix_Xnew.RData")
load(file = "/n/home11/nbiswas/datasets/maize/response_ynew.RData")

X_transpose <- t(X)
n <- length(y)
p <-dim(X)[2]

max_iterations <- 2000

t_dist_df <- 3
iterations <- 1
lag <- 1
# nrepeats_eta <- 1
print('Here we go')
gwas_meet_df <- 
  foreach(nrepeats_eta = seq(1,9,2), .combine = rbind) %:% 
  foreach(i = 1:iterations, .combine = rbind) %dopar% {
  meeting_time_2 <- 
    meetingtime_half_t(X, X_transpose, y, a0=1, b0=1, std_MH=0.8,
                       rinit=NULL, epsilon_eta = 0.5, 
                       max_iterations=max_iterations, verbose = TRUE,
                       nrepeats_eta = nrepeats_eta, t_dist_df=t_dist_df, lag=lag)
  return(data.frame(seed=seed, iteration = i, lag, meeting_time_2, t_dist_df, nrepeats_eta))
}
write.table(gwas_meet_df, "gwas_simulations/GWASLlag_new.csv", sep = ",", 
            col.names = !file.exists("gwas_simulations/GWASLlag_new.csv"), 
            append = TRUE, row.names = FALSE)

