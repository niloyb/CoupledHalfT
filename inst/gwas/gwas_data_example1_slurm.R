## R script for parallel computing on a cluster using Slurm

# Set working directory
setwd("/n/home11/nbiswas/coupledHalfT/")

# Set library path
.libPaths('/n/home11/nbiswas/r_packages_installed/')

# Libraries
rm(list = ls())
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

max_iterations <- 2500

t_dist_df <- 2
iterations <- 20
lag <- 1250

print('Here we go')
gwas_meet_df <- foreach(i = 1:iterations, .combine = rbind) %dopar% {
  meeting_time_2 <- 
    meetingtime_half_t(X, X_transpose, y, max_iterations=max_iterations, 
                       verbose = TRUE, t_dist_df=t_dist_df, lag=lag)
  return(data.frame(seed=seed, iteration = i, lag, meeting_time_2, t_dist_df))
}
write.table(gwas_meet_df, "gwas_simulations/GWASLlag_new.csv", sep = ",", 
            col.names = !file.exists("gwas_simulations/GWASLlag_new.csv"), 
            append = TRUE, row.names = FALSE)

