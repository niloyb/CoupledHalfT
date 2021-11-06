# Varying d_threshold plot
# Note: depending on configuration %do% may be faster than %dopar%

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
registerDoParallel(cores = detectCores()-2)


# Generate Data
n <- 100
p <- 1000
s <- 20
true_beta <- matrix(0,p,1)
true_beta[1:s] = 2^(-(seq(s)/4-9/4))
X <- matrix(rnorm(n*p), nrow = n, ncol = p)
X_transpose <- t(X)
#Error terms
error_std <- 2
error_terms = error_std*rnorm(n, mean = 0, sd = 1)
y = X%*%true_beta + error_terms

# Simulation parameters
iterations <- 1
mc_chain_size <- 1
max_iterations <- 2e3

threshold <- 0
full_couple_prob_trajectories_crn <-
  foreach(t_dist_df = seq(1,2,0.2), .combine = rbind) %:% 
  foreach(i = 1:iterations, .combine = rbind) %dopar% {
    chain_two_scale_1 <- 
      coupled_half_t_mcmc(mc_chain_size, X, X_transpose, y, 
                          a0=1, b0=1, std_MH=0.8, t_dist_df=t_dist_df,
                          max_iterations=max_iterations, epsilon_eta=threshold)
    # Plot of tvUB
    tvUB_full <- (chain_two_scale_1$metric_d)[2:length(chain_two_scale_1$metric_d)]
    
    print(c(threshold,t_dist_df, i)) # Will not print with %dopar%
    
    # Full couple prob plot
    return(data.frame(n, p, iteration = i, t_dist_df=t_dist_df,
               t= c(1:length(tvUB_full)), metric_d=tvUB_full))
  }

# save(full_couple_prob_trajectories_crn, file = "inst/half_t_degree_of_freedom/t_dist_df_plot_crn.RData")
# load("inst/half_t_degree_of_freedom/t_dist_df_plot_crn.RData")

plots_t_dist_df <- full_couple_prob_trajectories_crn %>%
  dplyr::group_by(n, p, t_dist_df, t) %>%
  dplyr::summarise(metric_d=mean(metric_d), count=n())
t_dist_df_plot_crn <-
  ggplot(plots_t_dist_df, 
         aes(y=metric_d, x=t, color=as.factor(t_dist_df)),
         main= 'Full eta vector couplings prob') + 
  geom_line(size=1) + 
  ylab(TeX('$d(C_{t}^{(1)}, C_{t}^{(2)})$ for $d_{threshold}=0$')) +
  xlab(TeX('iteration')) +
  scale_x_continuous(breaks=c(0,1e3,2e3)) +
  scale_y_log10() + 
  scale_colour_viridis_d(name=TeX('$\\nu$')) +
  theme_classic(base_size = 16) + theme(legend.position = 'right')
t_dist_df_plot_crn
# ggsave(filename = "inst/half_t_degree_of_freedom/t_dist_df_plot_crn.pdf", plot = t_dist_df_plot_crn, width = 4, height = 3)
