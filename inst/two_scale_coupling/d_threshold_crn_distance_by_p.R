# Varying d_threshold plot

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

iterations <- 10
metric_trajectories_crn_by_p <- data.frame()
n <- 100
threshold <- 0
for (p in seq(500,1000,100)){
  cat("p =", p, "\n")
  s <- 20
  true_beta <- matrix(0,p,1)
  true_beta[1:s] = 2^(-(seq(s)/4-9/4))
  
  X <- matrix(rnorm(n*p), nrow = n, ncol = p)
  X_transpose <- t(X)
  #Error terms
  error_std <- 2
  error_terms = error_std*rnorm(n, mean = 0, sd = 1)
  y = X%*%true_beta + error_terms
  
  mc_chain_size <- 1
  max_iterations <- 2e3
  
  results_ <- foreach (i = 1:iterations, .combine = rbind) %dopar% {
    chain_two_scale_1 <- 
      coupled_half_t_mcmc(mc_chain_size, X, X_transpose, y, 
                          a0=1, b0=1, std_MH=0.8,
                          t_dist_df=1, max_iterations=max_iterations,
                          epsilon_eta=threshold)
    # Plot of tvUB
    tvUB_full <- (chain_two_scale_1$metric_d)[2:length(chain_two_scale_1$metric_d)]
    data.frame(n, p, iteration = i, two_scale_threshold=threshold,
               t= c(1:length(tvUB_full)), metric_d=tvUB_full)
  }
  metric_trajectories_crn_by_p <- rbind(metric_trajectories_crn_by_p, results_)
}

# Saving/ loading data
# save(metric_trajectories_crn_by_p, file = "examples/two_scale_coupling/metric_trajectories_crn_by_p.RData")
# load("examples/two_scale_coupling/metric_trajectories_crn_by_p.RData")

plots_df_threshold <- metric_trajectories_crn_by_p %>%
  dplyr::group_by(n, p, two_scale_threshold, t) %>%
  dplyr::summarise(metric_d=mean(metric_d), count=n())
metric_d_hist_crn <-
  ggplot(data = metric_trajectories_crn_by_p, aes(y = p)) +
  geom_density_ridges(alpha=0.8, aes(x = metric_d, group = p, fill = as.factor(two_scale_threshold)))+
  # geom_density_ridges(alpha=0.8, aes(x = metric_d, group = p, fill = as.factor(p)))+
  ylab(TeX('Dimension $p$')) +
  xlab(TeX('$(d(C_t^{(1)}, C_t^{(2)}))_{t=0}^{2000}$ histogram')) +
  scale_y_continuous(breaks=seq(600,1e3,2e2)) +
  scale_x_continuous(breaks=seq(0,1,0.5)) +
  scale_fill_manual(name=TeX('$d_{threshold}$'),
                     values = c("grey")) +
  # scale_fill_viridis_d(name=TeX('$p$')) +
  theme_classic(base_size = 16) # + theme(legend.position = 'none')
metric_d_hist_crn
#ggsave(filename = "examples/two_scale_coupling/metric_d_hist_crn.pdf", plot = metric_d_hist_crn, width = 4, height = 3)


