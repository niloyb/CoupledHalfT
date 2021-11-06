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
full_couple_prob_trajectories_threshold <- data.frame()
n <- 100
p <- 500
for (threshold in c(0,seq(0.1,0.9,0.2),1)){
  print(threshold)
  results_ <- foreach (i = 1:iterations, .combine = rbind) %dopar% {
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
    max_iterations <- 1e3
    
    chain_two_scale_1 <- 
      coupled_half_t_mcmc(mc_chain_size, X, X_transpose, y, 
                          a0=1, b0=1, std_MH=0.8, 
                          t_dist_df=1, max_iterations=max_iterations, 
                          epsilon_eta=threshold)
    
    # Plot of tvUB
    tvUB_full <- (chain_two_scale_1$metric_d)[2:length(chain_two_scale_1$metric_d)]
    T <- length(tvUB_full)
    if(T<(max_iterations-1)) tvUB_full <- c(tvUB_full, rep(0,(max_iterations-1-T))) 
    
    # Full couple prob plot
    # full_couple_prob_trajectories_threshold <- 
      # rbind(full_couple_prob_trajectories_threshold,
    data.frame(n, p, iteration = i, two_scale_threshold=threshold, 
               t= 1:length(tvUB_full), metric_d=tvUB_full)
  }
  full_couple_prob_trajectories_threshold <-
    rbind(full_couple_prob_trajectories_threshold, results_)
}

# save(full_couple_prob_trajectories_threshold, file = "examples/two_scale_coupling/full_couple_prob_trajectories_threshold.RData")
# load("examples/two_scale_coupling/full_couple_prob_trajectories_threshold.RData")

plots_df_threshold <- full_couple_prob_trajectories_threshold %>%
  dplyr::group_by(n, p, two_scale_threshold, t) %>%
  dplyr::summarise(metric_d=mean(metric_d),
                   count=n())

metric_d_plot_p <-
  ggplot(plots_df_threshold, 
         aes(y=metric_d, x=t, color=as.factor(two_scale_threshold)),
         main= 'Full eta vector couplings prob') + 
  geom_line(size=1) + 
  ylab(TeX('$d(C_{t}^{(1)}, C_{t}^{(2)})$')) +
  xlab(TeX('iteration')) +
  scale_color_manual(name=TeX('$d_{threshold}$'),
                     values = c("grey", viridis::viridis(5), "black")) +
  # scale_colour_viridis_d(name=TeX('$d_{threshold}$')) +
  scale_x_continuous(breaks=c(0,5e2,1e3)) +
  scale_y_log10(breaks=c(1e-8,1e-4,1)) + 
  theme_classic(base_size = 16) + theme(legend.position = 'right')
metric_d_plot_p
#ggsave(filename = "/Users/niloybiswas/Dropbox/horseshoe_coupling/Drafts/images/two_scale_coupling_plot/metric_d_p_500.pdf", plot = metric_d_plot_p, width = 4, height = 3)




