# Metric d vs. Euclidean distance plot

# Libraries
rm(list = ls())
set.seed(1)
library(CoupledHalfT)
library(dplyr)
library(ggplot2)
library(latex2exp)
library(gridExtra)
library(ggridges)

iterations <- 1 # A single trajectory considered
metric_trajectories_crn_by_p <- data.frame()
n <- 100
threshold <- 0 # Common random numbers coupling
for (p in seq(1000,1000,100)){
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
  max_iterations <- 1e3
  
  for (i in c(1:iterations)){ # Indep trajectories for a fixed dataset
    chain_two_scale_1 <- 
      coupled_half_t_mcmc(mc_chain_size, X, X_transpose, y, 
                          a0=1, b0=1, std_MH=0.8,
                          t_dist_df=1, max_iterations=max_iterations, 
                          epsilon_eta=threshold)
    
    # Plot of tvUB
    tvUB_full <- (chain_two_scale_1$metric_d)[2:length(chain_two_scale_1$metric_d)]
    rates_1 <- chain_two_scale_1$beta_samples1[c(2:max_iterations),]^2*chain_two_scale_1$xi_samples1[c(2:max_iterations)]/(2*chain_two_scale_1$sigma2_samples1[c(2:max_iterations)])
    rates_2 <- chain_two_scale_1$beta_samples2[c(1:(max_iterations-1)),]^2*chain_two_scale_1$xi_samples2[c(1:(max_iterations-1))]/(2*chain_two_scale_1$sigma2_samples2[c(1:(max_iterations-1))])
    l1_distance <- rowSums(abs(rates_1-rates_2))
    l1_log_distance <- rowSums(abs(log(rates_1/rates_2)))
    
    # Full couple prob plot
    metric_trajectories_crn_by_p <- 
      rbind(metric_trajectories_crn_by_p,
            data.frame(n, p, iteration = i, two_scale_threshold=threshold,
                       t= c(1:length(tvUB_full)), metric_d=tvUB_full, 
                       metric_l1=l1_distance,
                       metric_l1_log=l1_log_distance))
    print(c(threshold,i))
  }
}

# Saving/ loading data
# save(metric_trajectories_crn_by_p, file = "examples/two_scale_coupling/metric_trajectories_crn_by_p.RData")
# load("examples/two_scale_coupling/metric_trajectories_crn_by_p.RData")


# Plots
plots_df_metric <- metric_trajectories_crn_by_p %>%
  dplyr::group_by(n, p, two_scale_threshold, t) %>%
  dplyr::summarise(metric_d=mean(metric_d), metric_l1=mean(metric_l1),
                   metric_l1_log=mean(metric_l1_log), count=n()) %>%
  tidyr::gather(metric, value, metric_d:metric_l1_log)

metric_d_plot_crn <-
  ggplot(plots_df_metric %>% filter(t <= 1000), 
         aes(y=value, x=t, color=as.factor(metric)),
         main= 'Full eta vector couplings prob') + 
  geom_line(size=1) + 
  ylab(TeX('Distance')) +
  xlab(TeX('iteration')) +
  #scale_x_continuous(breaks=c(0,1e3,2e3)) +
  #scale_y_continuous(limits = c(0,1)) +
  scale_y_log10() + 
  scale_colour_viridis_d(name=TeX('$Metric$'), 
                         labels = unname(TeX(c("$d(C, \\tilde{C})", 
                                               "$\\sum_j | m_j- \\tilde{m}_j  |",
                                               "$\\sum_j | \\log(m_j / \\tilde{m}_j) |")))) +
  theme_classic(base_size = 18) + theme(legend.position = 'right')
metric_d_plot_crn
#ggsave(filename = "/Users/niloybiswas/Dropbox/horseshoe_coupling/Drafts/images/two_scale_coupling_plot/choice_of_metric_d_plot.pdf", plot = metric_d_plot_crn, width = 8, height = 3)


