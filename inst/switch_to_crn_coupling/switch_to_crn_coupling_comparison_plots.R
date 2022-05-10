# Plots of meeting times for two-scale coupling and switch to CRN coupling
# Note: depending on configuration %do% may be faster than %dopar%

# Libraries
rm(list = ls())
set.seed(1)
library(CoupledHalfT)
library(doParallel)
registerDoParallel(cores = detectCores()-2)

library(dplyr)
library(ggplot2)
library(ggridges)
library(latex2exp)

## ## ## ## ## 
## n<=p plot ## 
## ## ## ## ##

# Simulation setup
iterations <- 100
meetingtimes_df <- data.frame()

# Horseshoe prior
t_dist_df <- 1

# Two scaling coupling
epsilon_eta <- 0.5

# Fixed n varying p simulations
n <- 100
meetingtimes_df_p <- 
  foreach(p = seq(200,1000,200), .combine = rbind) %:% 
  foreach(two_scale = c(TRUE, FALSE), .combine = rbind) %:% 
  foreach(i = c(1:iterations), .combine = rbind) %do% {
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
                         epsilon_eta = epsilon_eta, max_iterations = max_iterations, 
                         t_dist_df=t_dist_df, two_scale=two_scale)
    print(c(n,p,s,error_std,i, meetingtime$meetingtime))  # Will not print with %dopar%
    
    data.frame(n = n, p =p, meetingtime = meetingtime, t_dist_df = t_dist_df, 
               sparsity = s, error_std = error_std, two_scale=two_scale)
  }
meetingtimes_df <- rbind(meetingtimes_df, meetingtimes_df_p)
# save(meetingtimes_df, file = "inst/switch_to_crn_coupling/switch_to_crn__two_scale_coupling_plot.RData")
# load("inst/switch_to_crn_coupling/switch_to_crn__two_scale_coupling_plot.RData") # meetingtimes_df_two_scale

# One scale, two scale and switch to CRN comparison
plot_vary_p <-
  ggplot(data = meetingtimes_df, aes(y = p)) +
  geom_density_ridges(alpha=0.85, scale = 1, 
                      aes(x = meetingtime.meetingtime, group = interaction(p, two_scale), fill = two_scale)) +
  scale_fill_manual(name="Coupling", breaks = c("TRUE", "FALSE"),
                    labels = c("Two-scale", "Switch-to-CRN"),
                    values = c("lightgrey", "white"), guide = "legend") +
  scale_x_log10() + 
  xlab(TeX('Meeting time $\\tau$')) + ylab(TeX('Dimension $p$')) +
  scale_y_continuous(breaks = seq(200,1e3,200)) +
  theme_classic(base_size = 14) + coord_flip() +
  theme(legend.position = 'right')
plot_vary_p
# ggsave(filename = "/Users/niloybiswas/Dropbox/horseshoe_coupling/Drafts/images/switch_to_crn_coupling_plot/switch_to_crn_coupling_two_scale_comparison_plot_vary_p.pdf", plot = plot_vary_p, width = 5, height = 3) # base_size = 14
ggsave(filename = "/Users/niloybiswas/Dropbox/horseshoe_coupling/New_plots_May_2022/switch_to_crn_coupling_two_scale_comparison_plot_vary_p.pdf", plot = plot_vary_p, width = 5, height = 3) # base_size = 14


