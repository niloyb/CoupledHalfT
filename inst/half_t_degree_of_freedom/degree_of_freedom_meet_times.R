# Plots showing coupling times for difference degrees of freedom 

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

## ## ## ## ## 
## n<=p plot ## 
## ## ## ## ##

# Two scaling coupling simulation setup
iterations <- 100

epsilon_eta <- 0.5

# Varying p, nu simulations
n <- 100

tdist_meetingtimes_df2 <-
  foreach(t_dist_df = seq(2,1,-0.2), .combine = rbind) %:% 
  foreach(p = seq(100,300,100), .combine = rbind) %:% 
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
    
    max_iterations <- Inf
    meetingtime <- 
      meetingtime_half_t(X, X_transpose, y, a0 = 1, b0 = 1, std_MH = 0.8, 
                         epsilon_eta = epsilon_eta, 
                         max_iterations = max_iterations, t_dist_df=t_dist_df,
                         verbose = FALSE)
    print(c(t_dist_df, p, i, meetingtime$meetingtime))  # Will not print with %dopar%
    return(data.frame(n = n, p =p, sparsity = s, error_std = error_std, 
                      meetingtime = meetingtime, t_dist_df = t_dist_df))
  }

# # nu=2, varying p simulations
# tdist_meetingtimes_dfv2 <-
#   foreach(t_dist_df = seq(2,2,0), .combine = rbind) %:% 
#   foreach(p = seq(2e3,1e4,2e3), .combine = rbind) %:% 
#   foreach(i = 1:iterations, .combine = rbind) %dopar% {
#     ## Generate data
#     s <- 20
#     true_beta <- matrix(0,p,1)
#     true_beta[1:s] = 2^(-(seq(s)/4-9/4))
#     X <- matrix(rnorm(n*p), nrow = n, ncol = p)
#     X_transpose <- t(X)
#     #Error terms
#     error_std <- 2
#     error_terms = error_std*rnorm(n, mean = 0, sd = 1)
#     y = X%*%true_beta + error_terms
#     
#     max_iterations <- Inf
#     meetingtime <- 
#       meetingtime_half_t(X, X_transpose, y, a0 = 1, b0 = 1, std_MH = 0.8, 
#                          epsilon_eta = epsilon_eta, 
#                          max_iterations = max_iterations, t_dist_df=t_dist_df)
#     print(c(t_dist_df, p, i, meetingtime$meetingtime))  # Will not print with %dopar%
#     return(data.frame(n = n, p =p, sparsity = s, error_std = error_std,
#                       meetingtime = meetingtime, t_dist_df = t_dist_df))
#   }
# tdist_meetingtimes_df <- rbind(tdist_meetingtimes_df, tdist_meetingtimes_dfv3)

# save(tdist_meetingtimes_df, file = "inst/half_t_degree_of_freedom/degree_of_freedom_meeting_times_new.RData") # meetingtimes_df
# load("inst/half_t_degree_of_freedom/degree_of_freedom_meeting_times_new.RData") # meetingtimes_df_two_scale

# Plots
plot_vary_p <-
  ggplot(data = tdist_meetingtimes_df %>% dplyr::filter(p <= 500),
         aes(y = p)) +
  geom_density_ridges(alpha=0.8, scale = 0.8,
                      aes(x = meetingtime.meetingtime, 
                          group = interaction(p, t_dist_df), 
                          fill = as.factor(t_dist_df))) +
  scale_fill_viridis_d(name=TeX('$\\nu$')) +
  scale_x_log10() +
  scale_y_continuous(breaks = seq(0,1e3,100)) +
  xlab(TeX('Meeting time $\\tau$')) +
  ylab(TeX('Dimension $p$')) +
  theme_classic(base_size = 16) + coord_flip() +
  theme(legend.position = 'right')
plot_vary_p
# ggsave(filename = "examples/half_t_degree_of_freedom/half_t_degree_of_freedom_vary_p.pdf", plot = plot_vary_p, width = 4, height = 3)
# ggsave(filename = "/Users/niloybiswas/Dropbox/horseshoe_coupling/New_plots_May_2022/half_t_degree_of_freedom_vary_p.pdf", plot = plot_vary_p, width = 4, height = 3)

plot_nu_2_vary_p <-
  ggplot(data = tdist_meetingtimes_df %>% 
           dplyr::filter(t_dist_df ==2, p %in% seq(2e3,2e4,2e3)),
         aes(y = p)) +
  geom_density_ridges(alpha=0.8, scale = 0.8,
                      aes(x = meetingtime.meetingtime, 
                          group = interaction(p, t_dist_df), 
                          fill = as.factor(t_dist_df))) +
  scale_fill_viridis_d(name=TeX('$\\nu$'),direction = -1) +
  scale_x_log10() + 
  xlab(TeX('Meeting time $\\tau$')) +
  ylab(TeX('Dimension $p$')) +
  scale_y_continuous(breaks = seq(0,2e4,1e4)) +
  theme_classic(base_size = 16) + coord_flip() +
  theme(legend.position = 'right')
plot_nu_2_vary_p
# ggsave(filename = "examples/half_t_degree_of_freedom/half_t_degree_of_freedom_nu_2_vary_p.pdf", plot = plot_nu_2_vary_p, width = 4, height = 3)

