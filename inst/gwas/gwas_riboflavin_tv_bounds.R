# Libraries
rm(list = ls())
set.seed(1)

library(doParallel)
library(doRNG)
registerDoParallel(cores = detectCores())

library(CoupledHalfT)
library(ggplot2)
library(latex2exp)
library(dplyr)

################################################################################
# Loading riboflavin dataset
data(riboflavin)

y <- as.vector(riboflavin$y)
X <- as.matrix(riboflavin$x)
colnames(X) <- NULL
rownames(X) <- NULL
X_transpose <- t(X)
n <- length(y)
p <-dim(X)[2]

mc_chain_size <- 1
burnin <- 0
max_iterations <- 2000
t_dist_df <- 2
lag <- 750
iterations <- 100

# gwas_meet_nrepeats_eta <-
#   foreach(nrepeats_eta = seq(1,5,1), .combine = rbind) %:%
#   foreach(i = 1:iterations, .combine = rbind) %dopar% {
#   meeting_time <-
#     meetingtime_half_t(X, X_transpose, y, a0=1, b0=1, std_MH=0.8,
#                        rinit=NULL, epsilon_eta = 0.5, max_iterations=max_iterations,
#                        nrepeats_eta=nrepeats_eta, verbose = TRUE, t_dist_df=t_dist_df, lag=1)
#   return(data.frame(iteration = i, lag, meeting_time, t_dist_df, nrepeats_eta))
# }
# gwas_meet_nrepeats_eta %>%
#   group_by(nrepeats_eta) %>%
#   summarise(meet_time_mean = mean(meetingtime), meet_time_max = max(meetingtime))

gwas_meet_df <- 
  foreach(i = 1:iterations, .combine = rbind) %dopar% {
    nrepeats_eta=1
    meeting_time <- 
      meetingtime_half_t(X, X_transpose, y, a0=1, b0=1, std_MH=0.8,
                         rinit=NULL, epsilon_eta = 0.5, max_iterations=max_iterations, 
                         nrepeats_eta=nrepeats_eta, verbose = TRUE, t_dist_df=t_dist_df, 
                         lag=lag)
    return(data.frame(iteration = i, lag, meeting_time, t_dist_df, nrepeats_eta))
  }
# save(gwas_meet_df, file = 'inst/gwas/gwas_riboflavin.RData')

# load(file = 'gwas/gwas_riboflavin.RData')
# Plotting L-Lag TV bounds
tv_upper_bound_estimates <- 
  function(coupling_times, L, t){return(pmax(0,ceiling((coupling_times-L-t)/L)))}
time <- c(0:1e3)
tv_ub <- colMeans(sapply(time, function(x) tv_upper_bound_estimates(coupling_times = gwas_meet_df$meetingtime, L=unique(gwas_meet_df$lag), x)))
tv_ub_df <- data.frame(t=time, t_dist_df=t_dist_df, tv_ub=tv_ub)

tv_ub_plot_riboflavin <- 
  ggplot(tv_ub_df %>% dplyr::filter(t <= 1e3), aes(x=t, y=tv_ub)) +
  geom_line() + # scale_x_continuous(trans='log10') +
  # geom_errorbar(aes(ymin=tv_ub-tv_ub_sd/sqrt(indep_runs), ymax=tv_ub+tv_ub_sd/sqrt(indep_runs)), width=.2, position=position_dodge(.9)) +
  ylab(TeX('Total variation distance')) + xlab(TeX('iteration')) +
  #scale_linetype_manual(name='Dataset', values = c('solid','dashed')) +
  scale_y_continuous(breaks = seq(0,1,0.25)) +
  guides(linetype=guide_legend(nrow=2,byrow=TRUE)) +
  theme_classic(base_size = 14) + 
  theme(legend.position = 'right', legend.key.width=unit(1,"cm"))
tv_ub_plot_riboflavin
# ggsave(filename = "images/dataset_examples/gwas_LlagTVUB_riboflavin.pdf", plot = tv_ub_plot_riboflavin, width = 4, height = 2.75)
# ggsave(filename = "/Users/niloybiswas/Dropbox/horseshoe_coupling/New_plots_May_2022/gwas_LlagTVUB_riboflavin.pdf", plot = tv_ub_plot_riboflavin, width = 4, height = 2.75)
# ggsave(filename = "/Users/niloybiswas/Dropbox/horseshoe_coupling/New_plots_May_2022/gwas_LlagTVUB_riboflavin_nrepeats_eta=5.pdf", plot = tv_ub_plot_riboflavin, width = 4, height = 2.75)



