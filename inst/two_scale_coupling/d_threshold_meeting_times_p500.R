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

n <- 100
p <- 500
iterations <- 10

meetingtimes_df <- 
  foreach(threshold = seq(0.9,0.1,-0.2), .combine = rbind) %:% 
  foreach(i = 1:iterations, .combine = rbind) %dopar% {
    s <- 20
    true_beta <- matrix(0,p,1)
    true_beta[1:s] = 2^(-(seq(s)/4-9/4))
    X <- matrix(rnorm(n*p), nrow = n, ncol = p)
    X_transpose <- t(X)
    #Error terms
    error_std <- 2
    error_terms = error_std*rnorm(n, mean = 0, sd = 1)
    y = X%*%true_beta + error_terms
    
    t_dist_df <- 1
    max_iterations <- 1e5
    lag <- 1
    
    meetingtime <- 
      meetingtime_half_t(X, X_transpose, y, a0=1, b0=1, std_MH = 0.8, 
                         epsilon_eta = threshold, max_iterations = max_iterations, 
                         t_dist_df=t_dist_df, lag=lag)
    print(c(n,p,s,error_std,i, meetingtime$meetingtime))  # Will not print with %dopar%
    return(data.frame(n=n, p=p, meetingtime=meetingtime,
                      threshold=threshold, t_dist_df = t_dist_df, lag=lag))
  }

# save(meetingtimes_df, file = "examples/two_scale_coupling/meeting_times_p500.RData")
# load("examples/two_scale_coupling/meeting_times_p500.RData")

plot_meeting_times_p <-
  ggplot(data = meetingtimes_df, aes(y = threshold)) +
  geom_density_ridges(alpha=0.8, aes(x = meetingtime.meetingtime, 
                                     group = threshold,
                                     fill = as.factor(threshold))) +
  scale_fill_viridis_d(name=TeX('$d_{threshold}$')) +
  xlab(TeX('Meeting time $\\tau$')) +
  ylab(TeX('Metric threshold $d_{threshold}$')) +
  scale_y_continuous(breaks=c(0,0.5,1)) +
  #scale_y_log10() +
  theme_classic(base_size = 16) + coord_flip() # + theme(legend.position = 'none')
plot_meeting_times_p
#ggsave(filename = "examples/two_scale_coupling/meeting_times_p500.pdf", plot = plot_meeting_times_p, width = 4, height = 3)

