# Plots showing coupling times for difference degrees of freedom 
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
registerDoParallel(cores = detectCores())

## ## ## ## ## 
## n<=p plot ## 
## ## ## ## ##

# Two scaling coupling simulation setup
iterations <- 10
epsilon_eta <- 0.5

# Varying p simulations
n <- 100

tdist_meetingtimes_df <-
  foreach(t_dist_df = seq(2,1,-0.2), .combine = rbind) %:% 
  foreach(p = seq(100,500,100), .combine = rbind) %:% 
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
    
    lag <- 500
    
    max_iterations <- Inf
    meetingtime <- 
      meetingtime_half_t(X, X_transpose, y, a0 = 1, b0 = 1, std_MH = 0.8, 
                         epsilon_eta = epsilon_eta, max_iterations = max_iterations, 
                         t_dist_df=t_dist_df, lag=lag)
    print(c(t_dist_df, p, i, meetingtime$meetingtime)) # Will not print with %dopar%
    return(data.frame(n = n, p =p, sparsity = s, error_std = error_std,
                      meetingtime = meetingtime, t_dist_df = t_dist_df,
                      lag=lag))
  }

# t_dist_df=2 simulations for larger p  
tdist_meetingtimes_df_nu2 <-
  foreach(t_dist_df = seq(2,2,0), .combine = rbind) %:% 
  foreach(p = seq(2e3,2e4,2e3), .combine = rbind) %:% 
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
    
    lag <- 500
    max_iterations <- Inf
    meetingtime <- 
      meetingtime_half_t(X, X_transpose, y, a0 = 1, b0 = 1, std_MH = 0.8, 
                         epsilon_eta = epsilon_eta, max_iterations = max_iterations, 
                         t_dist_df=t_dist_df, lag=lag)
    print(c(t_dist_df, p, i, meetingtime$meetingtime)) # Will not print with %dopar%
    return(data.frame(n = n, p =p, sparsity = s, error_std = error_std,
                      meetingtime = meetingtime, t_dist_df = t_dist_df,
                      lag=lag))
  }
tdist_meetingtimes_df <- rbind(tdist_meetingtimes_df, tdist_meetingtimes_df_nu2)

# save(tdist_meetingtimes_df, file = "examples/half_t_degree_of_freedom/lagged_meetingtimes_df.RData")
# load(file = "examples/half_t_degree_of_freedom/lagged_meetingtimes_df.RData")

# Plotting L-Lag TV bounds
tv_upper_bound_estimates <- 
  function(coupling_times, L, t){return(pmax(0,ceiling((coupling_times-L-t)/L)))}

# TV UB from a dataframe
tv_upper_bound_estimates_df <- function(df, time){
  Llagmeetingtimes <- df$meetingtime.meetingtime
  lag <- (df$lag) # df$lag should have one unique value
  tv_ub <- colMeans(sapply(time, function(x) tv_upper_bound_estimates(coupling_times = Llagmeetingtimes, L=lag, x)))
  return(data.frame(t=time, tv_ub=tv_ub))
}

# TV UB Plot
time <- c(0:600)
tv_ub_df <- tdist_meetingtimes_df %>% 
  dplyr::group_by(n,p,t_dist_df) %>%
  dplyr::do(tv_upper_bound_estimates_df(df=., time=time))
tv_ub_plot <- ggplot(tv_ub_df %>% dplyr::filter(p==500),
                     aes(x=t, y=tv_ub, col=as.factor(t_dist_df))) + 
  geom_line(size=1) + 
  scale_color_viridis_d(name=TeX('$\\nu$')) +
  # scale_x_log10() + 
  scale_y_continuous(breaks = seq(0,1,0.5)) +
  ylab(TeX('Total variation distance')) +
  xlab(TeX('iteration')) +
  theme_classic(base_size = 16) + theme(legend.position = 'right')
tv_ub_plot
#ggsave(filename = "examples/half_t_degree_of_freedom/half_t_llag_tv_bounds.pdf", plot = tv_ub_plot, width = 4, height = 3)

mixing_epsilon <- 0.25
mixing_times_df <- tv_ub_df %>%
  dplyr::group_by(n,p,t_dist_df) %>%
  dplyr::filter(tv_ub<mixing_epsilon) %>%
  dplyr::filter(row_number()==1) %>%
  dplyr::rename(mixing_time=t)
# Mixing time plot
mixing_times_plot <- ggplot(mixing_times_df  %>% dplyr::filter(p<=500),
                     aes(x=p, y=mixing_time, col=as.factor(t_dist_df))) + 
  geom_line(size=1) +
  scale_color_viridis_d(name=TeX('$\\nu$')) +
  # scale_y_continuous(limits = c(0,1e4)) + 
  # scale_y_log10() + 
  ylab(TeX('Mixing time $t_{mix}(0.25)$')) +
  xlab(TeX('Dimension $p$')) +
  theme_classic(base_size = 16) + theme(legend.position = 'right')
mixing_times_plot
#ggsave(filename = "examples/half_t_degree_of_freedom/half_t_mixing_times.pdf", plot = mixing_times_plot, width = 4, height = 3)

# Mixing time plot 2: t_dist_df==2
mixing_times_plot_degree_2 <- 
  ggplot(mixing_times_df  %>% dplyr::filter(t_dist_df==2, p>=2e3),
         aes(x=p, y=mixing_time, col=as.factor(t_dist_df))) + 
  geom_line(size=1) +
  scale_color_viridis_d(name=TeX('$\\nu$'), direction = -1) +
  #scale_y_log10() + 
  scale_y_continuous(breaks = seq(0,200,50), limits = c(0,150)) +
  scale_x_continuous(breaks = seq(0,2e4,1e4)) +
  ylab(TeX('Mixing time $t_{mix}(0.25)$')) +
  xlab(TeX('Dimension $p$')) +
  theme_classic(base_size = 16) + theme(legend.position = 'right')
mixing_times_plot_degree_2
#ggsave(filename = "examples/half_t_degree_of_freedom/half_t_mixing_times_degree_2.pdf", plot = mixing_times_plot_degree_2, width = 4, height = 3)

