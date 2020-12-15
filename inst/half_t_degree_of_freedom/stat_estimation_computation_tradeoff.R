# Plot showing the tradeoff between statistical estimation and computation

# Libraries
rm(list = ls())
set.seed(1)
library(dplyr)
library(ggplot2)
library(latex2exp)
library(gridExtra)
library(ggridges)

# load("examples/half_t_degree_of_freedom/degree_of_freedom_meeting_times.RData") # meetingtimes_df_two_scale
# Plots
meetings_by_t_dist_df <-
  ggplot(data = tdist_meetingtimes_df %>% 
           dplyr::filter(p==1000),
         aes(y = meetingtime.meetingtime,x=t_dist_df, group=t_dist_df)) +
  geom_boxplot() +
  #scale_fill_viridis_d(name=TeX('$\\nu$')) +
  scale_y_log10() + 
  scale_x_continuous(breaks = seq(0,5,0.2)) +
  ylab(TeX('Meeting time $\\tau$')) +
  xlab(TeX('Degree of freedom $\\nu$')) +
  theme_classic(base_size = 12) +
  theme(legend.position = 'right')
meetings_by_t_dist_df
#ggsave(filename = "examples/half_t_degree_of_freedom/meetings_by_t_dist_df.pdf", plot = meetings_by_t_dist_df, width = 4, height = 3)

# load(file = "examples/half_t_degree_of_freedom/stat_performance_t_dist_prior.RData")
# Beta MSE
beta_mse_by_t_dist_df <- ggplot(stat_perform_df %>% dplyr::filter(t_dist_df <= 2),
                   aes(x=as.factor(t_dist_df), y=BetaMSE, group=t_dist_df)) +
  geom_boxplot() + xlab(TeX('Degree of freedom $\\nu$')) +
  ylab(TeX('Mean Squared Error')) +
  scale_y_continuous(breaks = seq(0,1,0.005)) +
  theme_classic(base_size = 12)
beta_mse_by_t_dist_df
#ggsave(filename = "examples/half_t_degree_of_freedom/beta_mse_by_t_dist_df.pdf", plot = beta_mse_by_t_dist_df, width = 4, height = 3)



