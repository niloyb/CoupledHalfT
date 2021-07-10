# GWAS plots

rm(list = ls())
library(dplyr)
library(tidyr)
library(ggplot2)
library(latex2exp)
library(ggridges)

# Loading coupling times
gwas_meet_df <- read.csv("GWASLlag.csv")

#Plotting L-Lag TV bounds
gwas_llag <- gwas_meet_df %>% dplyr::filter(lag==750)
tv_upper_bound_estimates <- function(coupling_times, L, t){return(pmax(0,ceiling((coupling_times-L-t)/L)))}
time <- c(0:1000)
tv_ub <- colMeans(sapply(time, function(x) tv_upper_bound_estimates(coupling_times = gwas_llag$meetingtime, L=gwas_llag$lag, x)))
tv_ub_df <- data.frame(t=time, t_dist_df=unique(gwas_llag$t_dist_df), tv_ub=tv_ub)

tv_ub_plot_maize <- 
  ggplot(tv_ub_df, aes(x=t, y=tv_ub)) +
  geom_line() + # scale_x_continuous(trans='log10') +
  # geom_errorbar(aes(ymin=tv_ub-tv_ub_sd/sqrt(indep_runs), ymax=tv_ub+tv_ub_sd/sqrt(indep_runs)), width=.2, position=position_dodge(.9)) +
  ylab(TeX('Total variation distance')) + xlab(TeX('iteration')) +
  #scale_linetype_manual(name='Dataset', values = c('solid','dashed')) +
  scale_y_continuous(breaks = seq(0,1,0.25)) +
  guides(linetype=guide_legend(nrow=2,byrow=TRUE)) +
  theme_classic(base_size = 14) + 
  theme(legend.position = 'right', legend.key.width=unit(1,"cm"))
tv_ub_plot_maize
# ggsave(filename = "images/dataset_examples/gwas_LlagTVUB_maize.pdf", plot = tv_ub_plot_maize, width = 4, height = 2.75)


