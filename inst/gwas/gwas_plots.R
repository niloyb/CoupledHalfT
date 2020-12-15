# GWAS plots

rm(list = ls())
library(dplyr)
library(tidyr)
library(ggplot2)
library(latex2exp)
library(ggridges)

# Loading coupling times
setwd('/Users/niloybiswas/Google Drive/Niloy_Files/github/CoupledHalfT/')
gwas_meet_df <- read.csv("GWASLlag_new.csv")

#Plotting Coupling times
one_lag_meeting_plot <- 
  ggplot(gwas_meet_df %>% dplyr::filter(lag==1), 
                            aes(x=meetingtime)) + 
  geom_histogram(aes(y=..density..), color="black", fill="white") +
  geom_density(alpha=.2, color="black",size=1) +
  xlab(TeX('Meeting time $\\tau$')) +
  # scale_x_continuous(breaks = seq(0,1000,2000)) +
  theme_classic(base_size = 16)
one_lag_meeting_plot
# ggsave(filename = "/Users/niloybiswas/Dropbox/horseshoe_coupling/Drafts/images/dataset_examples/gwas_meetings.pdf", plot = one_lag_meeting_plot, width = 4, height = 3)


#Plotting L-Lag TV bounds
gwas_llag <- gwas_meet_df %>% dplyr::filter(lag==750)
tv_upper_bound_estimates <- function(coupling_times, L, t){return(pmax(0,ceiling((coupling_times-L-t)/L)))}
time <- c(0:1000)
tv_ub <- colMeans(sapply(time, function(x) tv_upper_bound_estimates(coupling_times = gwas_llag$meetingtime, L=gwas_llag$lag, x)))
tv_ub_df <- data.frame(t=time, t_dist_df=unique(gwas_llag$t_dist_df), tv_ub=tv_ub)
# TV UB plots
tv_ub_plot <- ggplot(tv_ub_df, aes(x=t, y=tv_ub)) +
  geom_line() + # scale_x_continuous(trans='log10') +
  ylab(TeX('Total variation distance')) +
  xlab(TeX('iteration')) +
  scale_y_continuous(breaks = seq(0,1,0.25)) +
  theme_classic(base_size = 16)
tv_ub_plot
# ggsave(filename = "/Users/niloybiswas/Dropbox/horseshoe_coupling/Drafts/images/dataset_examples/gwas_LlagTVUB.pdf", plot = tv_ub_plot, width = 4, height = 3) # base_size 16
# ggsave(filename = "/Users/niloybiswas/Dropbox/horseshoe_coupling/Drafts/images/dataset_examples/gwas_LlagTVUB_highlight.pdf", plot = tv_ub_plot, width = 4, height = 3) # base_size 12

# Single trajectory traceplots
load(file = "/Users/niloybiswas/Google Drive/Niloy_Files/Harvard/PhD/Research/Pierre Jacob Lab/CoupledHalfT/simulation_data/gwas_single_chain.RData")

beta_means <- colMeans(single_chain$beta_samples)

component_null <- 1 #
component_unimodal <- which(beta_means==min(beta_means)) #
component_bimodal <- which(beta_means==max(beta_means))  #

traceplots_df <- 
  data.frame(t=c(1:dim(single_chain$beta_samples)[1]),
             Null= single_chain$beta_samples[,component_null],
             Unimodal= single_chain$beta_samples[,component_unimodal],
             Bimodal= single_chain$beta_samples[,component_bimodal])
traceplots_df <- traceplots_df %>% tidyr::gather("Component", "Value", -t)
  
traceplots <- ggplot(traceplots_df, aes(y = Component)) +
  geom_density_ridges(alpha=0.8, scale = 1,
                      aes(x=Value, fill=Component)) +
  scale_fill_viridis_d(name="Component",
                       breaks=c("Unimodal", "Null", "Bimodal")) +
  scale_x_continuous(breaks=seq(-10,10,0.5), limits = c(-1,1)) +
  theme_classic(base_size = 16) +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
traceplots
#ggsave(filename = "/Users/niloybiswas/Dropbox/horseshoe_coupling/Drafts/images/dataset_examples/gwas_traceplots.pdf", plot = traceplots, width = 4, height = 3)

