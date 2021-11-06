rm(list = ls())
set.seed(1)

library(dplyr)
library(ggplot2)
library(latex2exp)
library(gridExtra)
library(ggridges)


# One-scale coupling vs. Two-scale coupling plots
load("../BackUpFiles/CoupledHalfT_old/examples/truncated_prior_xi/one_two_scale_coupling_meetings.RData") # meetingtimes_df_one_scale
meetingtimes_df

meetingtimes_df$epsilon_eta[meetingtimes_df$epsilon_eta==0.5] <- 'two'
meetingtimes_df$epsilon_eta[meetingtimes_df$epsilon_eta==Inf] <- 'one'
colnames(meetingtimes_df)[colnames(meetingtimes_df)=='epsilon_eta'] <- 'scale'

plot_vary_p <-
  ggplot(data = meetingtimes_df %>% 
           dplyr::filter(n==100, sparsity==10, p<=150, error_std==0.5),
         aes(y = p)) +
  geom_density_ridges(alpha=0.8, aes(x = meetingtime.meetingtime, group = interaction(p, scale),
                                     fill = scale)) +
  scale_fill_manual(name="Scale",
                    breaks = c("one", "two"),
                    labels = c("One", "Two"),
                    values = c("lightgrey", "white"),
                    guide = "legend") +
  scale_x_log10() + 
  xlab(TeX('Meeting time $\\tau$')) +
  ylab(TeX('dimension $p$')) +
  scale_y_continuous(breaks = seq(50,200,50)) +
  theme_classic(base_size = 16) + coord_flip() +
  theme(legend.justification='right',
        legend.direction='vertical',
        legend.title = element_text(size=13),
        legend.text = element_text(size=13))
plot_vary_p
# ggsave(filename = "../BackUpFiles/CoupledHalfT_old/examples/truncated_prior_xi/coupling_comparison_plot_vary_p.pdf", plot = plot_vary_p, width = 4, height = 3) # base_size = 18


# Meeting time plots
load("../BackUpFiles/CoupledHalfT_old/examples/truncated_prior_xi/degree_of_freedom_meeting_times.RData") # meetingtimes_df_two_scale
plot_vary_p <-
  ggplot(data = tdist_meetingtimes_df %>% dplyr::filter(p <=1000),
         aes(y = p)) +
  geom_density_ridges(alpha=0.8, scale = 0.8,
                      aes(x = meetingtime.meetingtime, 
                          group = interaction(p, t_dist_df), 
                          fill = as.factor(t_dist_df))) +
  scale_fill_viridis_d(name=TeX('$\\nu$')) +
  scale_x_log10() + 
  scale_y_continuous(breaks = seq(0,1e3,500)) +
  xlab(TeX('Meeting time $\\tau$')) +
  ylab(TeX('Dimension $p$')) +
  theme_classic(base_size = 16) + coord_flip() +
  theme(legend.position = 'right')
plot_vary_p
# ggsave(filename = "../BackUpFiles/CoupledHalfT_old/examples/truncated_prior_xi/half_t_degree_of_freedom_vary_p.pdf", plot = plot_vary_p, width = 4, height = 3)
# ggsave(filename = "/Users/niloybiswas/Dropbox/horseshoe_coupling/Drafts/images/truncated_prior_xi/half_t_degree_of_freedom_vary_p.pdf", plot = plot_vary_p, width = 4, height = 3)

# Statistical performance plots
load("../BackUpFiles/CoupledHalfT_old/examples/truncated_prior_xi/stat_performance_t_dist_prior.RData") # meetingtimes_df_two_scale

# Beta MSE
beta_mse <- ggplot(stat_perform_df,
                   aes(x=as.factor(t_dist_df), y=BetaMSE, fill=as.factor(t_dist_df))) +
  geom_boxplot() + xlab(TeX('$\\nu$')) +
  ylab(TeX('$ |\\beta_* - \\hat{\\beta} |_2^2 / p $')) +
  scale_fill_manual(name=TeX('$\\nu$'),
                    values = c(rep("white",6), rep("lightgrey",5))) +
  # scale_fill_viridis_d(name=TeX('$\\nu$')) +
  theme_classic(base_size = 16) +
  theme(axis.text.x = element_text(size=12)) +
  scale_x_discrete(breaks = c(seq(1,2,0.2),seq(2,8,2), 'lasso'),
                   labels = c(seq(1,2,0.2),seq(2,8,2), 'Lasso')) +
  # scale_y_continuous(breaks = seq(0,1,0.02)) +
  theme(legend.position = "none")
beta_mse
# ggsave(filename = "../BackUpFiles/CoupledHalfT_old/examples/truncated_prior_xi/beta_mse.pdf", plot = beta_mse, width = 4, height = 4)
# ggsave(filename = "/Users/niloybiswas/Dropbox/horseshoe_coupling/Drafts/images/truncated_prior_xi/beta_mse2.pdf", plot = beta_mse, width = 4, height = 3)
