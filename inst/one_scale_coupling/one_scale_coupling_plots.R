# Plot of meeting times for one scale coupling

# Libraries
library(dplyr)
library(ggplot2)
library(latex2exp)
library(gridExtra)
library(ggridges)

# Load Data
load("examples/one_scale_coupling/one_scale_coupling_meetings.RData")

# Plots
plot_vary_p <-
  ggplot(data = meetingtimes_df %>% dplyr::filter(n==100, sparsity==10, error_std==0.5),
         aes(x = meetingtime.meetingtime, y = p, group = p)) +
  geom_density_ridges(fill = "lightgrey") + 
  scale_x_log10() + 
  xlab(TeX('Meeting time $\\tau$')) +
  scale_y_continuous(breaks = seq(50,150,50)) +
  theme_classic(base_size = 18) + coord_flip()
plot_vary_p
# ggsave(filename = "examples/one_scale_coupling/one_scale_plot_vary_p.pdf", plot = plot_vary_p, width = 3, height = 3.5)

plot_vary_n <-
  ggplot(data = meetingtimes_df %>% dplyr::filter(p==150, sparsity==10, error_std==0.5), 
         aes(x = meetingtime.meetingtime, y = n, group = n)) +
  geom_density_ridges(fill = "lightgrey") +
  scale_x_log10() + 
  xlab(TeX('Meeting time $\\tau$')) +
  theme_classic(base_size = 18) +
  coord_flip()
plot_vary_n
#ggsave(filename = "examples/one_scale_coupling/one_scale_plot_vary_n.pdf", plot = plot_vary_n, width = 3, height = 3.5)

plot_vary_sparsity <-
  ggplot(data = meetingtimes_df %>% dplyr::filter(p==150,n==100, error_std==0.5), 
         aes(x = meetingtime.meetingtime, y = sparsity, group = sparsity)) +
  geom_density_ridges(fill = "lightgrey") +
  scale_x_log10() + 
  xlab(TeX('Meeting time $\\tau$')) +
  ylab(TeX('$s$')) +
  theme_classic(base_size = 18) +
  scale_y_continuous(breaks = seq(0,10,5)) +
  coord_flip()
plot_vary_sparsity
#ggsave(filename = "examples/one_scale_coupling/one_scale_plot_vary_sparsity.pdf", plot = plot_vary_sparsity, width = 3, height = 3.5)

plot_vary_error_std <-
  ggplot(data = meetingtimes_df %>% dplyr::filter(p==150,n==100,sparsity==10), 
         aes(x = meetingtime.meetingtime, y = error_std, group = error_std)) +
  geom_density_ridges(fill = "lightgrey") +
  scale_x_log10() + 
  xlab(TeX('Meeting time $\\tau$')) +
  ylab(TeX('$\\sigma_*$')) +
  theme_classic(base_size = 18) +
  scale_y_continuous(breaks = seq(0,0.5,0.25)) +
  coord_flip()
plot_vary_error_std
#ggsave(filename = "examples/one_scale_coupling/one_scale_plot_vary_error_std.pdf", plot = plot_vary_error_std, width = 3, height = 3.5)
