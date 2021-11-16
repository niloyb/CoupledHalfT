
library(dplyr)
library(ggplot2)

ksd_data <- read.csv('/Users/niloybiswas/Dropbox/horseshoe_coupling/Code/stein_discrepancy/ksd_half_t_prior.csv')
ksd_df <- data.frame('n'=ksd_data[,1], ksd=ksd_data[,2])
ksdplot1 <- ggplot(ksd_df) + 
  geom_line(aes(x=n, y=ksd)) + xlab('Number of samples') + 
  ylab('Kernel stein discrepancy') +
  theme_classic(base_size = 12) + scale_y_log10()
ksdplot1
# ggsave(filename = '/Users/niloybiswas/Dropbox/horseshoe_coupling/Drafts/images/stein_discrepancy/ksd_over_time_plot.pdf',ksdplot1, width = 4, height = 2.5)




std_normal_ksd_by_dimension <- read.csv('/Users/niloybiswas/Dropbox/horseshoe_coupling/Code/stein_discrepancy/std_normal_ksd_by_dimension.csv')
std_normal_ksd_df <- data.frame('dimension'=std_normal_ksd_by_dimension[,1], ksd=std_normal_ksd_by_dimension[,2])
ksdplot2 <- ggplot(std_normal_ksd_df) + 
  geom_line(aes(x=dimension, y=ksd)) + xlab('Dimension') + 
  ylab('Kernel stein discrepancy') + ylim(c(0,1.5)) +
  theme_classic(base_size = 12) 
ksdplot2
# ggsave(filename = '/Users/niloybiswas/Dropbox/horseshoe_coupling/Drafts/images/stein_discrepancy/std_normal_ksd_by_dimension.png', ksdplot2, width = 4, height = 2.5)










