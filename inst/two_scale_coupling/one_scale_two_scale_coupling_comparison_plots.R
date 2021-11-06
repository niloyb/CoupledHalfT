# One scale v two scale comparison plot

# Libraries
rm(list = ls())
set.seed(1)
library(dplyr)
library(ggplot2)
library(latex2exp)
library(gridExtra)
library(ggridges)

# Load data
load("two_scale_coupling_plot.RData") # Two scale meetings
meetingtimes_df_two_scale <- meetingtimes_df
load("one_scale_coupling_meetings.RData") # One scale meetings
meetingtimes_df_one_scale <- meetingtimes_df

meetingtimes_df_combined <- 
  meetingtimes_df_one_scale %>% 
  dplyr::mutate(scale='one') %>%
  dplyr::bind_rows(meetingtimes_df_two_scale %>% dplyr::mutate(scale='two'))

plot_vary_p <-
  ggplot(data = meetingtimes_df_combined %>% 
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
  theme_classic(base_size = 18) + coord_flip() +
  theme(legend.position='bottom', legend.justification='left',
        legend.direction='horizontal',
        legend.title = element_text(size=13),
        legend.text = element_text(size=13))
plot_vary_p
#ggsave(filename = "examples/two_scale_coupling/coupling_comparison_plot_vary_p.pdf", plot = plot_vary_p, width = 4, height = 3) # base_size = 18

# plot_vary_p_highlight <-
#   ggplot(data = meetingtimes_df_combined %>%
#            dplyr::filter(n==100, sparsity==10, p<=150, error_std==0.5),
#          aes(x=p, y=meetingtime.meetingtime, group=interaction(p, scale))) +
#   geom_boxplot(aes(fill = scale)) +
#   scale_y_log10() +
#   scale_fill_manual(name="Scale",
#                     breaks = c("one", "two"),
#                     labels = c("One", "Two"),
#                     values = c("lightgrey", "white"),
#                     guide = "legend") +
#   ylab(TeX('Meeting time $\\tau$')) +
#   xlab(TeX('dimension $p$')) +
#   theme_classic(base_size = 18) +
#   theme(legend.position='bottom', legend.justification='left',
#         legend.direction='horizontal',
#         legend.title = element_text(size=13),
#         legend.text = element_text(size=13))
# plot_vary_p_highlight
# #ggsave(filename = "examples/two_scale_coupling/coupling_comparison_plot_vary_p_highlight.pdf", plot = plot_vary_p_highlight, width = 4.25, height = 3) # base_size = 12
# #ggsave(filename = "/Users/niloybiswas/Dropbox/horseshoe_coupling/Drafts/jasa_submission/images/two_scale_coupling_plot/coupling_comparison_plot_vary_p_highlight.pdf", plot = plot_vary_p_highlight, width = 3, height = 4) # base_size = 12
# 
# plot_vary_n <-
#   ggplot(data = meetingtimes_df_combined %>% 
#            dplyr::filter(p==150, n>=100, sparsity==10, error_std==0.5),
#          aes(y = n)) +
#   geom_density_ridges(alpha=0.8, aes(x = meetingtime.meetingtime, 
#                                      group = interaction(n, scale),
#                                      fill = scale)) +
#   scale_fill_manual(name="Scale",
#                     breaks = c("one", "two"),
#                     labels = c("One", "Two"),
#                     values = c("lightgrey", "white"),
#                     guide = "legend") +
#   scale_x_log10() + 
#   xlab(TeX('Meeting time $\\tau$')) +
#   ylab(TeX('observations $n$')) +
#   scale_y_continuous(breaks = seq(100,200,50)) +
#   theme_classic(base_size = 18) + coord_flip() +
#   theme(legend.position='bottom', legend.justification='left',
#         legend.direction='horizontal',
#         legend.title = element_text(size=13),
#         legend.text = element_text(size=13))
# plot_vary_n
# #ggsave(filename = "examples/two_scale_coupling/coupling_comparison_plot_vary_n.pdf", plot = plot_vary_n, width = 4, height = 3)
# 
# plot_vary_sparsity <-
#   ggplot(data = meetingtimes_df_combined %>% 
#            dplyr::filter(p==150,n==100, sparsity<=10, error_std==0.5),
#          aes(y = sparsity)) +
#   geom_density_ridges(alpha=0.8, aes(x = meetingtime.meetingtime, 
#                                      group = interaction(sparsity, scale),
#                                      fill = scale)) +
#   scale_fill_manual(name="Scale",
#                     breaks = c("one", "two"),
#                     labels = c("One", "Two"),
#                     values = c("lightgrey", "white"),
#                     guide = "legend") +
#   scale_x_log10() + 
#   xlab(TeX('Meeting time $\\tau$')) +
#   ylab(TeX('sparsity $s$')) +
#   scale_y_continuous(breaks = seq(0,10,4)) +
#   theme_classic(base_size = 18) + coord_flip() +
#   theme(legend.position='bottom', legend.justification='left',
#         legend.direction='horizontal',
#         legend.title = element_text(size=13),
#         legend.text = element_text(size=13))
# plot_vary_sparsity
# #ggsave(filename = "examples/two_scale_coupling/coupling_comparison_plot_vary_sparsity.pdf", plot = plot_vary_sparsity, width = 4, height = 3)
# 
# plot_vary_error_std <-
#   ggplot(data = meetingtimes_df_combined %>% 
#            dplyr::filter(p==150,n==100,sparsity==10, error_std<=0.5),
#          aes(y = error_std)) +
#   geom_density_ridges(alpha=0.8, aes(x = meetingtime.meetingtime, 
#                                      group = interaction(error_std, scale),
#                                      fill = scale)) +
#   scale_fill_manual(name="Scale",
#                     breaks = c("one", "two"),
#                     labels = c("One", "Two"),
#                     values = c("lightgrey", "white"),
#                     guide = "legend") +
#   scale_x_log10() + 
#   xlab(TeX('Meeting time $\\tau$')) +
#   ylab(TeX('error std. dev. $\\sigma_*$')) +
#   scale_y_continuous(breaks = seq(0,0.5,0.2)) +
#   #scale_y_log10() +
#   theme_classic(base_size = 18) + coord_flip() +
#   theme(legend.position='bottom', legend.justification='left',
#         legend.direction='horizontal',
#         legend.title = element_text(size=13),
#         legend.text = element_text(size=13))
# plot_vary_error_std
# #ggsave(filename = "examples/two_scale_coupling/coupling_comparison_plot_vary_error_std.pdf", plot = plot_vary_error_std, width = 4, height = 3)

