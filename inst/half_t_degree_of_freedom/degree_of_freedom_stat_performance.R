### ## ## ## ## ## ## ## ## ## ## ## ## #
## Statistical Performance Simulations ##
### ## ## ## ## ## ## ## ## ## ## ## ## #

# Libraries
rm(list = ls())
set.seed(1)
library(CoupledHalfT)
library(dplyr)
library(ggplot2)
library(latex2exp)
library(gridExtra)
library(ggridges)

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
# Statistical Performance Functions
# MSE of Beta by Estimation of Time-Averages
BetaMSE <- function(true_beta, mc_beta_samples){
  return(mean((true_beta-colMeans(mc_beta_samples))^2))
}
# MSE of Beta by Estimation of Time-Averages
XBetaMSE <- function(true_beta, mc_beta_samples, X){
  return(mean((X%*%(true_beta-colMeans(mc_beta_samples)))^2))
}
ComponentwiseBetaMSE <- function(true_beta, mc_beta_samples){
  return((true_beta-colMeans(mc_beta_samples))^2)
}
# Covergae of credible intervals
BetaCoverageSymmetric <- function(true_beta, mc_beta_samples, qtl){
  componentwise_quantiles <- t(apply(mc_beta_samples, 2 , quantile , probs = c((1-qtl)/2,(1+qtl)/2)))
  coverage <- (componentwise_quantiles[,2]-true_beta)*(true_beta-componentwise_quantiles[,1])>0
  return(coverage)
}
# Covergae of credible intervals
BetaCoverageHD <- function(true_beta, mc_beta_samples, qtl){
  p <- length(true_beta)
  s <- max(which(true_beta!=0))
  coverage <- rep(NA,p)
  # Covering bimodalities w allowSplit
  for(i in 1:s){
    hdint <- HDInterval::hdi(density(mc_beta_samples[,i]), credMass = qtl, allowSplit=TRUE)
    coverage[i] <- any((hdint[,2]-true_beta[i])*(true_beta[i]-hdint[,1]) > 0)
  }
  componentwise_quantiles <- t(HDInterval::hdi(mc_beta_samples[,c((s+1):p)], credMass = qtl))
  coverage[(s+1):p] <- (componentwise_quantiles[,2]-true_beta[c((s+1):p)])*(true_beta[c((s+1):p)]-componentwise_quantiles[,1])>0
  return(coverage)
}


## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
# Run simulations

## Generate data
n <- 100
p <- 500
s <- 20
true_beta <- matrix(0,p,1)
true_beta[1:s] = 2^(-(seq(s)/4-9/4))

# X <- matrix(rnorm(n*p, mean = 0, sd = 1), nrow = n, ncol = p)
# #Error terms
# error_std <- 2
# error_terms = error_std*rnorm(n, mean = 0, sd = 1)
# y = X%*%true_beta + error_terms
# X_transpose <- t(X)

stat_perform_df <- data.frame()
stat_perform_componentwise_df <- data.frame()
traceplots_df <- data.frame()
iterations <- 100
for (i in 1:iterations){
  # Generate data
  X <- matrix(rnorm(n*p, mean = 0, sd = 1), nrow = n, ncol = p)
  #Error terms
  error_std <- 2
  error_terms = error_std*rnorm(n, mean = 0, sd = 1)
  y = X%*%true_beta + error_terms
  X_transpose <- t(X)
  
  # Lasso Estimates
  cv.lasso <- glmnet::cv.glmnet(X, y, alpha = 1, family = 'gaussian', intercept = FALSE)
  lasso.estimates <- coef(cv.lasso, cv.lasso$lambda.min)[2:(p+1)]
  lasso.estimates <- data.frame(matrix(lasso.estimates, nrow = 1))
  colnames(lasso.estimates) <- paste("X", 1:p, sep = '')
  
  chain_df <- lasso.estimates[,c(1:(s+10))] %>%
    tidyr::gather(component, value, X1:X30) %>%
    dplyr::mutate(n = n, p =p, t_dist_df = 'lasso')
  traceplots_df <- rbind(traceplots_df, chain_df)
  stat_perform_df <-
    rbind(stat_perform_df,
          data.frame(n = n, p =p, t_dist_df = 'lasso',
                     BetaMSE=BetaMSE(true_beta, lasso.estimates),
                     XBetaMSE=XBetaMSE(true_beta, lasso.estimates, X),
                     BetaHDCoverage0.95=NA))
  stat_perform_componentwise_df <-
    rbind(stat_perform_componentwise_df,
          data.frame(n=n, p=p, t_dist_df = 'lasso', component = c(1:p),
                     BetaHDCoverage0.95 = NA,
                     ComponentwiseBetaMSE =
                       ComponentwiseBetaMSE(true_beta, lasso.estimates)))
  for (t_dist_df in c(seq(1,2,0.2),4,6,8)){
      if (t_dist_df==1){burnin <- 600} else {burnin <- 300}
      mc_chain_size <- burnin + 1000

      chain <- half_t_mcmc(mc_chain_size, burnin=0, X, X_transpose,y, 
                           a0=1, b0=1, std_MH=0.8, t_dist_df=t_dist_df)
      mc_beta_samples <- chain$beta_samples[c(burnin:mc_chain_size),]

      chain_df <- data.frame(mc_beta_samples[,c(1:(s+10))]) %>%
        tidyr::gather(component, value, X1:X30) %>%
        dplyr::mutate(n = n, p =p, t_dist_df = t_dist_df)
      traceplots_df <- rbind(traceplots_df, chain_df)

      stat_perform_df <-
        rbind(stat_perform_df,
              data.frame(n = n, p =p, t_dist_df = t_dist_df,
                         BetaMSE=BetaMSE(true_beta, mc_beta_samples),
                         XBetaMSE=XBetaMSE(true_beta, mc_beta_samples, X),
                         BetaHDCoverage0.95=mean(BetaCoverageHD(true_beta, mc_beta_samples, 0.95))))
      stat_perform_componentwise_df <-
        rbind(stat_perform_componentwise_df,
              data.frame(n=n, p=p, t_dist_df = t_dist_df, component = c(1:p),
                         BetaHDCoverage0.95 =
                           BetaCoverageHD(true_beta, mc_beta_samples, 0.95),
                         ComponentwiseBetaMSE =
                           ComponentwiseBetaMSE(true_beta, mc_beta_samples)))
      print(c(p, t_dist_df, i))
    }
  }

traceplots_df <- traceplots_df %>% 
  dplyr::filter(component %in% c('X1','X9','X13','X25'))

# save(stat_perform_df, file = "inst/half_t_degree_of_freedom/stat_performance_t_dist_prior.RData")
#load(file = "inst/half_t_degree_of_freedom/stat_performance_t_dist_prior.RData")
#save(stat_perform_componentwise_df, file = "inst/half_t_degree_of_freedom/stat_perform_componentwise_df.RData")
#load(file = "inst/half_t_degree_of_freedom//stat_perform_componentwise_df.RData")

#save(traceplots_df, file = "inst/half_t_degree_of_freedom/traceplots_t_dist_prior.RData")
#load(file = "inst/half_t_degree_of_freedom/traceplots_t_dist_prior.RData")

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
## Plots ##
# Lasso Estimates
lasso.estimates <- traceplots_df %>% 
  dplyr::filter(t_dist_df == 'lasso') %>%
  dplyr::group_by(n, p, component) %>%
  dplyr::summarise(lasso.mean = mean(value))
# Traceplots
traceplot_1 <- ggplot(traceplots_df %>% dplyr::group_by(n, p, t_dist_df) %>%
                        dplyr::filter(component == 'X1',t_dist_df != 'lasso'),
                      aes(y = as.factor(t_dist_df))) +
  geom_density_ridges(alpha=0.8, scale = 0.8,
                      aes(x = value, group = as.factor(t_dist_df), 
                          fill = as.factor(t_dist_df))) +
  scale_fill_manual(name=TeX('$\\nu$'),
                     values = c(viridis::viridis(6), rep("lightgrey",4))) +
  xlab(TeX('$\\beta_{1}$ density')) +
  scale_y_discrete(limits = rev(levels(as.factor(traceplots_df$t_dist_df)))) +
  scale_x_continuous(breaks = seq(-10,10,2)) +
  geom_vline(xintercept=true_beta[1], linetype="solid", color = "black", size=1) +
  theme_classic(base_size = 24) +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
traceplot_1
#ggsave(filename = "examples/half_t_degree_of_freedom/traceplot_1.pdf", plot = traceplot_1, width = 4, height = 4.5)
traceplot_2 <- ggplot(traceplots_df %>%
                        dplyr::group_by(n, p, t_dist_df) %>%
                        dplyr::filter(component == 'X9',t_dist_df != 'lasso'),
                      aes(y = as.factor(t_dist_df))) +
  geom_density_ridges(alpha=0.8, scale = 0.8,
                      aes(x = value,
                          group = as.factor(t_dist_df),
                          fill = as.factor(t_dist_df))) +
  scale_fill_manual(name=TeX('$\\nu$'),
                    values = c(viridis::viridis(6), rep("lightgrey",4))) +
  # scale_fill_viridis_d(name=TeX('$\\nu$')) +
  xlab(TeX('$\\beta_{9}$ density')) +
  # scale_y_reverse() +
  scale_y_discrete(limits = rev(levels(as.factor(traceplots_df$t_dist_df)))) +
  scale_x_continuous(breaks = seq(-10,10,1)) +
  geom_vline(xintercept=true_beta[9], linetype="solid", color = "black", size=1) +
  theme_classic(base_size = 24) +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
traceplot_2
#ggsave(filename = "examples/half_t_degree_of_freedom/traceplot_2.pdf", plot = traceplot_2, width = 4, height = 4.5)
traceplot_3 <- ggplot(traceplots_df %>%
                        dplyr::group_by(n, p, t_dist_df) %>%
                        dplyr::filter(component == 'X13', t_dist_df != 'lasso'),
                      aes(y = as.factor(t_dist_df))) +
  geom_density_ridges(alpha=0.8, scale = 0.8,
                      aes(x = value,
                          group = as.factor(t_dist_df),
                          fill = as.factor(t_dist_df))) +
  # scale_fill_viridis_d(name=TeX('$\\nu$')) +
  scale_fill_manual(name=TeX('$\\nu$'),
                    values = c(viridis::viridis(6), rep("lightgrey",4))) +
  xlab(TeX('$\\beta_{13}$ density')) +
  # scale_y_reverse() +
  scale_y_discrete(limits = rev(levels(as.factor(traceplots_df$t_dist_df)))) +
  scale_x_continuous(breaks = seq(-10,10,0.5), limits = c(-0.75,0.75)) +
  geom_vline(xintercept=true_beta[13], linetype="solid", color = "black", size=1) +
  #annotate(geom="text", x=(true_beta[13]+0.4), y=2.05, label=TeX('$\\beta_{*,13}$'), color="blue", size=5) +
  theme_classic(base_size = 24) +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
traceplot_3
#ggsave(filename = "examples/half_t_degree_of_freedom/traceplot_3.pdf", plot = traceplot_3, width = 4, height = 4.5)
traceplot_4 <- ggplot(traceplots_df %>%
                        dplyr::group_by(n, p, t_dist_df) %>%
                        dplyr::filter(component == 'X25', t_dist_df != 'lasso'),
                      aes(y = as.factor(t_dist_df))) +
  geom_density_ridges(alpha=0.8, scale = 0.8,
                      aes(x = value,
                          group = as.factor(t_dist_df),
                          fill = as.factor(t_dist_df))) +
  scale_fill_manual(name=TeX('$\\nu$'),
                    values = c(viridis::viridis(6), rep("lightgrey",4))) +
  # scale_fill_viridis_d(name=TeX('$\\nu$')) +
  xlab(TeX('$\\beta_{25}$ density')) +
  # scale_y_reverse() +
  scale_y_discrete(limits = rev(levels(as.factor(traceplots_df$t_dist_df)))) +
  scale_x_continuous(breaks = seq(-10,10,0.25), limits = c(-0.25,0.25)) +
  geom_vline(xintercept=true_beta[25], linetype="solid", color = "black", size=1) +
  #annotate(geom="text", x=(true_beta[25]+0.4), y=2.05, label=TeX('$\\beta_{*,25}$'), color="blue", size=5) +
  theme_classic(base_size = 24) +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
traceplot_4
#ggsave(filename = "examples/half_t_degree_of_freedom/traceplot_4.pdf", plot = traceplot_4, width = 4, height = 4.5)


# Beta MSE
beta_mse <- ggplot(stat_perform_df,
                   aes(x=as.factor(t_dist_df), y=BetaMSE, fill=as.factor(t_dist_df))) +
  geom_boxplot() + xlab(TeX('$\\nu$')) +
  ylab(TeX('$ |\\beta_* - \\bar{\\beta}_{MCMC} |_2^2 / p $')) +
  scale_fill_manual(name=TeX('$\\nu$'),
                    values = c(viridis::viridis(6), rep("lightgrey",4))) +
  # scale_fill_viridis_d(name=TeX('$\\nu$')) +
  theme_classic(base_size = 24) +
  # scale_x_discrete(breaks = c(seq(1,2,0.2),10)) +
  scale_y_continuous(breaks = seq(0,1,0.02)) +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) + 
  # scale_x_reverse() + 
  scale_x_discrete(limits = rev(levels(as.factor(stat_perform_df$t_dist_df)))) +
  coord_flip()
beta_mse
#ggsave(filename = "examples/half_t_degree_of_freedom/beta_mse.pdf", plot = beta_mse, width = 4, height = 4)

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
# ggsave(filename = "examples/half_t_degree_of_freedom/beta_mse.pdf", plot = beta_mse, width = 4, height = 4)
# ggsave(filename = "/Users/niloybiswas/Dropbox/horseshoe_coupling/Drafts/images/half_t_degree_of_freedom_plot/beta_mse2.pdf", plot = beta_mse, width = 4, height = 3)
# ggsave(filename = "/Users/niloybiswas/Dropbox/horseshoe_coupling/New_plots_May_2022/beta_mse2.pdf", plot = beta_mse, width = 4, height = 3)

# XBeta MSE
# Xbeta_mse <- ggplot(stat_perform_df,
#                    aes(x=t_dist_df, y=XBetaMSE, fill=as.factor(t_dist_df))) +
#   geom_boxplot() + xlab(TeX('$\\nu$')) +
#   ylab(TeX('$ |X \\beta_* - X \\hat{\\beta}_{MCMC} |_2^2 / n $')) +
#   scale_fill_viridis_d(name=TeX('$\\nu$')) +
#   theme_classic(base_size = 24) +
#   theme(axis.title.y=element_blank(),
#         axis.text.y=element_blank(),
#         axis.ticks.y=element_blank()) + scale_x_reverse() + coord_flip()
# Xbeta_mse
#ggsave(filename = "examples/half_t_degree_of_freedom/Xbeta_mse.pdf", plot = Xbeta_mse, width = 4, height = 3)

# Coverage
Beta_HPDCoverage0.95 <-
  ggplot(stat_perform_df,
         aes(x=as.factor(t_dist_df), y=BetaHDCoverage0.95, fill=as.factor(t_dist_df))) +
  geom_boxplot() + xlab(TeX('$\\nu$')) +
  ylab(TeX('$95%$ HPD Int.')) +
  scale_fill_manual(name=TeX('$\\nu$'),
                    values = c(viridis::viridis(6), rep("lightgrey",4))) +
  # scale_fill_viridis_d(name=TeX('$\\nu$')) +
  scale_y_continuous(breaks = seq(0.975,1,0.02),
                     limits = c(0.97, 1),
                     labels = scales::percent_format(accuracy = 0.1)) +
  # scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
  theme_classic(base_size = 24) +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) + 
  # scale_x_reverse() + 
  scale_x_discrete(limits = rev(levels(as.factor(stat_perform_df$t_dist_df)))) +
  coord_flip()
Beta_HPDCoverage0.95
#ggsave(filename = "examples/half_t_degree_of_freedom/Beta_HPDCoverage.pdf", plot = Beta_HPDCoverage0.95, width = 4, height = 4)

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
# Componentwise Plots
stat_perform_componentwise_summary <-
  stat_perform_componentwise_df %>%
  dplyr::group_by(p, t_dist_df, component) %>%
  dplyr::summarise(meanBetaHDCoverage0.95 = mean(BetaHDCoverage0.95),
                   meanComponentwiseBetaMSE = mean(ComponentwiseBetaMSE))

# Componentwise Coverage plots
Beta_HPDCoverage0.95_componentwise <-
  ggplot(stat_perform_componentwise_summary %>% dplyr::filter(component <= (s+10)),
       aes(x=component, y=meanBetaHDCoverage0.95, col=as.factor(t_dist_df))) +
  geom_line(size=1) +
  scale_color_manual(name=TeX('$\\nu$'),
                    values = c(viridis::viridis(6), rep("lightgrey",4))) +
  # scale_color_viridis_d(name=TeX('$\\nu$')) +
  xlab(TeX('Component $j$')) +
  # ylab(TeX('$95%$ HPD Int.')) +
  ylab(TeX('')) +
  scale_y_continuous(breaks = seq(0,1,0.25), labels = scales::percent_format(accuracy = 1)) +
  #scale_x_continuous(breaks = seq(1,2,0.5)) +
  theme_classic(base_size = 24)
Beta_HPDCoverage0.95_componentwise
#ggsave(filename = "examples/half_t_degree_of_freedom/Beta_HPDCoverage_componentwise.pdf", plot = Beta_HPDCoverage0.95_componentwise, width = 4, height = 4)

# Componentwise MSE plots
Beta_mse_componentwise <-
  ggplot(stat_perform_componentwise_summary %>% dplyr::filter(component <= (s+10)),
         aes(x=component, y=meanComponentwiseBetaMSE, col=as.factor(t_dist_df))) +
  geom_line(size=1) +
  scale_color_manual(name=TeX('$\\nu$'),
                     values = c(viridis::viridis(6), rep("lightgrey",4))) +
  # scale_color_viridis_d(name=TeX('$\\nu$')) +
  xlab(TeX('Component $j$')) +
  # ylab(TeX('$(\\beta_{*,j} - \\hat{\\beta}_{MCMC,j})^2$')) +
  ylab(TeX('')) +
  #scale_x_continuous(breaks = seq(1,2,0.5)) +
  theme_classic(base_size = 24)
Beta_mse_componentwise
#ggsave(filename = "examples/half_t_degree_of_freedom/Beta_mse_componentwise.pdf", plot = Beta_mse_componentwise, width = 4, height = 4)


