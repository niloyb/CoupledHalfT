# Plot of marginal density
rm(list = ls())
set.seed(1)

library(mvtnorm)
library(dplyr)
library(ggplot2)
library(latex2exp)
library(gridExtra)

# Numerical approximation of confluent hypergeometric function of the second kind
confHyp <- function(a, b, x){
  if (x==0) return(+Inf)
  integrand_ <- function(t) t^(a-1) * (1+t)^(b-a-1) * exp(-x*t)
  return(integrate(f = integrand_, lower = 0, upper = +Inf)$val / gamma(a))
}

# Marginal Beta posterior pdf
# xi, sigma2, X, y, t_dist_df globally defined
betaMarginalPostPdf <- function(beta, log=FALSE){
  m <- beta^2*xi/(2*sigma2)
  prior_density <- sapply(m, function(x)confHyp((t_dist_df+1)/1,1,x))
  #l <- mvtnorm::dmvnorm(as.vector(y/sqrt(sigma2)), mean = (X%*%(beta))/sqrt(sigma2))
  l <- dnorm(as.vector(y-X%*%(beta))/sqrt(sigma2))
  log_pdf <- sum(log(prior_density))+sum(log(l))
  if (log) return(log_pdf)
  return(exp(log_pdf))
}


# Generating data
X <- as.matrix(cbind(c(1,1),c(1,0), c(0,1)))
X_transpose <- t(X)
n <- dim(X)[1]
n <- dim(X)[2]
error_std <- 0.5
error_terms = error_std*rnorm(n, mean = 0, sd = 1)
true_beta <- c(1,0,0)
y = X%*%true_beta # + error_terms

# Sigma2, xi, t_dist_df fixed
sigma2 <- error_std^2
xi <- sigma2
t_dist_df <- 2

# Exclude zero from grid
beta1 <- seq(-1.5, 1.5, by=0.016)
beta2 <- seq(-1.5, 1.5, by=0.04)
beta3 <- seq(-1.5, 1.5, by=0.04)

beta_pdf_df <- data.frame(expand.grid(beta1, beta2, beta3))
beta_pdf_df <- beta_pdf_df %>%
  dplyr::group_by(Var1,  Var2, Var3) %>%
  dplyr::mutate(pdf = betaMarginalPostPdf(c(Var1,  Var2, Var3)))

#save(beta_pdf_df, file = '/Users/niloybiswas/Dropbox/horseshoe_coupling/Drafts/images/beta_marginal_plot.RData')
#load('/Users/niloybiswas/Dropbox/horseshoe_coupling/Drafts/images/beta_marginal_plot.RData')

beta_pdf_df12 <- beta_pdf_df %>% 
  dplyr::group_by(Var1, Var2) %>%
  dplyr::summarise(pdf = sum(pdf)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(pdf = pdf/sum(pdf))
beta_pdf_df13 <- beta_pdf_df %>% 
  dplyr::group_by(Var1, Var3) %>%
  dplyr::summarise(pdf = sum(pdf)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(pdf = pdf/sum(pdf))
beta_pdf_df23 <- beta_pdf_df %>% 
  dplyr::group_by(Var2, Var3) %>%
  dplyr::summarise(pdf = sum(pdf)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(pdf = pdf/sum(pdf))

# Generating ggplot plots

beta_12 <- ggplot(beta_pdf_df12, aes(Var1, Var2, z = pdf)) +
  stat_contour(geom="polygon", aes(fill=-..level..)) + theme_classic(base_size = 18) + 
  xlab(TeX('$\\beta_1$')) + ylab(TeX('$\\beta_2$')) +
  scale_fill_gradientn(colours = gray.colors(20)) +
  theme(legend.position = "none") +
  coord_cartesian(xlim=c(-0.5,1.5),
                  ylim=c(-0.5,1.25))
beta_13 <- ggplot(beta_pdf_df13, aes(Var1, Var3, z = pdf)) +
  stat_contour(geom="polygon", aes(fill=-..level..)) + theme_classic(base_size = 18) + 
  xlab(TeX('$\\beta_1$')) + ylab(TeX('$\\beta_3$')) +
  scale_fill_gradientn(colours = gray.colors(20)) +
  theme(legend.position = "none") +
  coord_cartesian(xlim=c(-0.5,1.5),
                  ylim=c(-0.5,1.25))
beta_23 <- ggplot(beta_pdf_df23, aes(Var2, Var3, z = pdf)) +
  stat_contour(geom="polygon", aes(fill=-..level..)) + theme_classic(base_size = 18) + 
  xlab(TeX('$\\beta_2$')) + ylab(TeX('$\\beta_3$')) +
  scale_fill_gradientn(colours = gray.colors(20)) +
  theme(legend.position = "none") +
  coord_cartesian(xlim=c(-0.5,1.5),
                  ylim=c(-0.5,1.25))
beta_plot_combined <- grid.arrange(beta_12, beta_13, beta_23, ncol=3)

# ggsave(filename = "examples/beta_marginal_density_plot/beta_marginal_plot.pdf", plot = beta_plot_combined, width = 10, height = 3)
