################################################################################
##################### Functions for unbiased estimation #####################
## Functions are from https://github.com/pierrejacob/unbiasedmcmc

#' sample_unbiasedestimator
#' @description Outputs unbiased estimates without storing the full chain
#' @param single_kernel kernel of marginal chain
#' @param coupled_kernel kernel of coupled chain
#' @param rinit initial distribution
#' @param h real-valued test function
#' @param k unbiased mcmc parameter corresponding to burn-in
#' @param m unbiased mcmc parameter corresponding to chain length
#' @param lag lag of the coupled chain
#' @param max_iterations maximum iterations
#' @return list containing mcmc estimator, bias correction, unbiased estimator, meetinging, number of iterations, time elapsed, cost
#' @export 
sample_unbiasedestimator <- 
  function(single_kernel, coupled_kernel, rinit, h = function(x) x, 
           k = 0, m = 1, lag = 1, max_iterations = Inf){
    if (k > m){
      print("error: k has to be less than m")
      return(NULL)
    }
    if (lag > m){
      print("error: lag has to be less than m")
      return(NULL)
    }
    starttime <- Sys.time()
    state1 <- rinit(); state2 <- rinit()
    # mcmcestimator computes the sum of h(X_t) for t=k,...,m
    mcmcestimator <- h(state1$chain_state)
    dimh <- length(mcmcestimator)
    if (k > 0){
      mcmcestimator <- rep(0, dimh)
    }
    # correction computes the sum of min(m-k+1, ceiling((t - k)/lag)) * (h(X_{t}) - h(Y_{t-lag})) for t=k+lag,..., tau - 1
    correction <- rep(0, dimh)
    time <- 0
    for (t in 1:lag){
      time <- time + 1
      state1 <- single_kernel(state1)
      if (time >= k){
        mcmcestimator <- mcmcestimator + h(state1$chain_state)
      }
    }
    if (time >= k + lag){
      correction <- correction + min(m-k+1, ceiling((time - k)/lag)) * (h(state1$chain_state) - h(state2$chain_state))
    }
    meetingtime <- Inf
    # time here is equal to lag; at this point we have X_lag,Y_0 and we are going to generate successively X_{t},Y_{t-lag} where time t is >= lag+1
    while ((time < max(meetingtime, m)) && (time < max_iterations)){
      time <- time + 1 # time is lag+1,lag+2,...
      if (is.finite(meetingtime)){
        state1 <- single_kernel(state1)
        state2 <- state1
        if (k <= time && time <= m){
          mcmcestimator <- mcmcestimator + h(state1$chain_state)
        }
      } else {
        res_coupled_kernel <- coupled_kernel(state1, state2)
        state1 <- res_coupled_kernel$state1
        state2 <- res_coupled_kernel$state2
        if (res_coupled_kernel$identical){
          meetingtime <- time
        }
        if (k <= time && time <= m){
          mcmcestimator <- mcmcestimator + h(state1$chain_state)
        }
        if (time >= k + lag){
          correction <- correction + min(m-k+1, ceiling((time - k)/lag)) * (h(state1$chain_state) - h(state2$chain_state))
        }
      }
    }
    uestimator <- mcmcestimator + correction
    cost <- lag + 2*(meetingtime - lag) + max(0, time - meetingtime)
    currentime <- Sys.time()
    elapsedtime <- as.numeric(lubridate::as.duration(lubridate::ymd_hms(currentime) - lubridate::ymd_hms(starttime)), "seconds")
    return(list(mcmcestimator = mcmcestimator / (m - k + 1), correction = correction / (m - k + 1), uestimator = uestimator / (m - k + 1),
                meetingtime = meetingtime, iteration = time, elapsedtime = elapsedtime, cost = cost))
  }

#' beta_unbiased_histogram
#' @description Outputs unbiased estimates for a histogram of a component
#' @param nchains no of chains
#' @param component index of component
#' @param breaks histogram breaks
#' @param nclass histogram number of classes
#' @return list containing histogram data with mids = mids, breaks = breaks, proportions = prop_mean, sd = prop_sd, width = width
#' @export 
beta_unbiased_histogram <- function(nchains,component,breaks=NULL,nclass=30){
  if(is.null(breaks)){
    # MCMC estimate
    mcmc_beta <- foreach(i = c(1:nchains), .combine = rbind)%dopar%{
      output <- half_t_mcmc(chain_length=m, burnin=k, X, X_transpose, y, t_dist_df=t_dist_df)
      return(output$beta_samples[,component, drop=FALSE])
    }
    breaks <- hist(mcmc_beta, plot = F, nclass = nclass)$breaks
  }
  mids <- (breaks[-1] + breaks[-length(breaks)])/2
  width <- diff(breaks)[1]
  
  breaks_function <- function(x){
    betas <- x$beta_samples[component, drop=FALSE]
    output <- (betas<breaks[-1])&(betas>=breaks[-length(breaks)])
    return(output)
  }
  
  unbiased_estimates <- foreach(i = c(1:nchains), .combine = rbind)%dopar%{
    output <- sample_unbiasedestimator(single_kernel, coupled_kernel, rinit, 
                                       h = breaks_function, 
                                       k = k, m = m, lag = lag)
    return(output$uestimator)
  }
  
  prop_mean <- apply(unbiased_estimates, 2, mean)
  prop_sd <- apply(unbiased_estimates, 2, sd)/sqrt(nchains)
  return(list(mids = mids, breaks = breaks, proportions = prop_mean, 
              sd = prop_sd, width = width))
}

#' beta_unbiased_histogram
#' @description Plots histogram given histogram data
#' @param with_bar boolean to include/exclude error bars
#' @return histrogram plot
#' @export 
plot_histogram <- 
  function (histogram, with_bar = TRUE) {
    df_ <- data.frame(xmin = histogram$mids - histogram$width/2, 
                      xmax = histogram$mids + histogram$width/2, x = histogram$mids, 
                      y = histogram$proportions/histogram$width, ymin = (histogram$proportions - 2 * histogram$sd)/histogram$width, 
                      ymax = (histogram$proportions + 2 * histogram$sd)/histogram$width)
    g <- ggplot(df_, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)) + geom_rect(alpha = 0.5)
    if (with_bar) {
      g <- g + geom_segment(aes(x = x, xend = x, y = 0, yend = y)) + ylab("density")
    }
    return(g)
  }

