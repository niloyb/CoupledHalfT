#'@rdname get_max_coupling
#'@title Sample from maximally coupled distributions p and q
#'@description Takes two distributions (sampling function and pdf),
#' and return a function to sample from a maximal coupling of these distributions
#'@export
get_max_coupling <- function(rp, dp, rq, dq){
  function(){
    x <- rp()
    logu <- log(runif(1))
    if (dp(x) + logu < dq(x)){
      return(c(x,x))
    } else {
      reject <- TRUE
      y <- NA
      while (reject){
        y <- rq()
        logu <- log(runif(1))
        reject <- (dq(y) + logu < dp(y))
      }
      return(c(x,y))
    }
  }
}

### sample from maximal coupling of truncated exponentials (truncated on [0,trunc])
### with rate rate1, rate2, and truncations trunc1, trunc2
rtruncexp_maxcoupling <- function(rate1, trunc1, rate2, trunc2, eta_lower_bound){
  rp <- function() gen_truncated_exp(rate1, trunc1, eta_lower_bound)
  rq <- function() gen_truncated_exp(rate2, trunc2, eta_lower_bound)
  dp <- function(x){
    if (x > trunc1) return(-Inf) 
    else return(dexp((x-eta_lower_bound), rate1, log = TRUE) - pexp((trunc1-eta_lower_bound), rate1, log.p = TRUE))
  }
  dq <- function(x){
    if (x > trunc2) return(-Inf)
    else return(dexp((x-eta_lower_bound), rate2, log = TRUE) - pexp((trunc2-eta_lower_bound), rate2, log.p = TRUE))
  }
  f <- get_max_coupling(rp, dp, rq, dq)
  return(f())
}

### Sample from maximal coupling of two uniform distribution on [0, s] and [0,t]
### For s<t, sample U~Unif(0,s) and set X=sU. Set Y=sU w.p. s/t, Y=(t-s)U+s w.p. 1-s/t
runif_maxcoupling <- function(s, t){
  trun_min <- exp((log(s)+log(t) - abs(log(s)-log(t)))/2) # calculate min(s,t)
  trun_min[(s+t)==0] <- 0
  trun_max <- exp((log(s)+log(t) + abs(log(s)-log(t)))/2) # calculate max(s,t)
  p <- length(s)
  unif_draw <- runif(p)
  
  x <- unif_draw*trun_min
  y <- x
  coins <- runif(p)
  y[coins>trun_min/trun_max] <-
    (trun_max*unif_draw-x+trun_min)[coins>trun_min/trun_max]

  # Swapping back for all components
  swap_bool <- s>t
  x[swap_bool] <- x[swap_bool] + y[swap_bool]
  y[swap_bool] <- x[swap_bool] - y[swap_bool]
  x[swap_bool] <- x[swap_bool] - y[swap_bool]
  return(cbind(x,y))
}
