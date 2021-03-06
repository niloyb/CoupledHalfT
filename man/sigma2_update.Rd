% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mcmc_half_t.R
\name{sigma2_update}
\alias{sigma2_update}
\title{sigma2_update}
\usage{
sigma2_update(xi, eta, n, ssr, a0, b0, sigma2_fixed_value = NULL)
}
\arguments{
\item{xi}{positive scalar}

\item{eta}{length p vector}

\item{n}{positive integer}

\item{ssr}{positive scalar}

\item{a0}{positive scalar}

\item{b0}{positive scalar}

\item{sigma2_fixed_value}{NULL for sigma2 varying, positive scalar when fixed}
}
\value{
sigma^2 update step given eta and xi
}
\description{
sigma^2 update step given eta and xi
}
