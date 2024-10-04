\name{log_gamma_mixture}
\alias{log_gamma_mixture}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
PDF and CDF of a Gamma mixture distribution
}
\description{
Computes the logarithm of the probability density function and cumulative density function per observation for each observation under a Gamma mixture model.
}
\usage{
log_gamma_mixture(y, a1, a2, p, c_under = 1e-09)
}
\arguments{
  \item{y}{
observed data
}
  \item{a1}{
vector containing the shape parameters of each Gamma mixture component
}
  \item{a2}{
vector containing the rate parameters of each Gamma mixture component
}
  \item{p}{
vector of mixing proportions
}
  \item{c_under}{
threshold for underflows.
}
}
\value{
A list containing the following entries
\item{log_f}{natural logarithm of the pdf, evaluated at each datapoint.}
\item{log_F }{natural logarithm of the CDF, evaluated at each datapoint.}
}
\author{
Panagiotis Papastamoulis
}
\examples{
y <- runif(10)
a1 <- c(1,2)
a2 <- c(1,1)
p <- c(0.9,0.1)
log_gamma_mixture(y, a1, a2, p)
}

