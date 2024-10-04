\name{log_user_mixture}
\alias{log_user_mixture}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Define a finite mixture of a given family of distributions. 
}
\description{
This function computes the logarithm of the probability density function and cumulative density function per observation for each observation under a user-defined mixture of a given family of distributions. The parameters of the given family of distributions should belong to (0, inf). 
}
\usage{
log_user_mixture(user_f, y, a, p, c_under = 1e-09)
}
\arguments{
  \item{user_f}{
a user defined function that returns the logarithm of a given probability density and the corresponding logarithm of the cumulative distribution function. These arguments should be returned in the form of a list with two entries: \code{log_f} and \code{log_F}, containing the logarithm of the pdf and cdf values of \code{y}, respectively, for a given set of parameter values.
}
  \item{y}{
observed data
}
  \item{a}{a matrix where each column corresponds to component specific parameters
and the columns to different components. All parameters should be positive. The number of columns should be the same with the number of mixture components.
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
# We will define a mixture of 2 exponentials distributions.
# First we pass the exponential distribution at user_f
user_f <- function(y, a){
	log_f <- dexp(y, rate = a, log = TRUE)
	log_F <- pexp(y, rate = a, log.p = TRUE)
	result <- vector('list', length = 2)
	names(result) <- c('log_f', 'log_F')
	result[["log_f"]] = log_f
	result[["log_F"]] = log_F
	return(result)
}
#	simulate some date
y <- runif(10)
# Now compute the log of pdf and cdf for a mixture of K=2 exponentials
p <- c(0.9,0.1)
a <- matrix(c(0.1, 1.5), nrow = 1, ncol = 2)
log_user_mixture(user_f = user_f, y = y, a = a, p = p)
}

