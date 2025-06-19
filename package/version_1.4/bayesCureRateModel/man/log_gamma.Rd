\name{log_gamma}
\alias{log_gamma}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
PDF and CDF of the Gamma distribution
}
\description{
Computes the pdf and cdf of the Gamma distribution. 
}
\usage{
log_gamma(y, a1, a2, c_under = 1e-09)
}
\arguments{
  \item{y}{
	observed data
}
  \item{a1}{
	shape parameter
}
  \item{a2}{
	rate parameter
}
  \item{c_under}{
	A small positive value corresponding to the underflow threshold, e.g. c_under = 1e-9.
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
\seealso{
\code{\link[stats]{dgamma}}
}
\examples{
log_gamma(y = 1:10, a1 = 1, a2 = 1, c_under = 1e-9)
}

