\name{log_weibull}
\alias{log_weibull}
\title{
PDF and CDF of the Weibull distribution
}
\description{
Computes the log pdf and cdf of the weibull distribution. 
}
\usage{
log_weibull(y, a1, a2, c_under)
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
\examples{
log_weibull(y = 1:10, a1 = 1, a2 = 1, c_under = 1e-9)

}

