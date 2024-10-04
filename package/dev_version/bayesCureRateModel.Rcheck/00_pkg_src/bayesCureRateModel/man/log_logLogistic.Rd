\name{log_logLogistic}
\alias{log_logLogistic}
\title{
PDF and CDF of the log-Logistic distribution.
}
\description{
The log-Logistic distribution as evaluated at the \pkg{flexsurv} package.
}
\usage{
log_logLogistic(y, a1, a2, c_under = 1e-09)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
	observed data
}
  \item{a1}{
	shape parameter
}
  \item{a2}{
	scale parameter
}
  \item{c_under}{
	A small positive value corresponding to the underflow threshold, e.g. c_under = 1e-9.
}
}
\details{
The log-logistic distribution is the probability distribution of a random variable whose logarithm has a logistic distribution. 
}
\value{
A list containing the following entries
\item{log_f}{natural logarithm of the pdf, evaluated at each datapoint.}
\item{log_F }{natural logarithm of the CDF, evaluated at each datapoint.}
}
\references{
Christopher Jackson (2016). flexsurv: A Platform for Parametric
  Survival Modeling in R. Journal of Statistical Software, 70(8), 1-33.
  doi:10.18637/jss.v070.i08
}
\author{
Panagiotis Papastamoulis
}
\seealso{
\code{\link[flexsurv]{dllogis}}
}
\examples{
log_logLogistic(y = 1:10, a1 = 1, a2 = 1, c_under = 1e-9)
}

