\name{log_dagum}
\alias{log_dagum}
\title{
PDF and CDF of the Dagum distribution
}
\description{
The Dagum distribution as evaluated at the \pkg{VGAM} package.
}
\usage{
log_dagum(y, a1, a2, a3, c_under = 1e-09)
}
\arguments{
  \item{y}{
	observed data
}
  \item{a1}{
	scale parameter
}
  \item{a2}{
	shape1.a parameter
}
\item{a3}{
	shape2.p parameter
}
  \item{c_under}{
A small positive value corresponding to the underflow threshold, e.g. c_under = 1e-9.
}
}
\details{
 The Dagum distribution is a special case of the 4-parameter
     generalized beta II distribution.
}
\value{
A list containing the following entries
\item{log_f}{natural logarithm of the pdf, evaluated at each datapoint.}
\item{log_F }{natural logarithm of the CDF, evaluated at each datapoint.}
}
\references{
Thomas W. Yee (2015). Vector Generalized Linear and Additive Models:
  With an Implementation in R. New York, USA: Springer.
}
\author{
Panagiotis Papastamoulis
}
\seealso{
\code{\link[VGAM]{ddagum}}
}
\examples{
log_dagum(y = 1:10, a1 = 1, a2 = 1, a3 = 1, c_under = 1e-9)
}

