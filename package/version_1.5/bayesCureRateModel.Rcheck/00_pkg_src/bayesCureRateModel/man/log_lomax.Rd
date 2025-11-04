\name{log_lomax}
\alias{log_lomax}
\title{
PDF and CDF of the Lomax distribution
}
\description{
The Lomax distribution as evaluated at the \pkg{VGAM} package.
}
\usage{
log_lomax(y, a1, a2, c_under = 1e-09)
}
\arguments{
  \item{y}{
	observed data
}
  \item{a1}{
	scale parameter
}
  \item{a2}{
	shape parameter
}
  \item{c_under}{
	A small positive value corresponding to the underflow threshold, e.g. c_under = 1e-9.
}
}
\details{
   The Lomax distribution is a special case of the 4-parameter
     generalized beta II distribution.}
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
\code{\link[VGAM]{dlomax}}
}
\examples{
log_lomax(y = 1:10, a1 = 1, a2 = 1, c_under = 1e-9)
}

