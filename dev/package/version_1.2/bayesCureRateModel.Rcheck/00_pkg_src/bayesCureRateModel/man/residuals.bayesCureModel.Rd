\name{residuals.bayesCureModel}
\alias{residuals.bayesCureModel}
\title{
Computation of residuals.
}
\description{
Methods for computing residuals for an object of class \code{bayesCureModel}. The Cox-Snell residuals are available for now. 
}
\usage{
\method{residuals}{bayesCureModel}(object, type = "cox-snell",...)
}
\arguments{
  \item{object}{
An object of class \code{bayesCureModel}
}
 \item{type}{
The type of residuals to be computed. 
}
\item{...}{
        ignored.
}
}
\value{
A vector of residuals.
}
\references{
Papastamoulis and Milienos (2024). Bayesian inference and cure rate modeling for event history data. TEST doi: 10.1007/s11749-024-00942-w.
}
\author{
Panagiotis Papastamoulis
}
\seealso{
\code{\link{cure_rate_MC3}}
}
\examples{
# simulate toy data just for cran-check purposes        
	set.seed(10)
        n = 4
        # censoring indicators
        stat = rbinom(n, size = 1, prob = 0.5)
        # covariates
        x <- matrix(rnorm(2*n), n, 2)
        # observed response variable 
        y <- rexp(n)
#	define a data frame with the response and the covariates        
        my_data_frame <- data.frame(y, stat, x1 = x[,1], x2 = x[,2])
# run a weibull model with default prior setup
# considering 2 heated chains 
	fit1 <- cure_rate_MC3(survival::Surv(y, stat) ~ x1 + x2, 
		data = my_data_frame, 
		promotion_time = list(distribution = 'exponential'),
		nChains = 2, 
		nCores = 1, 
		mcmc_cycles = 3, sweep=2)
	my_residuals <- residuals(fit1)

}

