\name{logLik.bayesCureModel}
\alias{logLik.bayesCureModel}
\title{
Extract the log-likelihood.
}
\description{
Method to extract the log-likelihood of a \code{bayesCureModel} object. 
}
\usage{
\method{logLik}{bayesCureModel}(object, ...)
}
\arguments{
  \item{object}{
An object of class \code{bayesCureModel}
}
\item{...}{
        ignored.
}
}
\value{
The maximum (observed) log-likelihood value obtained across the MCMC run.
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
		promotion_time = list(family = 'exponential'),
		nChains = 2, 
		nCores = 1, 
		mcmc_cycles = 3, sweep=2)
	ll <- logLik(fit1)

}

