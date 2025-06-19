\name{summary.predict_bayesCureModel}
\alias{summary.predict_bayesCureModel}
\title{
Summary method for predictions.
}
\description{
This function produces MCMC summaries for an object of class \code{predict_bayesCureModel}. 
}
\usage{
\method{summary}{predict_bayesCureModel}(object, ...)
}
\arguments{
  \item{object}{
An object of class \code{predict_bayesCureModel}.
}
\item{...}{
        Other options passed to the \code{summary.mcmc} method of the \code{coda} package.
}
}
\value{
A list with the following entries
\item{survival}{MCMC summaries (quantiles) for the survival function.}
\item{cured_probability}{MCMC summaries (quantiles) for the conditional cured probability.}
\item{cumulative_hazard}{MCMC summaries (quantiles) for the cumulative hazard function.}
\item{hazard_rate}{MCMC summaries (quantiles) for the hazard rate function.}
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
#       define a data frame with the response and the covariates        
	my_data_frame <- data.frame(y, stat, x1 = x[,1], x2 = x[,2])
# run a weibull model with default prior setup
# considering 2 heated chains 
	fit1 <- cure_rate_MC3(survival::Surv(y, stat) ~ x1 + x2, data = my_data_frame, 
	     promotion_time = list(family = 'exponential'),
	     nChains = 2, 
	     nCores = 1, 
	     mcmc_cycles = 3, sweep=2)
	newdata <- data.frame(x1 = c(0.2,-1), x2 = c(-1,0))
	# return predicted values at tau = c(0.5, 1)
	my_prediction <- predict(fit1, newdata = newdata, 
	     burn = 0, tau_values = c(0.5, 1))
	my_summary <- summary(my_prediction, quantiles = c(0.1,0.9))		
}

