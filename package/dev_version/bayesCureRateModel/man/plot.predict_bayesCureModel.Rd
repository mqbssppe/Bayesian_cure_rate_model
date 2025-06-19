\name{plot.predict_bayesCureModel}
\alias{plot.predict_bayesCureModel}
\title{
Plot method
}
\description{
	Plot the output of the predict method.
}
\usage{
\method{plot}{predict_bayesCureModel}(x, what = 'survival', draw_legend = TRUE,...)
}
\arguments{
  \item{x}{
	An object of class \code{predict_bayesCureModel}
}
  \item{what}{
  	Character with possible values: \code{'cured_prob'} or \code{'survival'}, corresponding to the estimated cured probability or survival function.
}
  \item{draw_legend}{
	Boolean. If TRUE (default), a legend is plotted in the case where \code{what = 'survival'} or \code{what = 'cured_prob'}.  
}
\item{...}{
        arguments passed by other methods.
}
}
\value{
No value, just a plot. 
}
\author{
Panagiotis Papastamoulis
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
	fit1 <- cure_rate_MC3(survival::Surv(y, stat) ~ x1 + x2, data = my_data_frame, 
		promotion_time = list(family = 'exponential'),
		nChains = 2, 
		nCores = 1, 
		mcmc_cycles = 3, sweep=2)
	#compute predictions for two individuals with 
	#	x1 = 0.2 and x2 = -1
	#	and 
	#	x1 = -1 and x2 = 0
	covariate_levels1 <- data.frame(x1 = c(0.2,-1), x2 = c(-1,0))
	predictions <- predict(fit1, newdata = covariate_levels1, burn = 0)
	# plot cured probabilities based on the previous output
	plot(predictions, what='cured_prob')


}
