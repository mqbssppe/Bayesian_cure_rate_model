\name{plot.bayesCureModel}
\alias{plot.bayesCureModel}
\title{
Plot method
}
\description{
	Plots and computes HDIs.
}
\usage{
\method{plot}{bayesCureModel}(x, burn = NULL, alpha0 = 0.05, gamma_mix = TRUE, 
	K_gamma = 5, plot_graphs = TRUE, bw = "nrd0", what = NULL, predict_output = NULL,  
	index_of_main_mode = NULL, draw_legend = TRUE,...)
}
\arguments{
  \item{x}{
	An object of class \code{bayesCureModel}
}
  \item{burn}{
	Number of iterations to discard as burn-in period.
}
  \item{alpha0}{
 A value between 0 and 1 in order to compute the 1-\eqn{\alpha0} Highest Posterior Density regions. 
}
  \item{gamma_mix}{
Boolean. If TRUE, the density of the marginal posterior distribution of the \eqn{\gamma} parameter is estimated from the sampled MCMC values by fitting a normal mixture  model.
}
  \item{K_gamma}{
Used only when \code{gamma_mix = TRUE} and corresponds to the number of normal mixture components used to estimate the marginal posterior density of the \eqn{\gamma} parameter. 
}
  \item{plot_graphs}{
Boolean, if FALSE only HDIs are computed.
}
  \item{bw}{
bandwidth
}
  \item{what}{
  	Integer or a character string with possible values equal to \code{'cured_prob'}, \code{'survival'} or \code{'residuals'}.
	An integer entry indicates which parameter should be plotted. If set to \code{NULL} (default), all parameters are plotted one by one. If set to \code{'cured_prob'} or \code{'survival'} the estimated cured probability or survival function is plotted, conditional on a set of covariates defined in the \code{p_cured_output} argument. In case where \code{what = 'residuals'} the residuals of the fitted model are plotted versus the quantity \eqn{-log(S)} where \eqn{S} denotes the estimated survival function arising from the Kaplan-Meier estimate based on the residuals and the censoring times.
}
  \item{predict_output}{
	Optional argument which is required only when \code{what = 'cured_prob'} or \code{what = 'survival'}. It is returned by a call to the \code{predict.bayesCureModel} function.  
}
  \item{index_of_main_mode}{
If NULL (default), the whole MCMC output is used for plotting. Otherwise, it is a subset of the retained MCMC iterations in order to identify the main mode of the posterior distribution, as returned by the \code{index_of_main_mode} value of the \code{summary.bayesCureRateModel} function. 
}
  \item{draw_legend}{
	Boolean. If TRUE (default), a legend is plotted in the case where \code{what = 'survival'} or \code{what = 'cured_prob'}.  
}
\item{...}{
        arguments passed by other methods.
}
}
\value{
The function plots graphic output on the plot device if \code{plot_graphs = TRUE}. Furthermore, a list of \eqn{100(1-\alpha)\%} Highest Density Intervals per parameter is returned. 
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
	mySummary <- summary(fit1, burn = 0)
	# plot the marginal posterior distribution of the first parameter in returned mcmc output
	plot(fit1, what = 1, burn = 0)


}
