\name{plot.bayesCureModel}
\alias{plot.bayesCureModel}
\title{
Plot method
}
\description{
	Plots and computes HDIs.
}
\usage{
\method{plot}{bayesCureModel}(x, burn = NULL, alpha = 0.05, gamma_mix = TRUE, 
	K_gamma = 5, plot_graphs = TRUE, bw = "nrd0", what = NULL, p_cured_output = NULL,  
	index_of_main_mode = NULL, draw_legend = TRUE,...)
}
\arguments{
  \item{x}{
	An object of class \code{bayesCureModel}
}
  \item{burn}{
	Number of iterations to discard as burn-in period.
}
  \item{alpha}{
 A value between 0 and 1 in order to compute the 1-\eqn{\alpha} Highest Posterior Density regions. 
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
  	Integer or a character string with possible values equal to \code{'cured_prob'} or \code{'survival'}.
	An integer entry indicates which parameter should be plotted. If set to \code{NULL} (default), all parameters are plotted one by one. If set to \code{'cured_prob'} or \code{'survival'} the estimated cured probability or survival function is plotted, conditional on a set of covariates defined in the \code{p_cured_output} argument. 
}
  \item{p_cured_output}{
	Optional argument (list) which is required only when \code{what = 'cured_prob'} or \code{what = 'survival'}. It is returned by the \code{summary.bayesCureRateModel} function.  
}
  \item{index_of_main_mode}{
If NULL (default), all modes are plotted. Otherwise, it is a subset of the retained MCMC iterations in order to identify the main mode of the posterior distribution, as returned by the \code{index_of_main_mode} value of the \code{summary.bayesCureRateModel} function. 
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
        set.seed(1)
        n = 4
        stat = rbinom(n, size = 1, prob = 0.5)
        # simulate design matrix
        #	first column consists of 1s (const)
        #	and second and third column contain
        #	the values of two covariates
        x <- cbind(1, matrix(rnorm(2*n), n, 2))
        colnames(x) <- c('const', 'x1', 'x2')
        y <- rexp(n)
	fit1 <- cure_rate_MC3(y = y, X = x, Censoring_status = stat, 
		promotion_time = list(distribution = 'exponential'),
		nChains = 2, nCores = 1,
		mcmc_cycles = 3, sweep = 2)
	# plot the marginal posterior distribution of the first parameter in returned mcmc output
	plot(fit1, what = 1, burn = 0)
# using 'cured_prob'
\donttest{
	#compute cured probability for two individuals with 
	#	x1 = 0.2 and x2 = -1
	#	and 
	#	x1 = -1 and x2 = 0
	covariate_levels1 <- rbind(c(1,0.2,-1), c(1,-1,0))
	summary1 <- summary(fit1, covariate_levels = covariate_levels1, burn = 0)
	plot(fit1, what='cured_prob', p_cured_output = summary1$p_cured_output, 
	  ylim = c(0,1))
	
}


}
