\name{summary.bayesCureModel}
\alias{summary.bayesCureModel}
\title{
Summary method.
}
\description{
This function produces all summaries after fitting a cure rate model. 
}
\usage{
\method{summary}{bayesCureModel}(object, burn = NULL, gamma_mix = TRUE, 
	K_gamma = 3,  K_max = 3, fdr = 0.1, 
	covariate_levels = NULL, yRange = NULL, alpha = 0.1, ...)
}
\arguments{
  \item{object}{
An object of class \code{bayesCureModel}
}
  \item{burn}{
Positive integer corresponding to the number of mcmc iterations to discard as burn-in period
}
  \item{gamma_mix}{
Boolean. If TRUE, the density of the marginal posterior distribution of the \eqn{\gamma} parameter is estimated from the sampled MCMC values by fitting a normal mixture  model.
}
  \item{K_gamma}{
Used only when \code{gamma_mix = TRUE} and corresponds to the number of normal mixture components used to estimate the marginal posterior density of the \eqn{\gamma} parameter. 
}
  \item{K_max}{
 Maximum number of components in order to cluster the (univariate) values of the joint posterior distribution across the MCMC run. Used to identify the main mode of the posterior distribution. See details. 
}
  \item{fdr}{
	The target value for controlling the False Discovery Rate when classifying subjects as cured or not. 
}
  \item{covariate_levels}{
	Optional levels for the covariates. It is only required when the user wishes to obtain a vector with the estimated posterior cured probabilities for a given combination of covariates. Include the value ``1'' in the case where the model contains constant term. 
}
  \item{yRange}{
	Optional range (a vector of two non-negative values) for computing the sequence of posterior probabilities for the given values in \code{covariate_levels}. 
}
  \item{alpha}{
	Scalar between 0 and 1 corresponding to 1 - confidencel level for computing Highest Density Intervals. If set to NULL, the confidence intervals are not computed.
}
\item{...}{
        ignored.
}
}
\value{
A list with the following entries
\item{map_estimate}{Maximum A Posteriori (MAP) estimate of the parameters of the model.}
\item{highest_density_intervals}{Highest Density Interval per parameter}
\item{latent_cured_status}{Estimated posterior probabilities of the latent cure status per censored subject.}
\item{cured_at_given_FDR}{Classification as cured or not, at given FDR level.}
\item{p_cured_output}{It is returned only in the case where the argument \code{covariate_values} is not \code{NULL}. See details.}
\item{main_mode_index}{The retained MCMC iterations which correspond to the main mode of the posterior distribution.}
}
\details{
The values of the posterior draws are clustered according to a (univariate) normal mixture model, and the main mode corresponds to the cluster with the largest mean. The maximum number of mixture components corresponds to the \code{K_max} argument. The \pkg{mclust} library is used for this purpose. The inference for the latent cure status of each (censored) observation is based on the MCMC draws corresponding to the main mode of the posterior distribution. The FDR is controlled according to the technique proposed in Papastamoulis and Rattray (2018). 

In case where \code{covariate_levels} is set to \code{TRUE}, the \code{summary} function also returns a list named \code{p_cured_output} with the following entries
\describe{
\item{mcmc}{It is returned only in the case where the argument \code{covariate_values} is not \code{NULL}. A vector of posterior cured probabilities for the given values in \code{covariate_values}, per retained MCMC draw.}
\item{map}{It is returned only in the case where the argument \code{covariate_values} is not \code{NULL}. The cured probabilities computed at the MAP estimate of the parameters, for the given values \code{covariate_values}.}
\item{tau_values}{tau values}
\item{covariate_levels}{covariate levels}
\item{index_of_main_mode}{the subset of MCMC draws allocated to the main mode of the posterior distribution.}
}
}
\references{
Papastamoulis and Milienos (2023). Bayesian inference and cure rate modeling for event history data. arXiv:2310.06926.

Papastamoulis and Rattray (2018). A Bayesian Model Selection Approach for Identifying Differentially Expressed Transcripts from RNA Sequencing Data, Journal of the Royal Statistical Society Series C: Applied Statistics, Volume 67, Issue 1.

Scrucca L, Fraley C, Murphy TB, Raftery AE (2023). Model-Based Clustering, Classification, and Density Estimation Using mclust in R. Chapman and Hall/CRC. ISBN 978-1032234953
}
\author{
Panagiotis Papastamoulis
}
\seealso{
\code{\link{cure_rate_MC3}}
}
\examples{
# simulate toy data just for cran-check purposes        
        set.seed(1)
        n = 4
        stat = rbinom(n, size = 1, prob = 0.5)
        x <- cbind(1, matrix(rnorm(2*n), n, 2))
        y <- rexp(n)
	fit1 <- cure_rate_MC3(y = y, X = x, Censoring_status = stat, 
		promotion_time = list(distribution = 'exponential'),
		nChains = 2, nCores = 1, 
		mcmc_cycles = 3, sweep = 2)
	mySummary <- summary(fit1, burn = 0)
}

