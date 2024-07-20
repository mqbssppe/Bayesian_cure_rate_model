\name{cure_rate_MC3}
\alias{cure_rate_MC3}
\title{
Main function of the package
}
\description{
Runs a Metropolis Coupled MCMC (MC\eqn{^3}) sampler in order to estimate the joint posterior distribution of the model. 
}
\usage{
cure_rate_MC3(y, X, Censoring_status, nChains = 12, mcmc_cycles = 15000, 
	alpha = NULL,nCores = 8, sweep = 5, mu_g = 1, s2_g = 1, 
	a_l = 2.1, b_l = 1.1, mu_b = rep(0, dim(X)[2]), 
	Sigma = 100 * diag(dim(X)[2]), g_prop_sd = 0.045, 
	lambda_prop_scale = 0.03, b_prop_sd = rep(0.022, dim(X)[2]), 
	initialValues = NULL, plot = TRUE, adjust_scales = FALSE, 
	verbose = FALSE, tau_mala = 1.5e-05, mala = 0.15, 
	promotion_time = list(distribution = "weibull", 
	prior_parameters = matrix(rep(c(2.1, 1.1), 2), byrow = TRUE, 2, 2), 
	prop_scale = c(0.1, 0.2)), single_MH_in_f = 0.2)
}
\arguments{
  \item{y}{
Observed data, that is, a vector of length \eqn{n} with positive entries.
}
  \item{X}{
Design matrix with \eqn{k > 1} columns. Should contain a column of 1's if the model has a constant term. Should be a matrix with dimension \eqn{n\times k}.
}
  \item{Censoring_status}{
A vector \eqn{\boldsymbol\delta = (\delta_1,\ldots,\delta_n)} of binary variables corresponding to censoring indicators. The \eqn{i}-th observation is treated as a time-to-event if \eqn{\delta_i = 1} or as a censoring time otherwise (\eqn{\delta_i = 0}).
}
  \item{nChains}{
Positive integer corresponding to the number of heated chains in the MC\eqn{^3} scheme.
}
  \item{mcmc_cycles}{
Length of the generated MCMC sample. Default value: 15000. Note that each MCMC cycle consists of \code{sweep} (see below) usual MCMC iterations. 
}
  \item{alpha}{
A decreasing sequence in \eqn{[1,0)} of \code{nChains} temperatures (or heat values). The first value should always be equal to 1, which corresponds to the target posterior distribution (that is, the first chain). 
}
  \item{nCores}{
The number of cores used for parallel processing.
}
  \item{sweep}{
The number of usual MCMC iterations per MCMC cycle. Default value: 10.
}
  \item{mu_g}{
Parameter \eqn{a_{\gamma}} of the prior distribution of \eqn{\gamma}.
}
  \item{s2_g}{
Parameter \eqn{b_{\gamma}} of the prior distribution of \eqn{\gamma}.
}
  \item{a_l}{
Shape parameter \eqn{a_{\lambda}} of the Inverse Gamma prior distribution of \eqn{\lambda}.
}
  \item{b_l}{
Scale parameter \eqn{b_{\lambda}} of the Inverse Gamma prior distribution of \eqn{\lambda}.
}
  \item{mu_b}{
Mean (\eqn{\boldsymbol\mu}) of the multivariate normal prior distribution of regression coefficients. Should be a vector whose length is equal to \eqn{k}, i.e. the number of columns of the design matrix \code{X}. Default value: rep(0, k). 
}
  \item{Sigma}{
Covariance matrix of the multivariate normal prior distribution of regression coefficients.
}
  \item{g_prop_sd}{
The scale of the proposal distribution for single-site updates of the \eqn{\gamma} parameter.
}
  \item{lambda_prop_scale}{
The scale of the proposal distribution for single-site updates of the \eqn{\lambda} parameter.
}
  \item{b_prop_sd}{
The scale of the proposal distribution for the update of the \eqn{\beta} parameter (regression coefficients).
}
  \item{initialValues}{
A list of initial values for each parameter (optional).
}
  \item{plot}{
Plot MCMC sample on the run. Default: TRUE.
}
  \item{adjust_scales}{
Boolean. If TRUE the MCMC sampler runs an initial phase of a small number of iterations in order to tune the scale of the proposal distributions in the Metropolis-Hastings steps.
}
  \item{verbose}{
Print progress on the terminal if TRUE.
}
  \item{tau_mala}{
Scale of the Metropolis adjusted Langevin diffussion proposal distribution.
}
  \item{mala}{
A number between \eqn{[0,1]} indicating the proportion of times the sampler attempts a MALA proposal. Thus, the probability of attempting a typical Metropolis-Hastings move is equal to 1 - \code{mala}.
}
  \item{promotion_time}{A list with details indicating the parametric family of distribution describing the promotion time and corresponding prior distributions. See `details`.
}
  \item{single_MH_in_f}{
The probability for attempting a series of single site updates in the typical Metropolis-Hastings move. Otherwise,  with probability 1 - \code{single_MH_in_f} a Metropolis-Hastings move will attempt to update all continuous parameters simultaneously. It only makes sense when \code{mala < 1}.
}
}
\details{
It is advised to scale all continuous explanatory variables in the design matrix, so their sample mean and standard deviations are equal to 0 and 1, respectively. 
The \code{promotion_time} argument should be a list containing the following entries
\describe{
\item{\code{distribution} }{  Character string specifying the family of distributions \eqn{\{F(\cdot;\boldsymbol\alpha);\boldsymbol\alpha\in\mathcal A\}} describing the promotion time.}
\item{\code{prior_parameters} }{  Values of hyper-parameters in the prior distribution of the parameters \eqn{\boldsymbol\alpha}.}
\item{\code{prop_scale} }{  The scale of the proposal distributions for each parameter in \eqn{\boldsymbol\alpha}.}
\item{\code{dirichlet_concentration_parameter} }{  Relevant only in the case of the \code{'gamma_mixture'}. Positive scalar (typically, set to 1) determining the (common) concentration parameter of the Dirichlet prior distribution of mixing proportions. }
}
The \code{distribution} entry should be one of the following: \code{'exponential'}, \code{'weibull'}, \code{'gamma'}, \code{'logLogistic'}, \code{'gompertz'}, \code{'lomax'}, \code{'dagum'}, \code{'gamma_mixture'}.

The joint prior distribution of \eqn{\boldsymbol\alpha = (\alpha_1,\ldots,\alpha_d)} factorizes into products of inverse Gamma distributions for all (positive) parameters of \eqn{F}. Moreover, in the case of \code{'gamma_mixture'}, the joint prior also consists of another term to the Dirichlet prior distribution  on the mixing proportions. 

The \code{prop_scale} argument should be a vector with length equal to the length of vector \eqn{d} (number of elements in \eqn{\boldsymbol\alpha}), containing (positive) numbers which correspond to the scale of the proposal distribution. Note that these scale parameters are used only as initial values in case where \code{adjust_scales = TRUE}.
}
\value{
An object of class \code{bayesCureModel}, i.e. a list with the following entries
\item{mcmc_sample }{Object of class \code{mcmc} (see the \pkg{coda} package), containing the generated MCMC sample for the target chain. The column names correspond to 
\describe{
\item{\code{g_mcmc} }{ Sampled values for parameter \eqn{\gamma}}
\item{\code{lambda_mcmc} }{ Sampled values for parameter \eqn{\lambda}}
\item{\code{alpha1_mcmc}... }{ Sampled values for parameter \eqn{\alpha_1}... of the promotion time distribution \eqn{F(\cdot;\alpha_1,\ldots,\alpha_d)}. The subsequent \eqn{d-1} columns contain the sampled values for all remaining parameters, \eqn{\alpha_2,...,\ldots,\alpha_d}, where \eqn{d} depens on the family used in \code{promotion_time}.}
\item{\code{b0_mcmc} }{ Sampled values for the constant term of the regression (present only in the case where the design matrix \code{X} contains a column of 1s).}
\item{\code{b1_mcmc}... }{ Sampled values for the regression coefficient for the first explanatory variable, and similar for all subsequent columns.}
}
}
\item{complete_log_likelihood }{The complete log-likelihood for the target chain.}
\item{all_cll_values }{The complete log-likelihood for all chains}
\item{latent_status_censored }{ A data frame with the simulated binary latent status for each censored item. }
\item{swap_accept_rate}{ the acceptance rate of proposed swappings between adjacent MCMC chains.}
\item{input_data_and_model_prior}{the input data and specification of the prior parameters.}
\item{log_posterior}{the logarithm of the posterior distribution, up to a normalizing constant.}
\item{map_estimate}{The Maximum A Posterior estimate of parameters}
\item{BIC}{Bayesian Information Criterion.}
\item{AIC}{Akaike Information Criterion.}
}
\references{
Papastamoulis and Milienos (2023). Bayesian inference and cure rate modeling for event history data. arXiv:2310.06926 
}
\author{
Panagiotis Papastamoulis
}
\note{
The core function is \code{\link{cure_rate_mcmc}}.
}
\seealso{
\code{\link{cure_rate_mcmc}}
}
\examples{
# simulate toy data just for cran-check purposes        
	set.seed(10)
        n = 4
        stat = rbinom(n, size = 1, prob = 0.5)
        x <- cbind(1, matrix(rnorm(2*n), n, 2))
        y <- rexp(n)
	fit1 <- cure_rate_MC3(y = y, X = x, Censoring_status = stat, 
		promotion_time = list(distribution = 'weibull'),
		nChains = 2, nCores = 1, 
		mcmc_cycles = 3, sweep = 2)

}
