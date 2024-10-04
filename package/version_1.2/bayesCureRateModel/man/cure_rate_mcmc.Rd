\name{cure_rate_mcmc}
\alias{cure_rate_mcmc}
\title{
The basic MCMC scheme. 
}
\description{
This is core MCMC function. The continuous parameters of the model are updated using (a) single-site Metropolis-Hastings steps and (b) a Metropolis adjusted Langevin diffusion step. The binary latent variables of the model (cured status per censored observation) are updated according to a Gibbs step. This function is embedded to the main function of the package \code{\link{cure_rate_MC3}} which runs parallel tempered MCMC chains. 
}
\usage{
cure_rate_mcmc(y, X, Censoring_status, m, alpha = 1, mu_g = 1, s2_g = 1, 
	a_l = 2.1, b_l = 1.1, promotion_time = list(distribution = "weibull", 
	prior_parameters = matrix(rep(c(2.1, 1.1), 2), byrow = TRUE, 2, 2), 
	prop_scale = c(0.2, 0.03)), mu_b = NULL, Sigma = NULL, g_prop_sd = 0.045, 
	lambda_prop_scale = 0.03, b_prop_sd = NULL, initialValues = NULL, 
	plot = FALSE, verbose = FALSE, tau_mala = 1.5e-05, mala = 0.15, 
	single_MH_in_f = 0.5, c_under = 1e-9)
}
\arguments{
  \item{y}{
observed data (time-to-event or censored time)
}
  \item{X}{
design matrix. Should contain a column of 1's if the model has a constant term. 
}
  \item{Censoring_status}{
binary variables corresponding to time-to-event and censoring.
}
  \item{m}{
number of MCMC iterations.
}
  \item{alpha}{
A value between 0 and 1, corresponding to the temperature of the complete posterior distribution. The target posterior distribution corresponds to \code{alpha = 1}.
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
  \item{promotion_time}{
A list containing the specification of the promotion time distribution. See `details`.
}
  \item{mu_b}{
Mean \eqn{\mu} of the multivariate normal prior distribution of regression coefficients. Should be a vector whose length is equal to the number of columns of the design matrix \code{X}.
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
Boolean for plotting on the run.
}
  \item{verbose}{
Boolean for printing progress on the run.
}
  \item{tau_mala}{
scale of the MALA proposal.
}
  \item{mala}{
Propability of attempting a MALA step. Otherwise, a simple MH move is attempted. 
}
  \item{single_MH_in_f}{
Probability of attempting a single-site MH move in the basic Metropolis-Hastings step. Otherwise, a joint update is attempted.
}
  \item{c_under}{
A small positive number (much smaller than 1) which is used as a threshold in the CDF of the promotion time for avoiding underflows in the computation of the log-likelihood function. Default value: 1e-9.
}
}
\value{
A list containing the following entries
\item{mcmc_sample}{The sampled MCMC values per parameter. See `note`.}
\item{complete_log_likelihood}{ Logarithm of the complete likelihood per MCMC iteration.}
\item{acceptance_rates}{The acceptance rate per move.}
\item{latent_status_censored}{ The MCMC sample of the latent status per censored observation.}
\item{log_prior_density}{Logarithm of the prior density per MCMC iteration.}
}
\note{
In the case where the promotion time distribution is a mixture model, the mixing proportions \eqn{w_1,\ldots,w_K} are reparameterized according to the following transformation
\deqn{w_j = \frac{\rho_j}{\sum_{i=1}^{K}\rho_i}, j = 1,\ldots,K}
where \eqn{\rho_i > 0} for \eqn{i=1,\ldots,K-1} and \eqn{\rho_{K}=1}. The sampler returns the parameters \eqn{\rho_1,\ldots,\rho_{K-1}}, not the mixing proportions. 
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
        set.seed(1)
        n = 10
        stat = rbinom(n, size = 1, prob = 0.5)
        x <- cbind(1, matrix(rnorm(2*n), n, 2))
        y <- rexp(n)
# run a weibull model (including const. term) 
#	for m = 10 mcmc iterations 
        fit1 <- cure_rate_mcmc(y = y, X = x, Censoring_status = stat, 
              	plot = FALSE,
                promotion_time = list(distribution = 'weibull', 
                        prior_parameters = matrix(rep(c(2.1, 1.1), 2), 
                                                byrow = TRUE, 2, 2),
                        prop_scale = c(0.1, 0.1)
                ),
                m = 10)
#	the generated mcmc sampled values         
	fit1$mcmc_sample
}

