\name{complete_log_likelihood_general}
\alias{complete_log_likelihood_general}
\title{
Logarithm of the complete log-likelihood for the general cure rate model. 
}
\description{
Compute the logarithm of the complete likelihood, given a realization of the latent binary vector of cure indicators \code{I_sim} and current values of the model parameters \code{g}, \code{lambda}, \code{b} and promotion time parameters (\eqn{\boldsymbol\alpha}) which yield log-density values (one per observation) stored to the vector \code{log_f} and log-cdf values stored to the vector \code{log_F}. 
}
\usage{
complete_log_likelihood_general(y, X, Censoring_status, 
	g, lambda, log_f, log_F, b, I_sim, alpha)
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
  \item{g}{The \eqn{\gamma} parameter of the model (real). }
  \item{lambda}{The \eqn{\lambda} parameter of the model (positive). }  
  \item{log_f}{A vector containing the natural logarithm of the density function of the promotion time distribution per observation, for the current set of parameters. Its length should be equal to the sample size.}
  \item{log_F}{A vector containing the natural logarithm of the cumulative density function of the promotion time distribution per observation, for the current set of parameters. Its length should be equal to the sample size.}
  \item{b}{Vector of regression coefficients. Its length should be equal to the number of columns of the design matrix.}
  \item{I_sim}{Binary vector of the current value of the latent cured status per observation. Its length should be equal to the sample size. A time-to-event cannot be cured.}
  \item{alpha}{A parameter between 0 and 1, corresponding to the temperature of the complete posterior distribution.}
  
}
\details{
The complete likelihood of the model is 
\deqn{L_c(\boldsymbol{\theta};\boldsymbol{y}, \boldsymbol{I}) = \prod_{i\in\Delta_1}(1-p_0(\boldsymbol{x}_i,\boldsymbol\theta))f_U(y_i;\boldsymbol\theta,\boldsymbol{x}_i)\\
\prod_{i\in\Delta_0}p_0(\boldsymbol{x}_i,\boldsymbol\theta)^{1-I_i}\{(1-p_0(\boldsymbol{x}_i,\boldsymbol\theta))S_U(y_i;\boldsymbol\theta,\boldsymbol{x}_i)\}^{I_i}.}
\eqn{f_U} and \eqn{S_U} denote the probability density and survival function of the susceptibles, respectively, that is
\deqn{
S_U(y_i;\boldsymbol\theta,{\boldsymbol x}_i)=\frac{S_P(y_i;\boldsymbol{\theta},{\boldsymbol x}_i)-p_0({\boldsymbol x}_i;\boldsymbol\theta)}{1-p_0({\boldsymbol x}_i;\boldsymbol\theta)}, f_U(y_i;\boldsymbol\theta,{\boldsymbol x}_i)=\frac{f_P(y_i;\boldsymbol\theta,{\boldsymbol x}_i)}{1-p_0({\boldsymbol x}_i;\boldsymbol\theta)}.}
}
\value{
A list with the following entries
\item{cll }{the complete log-likelihood for the current parameter values.}
\item{logS }{Vector of logS values (one for each observation).}
\item{logP0 }{Vector of logP0 values (one for each observation).}
}
\references{
Papastamoulis and Milienos (2023). Bayesian inference and cure rate modeling for event history data. arXiv:2310.06926.
}
\author{
Panagiotis Papastamoulis
}
\examples{
# simulate toy data 
	set.seed(1)
	n = 4
	stat = rbinom(n, size = 1, prob = 0.5)
	x <- cbind(1, matrix(rnorm(n), n, 1))
	y <- rexp(n)
	lw <- log_weibull(y, a1 = 1, a2 = 1, c_under = 1e-9)
# compute complete log-likelihood
complete_log_likelihood_general(y = y, X = x, 
	Censoring_status = stat, 
	g = 1, lambda = 1, 
	log_f = lw$log_f, log_F = lw$log_F, 
	b = c(-0.5,0.5), 
	I_sim = stat, alpha = 1)

}
