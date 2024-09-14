\name{bayesCureRateModel-package}
\alias{bayesCureRateModel-package}
\alias{bayesCureRateModel}
\docType{package}
\title{
\packageTitle{bayesCureRateModel}
}
\description{
\packageDescription{bayesCureRateModel}

The main function of the package is \code{\link{cure_rate_MC3}}. See details for a brief description of the model.
}
\details{

Let \eqn{\boldsymbol{y} = (y_1,\ldots,y_n)} denote the observed data, which correspond to time-to-event data or censoring times. Let also \eqn{\boldsymbol{x}_i = (x_{i1},\ldots,x_{x_{ip}})'} denote the covariates for subject \eqn{i}, \eqn{i=1,\ldots,n}. 

Assuming that the \eqn{n} observations are independent, the observed likelihood is defined as
\deqn{
L=L({\boldsymbol \theta}; {\boldsymbol y}, {\boldsymbol x})=\prod_{i=1}^{n}f_P(y_i;{\boldsymbol\theta},{\boldsymbol x}_i)^{\delta_i}S_P(y_i;{\boldsymbol \theta},{\boldsymbol x}_i)^{1-\delta_i},
}
where \eqn{\delta_i=1} if the \eqn{i}-th observation corresponds to time-to-event while \eqn{\delta_i=0} indicates censoring time. The parameter vector \eqn{\boldsymbol\theta} is decomposed as 
\deqn{
\boldsymbol\theta = (\boldsymbol\alpha', \boldsymbol\beta', \gamma,\lambda)
}
where 
\itemize{
\item \eqn{\boldsymbol\alpha = (\alpha_1,\ldots,\alpha_d)'\in\mathcal A} are the parameters of the promotion time distribution whose cumulative distribution and density functions are denoted as \eqn{F(\cdot,\boldsymbol\alpha)} and \eqn{f(\cdot,\boldsymbol\alpha)}, respectively. 
\item \eqn{\boldsymbol\beta\in\mathbf R^{k}} are the regression coefficients with \eqn{k} denoting the number of columns in the design matrix (it may include a constant term or not).
\item \eqn{\gamma\in\mathbf R}
\item \eqn{\lambda > 0}. 
}

The population survival and density functions are defined as
\deqn{S_P(y;\boldsymbol\theta) = \left(1 + \gamma\exp\{\boldsymbol{x}_i\boldsymbol{\beta}'\}c^{\gamma\exp\{\boldsymbol{x}_i\boldsymbol{\beta}'\}}F(y;\boldsymbol\alpha)^\lambda\right)^{-1/\gamma}}
whereas,
\deqn{f_P(y;\boldsymbol\theta)=-\frac{\partial S_P(y;\boldsymbol\theta)}{\partial y}.}
Finally, the cure rate is affected through covariates and parameters as follows
\deqn{p_0(\boldsymbol{x}_i;\boldsymbol{\theta}) = \left(1 + \gamma\exp\{\boldsymbol{x}_i\boldsymbol{\beta}'\}c^{\gamma\exp\{\boldsymbol{x}_i\boldsymbol{\beta}'\}}\right)^{-1/\gamma}}
where \eqn{c = e^{e^{-1}}}.

The promotion time distribution can be a member of standard families (currently available are the following: Exponential, Weibull, Gamma, Lomax, Gompertz, log-Logistic) and in this case \eqn{\alpha = (\alpha_1,\alpha_2)\in (0,\infty)^2}. Also considered is the Dagum distribution, which has three parameters \eqn{(\alpha_1,\alpha_2,\alpha_3)\in(0,\infty)^3}. In case that the previous parametric assumptions are not justified, the promotion time can belong to the more flexible family of finite mixtures of Gamma distributions. For example, assume a mixture of two Gamma distributions of the form
\deqn{
f(y;\boldsymbol \alpha) = \alpha_5 f_{\mathcal G}(y;\alpha_1,\alpha_3) + (1-\alpha_5) f_{\mathcal G}(y;\alpha_2,\alpha_4),
}
where \deqn{f_\mathcal{G}(y;\alpha,\beta)=\frac{\beta^{\alpha}}{\Gamma(\alpha)}y^{\alpha-1}\exp\{-\beta y\}, y>0} denotes the density of the Gamma distribution with parameters \eqn{\alpha > 0} (shape) and  \eqn{\beta > 0} (rate).
For the previous model, the parameter vector is \deqn{\boldsymbol\alpha = (\alpha_1,\alpha_2,\alpha_3,\alpha_4,\alpha_5)'\in\mathcal A}
where \eqn{\mathcal A = (0,\infty)^4\times (0,1)}. 

More generally, one can fit a mixture of \eqn{K>2} Gamma distributions. The appropriate model can be selected according to information criteria such as the BIC.


The binary vector \eqn{\boldsymbol{I} = (I_1,\ldots,I_n)} contains the (latent) cure indicators, that is, \eqn{I_i = 1} if the \eqn{i}-th subject is susceptible and \eqn{I_i = 0} if the \eqn{i}-th subject is cured.  \eqn{\Delta_0} denotes the subset of \eqn{\{1,\ldots,n\}} containing the censored subjects, whereas \eqn{\Delta_1 = \Delta_0^c} is the (complementary) subset of uncensored subjects. The complete likelihood of the model is 
\deqn{L_c(\boldsymbol{\theta};\boldsymbol{y}, \boldsymbol{I}) = \prod_{i\in\Delta_1}(1-p_0(\boldsymbol{x}_i,\boldsymbol\theta))f_U(y_i;\boldsymbol\theta,\boldsymbol{x}_i)\\
\prod_{i\in\Delta_0}p_0(\boldsymbol{x}_i,\boldsymbol\theta)^{1-I_i}\{(1-p_0(\boldsymbol{x}_i,\boldsymbol\theta))S_U(y_i;\boldsymbol\theta,\boldsymbol{x}_i)\}^{I_i}.}
\eqn{f_U} and \eqn{S_U} denote the probability density and survival function of the susceptibles, respectively, that is
\deqn{
S_U(y_i;\boldsymbol\theta,{\boldsymbol x}_i)=\frac{S_P(y_i;\boldsymbol{\theta},{\boldsymbol x}_i)-p_0({\boldsymbol x}_i;\boldsymbol\theta)}{1-p_0({\boldsymbol x}_i;\boldsymbol\theta)}, f_U(y_i;\boldsymbol\theta,{\boldsymbol x}_i)=\frac{f_P(y_i;\boldsymbol\theta,{\boldsymbol x}_i)}{1-p_0({\boldsymbol x}_i;\boldsymbol\theta)}.}

\packageIndices{bayesCureRateModel}



}
\author{
Panagiotis Papastamoulis and Fotios S. Milienos

Maintainer: \packageMaintainer{bayesCureRateModel}
}
\references{
Papastamoulis and Milienos (2024). Bayesian inference and cure rate modeling for event history data. TEST doi: 10.1007/s11749-024-00942-w.
}
\keyword{ package }
\seealso{
\code{\link{cure_rate_MC3}}
}
\examples{
# TOY EXAMPLE (very small numbers... only for CRAN check purposes)
# simulate toy data 
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
		promotion_time = list(distribution = 'weibull'),
		nChains = 2, 
		nCores = 1, 
		mcmc_cycles = 3, sweep=2)
#	print method
	fit1	
# 	summary method	
	summary1 <- summary(fit1)
	
# WARNING: the following parameters
#  mcmc_cycles, nChains
#        should take _larger_ values. E.g. a typical implementation consists of:
#        mcmc_cycles = 15000, nChains = 12
	
\donttest{
# run a Gamma mixture model with K = 2 components and default prior setup
	fit2 <- cure_rate_MC3(survival::Surv(y, stat) ~ x1 + x2, data = my_data_frame, 
		promotion_time = list(
			distribution = 'gamma_mixture',
		        K = 2),
		nChains = 8, nCores = 2, 
		mcmc_cycles = 10)
	summary2 <- summary(fit2)
}

}		
