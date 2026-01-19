pkgname <- "bayesCureRateModel"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('bayesCureRateModel')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("Surv")
### * Surv

flush(stderr()); flush(stdout())

### Name: Surv
### Title: Create a Survival Object
### Aliases: Surv
### Keywords: survival

### ** Examples

# Right-censored survival data
Surv(5, 1)
Surv(c(5, 10), c(1, 0))



cleanEx()
nameEx("bayesCureRateModel-package")
### * bayesCureRateModel-package

flush(stderr()); flush(stdout())

### Name: bayesCureRateModel-package
### Title: Bayesian Cure Rate Modeling for Time-to-Event Data
### Aliases: bayesCureRateModel-package bayesCureRateModel
### Keywords: package

### ** Examples

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
		promotion_time = list(family = 'weibull'),
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
	




cleanEx()
nameEx("complete_log_likelihood_general")
### * complete_log_likelihood_general

flush(stderr()); flush(stdout())

### Name: complete_log_likelihood_general
### Title: Logarithm of the complete log-likelihood for the general cure
###   rate model.
### Aliases: complete_log_likelihood_general

### ** Examples

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




cleanEx()
nameEx("compute_fdr_tpr")
### * compute_fdr_tpr

flush(stderr()); flush(stdout())

### Name: compute_fdr_tpr
### Title: Compute the achieved FDR and TPR
### Aliases: compute_fdr_tpr

### ** Examples

set.seed(1)
v1 <- sample(0:1, size = 100, replace=TRUE, prob=c(0.8,0.2) )
v2 <- runif(100)
compute_fdr_tpr(true_latent_status = v1, posterior_probs = v2)



cleanEx()
nameEx("cure_rate_MC3")
### * cure_rate_MC3

flush(stderr()); flush(stdout())

### Name: cure_rate_MC3
### Title: Main function of the package
### Aliases: cure_rate_MC3

### ** Examples

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
	fit1 <- cure_rate_MC3(survival::Surv(y, stat) ~ x1 + x2, 
		data = my_data_frame,
		promotion_time = list(family = 'weibull'),
		nChains = 2, nCores = 1, 
		mcmc_cycles = 3, sweep = 2)




cleanEx()
nameEx("cure_rate_mcmc")
### * cure_rate_mcmc

flush(stderr()); flush(stdout())

### Name: cure_rate_mcmc
### Title: The basic MCMC scheme.
### Aliases: cure_rate_mcmc

### ** Examples

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
                promotion_time = list(family = 'weibull', 
                        prior_parameters = matrix(rep(c(2.1, 1.1), 2), 
                                                byrow = TRUE, 2, 2),
                        prop_scale = c(0.1, 0.1)
                ),
                m = 10)
#	the generated mcmc sampled values         
	fit1$mcmc_sample



cleanEx()
nameEx("logLik.bayesCureModel")
### * logLik.bayesCureModel

flush(stderr()); flush(stdout())

### Name: logLik.bayesCureModel
### Title: Extract the log-likelihood.
### Aliases: logLik.bayesCureModel

### ** Examples

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




cleanEx()
nameEx("log_dagum")
### * log_dagum

flush(stderr()); flush(stdout())

### Name: log_dagum
### Title: PDF and CDF of the Dagum distribution
### Aliases: log_dagum

### ** Examples

log_dagum(y = 1:10, a1 = 1, a2 = 1, a3 = 1, c_under = 1e-9)



cleanEx()
nameEx("log_gamma")
### * log_gamma

flush(stderr()); flush(stdout())

### Name: log_gamma
### Title: PDF and CDF of the Gamma distribution
### Aliases: log_gamma

### ** Examples

log_gamma(y = 1:10, a1 = 1, a2 = 1, c_under = 1e-9)



cleanEx()
nameEx("log_gamma_mixture")
### * log_gamma_mixture

flush(stderr()); flush(stdout())

### Name: log_gamma_mixture
### Title: PDF and CDF of a Gamma mixture distribution
### Aliases: log_gamma_mixture

### ** Examples

y <- runif(10)
a1 <- c(1,2)
a2 <- c(1,1)
p <- c(0.9,0.1)
log_gamma_mixture(y, a1, a2, p)



cleanEx()
nameEx("log_gompertz")
### * log_gompertz

flush(stderr()); flush(stdout())

### Name: log_gompertz
### Title: PDF and CDF of the Gompertz distribution
### Aliases: log_gompertz

### ** Examples

log_gompertz(y = 1:10, a1 = 1, a2 = 1, c_under = 1e-9)



cleanEx()
nameEx("log_logLogistic")
### * log_logLogistic

flush(stderr()); flush(stdout())

### Name: log_logLogistic
### Title: PDF and CDF of the log-Logistic distribution.
### Aliases: log_logLogistic

### ** Examples

log_logLogistic(y = 1:10, a1 = 1, a2 = 1, c_under = 1e-9)



cleanEx()
nameEx("log_lomax")
### * log_lomax

flush(stderr()); flush(stdout())

### Name: log_lomax
### Title: PDF and CDF of the Lomax distribution
### Aliases: log_lomax

### ** Examples

log_lomax(y = 1:10, a1 = 1, a2 = 1, c_under = 1e-9)



cleanEx()
nameEx("log_user_mixture")
### * log_user_mixture

flush(stderr()); flush(stdout())

### Name: log_user_mixture
### Title: Define a finite mixture of a given family of distributions.
### Aliases: log_user_mixture

### ** Examples

# We will define a mixture of 2 exponentials distributions.
# First we pass the exponential distribution at user_f
user_f <- function(y, a){
	log_f <- dexp(y, rate = a, log = TRUE)
	log_F <- pexp(y, rate = a, log.p = TRUE)
	result <- vector('list', length = 2)
	names(result) <- c('log_f', 'log_F')
	result[["log_f"]] = log_f
	result[["log_F"]] = log_F
	return(result)
}
#	simulate some date
y <- runif(10)
# Now compute the log of pdf and cdf for a mixture of K=2 exponentials
p <- c(0.9,0.1)
a <- matrix(c(0.1, 1.5), nrow = 1, ncol = 2)
log_user_mixture(user_f = user_f, y = y, a = a, p = p)



cleanEx()
nameEx("log_weibull")
### * log_weibull

flush(stderr()); flush(stdout())

### Name: log_weibull
### Title: PDF and CDF of the Weibull distribution
### Aliases: log_weibull

### ** Examples

log_weibull(y = 1:10, a1 = 1, a2 = 1, c_under = 1e-9)




cleanEx()
nameEx("plot.bayesCureModel")
### * plot.bayesCureModel

flush(stderr()); flush(stdout())

### Name: plot.bayesCureModel
### Title: Plot method
### Aliases: plot.bayesCureModel

### ** Examples

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





cleanEx()
nameEx("plot.predict_bayesCureModel")
### * plot.predict_bayesCureModel

flush(stderr()); flush(stdout())

### Name: plot.predict_bayesCureModel
### Title: Plot method
### Aliases: plot.predict_bayesCureModel

### ** Examples

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





cleanEx()
nameEx("predict.bayesCureModel")
### * predict.bayesCureModel

flush(stderr()); flush(stdout())

### Name: predict.bayesCureModel
### Title: Predict method.
### Aliases: predict.bayesCureModel

### ** Examples

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
	newdata <- data.frame(x1 = c(0.2,-1), x2 = c(-1,0))
	# return predicted values at tau = c(0.5, 1)
	my_prediction <- predict(fit1, newdata = newdata, 
		burn = 0, tau_values = c(0.5, 1))




cleanEx()
nameEx("residuals.bayesCureModel")
### * residuals.bayesCureModel

flush(stderr()); flush(stdout())

### Name: residuals.bayesCureModel
### Title: Computation of residuals.
### Aliases: residuals.bayesCureModel

### ** Examples

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
	my_residuals <- residuals(fit1)




cleanEx()
nameEx("summary.bayesCureModel")
### * summary.bayesCureModel

flush(stderr()); flush(stdout())

### Name: summary.bayesCureModel
### Title: Summary method.
### Aliases: summary.bayesCureModel

### ** Examples

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
	mySummary <- summary(fit1, burn = 0)




cleanEx()
nameEx("summary.predict_bayesCureModel")
### * summary.predict_bayesCureModel

flush(stderr()); flush(stdout())

### Name: summary.predict_bayesCureModel
### Title: Summary method for predictions.
### Aliases: summary.predict_bayesCureModel

### ** Examples

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



### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
