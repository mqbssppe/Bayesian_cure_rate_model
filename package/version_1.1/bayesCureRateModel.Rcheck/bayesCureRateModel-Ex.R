pkgname <- "bayesCureRateModel"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
base::assign(".ExTimings", "bayesCureRateModel-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('bayesCureRateModel')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("bayesCureRateModel-package")
### * bayesCureRateModel-package

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: bayesCureRateModel-package
### Title: Bayesian Cure Rate Modeling for Time-to-Event Data
### Aliases: bayesCureRateModel-package bayesCureRateModel
### Keywords: package

### ** Examples

# TOY EXAMPLE (very small numbers... only for CRAN check purposes)
# simulate toy data 
	set.seed(10)
	n = 4
	stat = rbinom(n, size = 1, prob = 0.5)
	x <- cbind(1, matrix(rnorm(n), n, 1))
	y <- rexp(n)
# run a weibull model with default prior setup
# considering 2 heated chains 
	fit1 <- cure_rate_MC3(y = y, X = x, Censoring_status = stat, 
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
	
## No test: 
# run a Gamma mixture model with K = 2 components and default prior setup
	fit2 <- cure_rate_MC3(y = y, X = x, Censoring_status = stat, 
		promotion_time = list(
			distribution = 'gamma_mixture',
		        K = 2),
		nChains = 8, nCores = 2, 
		mcmc_cycles = 10)
	summary2 <- summary(fit2)
## End(No test)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("bayesCureRateModel-package", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("complete_log_likelihood_general")
### * complete_log_likelihood_general

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
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




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("complete_log_likelihood_general", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("cure_rate_MC3")
### * cure_rate_MC3

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: cure_rate_MC3
### Title: Main function of the package
### Aliases: cure_rate_MC3

### ** Examples

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




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("cure_rate_MC3", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("cure_rate_mcmc")
### * cure_rate_mcmc

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
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
                promotion_time = list(distribution = 'weibull', 
                        prior_parameters = matrix(rep(c(2.1, 1.1), 2), 
                                                byrow = TRUE, 2, 2),
                        prop_scale = c(0.1, 0.1)
                ),
                m = 10)
#	the generated mcmc sampled values         
	fit1$mcmc_sample



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("cure_rate_mcmc", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("log_dagum")
### * log_dagum

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: log_dagum
### Title: PDF and CDF of the Dagum distribution
### Aliases: log_dagum

### ** Examples

log_dagum(y = 1:10, a1 = 1, a2 = 1, a3 = 1, c_under = 1e-9)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("log_dagum", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("log_gamma")
### * log_gamma

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: log_gamma
### Title: PDF and CDF of the Gamma distribution
### Aliases: log_gamma

### ** Examples

log_gamma(y = 1:10, a1 = 1, a2 = 1, c_under = 1e-9)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("log_gamma", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("log_gamma_mixture")
### * log_gamma_mixture

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: log_gamma_mixture
### Title: PDF and CDF of a Gamma mixture distribution
### Aliases: log_gamma_mixture

### ** Examples

y <- runif(10)
a1 <- c(1,2)
a2 <- c(1,1)
p <- c(0.9,0.1)
log_gamma_mixture(y, a1, a2, p)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("log_gamma_mixture", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("log_gompertz")
### * log_gompertz

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: log_gompertz
### Title: PDF and CDF of the Gompertz distribution
### Aliases: log_gompertz

### ** Examples

log_gompertz(y = 1:10, a1 = 1, a2 = 1, c_under = 1e-9)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("log_gompertz", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("log_logLogistic")
### * log_logLogistic

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: log_logLogistic
### Title: PDF and CDF of the log-Logistic distribution.
### Aliases: log_logLogistic

### ** Examples

log_logLogistic(y = 1:10, a1 = 1, a2 = 1, c_under = 1e-9)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("log_logLogistic", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("log_lomax")
### * log_lomax

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: log_lomax
### Title: PDF and CDF of the Lomax distribution
### Aliases: log_lomax

### ** Examples

log_lomax(y = 1:10, a1 = 1, a2 = 1, c_under = 1e-9)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("log_lomax", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("log_weibull")
### * log_weibull

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: log_weibull
### Title: PDF and CDF of the Weibull distribution
### Aliases: log_weibull

### ** Examples

log_weibull(y = 1:10, a1 = 1, a2 = 1, c_under = 1e-9)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("log_weibull", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("plot.bayesCureModel")
### * plot.bayesCureModel

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: plot.bayesCureModel
### Title: Plot method
### Aliases: plot.bayesCureModel

### ** Examples

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
## No test: 
	#compute cured probability for two individuals with 
	#	x1 = 0.2 and x2 = -1
	#	and 
	#	x1 = -1 and x2 = 0
	covariate_levels1 <- rbind(c(1,0.2,-1), c(1,-1,0))
	summary1 <- summary(fit1, covariate_levels = covariate_levels1, burn = 0)
	plot(fit1, what='cured_prob', p_cured_output = summary1$p_cured_output, 
	  ylim = c(0,1))
	
## End(No test)





base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("plot.bayesCureModel", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("summary.bayesCureModel")
### * summary.bayesCureModel

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: summary.bayesCureModel
### Title: Summary method.
### Aliases: summary.bayesCureModel

### ** Examples

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




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("summary.bayesCureModel", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
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
