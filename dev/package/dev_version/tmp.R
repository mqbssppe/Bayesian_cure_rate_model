library("bayesCureRateModel")
library("VGAM")
library("calculus")
library("coda")
library("doParallel")
source('bayesian_cure_rate_model.R')


user_promotion_time <- function(y, a, c_under = 1e-9){
	c_under <- log(c_under)
	log_f <- ddagum(y, scale = a[1], shape1.a = a[2], shape2.p = a[3], log = TRUE)
	log_F <- pdagum(y, scale = a[1], shape1.a = a[2], shape2.p = a[3], log.p = TRUE)
	ind <- (log_F < c_under)
	if(length(ind) > 0){
		log_F[ind] <- c_under
	}
	result <- vector('list', length = 2)
	names(result) <- c('log_f', 'log_F')
	result[["log_f"]] = log_f
	result[["log_F"]] = log_F
	return(result)

}
promotion_time = list(distribution = 'user', 
		define = user_promotion_time,
		prior_parameters = matrix(rep(c(2.1, 1.1), 3), byrow = TRUE, 3, 2),
		prop_scale = c(0.1, 0.1, 0.1)
	)




# lognormal with parameter μ = log a1, a1 > 0 and sigma2 = a2, a2 > 0
user_promotion_time <- function(y, a){
	log_f <- dnorm(log( y ), log(a[1]), a[2], log = TRUE) - log( y )
	log_F <- pnorm(log( y ), log(a[1]), a[2], log = TRUE)
	result <- vector('list', length = 2)
	names(result) <- c('log_f', 'log_F')
	result[["log_f"]] = log_f
	result[["log_F"]] = log_F
	return(result)

}

# option a: 
promotion_time = list(distribution = 'user', 
		define = user_promotion_time,
		prior_parameters = matrix(rep(c(2.1, 1.1), 2), byrow = TRUE, 2, 2),
		prop_scale = c(0.1, 0.1)
	)
############################################################################################################
############################################################################################################
############################################################################################################
# option b: user_mixture
#	---log-normal
user_promotion_time <- function(y, a){
	log_f <- dnorm(log( y ), log(a[1]), a[2], log = TRUE) - log( y )
	log_F <- pnorm(log( y ), log(a[1]), a[2], log = TRUE)
	result <- vector('list', length = 2)
	names(result) <- c('log_f', 'log_F')
	result[["log_f"]] = log_f
	result[["log_F"]] = log_F
	return(result)
}

K = 2
n_f = 2
prior_parameters = array(data = NA, dim = c(n_f,2,K))
for(k in 1:K){
	prior_parameters[,,k] = matrix(rep(c(2.1, 1.1), 2), byrow = TRUE, n_f, 2)			
}
#	---exponential
user_promotion_time <- function(y, a){
	log_f <- dexp(y, rate = a, log = TRUE)
	log_F <- pexp(y, rate = a, log.p = TRUE)
	result <- vector('list', length = 2)
	names(result) <- c('log_f', 'log_F')
	result[["log_f"]] = log_f
	result[["log_F"]] = log_F
	return(result)
}

K = 2
n_f = 1
prior_parameters = array(data = NA, dim = c(n_f,2,K))
for(k in 1:K){
	prior_parameters[,,k] = matrix(rep(c(2.1, 1.1), n_f), byrow = TRUE, n_f, 2)			
}


promotion_time = list(distribution = 'user_mixture', 
		define = user_promotion_time,
		prior_parameters = prior_parameters,
		prop_scale =  rep(0.1, K*n_f + K - 1), 		
		K = K,
		dirichlet_concentration_parameter = 1
	)


############################################################################################################
############################################################################################################
############################################################################################################

data(marriage_dataset)
age_mean <- mean(marriage_dataset$age)
age_sd <- sd(marriage_dataset$age)
marriage_dataset$age <- as.numeric(scale(marriage_dataset$age))
y <- as.numeric(marriage_dataset$time)
X1 <- as.matrix(cbind(1, marriage_dataset$age, marriage_dataset$kids, 
  1*(marriage_dataset$race == 2), 1*(marriage_dataset$race == 4)))
colnames(X1) = c('const', 'age', 'kids', 'hispanic', 'other')
stat <- as.numeric(marriage_dataset$censoring)


y = y; X = X1; Censoring_status = stat
m = 5000 
alpha = 1
mu_g = 1
s2_g = 1
a_l = 2.1
b_l = 1.1
mu_b = NULL
Sigma = NULL
g_prop_sd = 0.045
lambda_prop_scale = 0.03
b_prop_sd = NULL
initialValues = NULL
plot = TRUE
verbose = TRUE
tau_mala = 0.000015
mala = 0.15
single_MH_in_f = 0.5
			
						
fit <- cure_rate_mcmc( y = y, X = X1, Censoring_status = stat,  m = m, alpha = 1,
				mu_g = 1, s2_g = 1, a_l = 2.1, b_l = 1.1, 
				promotion_time = promotion_time,
				mu_b = NULL, Sigma = NULL,
				g_prop_sd = 0.045, 
				lambda_prop_scale = 0.03, 
				b_prop_sd = NULL, 
						initialValues = NULL, 
						plot = TRUE,
						verbose = TRUE,
						tau_mala = 0.000015, mala = 0.15, single_MH_in_f = 0.5
					)
					
mcmc_cycles = 15000; nChains = 8; nCores = 4
set.seed(10)
run_exp <- cure_rate_MC3(y = y, X = X1, Censoring_status = stat,                
  nChains = nChains, mcmc_cycles = mcmc_cycles,  nCores = nCores, 
  promotion_time = promotion_time, verbose = TRUE)
					
head(model.matrix(time~age+as.factor(kids) + as.factor(race), data = marriage_dataset))
library(survival)	
surv_object <- with(marriage_dataset, Surv(time, censoring))
head(model.matrix(surv_object~age+as.factor(kids) + as.factor(race), data = marriage_dataset))	


formu <- formula(surv_object ~ age + kids + race)

model.matrix(formula(surv_object ~ age + kids + race), data = my_data_frame)
	
	
mcmc_cycles = 2000; nChains = 8; nCores = 4
set.seed(10)
my_data_frame <- data.frame(surv_object, age = marriage_dataset$age, kids = as.factor(marriage_dataset$kids), race = as.factor(marriage_dataset$race))
formu <- formula(surv_object ~ age + kids + race)
fit <- cure_rate_MC3(formula = formu, data = my_data_frame,                
  nChains = nChains, mcmc_cycles = mcmc_cycles,  nCores = nCores, 
  promotion_time = promotion_time, verbose = TRUE)


mcmc_cycles = 10000; nChains = 8; nCores = 4
set.seed(10, kind = "L'Ecuyer-CMRG")
fit <- cure_rate_MC3(formula = formu, data = my_data_frame,                
  nChains = nChains, mcmc_cycles = mcmc_cycles,  nCores = nCores, 
  promotion_time = list(distribution = "exponential"), verbose = TRUE)

set.seed(10, kind = "L'Ecuyer-CMRG")
fit2 <- cure_rate_MC3(formula = formu, data = my_data_frame,                
  nChains = nChains, mcmc_cycles = mcmc_cycles,  nCores = nCores, 
  promotion_time = list(distribution = "exponential"), verbose = TRUE)

set.seed(10, kind = "L'Ecuyer-CMRG")
fit3 <- cure_rate_MC3(formula = formu, data = my_data_frame,                
  nChains = nChains, mcmc_cycles = mcmc_cycles,  nCores = nCores, 
  promotion_time = list(distribution = "weibull"), verbose = TRUE)



covariate_levels = data.frame(age = 0, kids = 1, race = 1)
	
x1 <- (20 - age_mean)/age_sd
x2 <- (30 - age_mean)/age_sd
x3 <- (40 - age_mean)/age_sd
covariate_levels1 <- data.frame(age = c(x1, x2, x3), kids = c(0,0,0), race = c(1,1,1))
covariate_levels2 <- data.frame(age = c(x1, x2, x3), kids = c(1,1,1), race = c(1,1,1))
covariate_levels3 <- data.frame(age = c(x1, x2, x3), kids = c(0,0,0), race = c(2,2,2))
covariate_levels4 <- data.frame(age = c(x1, x2, x3), kids = c(1,1,1), race = c(2,2,2))
covariate_levels5 <- data.frame(age = c(x1, x2, x3), kids = c(0,0,0), race = c(4,4,4))
covariate_levels6 <- data.frame(age = c(x1, x2, x3), kids = c(1,1,1), race = c(4,4,4))

ss_exp1 <- summary(fit, covariate_levels = covariate_levels1)
ss_exp2 <- summary(fit, covariate_levels = covariate_levels2)
ss_exp3 <- summary(fit, covariate_levels = covariate_levels3)
ss_exp4 <- summary(fit, covariate_levels = covariate_levels4)
ss_exp5 <- summary(fit, covariate_levels = covariate_levels5)
ss_exp6 <- summary(fit, covariate_levels = covariate_levels6)


covariate_levels <- rbind(covariate_levels1, covariate_levels2, covariate_levels3, covariate_levels4, covariate_levels5, covariate_levels6)

my_predictions <- predict(fit, newdata = covariate_levels)

res <- residuals(fit)

km <- survfit(Surv(res, marriage_dataset$censoring)~1)
my_index <- numeric(length(km$time));for(i in 1:length(res)){my_index[i] <- which(res == km$time[i])[1]}
del <- which(is.na(my_index))

plot(res[my_index[-del]], -log(km$surv)); abline(0,1)


res2 <- -log(km$surv)
qqplot(x=qexp(ppoints(length(res2))), y=res2, main="Exponential Q-Q Plot",
       xlab="Theoretical Quantiles", ylab= "Your Data Quantiles")
abline(0,1)

###################################################
### code chunk number 9: visualization2
###################################################
par(mfrow = c(2,2), mar = c(4,6,1,1))
plot(fit, what='survival', p_cured_output = ss_exp1$p_cured_output, 
  ylim = c(0,1), cex.axis = 2.0, cex.lab = 2.5, draw_legend = FALSE, alpha = 0.1)
plot(fit, what='survival', p_cured_output = ss_exp2$p_cured_output, 
  ylim = c(0,1), cex.axis = 2.0, cex.lab = 2.5, alpha = 0.1, draw_legend = FALSE)
plot(fit, what='cured_prob', p_cured_output = ss_exp1$p_cured_output, 
  ylim = c(0,1), cex.axis = 2.0, cex.lab = 2.5, alpha = 0.1, draw_legend = TRUE)
plot(fit, what='cured_prob', p_cured_output = ss_exp2$p_cured_output, 
  ylim = c(0,1), cex.axis = 2.0, cex.lab = 2.5, alpha = 0.1)

	
	

#' @export
fitted.bayesCureModel <- function(object, fdr = 0.1, quant = 0.5, ...){
	y = object$input_data_and_model_prior$y
	X = object$input_data_and_model_prior$X
	Censoring_status = object$input_data_and_model_prior$Censoring_status
	promotion_time = object$input_data_and_model_prior$promotion_time
	n <- dim(X)[1]
	ss <- summary(object, fdr = fdr)	
	cured_status <- rep('susceptible', n)
	cured_status[Censoring_status == 0] <- ss$cured_at_given_FDR
	fitted_values <- rep(NA, n)
	fitted_values[cured_status == 'cured'] <- Inf
	
	burn = 0
	K_max = 3	
	if(is.null(burn)){
		burn = floor(dim(object$mcmc_sample)[1]/3)
	}else{
		if(burn > dim(object$mcmc_sample)[1] - 1){stop('burn in period not valid.')}
		if(burn < 0){stop('burn in period not valid.')}		
	}

	if(burn > 0){
	retained_mcmc = object$mcmc_sample[-(1:burn),]
	}else{retained_mcmc = object$mcmc_sample}
	mu_g = object$input_data_and_model_prior$mu_g
	s2_g = object$input_data_and_model_prior$s2_g	
	mu_b = object$input_data_and_model_prior$mu_b
	Sigma = object$input_data_and_model_prior$Sigma
	a_l = object$input_data_and_model_prior$a_l
	b_l = object$input_data_and_model_prior$b_l
	map_index <- 1
	retained_mcmc <- rbind(object$map_estimate, retained_mcmc)	

	c_under = 1e-9
	ct = exp(exp(-1))
	X <- as.matrix(X)
	x <- X
	nCov <- dim(x)[2]

	n_pars_f <- dim(promotion_time$prior_parameters)[1]

	if(promotion_time$distribution == 'gamma_mixture'){
		K = promotion_time$K
		if(dim(promotion_time$prior_parameters)[3] != K){stop("incosistent number of mixture components")}		
		n_pars_f = K * n_pars_f + K - 1# note that we include all mixing weights
	}
	if(promotion_time$distribution == 'user_mixture'){
		K = promotion_time$K
		n_pars_f_k = n_pars_f
		if(dim(promotion_time$prior_parameters)[3] != K){stop("incosistent number of mixture components")}		
		n_pars_f = K * n_pars_f + K - 1# note that we include all mixing weights
	}

	
	a <- numeric(n_pars_f)

	log_S_p <- function(g, lambda, log_F, b, x){
		theta <- exp(x %*% b)
		return(-log(1 + c(g * theta * ct^{g*theta}) * exp(log_F)^lambda)/g)

	}

	log_p0 <- function(g, b, x){
		theta <- exp(x %*% b)
		return(-(log(1 + g*theta*ct^(g*theta)))/g)
	}



	log_f_p <- function(g, lambda, log_f, log_F, b, logS){
		# logS = log_S_p(tau = tau, g = g, lambda = lambda, a1 = a1, a2 = a2, b0 = b0, b1 = b1, b2 = b2)
		log_theta <- x %*% b
		return(
		        (1 + g) * logS + log(lambda) + log_theta +
		        g*exp(log_theta)*log(ct) + 
		        (lambda - 1)*log_F + #### NOTE: this causes the log -> -inf, so now it regulated.
		        log_f
		)
	 }

	if(length(unique(X[,1])) == 1){
		bIndices <- paste0('b',1:nCov - 1,'_mcmc')	
	}else{
		bIndices <- paste0('b',1:nCov,'_mcmc')
	}


	F_u <- function(tau, x){
		# tau is the time
		# x denotes the covariates
		g = retained_mcmc[1,'g_mcmc']
		lambda = retained_mcmc[1,'lambda_mcmc']
		b = retained_mcmc[1, bIndices]
		if(promotion_time$distribution == 'exponential'){
			lw <- log_weibull(y = tau, a1 = retained_mcmc[iter,'a1_mcmc'], 
			a2 = 1,  c_under = c_under)
		}
		if(promotion_time$distribution == 'weibull'){
			lw <- log_weibull(y = tau, a1 = retained_mcmc[iter,'a1_mcmc'], 
				a2 = retained_mcmc[iter,'a2_mcmc'],  c_under = c_under)
		}
		if(promotion_time$distribution == 'gamma'){
			lw <- log_gamma(y = tau, a1 = retained_mcmc[iter,'a1_mcmc'], 
				a2 = retained_mcmc[iter,'a2_mcmc'],  c_under = c_under)
		}
		if(promotion_time$distribution == 'gompertz'){
			lw <- log_gompertz(y = tau, a1 = retained_mcmc[iter,'a1_mcmc'], 
				a2 = retained_mcmc[iter,'a2_mcmc'],  c_under = c_under)
		}
		if(promotion_time$distribution == 'logLogistic'){
			lw <- log_logLogistic(y = tau, a1 = retained_mcmc[iter,'a1_mcmc'], 
				a2 = retained_mcmc[iter,'a2_mcmc'],  c_under = c_under)
		}
		if(promotion_time$distribution == 'lomax'){
			lw <- log_lomax(y = tau, a1 = retained_mcmc[iter,'a1_mcmc'], 
				a2 = retained_mcmc[iter,'a2_mcmc'],  c_under = c_under)
		}
		if(promotion_time$distribution == 'dagum'){
			lw <- log_dagum(y = tau, a1 = retained_mcmc[iter,'a1_mcmc'], a2 = retained_mcmc[iter,'a2_mcmc'], 
				a3 = retained_mcmc[iter,'a3_mcmc'], c_under = c_under)
		}
		if(promotion_time$distribution == 'gamma_mixture'){
			K = promotion_time$K
			a = retained_mcmc[iter,3:(3+n_pars_f - 1)]
			w <- c(a[-(1:(2*K))], 1)
			p <- w/sum(w)
			a1_ind <- seq(1, 2*K, by = 2)
			a2_ind <- seq(2, 2*K, by = 2)		
			a1 = a[a1_ind]
			a2 = a[a2_ind]		
			lw <- log_gamma_mixture(y = tau, a1 = a1, a2 = a2, p = p, c_under = c_under)
		}
		if(promotion_time$distribution == 'user_mixture'){
			K = promotion_time$K
			a = retained_mcmc[iter,3:(3+n_pars_f - 1)]
			w <- c(a[-(1:(n_pars_f_k*K))], 1)
			p <- w/sum(w)
			a_matrix = matrix(a[1:(n_pars_f_k*K)], n_pars_f_k, K, byrow = TRUE)
			lw <- log_user_mixture(user_f = promotion_time$define, y = tau, a = a_matrix, p = p, c_under = c_under)
		}
		if(promotion_time$distribution == 'user'){
			lw <- promotion_time$define(y = tau, a = retained_mcmc[iter, 3:(2+n_pars_f)])	
	
		}

		log_F = lw$log_F
		Sp <- exp(log_S_p(g, lambda, log_F, b, x))
		p0 <- exp(log_p0(g, b, x))
		return((1 - Sp)/c((1-p0)))
	}
#	define equation to solve
	my_objective_function <- function(tau, x, a){
		if(a < 0){stop("a should be between 0 and 1.")}
		if(a > 1){stop("a should be between 0 and 1.")}
		return(-a + F_u(tau, x))
	}
#	define function for computing numerically the inverse of F_u
	tau_max <- 2 * max(object$input_data_and_model_prior$y)
	event_subjects <- which(Censoring_status == 1)
	for(i in event_subjects){
		invert_cdf <- function(a_quant)return(uniroot(my_objective_function, 
				c(0, tau_max), x = X[i, ], a = a_quant, tol = 0.00001)$root)
		fitted_values[i] <- invert_cdf(quant)
	
	}





	n <- dim(X)[1]
	n_parameters <- dim(x)[2] + 2 + n_pars_f
	map_estimate <- object$map_estimate

	p_cured_given_tau <-  NULL 
	covariate_levels <- X
	tau_values <- y
	cox_snell <- numeric(n)
	
	
	for(iter in 1:1){

		i <- 0
		for(tau in tau_values){
			i <- i + 1
			cox_snell[i] <- - log_S_p(g = g, lambda = lambda,
							log_F = log_F[i], 
							b = b,
							x = covariate_levels[i,]
							)
		}
	}
	return(cox_snell)
}

	
