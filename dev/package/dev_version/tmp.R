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


covariate_levels = data.frame(age = 0, kids = 1, race = 1)
	
x1 <- (20 - age_mean)/age_sd
x2 <- (30 - age_mean)/age_sd
x3 <- (40 - age_mean)/age_sd
covariate_levels1 <- data.frame(age = c(x1, x2, x3), kids = c(0,0,0), race = c(1,1,1))
covariate_levels2 <- data.frame(age = c(x1, x2, x3), kids = c(1,1,1), race = c(1,1,1))

ss_exp1 <- summary(fit, covariate_levels = covariate_levels1, burn = 700)
ss_exp2 <- summary(fit, covariate_levels = covariate_levels2, burn = 700)


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

	
	
	
