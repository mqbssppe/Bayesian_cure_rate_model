source('sim_data_function_2.R')
source('bayesian_cure_rate_model.R')

# sample size
n <- 500

#	Simulation scenarios
# true Parameter values (14 scenarios)

scenario <- vector("list", length = 14)
scenario[[1]] <- c(1,1.5,0.8,0.8,1.5,1.5,-0.8, 0.85)
scenario[[2]] <- c(1,1.5,0.8,0.8,1.5,1.5,-0.8, 2.40)
scenario[[3]] <- c(1,1,0.5,0.5,-0.8,1.5,1.5, 0.5)
scenario[[4]] <- c(1,1,0.5,0.5,-0.8,1.5,1.5, 2.5)
scenario[[5]] <- c(1, 1, 1, 1, -4, 1, 1, 0.8)
scenario[[6]] <- c(1, 1, 1, 1, -4, 1, 1, 4)
scenario[[7]] <- c(-0.05,1,0.8,1,2,-1,1, 0.45)
scenario[[8]] <- c(-0.05,1,0.8,1,2,-1,1, 1.2)
scenario[[9]] <- c(-0.5,1,0.8,1,2,-0.7,1, 0.2)
scenario[[10]] <- c(-0.5,1,0.8,1,2,-0.7,1, 0.5)
scenario[[11]] <- c(-0.9999999,0.5,0.5,0.5,1,0,0, 0.1)
scenario[[12]] <- c(-0.9999999,0.5,0.5,0.5,1,0,0, 0.3)
scenario[[13]] <- c(-0.9999999,0.5,0.5,0.5,1,0,0, 0.2)
scenario[[14]] <- c(-0.9999999,0.5,0.5,0.5,1,0,0, 0.42)

# example with scenario 1

scenario_id <- 1

truePars <- scenario[[scenario_id]][1:7]
ab = scenario[[scenario_id]][8]

if(scenario_id < 5){
        source("sim_data_function_2.R")
}else{
        source('sim_data_extra_scenario.R')               
}
sim <- sim_fotis_model(n = n, truePars = truePars, ab = ab, seed = 1)           

# the first four columns of the `myData` data.frame contains the simulated dataset. 
#	column 1 ('Y'): observed times
#	column 2 ('Censoring_status'): censoring status per subject
#	column 3 ('Covariate1'): first covariate
#	column 4 ('Covariate2'): second covariate
#	The 5th column contains the cured status (which is latent for each subject however - this is ignored in the algorithm)

y <- sim[,1]
X1 <- cbind(1, sim[,3:4])
stat <- sim[,2]

a1 = truePars[3]
a2 = truePars[4]

lw <- log_weibull(y, a1 = a1, a2 = a2,  c_under = 1e-9)

cll_new <- complete_log_likelihood_general(y = y, X = X1, Censoring_status = stat, g = truePars[1], lambda = truePars[2], 
	log_f = lw$log_f, log_F = lw$log_F, 
	b = truePars[5:7], I_sim = sim[,'cured_status'], 
	alpha = 1)

cll_old <- complete_log_likelihood(y = y, X = X1, Censoring_status = stat, g = truePars[1], lambda = truePars[2], 
	a1 = a1, a2 = a2, 
	b = truePars[5:7], I_sim = sim[,'cured_status'], 
	alpha = 1)





run <- cure_rate_metropolis_hastings_cpp_general( y = y, X = X1, 
				Censoring_status = stat,  m = 30000, alpha = 1,
				promotion_time = list(distribution = 'weibull', 
						prior_parameters = matrix(rep(c(2.1, 1.1), 2), byrow = TRUE, 2, 2),
						prop_scale = c(0.2, 0.03)
					)
				)

run2 <- cure_rate_metropolis_hastings_cpp_general( y = y, X = X1, 
				Censoring_status = stat,  m = 500000, alpha = 1,
				promotion_time = list(distribution = 'exponential', 
						prior_parameters = matrix(rep(c(2.1, 1.1), 1), byrow = TRUE, 1, 2),
						prop_scale = c(0.2)
					)
				)

par(mfrow = c(2,2))
plot(run2$mcmc_sample[seq(1,500000, by = 100),1], col = 2, type = 'l', ylab = bquote(gamma))
abline(h = truePars[1], col = 1, lwd = 2)
legend('topright', 'true value', col = 1, lty = 1)
plot(run2$mcmc_sample[seq(1,500000, by = 100),2], col = 2, type = 'l', ylab = bquote(lambda))
abline(h = truePars[2], col = 1, lwd = 2)
legend('topright', 'true values', col = 1, lty = 1)
matplot(run2$mcmc_sample[seq(1,500000, by = 100),4:6], col = 2:4, type = 'l', ylab = 'regression coefficients')
abline(h = truePars[5:7], col = 2:4, lwd = 2, lty = 1:3)
legend('bottomright', rep('true values', 3), col = 2:4, lty = 1:3)
matplot(run2$mcmc_sample[seq(1,500000, by = 100),3], col = 2, type = 'l', ylab = 'value', main = "Exponential distribution parameters")
legend('topright', paste0('a', 1), col = 2, lty = 1)


run3 <- cure_rate_metropolis_hastings_cpp_general( y = y, X = X1, 
				Censoring_status = stat,  m = 500000, alpha = 1,
				promotion_time = list(distribution = 'gamma', 
						prior_parameters = matrix(rep(c(2.1, 1.1), 2), byrow = TRUE, 2, 2),
						prop_scale = c(0.1, 0.2)
					)
				)
par(mfrow = c(2,2))
plot(run3$mcmc_sample[seq(1,500000, by = 100),1], col = 2, type = 'l', ylab = bquote(gamma))
abline(h = truePars[1], col = 1, lwd = 2)
legend('topright', 'true value', col = 1, lty = 1)
plot(run3$mcmc_sample[seq(1,500000, by = 100),2], col = 2, type = 'l', ylab = bquote(lambda))
abline(h = truePars[2], col = 1, lwd = 2)
legend('topright', 'true values', col = 1, lty = 1)
matplot(run3$mcmc_sample[seq(1,500000, by = 100),5:7], col = 2:4, type = 'l', ylab = 'regression coefficients')
abline(h = truePars[5:7], col = 2:4, lwd = 2, lty = 1:3)
legend('bottomright', rep('true values', 3), col = 2:4, lty = 1:3)
matplot(run3$mcmc_sample[seq(1,500000, by = 100),3:4], col = 2:3, type = 'l', ylab = 'value', main = "Gamma distribution parameters")
legend('topright', paste0('a', 1:2), col = 2:3, lty = 1)





run4 <- cure_rate_metropolis_hastings_cpp_general( y = y, X = X1, 
				Censoring_status = stat,  m = 500000, alpha = 1,
				promotion_time = list(distribution = 'dagum', 
						prior_parameters = matrix(rep(c(2.1, 1.1), 3), byrow = TRUE, 3, 2),
						prop_scale = c(0.3, 0.1, 0.1)
					)
				)

par(mfrow = c(2,2))
plot(run4$mcmc_sample[seq(1,500000, by = 100),1], col = 2, type = 'l', ylab = bquote(gamma))
abline(h = truePars[1], col = 1, lwd = 2)
legend('topright', 'true value', col = 1, lty = 1)
plot(run4$mcmc_sample[seq(1,500000, by = 100),2], col = 2, type = 'l', ylab = bquote(lambda))
abline(h = truePars[2], col = 1, lwd = 2)
legend('topright', 'true values', col = 1, lty = 1)
matplot(run4$mcmc_sample[seq(1,500000, by = 100),6:8], col = 2:4, type = 'l', ylab = 'regression coefficients')
abline(h = truePars[5:7], col = 2:4, lwd = 2, lty = 1:3)
legend('bottomright', rep('true values', 3), col = 2:4, lty = 1:3)
matplot(run4$mcmc_sample[seq(1,500000, by = 100),3:5], col = 2:4, type = 'l', ylab = 'value', main = "Dagum distribution parameters")
legend('topright', paste0('a', 1:3), col = 2:4, lty = 1)




nChains <- 8
mcmcIterations <- 10000
alpha <- 1/1.00001^{(1:nChains)^2.5 - 1}

# run the `true` model (including const. term)

run1 <- cure_rate_MC3( 
	y = y, 
	X = X1,
	Censoring_status = stat,		
	nChains = nChains, 		# number of heated chains (denoted as `C` in the paper)
	mcmc_cycles = mcmcIterations,	# number of MCMC cycles (denoted as `N` in the paper)
        nCores = 8, 			# number of cores to use for parallel processing
        sweep = 10,			# number of iterations per MCMC cycle (denoted  as `m_1` in the paper)
        alpha = alpha,			# a decreasing sequence of temperatures per chain 
#					        (starting from 1 - the target chain is the first one). 
#
#				PRIOR PARAMETERS        
#######################################################################################################
        mu_g = 1, s2_g = 1,	# this denotes a_g and b_g in the paper: gamma ~ f_g(a_g, b_g
        a_l = 2.1, b_l = 1.1,   # mean = 1, var = 10	lambda ~ IG(a_l, b_l)
        a_1 = 2.1, b_1 = 1.1,   # mean = 1, var = 10	alpha1 ~ IG(a_1, b_1)
        a_2 = 2.1, b_2 = 1.1,   # mean = 1, var = 10  	alpha2 ~ IG(a_2, b_2)
########################################################################################################	        
        plot = TRUE, 		# opens plot device for plotting on the run the MCMC output
        adjust_scales = TRUE, 	# uses an initial phase in order to adjust the MH scale per parameter
        adjust_alpha = FALSE, 	# always set to false, is not working otherwise
        tau_mala = 0.0001, 	# initial value of the MALA scale parameter - it will be adjusted as well.
        mala = 0.33		# probability of proposing a MALA update  - otherwise a single MH is attempted
)

# run a model without constant term plus an irrelevant explanatory variable
X2 <- cbind(sim[,3:4], rnorm(n))

run2 <- cure_rate_MC3( 
	y = y, 
	X = X2,
	Censoring_status = stat,		
	nChains = nChains, 		# number of heated chains (denoted as `C` in the paper)
	mcmc_cycles = mcmcIterations,	# number of MCMC cycles (denoted as `N` in the paper)
        nCores = 8, 			# number of cores to use for parallel processing
        sweep = 10,			# number of iterations per MCMC cycle (denoted  as `m_1` in the paper)
        alpha = alpha,			# a decreasing sequence of temperatures per chain 
#					        (starting from 1 - the target chain is the first one). 
#
#				PRIOR PARAMETERS        
#######################################################################################################
        mu_g = 1, s2_g = 1,	# this denotes a_g and b_g in the paper: gamma ~ f_g(a_g, b_g
        a_l = 2.1, b_l = 1.1,   # mean = 1, var = 10	lambda ~ IG(a_l, b_l)
        a_1 = 2.1, b_1 = 1.1,   # mean = 1, var = 10	alpha1 ~ IG(a_1, b_1)
        a_2 = 2.1, b_2 = 1.1,   # mean = 1, var = 10  	alpha2 ~ IG(a_2, b_2)
########################################################################################################	        
        plot = TRUE, 		# opens plot device for plotting on the run the MCMC output
        adjust_scales = TRUE, 	# uses an initial phase in order to adjust the MH scale per parameter
        adjust_alpha = FALSE, 	# always set to false, is not working otherwise
        tau_mala = 0.0001, 	# initial value of the MALA scale parameter - it will be adjusted as well.
        mala = 0.33		# probability of proposing a MALA update  - otherwise a single MH is attempted
)




burn = 2000
getSummary1 <- compute_map_and_hdis(y = y, X = X1, Censoring_status = stat, 
	retained_mcmc = run1$mcmc_sample[-(1:burn), ], prior_parameters = NULL)

getSummary2 <- compute_map_and_hdis(y = y, X = X2, Censoring_status = stat, 
	retained_mcmc = run2$mcmc_sample[-(1:burn), ], prior_parameters = NULL)

# print bic (smaller = better)
getSummary1$bic
getSummary2$bic

# plot estimated marginal posterior distribution along with highest density intervals
plot.bayesCureModel(retained_mcmc = run1$mcmc_sample[-(1:burn), ], map_estimate = getSummary1$map_estimate, alpha = 0.05)




