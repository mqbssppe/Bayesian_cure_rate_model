source('sim_data_mix.R')
source('bayesian_cure_rate_model.R')

#	Simulation scenarios
# true Parameter values (14 scenarios)


F_parameters <- list(
	K = 2,
	p = c(0.5,0.5),
	gamma_parameters = matrix(c(1,1.5,4,2), 2, 2, byrow =TRUE)
)


F_parameters <- list(
	K = 2,
	p = c(0.5,0.5),
	gamma_parameters = matrix(c(1,1,6,1), 2, 2, byrow =TRUE)
)


sim <- sim_fotis_model(n = 1000, F_parameters = F_parameters, ab = 0.5, ranunL=0, ranunU=1, seed = 1, discrcov = 5)

# the first four columns of the `myData` data.frame contains the simulated dataset. 
#	column 1 ('Y'): observed times
#	column 2 ('Censoring_status'): censoring status per subject
#	column 3 ('Covariate1'): first covariate
#	column 4 ('Covariate2'): second covariate
#	The 5th column contains the cured status (which is latent for each subject however - this is ignored in the algorithm)

y <- sim[,1]
X1 <- cbind(1, sim[,3:4])
stat <- sim[,2]




#prior_parameters will be a three-dimensional array where the last slice contains the (inv-gamma) prior parameters for each mixture component
K = 2
prior_parameters = array(data = NA, dim = c(2,2,K))
prior_parameters[,,1] = matrix(rep(c(2.1, 1.1), 2), byrow = TRUE, 2, 2)
prior_parameters[,,2] = matrix(rep(c(2.1, 1.1), 2), byrow = TRUE, 2, 2)



run <- cure_rate_metropolis_hastings_cpp_general( y = y, X = X1, 
                                Censoring_status = stat,  m = 30000, alpha = 1,
                                promotion_time = list(distribution = 'gamma_mixture',
                                                K = 2, 
                                                prior_parameters = prior_parameters,
                                                prop_scale = c(0.2, 0.03, 0.2, 0.03, 1 ),# a11, a21, ..., a1K, a2K, w_1, ...w_K
                                                dirichlet_concentration_parameter = 1
                                        ),
                                mala = 0.1, single_MH_in_f=0.5
                                )




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
				Censoring_status = stat,  m = 50000, alpha = 1,
				promotion_time = list(distribution = 'gamma', 
						prior_parameters = matrix(rep(c(2.1, 1.1), 2), byrow = TRUE, 2, 2),
						prop_scale = c(0.1, 0.2)
					), 
				mala = 0.1, single_MH_in_f=0.5
				)

run3 <- cure_rate_metropolis_hastings_cpp_general( y = y, X = X1, 
				Censoring_status = stat,  m = 50000, alpha = 1,
				promotion_time = list(distribution = 'exponential', 
						prior_parameters = matrix(rep(c(2.1, 1.1), 1), byrow = TRUE, 1, 2),
						prop_scale = c(0.1)
					), 
				mala = 0.1, single_MH_in_f=0.5
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
abline(h = truePars[3:4], col = 2:4, lwd = 2, lty = 1:3)
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
mcmcIterations <- 2000
alpha <- 1/1.00001^{(1:nChains)^2.5 - 1}


# run an exponential model
promotion_time0 = list(distribution = 'exponential', 
		prior_parameters = matrix(rep(c(2.1, 1.1), 1), byrow = TRUE, 1, 2),
		prop_scale = c(0.1)
	)


run0 <- cure_rate_MC3_general( 
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
	promotion_time = promotion_time0,
########################################################################################################	        
        plot = TRUE, 		# opens plot device for plotting on the run the MCMC output
        adjust_scales = FALSE, 	# uses an initial phase in order to adjust the MH scale per parameter
        tau_mala = 0.0001, 	# initial value of the MALA scale parameter - it will be adjusted as well.
        mala = 0.15,
        single_MH_in_f = 0.5
)




# run a gamma model (including const. term)


promotion_time1 = list(distribution = 'gamma', 
		prior_parameters = matrix(rep(c(2.1, 1.1), 2), byrow = TRUE, 2, 2),
		prop_scale = c(0.1, 0.2)
	)


run1 <- cure_rate_MC3_general( 
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
	promotion_time = promotion_time1,
########################################################################################################	        
        plot = TRUE, 		# opens plot device for plotting on the run the MCMC output
        adjust_scales = FALSE, 	# uses an initial phase in order to adjust the MH scale per parameter
        tau_mala = 0.0001, 	# initial value of the MALA scale parameter - it will be adjusted as well.
        mala = 0.15,
        single_MH_in_f = 0.5
)



# run a weibull model (including const. term)

promotion_time2 = list(distribution = 'weibull', 
		prior_parameters = matrix(rep(c(2.1, 1.1), 2), byrow = TRUE, 2, 2),
		prop_scale = c(0.1, 0.2)
	)


run2 <- cure_rate_MC3_general( 
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
	promotion_time = promotion_time2,
########################################################################################################	        
        plot = TRUE, 		# opens plot device for plotting on the run the MCMC output
        adjust_scales = FALSE, 	# uses an initial phase in order to adjust the MH scale per parameter
        tau_mala = 0.0001, 	# initial value of the MALA scale parameter - it will be adjusted as well.
        mala = 0.15, single_MH_in_f = 0.5
)

# run a gamma mixture model

F_parameters <- list(
	K = 2,
	p = c(0.5,0.5),
	gamma_parameters = matrix(c(1,1,6,1), 2, 2, byrow =TRUE)
)

#
sim <- sim_fotis_model(n = 1000, lambda1 = 6, F_parameters = F_parameters, ab = 0.05, ranunL=0, ranunU=1, seed = 1, discrcov = 1)

# zero cured
sim <- sim_fotis_model(n = 1000, gamma1 = -0.9999999, lambda1 = 0.5, F_parameters = F_parameters, truebetas = c(1,0,0), ab = 0.1, ranunL=0, ranunU=1, seed = 1, discrcov = 1)


# exponential model

F_parameters <- list(
	K = 2,
	p = c(0.999999,0.000001),
	gamma_parameters = matrix(c(1,1,6,1), 2, 2, byrow =TRUE)
)

#
sim <- sim_fotis_model(n = 1000, gamma1 = 1, lambda1 = 1, F_parameters = F_parameters, ab = 1, ranunL=0, ranunU=1, seed = 1, discrcov = 1)



y <- sim[,1]
X1 <- cbind(1, sim[,3:4])
stat <- sim[,2]




K = 2
prior_parameters = array(data = NA, dim = c(2,2,K))
prior_parameters[,,1] = matrix(rep(c(2.1, 1.1), 2), byrow = TRUE, 2, 2)
prior_parameters[,,2] = matrix(rep(c(2.1, 1.1), 2), byrow = TRUE, 2, 2)



promotion_time3 = list(distribution = 'gamma_mixture',
		K = 2, 
		prior_parameters = prior_parameters,
		prop_scale = c(0.1, 0.1, 0.1, 0.1, 0.1 ),# a11, a21, ..., a1K, a2K, w_1, ...w_K
		dirichlet_concentration_parameter = 1
	)




run3 <- cure_rate_MC3_general( 
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
	promotion_time = promotion_time3,
########################################################################################################	        
        plot = TRUE, 		# opens plot device for plotting on the run the MCMC output
        adjust_scales = FALSE, 	# uses an initial phase in order to adjust the MH scale per parameter
        tau_mala = 0.0001, 	# initial value of the MALA scale parameter - it will be adjusted as well.
        mala = 0.15, 
        single_MH_in_f = 0.5
)



K = 3
prior_parameters = array(data = NA, dim = c(2,2,K))
prior_parameters[,,1] = matrix(rep(c(2.1, 1.1), 2), byrow = TRUE, 2, 2)
prior_parameters[,,2] = matrix(rep(c(2.1, 1.1), 2), byrow = TRUE, 2, 2)
prior_parameters[,,3] = matrix(rep(c(2.1, 1.1), 2), byrow = TRUE, 2, 2)



promotion_time4 = list(distribution = 'gamma_mixture',
		K = 3, 
		prior_parameters = prior_parameters,
		prop_scale = rep(0.1, 8),# a11, a21, ..., a1K, a2K, w_1, ...w_K
		dirichlet_concentration_parameter = 1
	)




run4 <- cure_rate_MC3_general( 
	y = y, 
	X = X1,
	Censoring_status = stat,		
	nChains = nChains, 		# number of heated chains (denoted as `C` in the paper)
	mcmc_cycles = 5000,	# number of MCMC cycles (denoted as `N` in the paper)
        nCores = 8, 			# number of cores to use for parallel processing
        sweep = 10,			# number of iterations per MCMC cycle (denoted  as `m_1` in the paper)
        alpha = alpha,			# a decreasing sequence of temperatures per chain 
#					        (starting from 1 - the target chain is the first one). 
#
#				PRIOR PARAMETERS        
#######################################################################################################
        mu_g = 1, s2_g = 1,	# this denotes a_g and b_g in the paper: gamma ~ f_g(a_g, b_g
        a_l = 2.1, b_l = 1.1,   # mean = 1, var = 10	lambda ~ IG(a_l, b_l)
	promotion_time = promotion_time4,
########################################################################################################	        
        plot = TRUE, 		# opens plot device for plotting on the run the MCMC output
        adjust_scales = FALSE, 	# uses an initial phase in order to adjust the MH scale per parameter
        tau_mala = 0.0001, 	# initial value of the MALA scale parameter - it will be adjusted as well.
        mala = 0, 
        single_MH_in_f = 0.5
)




# run a lomax model (including const. term)

promotion_time5 = list(distribution = 'lomax', 
		prior_parameters = matrix(rep(c(2.1, 1.1), 2), byrow = TRUE, 2, 2),
		prop_scale = c(0.1, 0.2)
	)



run5 <- cure_rate_MC3_general( 
	y = y, 
	X = X1,
	Censoring_status = stat,		
	nChains = nChains, 		# number of heated chains (denoted as `C` in the paper)
	mcmc_cycles = 5000,	# number of MCMC cycles (denoted as `N` in the paper)
        nCores = 8, 			# number of cores to use for parallel processing
        sweep = 10,			# number of iterations per MCMC cycle (denoted  as `m_1` in the paper)
        alpha = alpha,			# a decreasing sequence of temperatures per chain 
#					        (starting from 1 - the target chain is the first one). 
#
#				PRIOR PARAMETERS        
#######################################################################################################
        mu_g = 1, s2_g = 1,	# this denotes a_g and b_g in the paper: gamma ~ f_g(a_g, b_g
        a_l = 2.1, b_l = 1.1,   # mean = 1, var = 10	lambda ~ IG(a_l, b_l)
	promotion_time = promotion_time5,
########################################################################################################	        
        plot = TRUE, 		# opens plot device for plotting on the run the MCMC output
        adjust_scales = FALSE, 	# uses an initial phase in order to adjust the MH scale per parameter
        tau_mala = 0.0001, 	# initial value of the MALA scale parameter - it will be adjusted as well.
        mala = 0.15, 
        single_MH_in_f = 0.5
)


# run a dagum model (including const. term)

promotion_time6 = list(distribution = 'dagum', 
		prior_parameters = matrix(rep(c(2.1, 1.1), 3), byrow = TRUE, 3, 2),
		prop_scale = c(0.1, 0.1, 0.1)
	)



run6 <- cure_rate_MC3_general( 
	y = y, 
	X = X1,
	Censoring_status = stat,		
	nChains = nChains, 		# number of heated chains (denoted as `C` in the paper)
	mcmc_cycles = 2000,	# number of MCMC cycles (denoted as `N` in the paper)
        nCores = 8, 			# number of cores to use for parallel processing
        sweep = 10,			# number of iterations per MCMC cycle (denoted  as `m_1` in the paper)
        alpha = alpha,			# a decreasing sequence of temperatures per chain 
#					        (starting from 1 - the target chain is the first one). 
#
#				PRIOR PARAMETERS        
#######################################################################################################
        mu_g = 1, s2_g = 1,	# this denotes a_g and b_g in the paper: gamma ~ f_g(a_g, b_g
        a_l = 2.1, b_l = 1.1,   # mean = 1, var = 10	lambda ~ IG(a_l, b_l)
	promotion_time = promotion_time6,
########################################################################################################	        
        plot = TRUE, 		# opens plot device for plotting on the run the MCMC output
        adjust_scales = FALSE, 	# uses an initial phase in order to adjust the MH scale per parameter
        tau_mala = 0.0001, 	# initial value of the MALA scale parameter - it will be adjusted as well.
        mala = 0.15, 
        single_MH_in_f = 0.5
)



burn = 2000
getSummary0 <- compute_map_and_hdis_general(y = y, X = X1, Censoring_status = stat, 
	retained_mcmc = run0$mcmc_sample[-(1:burn), ], prior_parameters = NULL, promotion_time = promotion_time0)

getSummary1 <- compute_map_and_hdis_general(y = y, X = X1, Censoring_status = stat, 
	retained_mcmc = run1$mcmc_sample[-(1:burn), ], prior_parameters = NULL, promotion_time = promotion_time1)

getSummary2 <- compute_map_and_hdis(y = y, X = X1, Censoring_status = stat, 
	retained_mcmc = run2$mcmc_sample[-(1:burn), ], prior_parameters = NULL)

getSummary3 <- compute_map_and_hdis_general(y = y, X = X1, Censoring_status = stat, 
	retained_mcmc = run3$mcmc_sample[-(1:burn), ], prior_parameters = NULL, promotion_time = promotion_time3)

getSummary4 <- compute_map_and_hdis_general(y = y, X = X1, Censoring_status = stat, 
	retained_mcmc = run4$mcmc_sample[-(1:burn), ], prior_parameters = NULL, promotion_time = promotion_time4)

getSummary5 <- compute_map_and_hdis_general(y = y, X = X1, Censoring_status = stat, 
	retained_mcmc = run5$mcmc_sample[-(1:burn), ], prior_parameters = NULL, promotion_time = promotion_time5)

getSummary6 <- compute_map_and_hdis_general(y = y, X = X1, Censoring_status = stat, 
	retained_mcmc = run6$mcmc_sample[-(1:burn), ], prior_parameters = NULL, promotion_time = promotion_time6)


# print bic (smaller = better)
bics <- numeric(6)
bics[1] <- getSummary0$bic
bics[2] <- getSummary2$bic
bics[3] <- getSummary1$bic
bics[4] <- getSummary3$bic
bics[5] <- getSummary5$bic
bics[6] <- getSummary6$bic
names(bics) = c('Exponential', 'Weibull', 'Gamma', 'Gamma mix. 2',  'Lomax', 'Dagum')
bics <- bics[c(1,3,4,5,2)]
pdf(file = '../../../img/bics_sim_gamma_mix.pdf', width = 12, height = 6)
par(mar = c(8, 4, 2, 1), mfrow = c(1,2))
hist(sim[,1], 200, xlim = c(0,20), xlab = 'y', main = 'Mixture of 2 Gamma distr. under exponential censoring', freq = F )
plot(bics, type = 'b', ylab = 'Bayesian Information Criterion', xaxt = 'n', main = 'Fitted models', xlab = '')
axis(1, at = 1:length(bics), labels = names(bics), las = 2)
dev.off()




# plot estimated marginal posterior distribution along with highest density intervals
plot.bayesCureModel(retained_mcmc = run3$mcmc_sample[-(1:burn), ], map_estimate = getSummary3$map_estimate, alpha = 0.05)




