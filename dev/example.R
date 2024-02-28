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
myData <- sim_fotis_model(n = n, truePars = truePars, ab = ab, seed = 1)           

# the first four columns of the `myData` data.frame contains the simulated dataset. 
#	column 1 ('Y'): observed times
#	column 2 ('Censoring_status'): censoring status per subject
#	column 3 ('Covariate1'): first covariate
#	column 4 ('Covariate2'): second covariate
#	The 5th column contains the cured status (which is latent for each subject however - this is ignored in the algorithm)
latent_statuses <- myData[,5]
myData <- myData[,1:4]

# define the sets D0 and D1 containing censored items and time-to-events, respectively
D0 <- which(myData[,"Censoring_status"]==0)	# censored items
D1 <- which(myData[,"Censoring_status"] == 1) 	# events

cll <- complete_log_likelihood(myData = myData, g = truePars[1], lambda = truePars[2], a1 = truePars[3], a2 = truePars[4], b =  truePars[5:7], I_sim = latent_statuses, alpha = 1.0)



m = 1000
alpha = 1
mu_g = 0.2
s2_g = 0.1
a_l = 2.001
b_l = 1
a_1 = 2.001
b_1 = 1
a_2 = 2.001
b_2 = 1
mu_b = rep(0,5)
Sigma = 100*diag(5)
g_prop_sd = 0.045
lambda_prop_scale = 0.03
a1_prop_scale = 0.2
a2_prop_scale = 0.03
b_prop_sd = rep(0.022, 5)
initialValues = NULL
plot = TRUE
verbose = TRUE
tau_mala = 0.000015
mala = 0.33


myData2 <- cbind(myData, matrix(rnorm(n*2), n, 2))

run <- cure_rate_metropolis_hastings_cpp( myData = myData2,  m = 50000, alpha = 1,
				mu_g = 0.2, s2_g = 0.1, a_l = 2.001, b_l = 1, 
				a_1 = 2.001, b_1 = 1, a_2 = 2.001, b_2 = 1, 
				mu_b = rep(0,5), Sigma = 100*diag(5),
				g_prop_sd = 0.045, 
				lambda_prop_scale = 0.03, 
				a1_prop_scale = 0.2, 
				a2_prop_scale = 0.03, 
				b_prop_sd = rep(0.022, 5), 
						initialValues = NULL, 
						plot = TRUE,
						verbose = TRUE,
						tau_mala = 0.000015, mala = 0.33)

nChains <- 8    
mcmcIterations <- 20000
alpha <- 1/1.00001^{(1:nChains)^2.5 - 1}


	myData = myData2 		# this is the input dataset - should be in the format described above
	nChains = nChains 		# number of heated chains (denoted as `C` in the paper)
	mcmc_cycles = mcmcIterations	# number of MCMC cycles (denoted as `N` in the paper)
        nCores = 8 			# number of cores to use for parallel processing
        sweep = 10			# number of iterations per MCMC cycle (denoted  as `m_1` in the paper)
        alpha = alpha			# a decreasing sequence of temperatures per chain 
#					        (starting from 1 - the target chain is the first one). 
#
#				PRIOR PARAMETERS        
	nCov <- dim(myData)[2] - 1	
#######################################################################################################
        mu_g = 1; s2_g = 1	# this denotes a_g and b_g in the paper: gamma ~ f_g(a_g, b_g
        a_l = 2.1; b_l = 1.1   # mean = 1, var = 10	lambda ~ IG(a_l, b_l)
        a_1 = 2.1; b_1 = 1.1   # mean = 1, var = 10	alpha1 ~ IG(a_1, b_1)
        a_2 = 2.1; b_2 = 1.1   # mean = 1, var = 10  	alpha2 ~ IG(a_2, b_2)
        Sigma = 10*diag(nCov)     # beta = (beta_0, beta_1, beta_2) ~ N_3 (0, Sigma)
########################################################################################################	        
        plot = TRUE 		# opens plot device for plotting on the run the MCMC output
        adjust_scales = TRUE 	# uses an initial phase in order to adjust the MH scale per parameter
        adjust_alpha = FALSE 	# always set to false, is not working otherwise
        tau_mala = 0.0001 	# initial value of the MALA scale parameter - it will be adjusted as well.
        mala = 0		# probability of proposing a MALA update  - otherwise a single MH is attempted

nCov <- dim(myData)[2] - 1

run <- cure_rate_MC3( 
	myData = myData2, 		# this is the input dataset - should be in the format described above
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




# plot mcmc traces with corresponding real values after discarding first 10000 cycles
varnames <- as.expression(c(
        bquote(gamma), 
        bquote(lambda), 
        bquote(alpha[1]), 
        bquote(alpha[2]), 
        bquote(beta[0]), 
        bquote(beta[1]), 
        bquote(beta[2]), 
        bquote(beta[3]), 
        bquote(beta[4])        
        ))

par(mfrow = c(3,3))
truePars2 <- c(truePars, c(0,0))
for(i in 1:9){
	plot(run$mcmc_sample[-(1:10000),i], ylab = varnames[i], xlab = 'iteration')
	abline(h = truePars2[i], col = 'yellow', lwd = 2)
}


