library(bayesCureRateModel)

# simulate a gamma mixture model
F_parameters <- list(
        K = 2,
        p = c(0.5,0.5),
        gamma_parameters = matrix(c(2,2,12,2), 2, 2, byrow =TRUE)
)

#
source('sim_data_mix.R')
n <- 500
sim <- sim_fotis_model(n = n, lambda1 = 6, F_parameters = F_parameters, ab = 0.05, ranunL=0, ranunU=1, seed = 1, discrcov = 2)

true_status = rep('susceptible', n)
true_status[sim[,'cured_status'] == 1] = 'cured'
sim_mix_data <- data.frame(
			time = sim[,1], 
			censoring = sim[,2], 
			x1 = sim[,3], 
			x2 = as.factor(sim[,4]),
			true_status = as.factor(true_status)
			)

mcmc_cycles = 1000
nChains = 4
nCores = 1

	
user_promotion_time <- function(y, a){
	log_f <- -0.5*log(2*pi) - log(y) - log(a[2]) - ((log(y) - log(a[1]))^2)/(2 * a[2]^2)
	log_F <- pnorm((log(y) - log(a[1]))/a[2], log.p = TRUE)
	result <- vector('list', length = 2)
	names(result) <- c('log_f', 'log_F')
	result[["log_f"]] = log_f
	result[["log_F"]] = log_F
	return(result)

}


promotion_time = list(distribution = 'user', 
		define = user_promotion_time,
		prior_parameters = matrix(rep(c(2.1, 1.1), 2), byrow = TRUE, 2, 2),
		prop_scale = c(0.1, 0.1)
	)

set.seed(1)
run_ln <- cure_rate_MC3(survival::Surv(time, censoring) ~ x1 + x2, data = sim_mix_data, 
	mcmc_cycles = mcmc_cycles, promotion_time = promotion_time, nChains = nChains, nCores = nCores, verbose = TRUE)


K = 2
n_f = 2
prior_parameters = array(data = NA, dim = c(n_f,2,K))
for(k in 1:K){
	prior_parameters[,,k] = matrix(rep(c(2.1, 1.1), 2), byrow = TRUE, n_f, 2)			
}

promotion_time = list(distribution = 'user_mixture', 
		define = user_promotion_time,
		prior_parameters = prior_parameters,
		prop_scale =  rep(0.1, K*n_f + K - 1), 		
		K = K,
		dirichlet_concentration_parameter = 1
	)

set.seed(1)
run_ln_mix <- cure_rate_MC3(survival::Surv(time, censoring) ~ x1 + x2, data = sim_mix_data, 
	mcmc_cycles = mcmc_cycles, promotion_time = promotion_time, nChains = nChains, nCores = nCores, verbose = TRUE)



	bics <- numeric(2)
	bics[1] <- run_ln$BIC
	bics[2] <- run_ln_mix$BIC
	names(bics) = c('Log-normal', 'Log-normal mixture')
	
	ss_ln <- summary(run_ln, burn = 300)
	ss_ln_mix <- summary(run_ln_mix, burn = 300)	
	latent_cured_status_ln <- ss_ln$latent_cured_status
	latent_cured_status_ln_mix <- ss_ln_mix$latent_cured_status
		
	labels <- sim_mix_data$true_status[sim_mix_data$censoring == 0]
	labels <- factor(labels, levels = c('susceptible', 'cured'), ordered = TRUE)
	pred_ln <- prediction(latent_cured_status_ln, labels)
	pred_ln_mix <- prediction(latent_cured_status_ln_mix, labels)
	
	perf_ln <- performance(pred_ln, "tpr", "fpr")
	perf_ln_mix <- performance(pred_ln_mix, "tpr", "fpr")
	
	plot(perf_ln, lwd = 2, xlim = c(0,0.5), ylim = c(0.5,1), col = 2, lty = 2)
	plot(perf_ln_mix, add = TRUE,  lwd = 2, col = 3, lty = 1)	

	auc_ln <- performance(pred_ln, measure = "auc")@y.values[[1]]
	auc_ln_mix <- performance(pred_ln_mix, measure = "auc")@y.values[[1]]
	aucs <- c(auc_ln, auc_ln_mix)
	legend("bottomright", paste0(names(bics), ' (AUC: ', round(aucs,3),')'), col = 2:4, lty = c(2,1,3), lwd = 2)

	myCut = c(1,2, 5,10)/100
	true_latent_status <- as.numeric(labels) - 1
	def <-  c(0,1,2,5)
	def2 <- c(15, 16, 17, 18)
	plot(c(0,0.20), c(0,1), type = 'n',xlab = "achieved FDR",ylab = "True positive rate")
	abline(v = myCut, lty = 2, col = 'gray')
	
	posterior_probs <- latent_cured_status_ln
	fdr_tpr <- compute_fdr_tpr(true_latent_status, posterior_probs, myCut = myCut)
	my_pch <- def
	contr <- apply(fdr_tpr[,c(1,3)], 1, diff) > 0	
	my_pch[contr] <- def2[contr]
	points(fdr_tpr, pch = my_pch, lwd = 2, cex = 2, type = 'b', col = 2, lty = 2)

	posterior_probs <- latent_cured_status_ln_mix
	fdr_tpr <- compute_fdr_tpr(true_latent_status, posterior_probs, myCut = myCut)
	my_pch <- def
	contr <- apply(fdr_tpr[,c(1,3)], 1, diff) > 0	
	my_pch[contr] <- def2[contr]
	points(fdr_tpr, pch = my_pch, lwd = 2, cex = 2, type = 'b', col = 3, lty = 1)

	legend('bottomright', c('nominal FDR', myCut), pch = c(NA, def2), col = 'gray')
	
