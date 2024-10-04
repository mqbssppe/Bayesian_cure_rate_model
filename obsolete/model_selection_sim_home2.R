library(ROCR)
source('bayesian_cure_rate_model.R')


compute_fdr_tpr <- function(true_latent_status, posterior_probs, myCut = c(0.01, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4)){
        l <- length(myCut)
        realDE <- array(1-true_latent_status, dim = c(length(true_latent_status),1))
        p <- matrix(1 - posterior_probs, ncol = 1)
        perm <- order(p,decreasing = TRUE)
        orderedP <- p[perm,]
        nDE <- length(which(realDE==1))
        aoua <- array(data = NA, dim =c(l,3))
        iter <- 0
        for (alpha in myCut){
                iter <- iter + 1 

                K <- dim(p)[1]
                myList <-  1 - orderedP[1]
                k <- 1 
                criterion <- myList
                while ((criterion < alpha) & (k < length(orderedP))){
                        k <- k + 1 
                        myList <- myList + 1 - orderedP[k]
                        criterion <- myList/k
                }
                if(k > 1){
                        ind <- perm[1:(k-1)]
                }else{
                        ind <- c()
                }
                if (dim(table(realDE[ind,1])) > 1){
                        point1 <- as.numeric(table(realDE[ind,1])[1]/length(ind))  #achieved fdr
                        point2 <- as.numeric(table(realDE[ind,1])[2]/nDE)  #achieved tpr
                        if(length(nDE) == 0){
                                point2 = 0
                        }
                }else{
                        point1 <- 0
                        if(nDE>0){
                        point2 <- as.numeric(length(ind)/nDE)  #achieved tpr
                        }else{
                        point2 = 0
                        }
                }
                if(sum(realDE) == 0){
                        point1 <- min(length(ind), 1)
                }
                aoua[iter,] <- c(point1,point2,alpha)

        }
        return(aoua)

}



scenario <- vector("list", length = 14)
scenario[[1]] <- c(1, 1.5, 0.8, 0.8, 1.5, 1.5, -0.8, 0.85)
scenario[[3]] <- c(1,1,0.5,0.5,-0.8,1.5,1.5, 0.5)
scenario[[4]] <- c(1,1,0.5,0.5,-0.8,1.5,1.5, 2.5)
scenario[[7]] <- c(-0.05,1,0.8,1,2,-1,1, 0.45)
scenario[[8]] <- c(-0.05,1,0.8,1,2,-1,1, 1.2)
scenario[[9]] <- c(-0.5,1,0.8,1,2,-0.7,1, 0.2)
scenario[[10]] <- c(-0.5,1,0.8,1,2,-0.7,1, 0.5)

n <- 1000

# simulate a gamma model

scenario_id <- 9
truePars <- scenario[[scenario_id]][1:7]
ab = scenario[[scenario_id]][8]


source('sim_data_function_new_weibull_gamma.R')
sim <- sim_fotis_model(n = n, truePars = truePars, ab = ab, seed = 1, discrcov = 1, promodis = 2)

# simulate a Weibull model
source('sim_data_function_new_weibull_gamma.R')
sim <- sim_fotis_model(n = n, truePars = truePars, ab = ab, seed = 1, discrcov = 1, promodis = 1)


> table(sim[,2])

  0   1 
322 678 
> table(sim[,5])

  0   1 
818 182 

# simulate a gamma mixture model
F_parameters <- list(
        K = 2,
        p = c(0.5,0.5),
        gamma_parameters = matrix(c(1,1,6,1), 2, 2, byrow =TRUE)
)

#
source('sim_data_mix.R')
sim <- sim_fotis_model(n = n, lambda1 = 6, F_parameters = F_parameters, ab = 0.05, ranunL=0, ranunU=1, seed = 1, discrcov = 2)

# simulate an exponential model
F_parameters <- list(
        K = 2,
        p = c(0.999999,0.000001),
        gamma_parameters = matrix(c(1,1,6,1), 2, 2, byrow =TRUE)
)

#
source('sim_data_mix.R')
sim <- sim_fotis_model(n = n, gamma1 = 1, lambda1 = 1, F_parameters = F_parameters, ab = 1, ranunL=0, ranunU=1, seed = 1, discrcov = 2)

###################################################
#simulate from a Weibull  model			!!!
###################################################
y <- sim[,1]
X1 <- cbind(1, sim[,3:4])
stat <- sim[,2]



nChains <- 8
mcmcIterations <- 3000
alpha <- 1/1.00001^{(1:nChains)^2.5 - 1}


# run an exponential model
promotion_time = list(distribution = 'exponential', 
                prior_parameters = matrix(rep(c(2.1, 1.1), 1), byrow = TRUE, 1, 2),
                prop_scale = c(0.1)
        )


run_exp <- cure_rate_MC3_general( 
        y = y, 
        X = X1,
        Censoring_status = stat,                
        nChains = nChains,              # number of heated chains (denoted as `C` in the paper)
        mcmc_cycles = mcmcIterations,   # number of MCMC cycles (denoted as `N` in the paper)
        nCores = 8,                     # number of cores to use for parallel processing
        sweep = 10,                     # number of iterations per MCMC cycle (denoted  as `m_1` in the paper)
        alpha = alpha,                  # a decreasing sequence of temperatures per chain 
#                                               (starting from 1 - the target chain is the first one). 
#
#                               PRIOR PARAMETERS        
#######################################################################################################
        mu_g = 1, s2_g = 1,     # this denotes a_g and b_g in the paper: gamma ~ f_g(a_g, b_g
        a_l = 2.1, b_l = 1.1,   # mean = 1, var = 10    lambda ~ IG(a_l, b_l)
        promotion_time = promotion_time,
########################################################################################################                
        plot = TRUE,            # opens plot device for plotting on the run the MCMC output
        adjust_scales = TRUE,   # uses an initial phase in order to adjust the MH scale per parameter
        tau_mala = 0.0001,      # initial value of the MALA scale parameter - it will be adjusted as well.
        mala = 0,
        single_MH_in_f = 0.5
)
burn = 200
getSummary_exp <- compute_map_and_hdis_general(y = y, X = X1, Censoring_status = stat, 
        retained_mcmc = run_exp$mcmc_sample[-(1:burn), ], prior_parameters = NULL, promotion_time = promotion_time)





# run a gamma model (including const. term)
promotion_time = list(distribution = 'gamma', 
                prior_parameters = matrix(rep(c(2.1, 1.1), 2), byrow = TRUE, 2, 2),
                prop_scale = c(0.1, 0.2)
        )

run_gamma <- cure_rate_MC3_general( 
        y = y, 
        X = X1,
        Censoring_status = stat,                
        nChains = nChains,              # number of heated chains (denoted as `C` in the paper)
        mcmc_cycles = mcmcIterations,   # number of MCMC cycles (denoted as `N` in the paper)
        nCores = 8,                     # number of cores to use for parallel processing
        sweep = 10,                     # number of iterations per MCMC cycle (denoted  as `m_1` in the paper)
        alpha = alpha,                  # a decreasing sequence of temperatures per chain 
#                                               (starting from 1 - the target chain is the first one). 
#
#                               PRIOR PARAMETERS        
#######################################################################################################
        mu_g = 1, s2_g = 1,     # this denotes a_g and b_g in the paper: gamma ~ f_g(a_g, b_g
        a_l = 2.1, b_l = 1.1,   # mean = 1, var = 10    lambda ~ IG(a_l, b_l)
        promotion_time = promotion_time,
########################################################################################################                
        plot = TRUE,            # opens plot device for plotting on the run the MCMC output
        adjust_scales = TRUE,   # uses an initial phase in order to adjust the MH scale per parameter
        tau_mala = 0.0001,      # initial value of the MALA scale parameter - it will be adjusted as well.
        mala = 0,
        single_MH_in_f = 0.5
)
burn = 200
getSummary_gamma <- compute_map_and_hdis_general(y = y, X = X1, Censoring_status = stat, 
        retained_mcmc = run_gamma$mcmc_sample[-(1:burn), ], prior_parameters = NULL, promotion_time = promotion_time)


# run a weibull model (including const. term)

promotion_time = list(distribution = 'weibull', 
                prior_parameters = matrix(rep(c(2.1, 1.1), 2), byrow = TRUE, 2, 2),
                prop_scale = c(0.1, 0.2)
        )


run_weibull <- cure_rate_MC3_general( 
        y = y, 
        X = X1,
        Censoring_status = stat,                
        nChains = nChains,              # number of heated chains (denoted as `C` in the paper)
        mcmc_cycles = mcmcIterations,   # number of MCMC cycles (denoted as `N` in the paper)
        nCores = 8,                     # number of cores to use for parallel processing
        sweep = 10,                     # number of iterations per MCMC cycle (denoted  as `m_1` in the paper)
        alpha = alpha,                  # a decreasing sequence of temperatures per chain 
#                                               (starting from 1 - the target chain is the first one). 
#
#                               PRIOR PARAMETERS        
#######################################################################################################
        mu_g = 1, s2_g = 1,     # this denotes a_g and b_g in the paper: gamma ~ f_g(a_g, b_g
        a_l = 2.1, b_l = 1.1,   # mean = 1, var = 10    lambda ~ IG(a_l, b_l)
        promotion_time = promotion_time,
########################################################################################################                
        plot = TRUE,            # opens plot device for plotting on the run the MCMC output
        adjust_scales = TRUE,   # uses an initial phase in order to adjust the MH scale per parameter
        tau_mala = 0.0001,      # initial value of the MALA scale parameter - it will be adjusted as well.
        mala = 0, single_MH_in_f = 0.5
)

if(1 >2 ){
	run_weibull_0 <- cure_rate_MC3( 
		y = y, 
		X = X1,
		Censoring_status = stat,                
		nChains = nChains,              # number of heated chains (denoted as `C` in the paper)
		mcmc_cycles = mcmcIterations,   # number of MCMC cycles (denoted as `N` in the paper)
		nCores = 8,                     # number of cores to use for parallel processing
		sweep = 10,                     # number of iterations per MCMC cycle (denoted  as `m_1` in the paper)
		alpha = alpha,                  # a decreasing sequence of temperatures per chain 
	#                                               (starting from 1 - the target chain is the first one). 
	#
	#                               PRIOR PARAMETERS        
	#######################################################################################################
		mu_g = 1, s2_g = 1, a_l = 2.1, b_l = 1.1, 
		a_1 = 2.1, b_1 = 1.1, a_2 = 2.1, b_2 = 1.1, 
		mu_b = rep(0,dim(X1)[2]), Sigma = 100*diag(dim(X1)[2]),
		g_prop_sd = 0.045, 
		lambda_prop_scale = 0.03, 
		a1_prop_scale = 0.2, 
		a2_prop_scale = 0.03, 
		b_prop_sd = rep(0.022, dim(X1)[2]), 
	########################################################################################################                
		plot = TRUE,            # opens plot device for plotting on the run the MCMC output
		adjust_scales = TRUE,   # uses an initial phase in order to adjust the MH scale per parameter
		tau_mala = 0.000015,      # initial value of the MALA scale parameter - it will be adjusted as well.
		mala = 0
	)
}

burn = 200
getSummary_weibull <- compute_map_and_hdis(y = y, X = X1, Censoring_status = stat, 
        retained_mcmc = run_weibull$mcmc_sample[-(1:burn), ], prior_parameters = NULL)

#getSummary_weibull0 <- compute_map_and_hdis(y = y, X = X1, Censoring_status = stat, 
 #       retained_mcmc = run_weibull_0$mcmc_sample[-(1:burn), ], prior_parameters = NULL)


# run a gamma mixture model


K = 2
prior_parameters = array(data = NA, dim = c(2,2,K))
prior_parameters[,,1] = matrix(rep(c(2.1, 1.1), 2), byrow = TRUE, 2, 2)
prior_parameters[,,2] = matrix(rep(c(2.1, 1.1), 2), byrow = TRUE, 2, 2)

promotion_time = list(distribution = 'gamma_mixture',
                K = 2, 
                prior_parameters = prior_parameters,
                prop_scale = c(0.1, 0.1, 0.1, 0.1, 0.1 ),# a11, a21, ..., a1K, a2K, w_1, ...w_K
                dirichlet_concentration_parameter = 1
        )

run_gamma_mix2 <- cure_rate_MC3_general( 
        y = y, 
        X = X1,
        Censoring_status = stat,                
        nChains = nChains,              # number of heated chains (denoted as `C` in the paper)
        mcmc_cycles = mcmcIterations,   # number of MCMC cycles (denoted as `N` in the paper)
        nCores = 8,                     # number of cores to use for parallel processing
        sweep = 10,                     # number of iterations per MCMC cycle (denoted  as `m_1` in the paper)
        alpha = alpha,                  # a decreasing sequence of temperatures per chain 
#                                               (starting from 1 - the target chain is the first one). 
#
#                               PRIOR PARAMETERS        
#######################################################################################################
        mu_g = 1, s2_g = 1,     # this denotes a_g and b_g in the paper: gamma ~ f_g(a_g, b_g
        a_l = 2.1, b_l = 1.1,   # mean = 1, var = 10    lambda ~ IG(a_l, b_l)
        promotion_time = promotion_time,
########################################################################################################                
        plot = TRUE,            # opens plot device for plotting on the run the MCMC output
        adjust_scales = TRUE,  # uses an initial phase in order to adjust the MH scale per parameter
        tau_mala = 0.0001,      # initial value of the MALA scale parameter - it will be adjusted as well.
        mala = 0, 
        single_MH_in_f = 0.5
)

getSummary_gamma_mix2 <- compute_map_and_hdis_general(y = y, X = X1, Censoring_status = stat, 
        retained_mcmc = run_gamma_mix2$mcmc_sample[-(1:burn), ], prior_parameters = NULL, promotion_time = promotion_time)


# run a gamma mixture model

K = 3
prior_parameters = array(data = NA, dim = c(2,2,K))
prior_parameters[,,1] = matrix(rep(c(2.1, 1.1), 2), byrow = TRUE, 2, 2)
prior_parameters[,,2] = matrix(rep(c(2.1, 1.1), 2), byrow = TRUE, 2, 2)
prior_parameters[,,3] = matrix(rep(c(2.1, 1.1), 2), byrow = TRUE, 2, 2)



promotion_time = list(distribution = 'gamma_mixture',
                K = 3, 
                prior_parameters = prior_parameters,
                prop_scale = rep(0.1, 8),# a11, a21, ..., a1K, a2K, w_1, ...w_K
                dirichlet_concentration_parameter = 1
        )




run_gamma_mix3 <- cure_rate_MC3_general( 
        y = y, 
        X = X1,
        Censoring_status = stat,                
        nChains = nChains,              # number of heated chains (denoted as `C` in the paper)
        mcmc_cycles = mcmcIterations,     # number of MCMC cycles (denoted as `N` in the paper)
        nCores = 8,                     # number of cores to use for parallel processing
        sweep = 10,                     # number of iterations per MCMC cycle (denoted  as `m_1` in the paper)
        alpha = alpha,                  # a decreasing sequence of temperatures per chain 
#                                               (starting from 1 - the target chain is the first one). 
#
#                               PRIOR PARAMETERS        
#######################################################################################################
        mu_g = 1, s2_g = 1,     # this denotes a_g and b_g in the paper: gamma ~ f_g(a_g, b_g
        a_l = 2.1, b_l = 1.1,   # mean = 1, var = 10    lambda ~ IG(a_l, b_l)
        promotion_time = promotion_time,
########################################################################################################                
        plot = TRUE,            # opens plot device for plotting on the run the MCMC output
        adjust_scales = TRUE,  # uses an initial phase in order to adjust the MH scale per parameter
        tau_mala = 0.0001,      # initial value of the MALA scale parameter - it will be adjusted as well.
        mala = 0, 
        single_MH_in_f = 0.5
)

getSummary_gamma_mix3 <- compute_map_and_hdis_general(y = y, X = X1, Censoring_status = stat, 
        retained_mcmc = run_gamma_mix3$mcmc_sample[-(1:burn), ], prior_parameters = NULL, promotion_time = promotion_time)


	bics <- numeric(5)
	bics[1] <- getSummary_exp$bic
	bics[2] <- getSummary_weibull$bic
	bics[3] <- getSummary_gamma$bic
	bics[4] <- getSummary_gamma_mix2$bic
	bics[5] <- getSummary_gamma_mix3$bic
	names(bics) = c('Exponential', 'Weibull', 'Gamma', 'Gamma_mix_2',  'Gamma_mix_3')


	selected_model <- names(which.min(bics))

	if(selected_model == 'Exponential'){run <- run_exp}
	if(selected_model == 'Weibull'){run <- run_weibull}
	if(selected_model == 'Gamma'){run <- run_gamma}
	if(selected_model == 'Gamma_mix_2'){run <- run_gamma_mix2}
	if(selected_model == 'Gamma_mix_3'){run <- run_gamma_mix3}	

	latent_cured_status_exp <- 1-colMeans(run_exp$latent_status_censored)
	latent_cured_status_weibull <- 1-colMeans(run_weibull$latent_status_censored)
	latent_cured_status_gamma <- 1-colMeans(run_gamma$latent_status_censored)	
	latent_cured_status_gamma_mix2 <- 1-colMeans(run_gamma_mix2$latent_status_censored)
	latent_cured_status_gamma_mix3 <- 1-colMeans(run_gamma_mix3$latent_status_censored)	

	labels <- sim[sim[,2]==0,'cured_status']
	pred_exp <- prediction(latent_cured_status_exp, labels)
	pred_weibull <- prediction(latent_cured_status_weibull, labels)
	pred_gamma <- prediction(latent_cured_status_gamma, labels)
	pred_gamma_mix2 <- prediction(latent_cured_status_gamma_mix2, labels)			
	pred_gamma_mix3 <- prediction(latent_cured_status_gamma_mix3, labels)				
	perf_exp <- performance(pred_exp, "tpr", "fpr")
	perf_weibull <- performance(pred_weibull, "tpr", "fpr")
	perf_gamma <- performance(pred_gamma, "tpr", "fpr")				
	perf_gamma_mix2 <- performance(pred_gamma_mix2, "tpr", "fpr")
	perf_gamma_mix3 <- performance(pred_gamma_mix3, "tpr", "fpr")


	pdf(file = 'img/model_selection_weibull_n1000.pdf', width = 12, height = 5)
	par(mfrow = c(1,3), mar = c(5,4,2,1))
	plot(bics, ylab = 'BIC', xaxt = 'n', xlab = '')
	axis(side = 1, labels = FALSE)
	text(x = 1:length(bics),
	     y = par("usr")[3] - 0.45,
	     labels = names(bics),
	     xpd = NA,
	     ## Rotate the labels by 35 degrees.
	     srt = 35,
	     cex = 0.8, adj = 1.2)

	plot(perf_exp, lwd = 2, xlim = c(0,0.5), ylim = c(0.5,1))
	plot(perf_weibull, add = TRUE, col = 2, lty = 2, lwd = 2)	
	plot(perf_gamma, add = TRUE, col = 3, lty = 3, lwd = 2)	
	plot(perf_gamma_mix2, add = TRUE, col = 4, lty = 4, lwd = 2)	
	plot(perf_gamma_mix3, add = TRUE, col = 5, lty = 5, lwd = 2)				

	auc_exp <- performance(pred_exp, measure = "auc")@y.values[[1]]
	auc_weibull <- performance(pred_weibull, measure = "auc")@y.values[[1]]
	auc_gamma <- performance(pred_gamma, measure = "auc")@y.values[[1]]	
	auc_gamma_mix2 <- performance(pred_gamma_mix2, measure = "auc")@y.values[[1]]
	auc_gamma_mix3 <- performance(pred_gamma_mix3, measure = "auc")@y.values[[1]]	
	aucs <- c(auc_exp, auc_weibull, auc_gamma, auc_gamma_mix2, auc_gamma_mix3)

	legend("bottomright", paste0(names(bics), ' (AUC: ', round(aucs,3),')'), col = 1:length(bics), lty = 1:length(bics), lwd = 2)




	myCut = c(1,2.5,5,10)/100
	true_latent_status <- 1 - sim[sim[,2]==0,'cured_status']
	def <-  c(0,1,2,5)
	def2 <- c(15, 16, 17, 18)
	plot(c(0,0.15), c(0,1), type = 'n',xlab = "achieved FDR",ylab = "True positive rate")
	abline(v = myCut, lty = 2, col = 'gray')
	
	posterior_probs <- colMeans(run_gamma_mix2$latent_status_censored)
	fdr_tpr <- compute_fdr_tpr(true_latent_status, posterior_probs, myCut = myCut)
	my_pch <- def
	contr <- apply(fdr_tpr[,c(1,3)], 1, diff) > 0	
	my_pch[contr] <- def2[contr]
	points(fdr_tpr, pch = my_pch, lwd = 2, cex = 2, type = 'b', col = 4, lty = 4)

	posterior_probs <- colMeans(run_exp$latent_status_censored)
	fdr_tpr <- compute_fdr_tpr(true_latent_status, posterior_probs, myCut = myCut)
	my_pch <- def
	contr <- apply(fdr_tpr[,c(1,3)], 1, diff) > 0	
	my_pch[contr] <- def2[contr]
	points(fdr_tpr, pch = my_pch, lwd = 2, cex = 2, type = 'b', col = 1, lty = 1)

	posterior_probs <- colMeans(run_weibull$latent_status_censored)
	fdr_tpr <- compute_fdr_tpr(true_latent_status, posterior_probs, myCut = myCut)
	my_pch <- def
	contr <- apply(fdr_tpr[,c(1,3)], 1, diff) > 0	
	my_pch[contr] <- def2[contr]
	points(fdr_tpr, pch = my_pch, lwd = 2, cex = 2, type = 'b', col = 2, lty = 2)

	posterior_probs <- colMeans(run_gamma$latent_status_censored)
	fdr_tpr <- compute_fdr_tpr(true_latent_status, posterior_probs, myCut = myCut)
	my_pch <- def
	contr <- apply(fdr_tpr[,c(1,3)], 1, diff) > 0	
	my_pch[contr] <- def2[contr]
	points(fdr_tpr, pch = my_pch, lwd = 2, cex = 2, type = 'b', col = 3, lty = 3)

	posterior_probs <- colMeans(run_gamma_mix3$latent_status_censored)
	fdr_tpr <- compute_fdr_tpr(true_latent_status, posterior_probs, myCut = myCut)
	my_pch <- def
	contr <- apply(fdr_tpr[,c(1,3)], 1, diff) > 0	
	my_pch[contr] <- def2[contr]
	points(fdr_tpr, pch = my_pch, lwd = 2, cex = 2, type = 'b', col = 5, lty = 5)
	legend('bottomright', 'nominal FDR levels', lty = 2, col = 'gray')
	dev.off()
	
	save.image('model_selection_weibull_n1000.RData')

