
log_gamma_mixture <- function(y, a1, a2, p, c_under = 1e-9){
	c_under <- log(c_under)
	K <- length(a1)
	n <- length(y)
	logp <- log(p)
	if(K != length(p)){stop('length a not equal to length p')}
	if(min(p) < 0){stop('mixing weights should be positive')}
	log_f <- log_F <- numeric(n)
	log_f_k <- array(data = NA, dim = c(n,K))
	log_F_k <- array(data = NA, dim = c(n,K))	
	for (k in 1:K){
		log_f_k[,k] <- logp[k] + dgamma(y, shape = a1[k], rate = a2[k], log = TRUE)
		log_F_k[,k] <- logp[k] + pgamma(y, shape = a1[k], rate = a2[k], log.p = TRUE)
	}
	log_f_max <- apply(log_f_k, 1, max)
	log_F_max <- apply(log_F_k, 1, max)	
	log_f_k <- log_f_k - log_f_max
	log_F_k <- log_F_k - log_F_max	
	log_f <- log_f_max + log(rowSums(exp(log_f_k)))
	log_F <- log_F_max + log(rowSums(exp(log_F_k)))
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





log_logLogistic <- function(y, a1, a2, c_under = 1e-9){
	c_under <- log(c_under)
	log_f <- dllogis(y, shape = a1, scale = a2, log = TRUE)
	log_F <- pllogis(y, shape = a1, scale = a2, log.p = TRUE)
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



log_gompertz <- function(y, a1, a2, c_under = 1e-9){
	c_under <- log(c_under)
	log_f <- dgompertz(y, shape = a1, rate = a2, log = TRUE)
	log_F <- pgompertz(y, shape = a1, rate = a2, log.p = TRUE)
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


log_gamma <- function(y, a1, a2, c_under = 1e-9){
	c_under <- log(c_under)
	log_f <- dgamma(y, shape = a1, rate = a2, log = TRUE)
	log_F <- pgamma(y, shape = a1, rate = a2, log.p = TRUE)
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

log_lomax <- function(y, a1, a2, c_under = 1e-9){
	c_under <- log(c_under)
	log_f <- dlomax(y, scale = a1, shape3.q = a2, log = TRUE)
	log_F <- plomax(y, scale = a1, shape3.q = a2, log.p = TRUE)
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


log_dagum <- function(y, a1, a2, a3, c_under = 1e-9){
	c_under <- log(c_under)
	log_f <- ddagum(y, scale = a1, shape1.a = a2, shape2.p = a3, log = TRUE)
	log_F <- pdagum(y, scale = a1, shape1.a = a2, shape2.p = a3, log.p = TRUE)
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


cure_rate_mcmc <- function( y, X, Censoring_status,  m, alpha = 1,
				mu_g = 1, s2_g = 1, a_l = 2.1, b_l = 1.1, 
				promotion_time = list(distribution = 'weibull', 
						prior_parameters = matrix(rep(c(2.1, 1.1), 2), byrow = TRUE, 2, 2),
						prop_scale = c(0.2, 0.03)
					),
				mu_b = NULL, Sigma = NULL,
				g_prop_sd = 0.045, 
				lambda_prop_scale = 0.03, 
				b_prop_sd = NULL, 
						initialValues = NULL, 
						plot = FALSE,
						verbose = FALSE,
						tau_mala = 0.000015, mala = 0.15, single_MH_in_f = 0.5
					){
# 	prior_parameters should be a matrix with as many rows as the number of parameters in f. 
#		* positive parameters are assigned independent IG priors. 
#		* each row in prior_parameters corresponds to the IG(.,.) prior parameters
#		* only positive prior parameters are considered at the moment. 

	if(plot){
	oldpar <- par(no.readonly = TRUE)
	}

	n_pars_f <- dim(promotion_time$prior_parameters)[1]

	if(promotion_time$distribution == 'exponential'){
		if(n_pars_f != 1){stop("inconsistent number of parameters!")}
	}

	if(promotion_time$distribution == 'weibull'){
		if(n_pars_f != 2){stop("inconsistent number of parameters!")}
	}

	if(promotion_time$distribution == 'gamma'){
		if(n_pars_f != 2){stop("inconsistent number of parameters!")}
	}

	if(promotion_time$distribution == 'logLogistic'){
		if(n_pars_f != 2){stop("inconsistent number of parameters!")}
	}

	if(promotion_time$distribution == 'gompertz'){
		if(n_pars_f != 2){stop("inconsistent number of parameters!")}
	}

	if(promotion_time$distribution == 'lomax'){
		if(n_pars_f != 2){stop("inconsistent number of parameters!")}
	}
	if(promotion_time$distribution == 'dagum'){
		if(n_pars_f != 3){stop("inconsistent number of parameters!")}
	}
	if(promotion_time$distribution == 'gamma_mixture'){
		K = promotion_time$K
		if(dim(promotion_time$prior_parameters)[3] != K){stop("incosistent number of mixture components")}		
		n_pars_f = K * n_pars_f + K - 1 # note that we include all mixing weights
	}
	if(length(promotion_time$prop_scale) != n_pars_f){stop("inconsistent number of proposal scale parameters!")}


	f <- function(par_vec, family, K, I_sim){
	#			par_vec should be gamma, lambda, a1, a2, p, b
	#			a11, a12,..., a1K denotes gamma-shape parameter per component (1, 2, ... ,K)
	#			a21, a22,..., a2K denotes gamma-rate parameter per component (1, 2, ... ,K)
		g = par_vec[1]
		lambda = par_vec[2]	
		b = par_vec[-(1:(2+n_pars_f))]
		if(family == 'weibull'){
			lw <- log_weibull(y, a1 = par_vec[3], a2 = par_vec[4],  c_under = 1e-9)	
		}

		if(family == 'exponential'){
			lw <- log_weibull(y, a1 = par_vec[3], a2 = 1,  c_under = 1e-9)	
		}


		if(family == 'gamma'){
			lw <- log_gamma(y = y, a1 = par_vec[3], a2 = par_vec[4], c_under = 1e-9)	
		}

		if(family == 'logLogistic'){
			lw <- log_logLogistic(y = y, a1 = par_vec[3], a2 = par_vec[4], c_under = 1e-9)	
		}

		if(family == 'gompertz'){
			lw <- log_gompertz(y = y, a1 = par_vec[3], a2 = par_vec[4], c_under = 1e-9)	
		}
		
		if(family == 'gamma_mixture'){
			w = par_vec[(2+K*2 + 1):(2+K*2 + K-1)]
			p = c(w, 1)
			p = p/sum(p)
			lw <- log_gamma_mixture(y, a1 = par_vec[3:(2+K)], a2 = par_vec[(K+3):(2*K+2)], p = p, c_under = 1e-9)
		}

		
		if(family == 'lomax'){
			lw <- log_lomax(y, a1 = par_vec[3], a2 = par_vec[4], c_under = 1e-9)
		}

		
		if(family == 'dagum'){
			lw <- log_dagum(y, a1 = par_vec[3], a2 = par_vec[4], a3 = par_vec[5], c_under = 1e-9)
		}

		
		return(complete_log_likelihood_general(y = y, X = X, Censoring_status = Censoring_status, g = g, lambda = lambda, 
			log_f = lw$log_f, log_F = lw$log_F, b = b, I_sim = I_sim, alpha = alpha)$cll
			)
	}



	# the current value of the parameters of f will be stored to vector a
	a <- numeric(n_pars_f)
	c_max = 1e+9

#	nCov <- dim(myData)[2] - 1
	X <- as.matrix(X)
	nCov <- dim(X)[2]
	nPars <- nCov + 2 + n_pars_f
	if(is.null(mu_b)){mu_b <- rep(0, nCov)}
	if(is.null(Sigma)){Sigma <- 100 * diag(nCov)}
	if(is.null(b_prop_sd)){b_prop_sd = rep(0.022, nCov)}
#	input data: myData 
#	this should be a matrix with the following columns:
#	"Y", "Censoring_status", "Covariate1", "Covariate2"
#	mala: probability of performing a mala move
#		otherwise, usual MH update is performed

#	m: denotes the total number of MCMC iterations (arbitrary for now)

	# prior parameters
	#	gamma ~ 0.5G(a_g, b_g) + 0.5NG(a_g, b_g)
	#		NOTE mu_g = a_g and s2_g = b_g
	#	lambda ~ IG(a_l, b_l)
	# b = (b0, b1, b2) ~ N(mu_b, Sigma) 

	# proposal parameters
#	g_prop_sd <- 0.4
#	lambda_prop_scale <- 0.1
#	a1_prop_scale <- 0.2
#	a2_prop_scale <- 0.1
#	b_prop_sd = rep(0.05, 3)
#	initialValues: should be a list with entries named as g_mcmc, b0_mcmc, b1_mcmc, b2_mcmc, lambda_mcmc, a1_mcmc, a2_mcmc
#			and if it is NULL then random starting values are generated

#	n <- dim(myData)[1]
	n <- dim(X)[1]
#	X <- myData[,-(1:2)]
	ct = exp(exp(-1))

#	D1 <- which(myData[,"Censoring_status"] == 1)	# fixed
#	D0 <- which(myData[,"Censoring_status"] == 0)	# fixed

	D1 <- which(Censoring_status == 1)	# fixed
	D0 <- which(Censoring_status == 0)	# fixed


	g_mcmc <- lambda_mcmc <- numeric(m)
	a_mcmc <- array(data = NA, dim = c(m, n_pars_f))
	b_mcmc <- array(data = NA, dim = c(m, nCov))
	cllValues <- numeric(m)
	lpd <- numeric(m)
	I_sim_values_D0 <- matrix(data = NA, ncol = length(D0), nrow = m)

	if(length(unique(X[,1])) == 1){
		bRange <- 1:nCov - 1
	}else{
		bRange <- 1:nCov 
	}


	# initial values (random)
	if(is.null(initialValues)){
		g = rnorm(1, sd = 2)
		b = rnorm(nCov, sd = 2)
		lambda = rgamma(1, shape = 1, rate = 1) 
		a = rgamma(n_pars_f, shape = 1, rate = 1)		
	}else{
		g = initialValues[['g_mcmc']]
		b <- numeric(nCov)
		for(j in 1:nCov){
			b[j] <- initialValues[[2 + n_pars_f + j]]
		}
		lambda = initialValues[['lambda_mcmc']]								
		a = unlist(initialValues[3:(3+n_pars_f-1)])
	}
	


	g_mcmc[1] <- g
	b_mcmc[1,] <- b
	a_mcmc[1,] <- a
	lambda_mcmc[1] <- lambda


	# latent binary variables (I = 1 => susceptible, I = 0 => cured)
	I_sim <- rep(1, n)
	# for i in D1: I_i = 1
	if(is.null(initialValues)){
		susceptible_prob <- runif(n)
		I_sim[D0] <- rbinom(n = length(D0), size = 1, prob = susceptible_prob[D0])
	}else{
		I_sim[D0] <- initialValues[["I_sim_D0"]]
	}
	I_sim_values_D0[1, ] <- I_sim[D0]

#	compute log dens and cdf of promotion time density
	if(promotion_time$distribution == 'exponential'){
		lw <- log_weibull(y, a1 = a[1], a2 = 1,  c_under = 1e-9)
	}

	if(promotion_time$distribution == 'weibull'){
		lw <- log_weibull(y, a1 = a[1], a2 = a[2],  c_under = 1e-9)
	}

	if(promotion_time$distribution == 'gamma'){
		lw <- log_gamma(y, a1 = a[1], a2 = a[2],  c_under = 1e-9)
	}

	if(promotion_time$distribution == 'logLogistic'){
		lw <- log_logLogistic(y, a1 = a[1], a2 = a[2],  c_under = 1e-9)
	}

	if(promotion_time$distribution == 'gompertz'){
		lw <- log_gompertz(y, a1 = a[1], a2 = a[2],  c_under = 1e-9)
	}


	if(promotion_time$distribution == 'lomax'){
		lw <- log_lomax(y, a1 = a[1], a2 = a[2],  c_under = 1e-9)
	}
	if(promotion_time$distribution == 'dagum'){
		lw <- log_dagum(y, a1 = a[1], a2 = a[2], a3 = a[3], c_under = 1e-9)
	}
	if(promotion_time$distribution == 'gamma_mixture'){
		w <- a[-(1:(2*K))]
		w2 <- c(w, 1)
		p <- w2/sum(w2)
		a1_ind <- seq(1, 2*K, by = 2)
		a2_ind <- seq(2, 2*K, by = 2)		
		a1 = a[a1_ind]
		a2 = a[a2_ind]		
		lw <- log_gamma_mixture(y, a1 = a1, a2 = a2, p = p, c_under = 1e-9)
	}




	cll <- complete_log_likelihood_general(y = y, X = X, Censoring_status = Censoring_status, g = g, lambda = lambda, log_f = lw$log_f, log_F = lw$log_F, b = b, I_sim = I_sim, alpha = alpha)
	cllValues[1] <- cll$cll
	##########################################################edw stamatisa################################################3
#	print("0")
#	print(cllValues[1])
#	print(c(g, lambda, a1, a2, b0, b1, b2))

	log_dirichlet_pdf <- function (alpha, weights){
		min_weight <- 1e-300
		weights[weights < min_weight] <- min_weight
		weights <- weights/sum(weights)
		normConstant <- sum(lgamma(alpha)) - lgamma(sum(alpha))
		pdf <- sum((alpha - 1) * log(weights)) - normConstant
		return(pdf)
	}


	# a: shape, b: scale
	log_inv_gamma_kernel <- function(x, a, b){
		if(min(c(x,a,b)) < 0){
		stop("input should be positive")
		}
		return(-b/x - (a+1) * log(x))
	}
#	print(c(g, lambda, a1, a2, b0, b1, b2))
	
	
#	log_0.5 <- log(0.5)
	log_gamma_mix_kernel <- function(g, a, b){
		return((a-1)*log(abs(g)) - b*abs(g))
#		if(g < 0){
#			return((a-1)*log(-g) + b*g)
#		}else{
#			return((a-1)*log(g) - b*g)
#		}
	}
	
	log_prior_density <- -alpha*log_gamma_mix_kernel(g, mu_g, s2_g) + 
				alpha*log_inv_gamma_kernel(lambda, a_l, b_l) -
				alpha*0.5 * mahalanobis(b, mu_b, Sigma)^2

	if(promotion_time$distribution == 'gamma_mixture'){
		for(k in 1:K){
			log_prior_density <- log_prior_density + sum(alpha*log_inv_gamma_kernel(c(a1[k], a2[k]), 
						promotion_time$prior_parameters[,1,k], promotion_time$prior_parameters[,2,k]))
		}
		a_vec <- rep(promotion_time$dirichlet_concentration_parameter, K)
		log_prior_density <- log_prior_density + log_dirichlet_pdf(a_vec, p)
		
	}else{
		for(i in 1:n_pars_f){
			log_prior_density <- log_prior_density + alpha*log_inv_gamma_kernel(a[i], 
					promotion_time$prior_parameters[i,1], promotion_time$prior_parameters[i,2])
		}
	}
	lpd[1] <- log_prior_density

	mala_accept_rate <- 0
	g_accept_rate <- lambda_accept_rate <- b_accept_rate <- a_all_accept_rate <- 0
	a_accept_rate <- numeric(n_pars_f)



	if(plot){
	on.exit(par(oldpar)) 
	dev.new(width=18, height=8, unit="in")	
	}
	######################################
	# MCMC sampler
	######################################
	if(m > 1){
	mh_iter <- 0
	mala_iter <- 0
	mh_single <- 0
	mh_all <- 0
	for(iter in 2:m){

		joe <- runif(1)
		if(joe < 1 - mala){
	   #------------------------------------------------
		#       Metropolis-Hastings proposal for gamma: normal proposal distribution
		#------------------------------------------------
		        g_prop <- rnorm(1, mean = g, sd = g_prop_sd)

		        cll_prop <- complete_log_likelihood_general(y = y, X = X, Censoring_status = Censoring_status, 
		        		g = g_prop, lambda = lambda, log_f = lw$log_f, log_F = lw$log_F, 
		        		b = b, I_sim = I_sim, alpha = alpha)
		        cll_diff <- cll_prop$cll - cll$cll
	#               print("1")
	#               print(cll_diff)
		        log_prior_diff <- alpha*diff(log_gamma_mix_kernel(c(g, g_prop), mu_g, s2_g))
		        log_proposal_ratio <- 0
		        mh_ar <- cll_diff + log_prior_diff + log_proposal_ratio
		#       print("gamma")
		#       print(cll_prop)
		        if(is.na(mh_ar)==FALSE){
		        if(log(runif(1)) < mh_ar){
		                g = g_prop
		                cll <- cll_prop
		                g_accept_rate <- g_accept_rate + 1
		        }
		        }
		        g_mcmc[iter] <- g

		#------------------------------------------------
		#       Metropolis-Hastings proposal for lambda: lognormal proposal distribution
		#------------------------------------------------
		        lambda_prop = rlnorm(1, meanlog = log(lambda), sdlog = lambda_prop_scale)
		        cll_prop <- complete_log_likelihood_general(y = y, X = X, Censoring_status = Censoring_status, 
		        	g = g, lambda = lambda_prop, log_f = lw$log_f, log_F = lw$log_F,
		        	b = b, I_sim = I_sim, alpha = alpha)
		        cll_diff <- cll_prop$cll - cll$cll
	#               print("2")
	#               print(cll_diff)         
		        log_prior_diff <- alpha*diff(log_inv_gamma_kernel(c(lambda, lambda_prop), a_l, b_l))
		        log_proposal_ratio <- log(lambda_prop) - log(lambda)
		        mh_ar <- cll_diff + log_prior_diff + log_proposal_ratio
		#       print("lambda")
		#       print(cll_prop)

		        if(is.na(mh_ar)==FALSE){
		        if(log(runif(1)) < mh_ar){
		                lambda = lambda_prop
		                cll <- cll_prop
		                lambda_accept_rate <- lambda_accept_rate + 1
		        }
		        }
		        lambda_mcmc[iter] <- lambda
		

		u1 <- runif(1)
		if(u1 < single_MH_in_f){
		mh_single <- mh_single + 1
	
		for(i in 1:n_pars_f){
		#------------------------------------------------
		#       Metropolis-Hastings proposal for a[i]: lognormal proposal distribution
		#------------------------------------------------
		
			a_prop <- a
			a_prop[i] = rlnorm(1, meanlog = log(a[i]), sdlog = promotion_time$prop_scale[i])
			a_prop[i] = min(c(a_prop[i],c_max))
			
			if(promotion_time$distribution == 'exponential'){
				lw_prop <- log_weibull(y, a1 = a_prop[1], a2 = 1,  c_under = 1e-9)
				log_prior_diff <- alpha*diff(log_inv_gamma_kernel(c(a[i], a_prop[i]), 
					promotion_time$prior_parameters[i,1], promotion_time$prior_parameters[i,2]))
			}


			if(promotion_time$distribution == 'weibull'){
				lw_prop <- log_weibull(y, a1 = a_prop[1], a2 = a_prop[2],  c_under = 1e-9)
				log_prior_diff <- alpha*diff(log_inv_gamma_kernel(c(a[i], a_prop[i]), 
					promotion_time$prior_parameters[i,1], promotion_time$prior_parameters[i,2]))				
			}

			if(promotion_time$distribution == 'gamma'){
				lw_prop <- log_gamma(y, a1 = a_prop[1], a2 = a_prop[2],  c_under = 1e-9)
				log_prior_diff <-  alpha*diff(log_inv_gamma_kernel(c(a[i], a_prop[i]), 
					promotion_time$prior_parameters[i,1], promotion_time$prior_parameters[i,2]))				
			}

			if(promotion_time$distribution == 'gompertz'){
				lw_prop <- log_gompertz(y, a1 = a_prop[1], a2 = a_prop[2],  c_under = 1e-9)
				log_prior_diff <-  alpha*diff(log_inv_gamma_kernel(c(a[i], a_prop[i]), 
					promotion_time$prior_parameters[i,1], promotion_time$prior_parameters[i,2]))				
			}


			if(promotion_time$distribution == 'logLogistic'){
				lw_prop <- log_logLogistic(y, a1 = a_prop[1], a2 = a_prop[2],  c_under = 1e-9)
				log_prior_diff <-  alpha*diff(log_inv_gamma_kernel(c(a[i], a_prop[i]), 
					promotion_time$prior_parameters[i,1], promotion_time$prior_parameters[i,2]))				
			}



			if(promotion_time$distribution == 'lomax'){
				lw_prop <- log_lomax(y, a1 = a_prop[1], a2 = a_prop[2],  c_under = 1e-9)
				log_prior_diff <-  alpha*diff(log_inv_gamma_kernel(c(a[i], a_prop[i]), 
					promotion_time$prior_parameters[i,1], promotion_time$prior_parameters[i,2]))				
			}
			if(promotion_time$distribution == 'dagum'){
				lw_prop <- log_dagum(y, a1 = a_prop[1], a2 = a_prop[2],  a3 = a_prop[3], c_under = 1e-9)
				log_prior_diff <-  alpha*diff(log_inv_gamma_kernel(c(a[i], a_prop[i]), 
					promotion_time$prior_parameters[i,1], promotion_time$prior_parameters[i,2]))				
			}
			if(promotion_time$distribution == 'gamma_mixture'){
				w_prop <- a_prop[-(1:(2*K))]
				w_prop2 <- c(w_prop, 1)
				p_prop <- w_prop2/sum(w_prop2)
				a1_prop = a_prop[a1_ind]
				a2_prop = a_prop[a2_ind]		
				lw_prop <- log_gamma_mixture(y, a1 = a1_prop, a2 = a2_prop, p = p_prop, c_under = 1e-9)
				if(i < 2*K + 1){
					k <- which(a1_ind == i)
					if(length(k) == 0){
						k <- which(a2_ind == i)
					}
#					k <- floor((1:(2*K))/K + 0.5)[i]
					j = (i+1)%%2 + 1
					log_prior_diff <- alpha*diff(log_inv_gamma_kernel(c(a[i], a_prop[i]), 
						promotion_time$prior_parameters[j,1,k], promotion_time$prior_parameters[j,2,k]))
				}else{
					log_prior_diff <- log_dirichlet_pdf(a_vec, p_prop) - log_dirichlet_pdf(a_vec, p)
				}

			}

		
			cll_prop <- complete_log_likelihood_general(y = y, X = X, Censoring_status = Censoring_status, 
				g = g, lambda = lambda, log_f = lw_prop$log_f,  log_F = lw_prop$log_F,
				b = b, I_sim = I_sim, alpha = alpha)
			cll_diff <- cll_prop$cll - cll$cll
			#               print("3")
			#               print(cll_diff)

			log_proposal_ratio <- log(a_prop[i]) - log(a[i])
			mh_ar <- cll_diff + log_prior_diff + log_proposal_ratio
			#       print("a1")
			#       print(cll_prop)
			if(is.na(mh_ar)==FALSE){
			if(log(runif(1)) < mh_ar){
				a[i] = a_prop[i]
				if(promotion_time$distribution == 'gamma_mixture'){
					w = w_prop
					p = p_prop
				}
				cll <- cll_prop
				a_accept_rate[i] <- a_accept_rate[i] + 1
				lw <- lw_prop
			}
			}
			a_mcmc[iter,i] <- a[i]
		}
		}else{
			mh_all <- mh_all + 1
############################ joint update
		
			a_prop = rlnorm(n_pars_f, meanlog = log(a), sdlog = promotion_time$prop_scale)
			for(i in 1:n_pars_f){
				a_prop[i] = min(c(a_prop[i],c_max))
			}
			if(promotion_time$distribution == 'exponential'){
				lw_prop <- log_weibull(y, a1 = a_prop[1], a2 = 1,  c_under = 1e-9)
				log_prior_diff <- alpha*diff(log_inv_gamma_kernel(c(a[1], a_prop[1]), 
					promotion_time$prior_parameters[1,1], promotion_time$prior_parameters[1,2]))
			}


			if(promotion_time$distribution == 'weibull'){
				lw_prop <- log_weibull(y, a1 = a_prop[1], a2 = a_prop[2],  c_under = 1e-9)
				log_prior_diff <- 0
				for(i in 1:n_pars_f){
				log_prior_diff <- log_prior_diff + alpha*diff(log_inv_gamma_kernel(c(a[i], a_prop[i]), 
					promotion_time$prior_parameters[i,1], promotion_time$prior_parameters[i,2]))				
				}
			}

			if(promotion_time$distribution == 'gamma'){
				lw_prop <- log_gamma(y, a1 = a_prop[1], a2 = a_prop[2],  c_under = 1e-9)
				log_prior_diff <- 0
				for(i in 1:n_pars_f){
				log_prior_diff <- log_prior_diff + alpha*diff(log_inv_gamma_kernel(c(a[i], a_prop[i]), 
					promotion_time$prior_parameters[i,1], promotion_time$prior_parameters[i,2]))
				}
			}

			if(promotion_time$distribution == 'gompertz'){
				lw_prop <- log_gompertz(y, a1 = a_prop[1], a2 = a_prop[2],  c_under = 1e-9)
				log_prior_diff <- 0
				for(i in 1:n_pars_f){
				log_prior_diff <- log_prior_diff + alpha*diff(log_inv_gamma_kernel(c(a[i], a_prop[i]), 
					promotion_time$prior_parameters[i,1], promotion_time$prior_parameters[i,2]))
				}
			}



			if(promotion_time$distribution == 'logLogistic'){
				lw_prop <- log_logLogistic(y, a1 = a_prop[1], a2 = a_prop[2],  c_under = 1e-9)
				log_prior_diff <- 0
				for(i in 1:n_pars_f){
				log_prior_diff <- log_prior_diff + alpha*diff(log_inv_gamma_kernel(c(a[i], a_prop[i]), 
					promotion_time$prior_parameters[i,1], promotion_time$prior_parameters[i,2]))
				}
			}



			if(promotion_time$distribution == 'lomax'){
				lw_prop <- log_lomax(y, a1 = a_prop[1], a2 = a_prop[2],  c_under = 1e-9)
				log_prior_diff <- 0
				for(i in 1:n_pars_f){
				log_prior_diff <- log_prior_diff + alpha*diff(log_inv_gamma_kernel(c(a[i], a_prop[i]), 
					promotion_time$prior_parameters[i,1], promotion_time$prior_parameters[i,2]))				
				}
			}
			if(promotion_time$distribution == 'dagum'){
				lw_prop <- log_dagum(y, a1 = a_prop[1], a2 = a_prop[2],  a3 = a_prop[3], c_under = 1e-9)
				log_prior_diff <- 0				
				for(i in 1:n_pars_f){
				log_prior_diff <- log_prior_diff +  alpha*diff(log_inv_gamma_kernel(c(a[i], a_prop[i]), 
					promotion_time$prior_parameters[i,1], promotion_time$prior_parameters[i,2]))	
				}			
			}
			if(promotion_time$distribution == 'gamma_mixture'){
				w_prop <- a_prop[-(1:(2*K))]
				w_prop2 <- c(w_prop, 1)
				p_prop <- w_prop2/sum(w_prop2)
				a1_prop = a_prop[a1_ind]
				a2_prop = a_prop[a2_ind]		
				lw_prop <- log_gamma_mixture(y, a1 = a1_prop, a2 = a2_prop, p = p_prop, c_under = 1e-9)
				log_prior_diff <- 0				
				for(i in 1:(2*K)){
					k <- which(a1_ind == i)
					if(length(k) == 0){
						k <- which(a2_ind == i)
					}
					j = (i+1)%%2 + 1
#					k <- floor((1:(2*K))/K + 0.5)[i]
#					j = i%%2 + 1
					log_prior_diff <- log_prior_diff + alpha*diff(log_inv_gamma_kernel(c(a[i], a_prop[i]), 
						promotion_time$prior_parameters[j,1,k], promotion_time$prior_parameters[j,2,k]))
				}
				log_prior_diff <-  log_prior_diff + alpha*log_dirichlet_pdf(a_vec, p_prop) - alpha*log_dirichlet_pdf(a_vec, p)

			}
			cll_prop <- complete_log_likelihood_general(y = y, X = X, Censoring_status = Censoring_status, 
				g = g, lambda = lambda, log_f = lw_prop$log_f,  log_F = lw_prop$log_F,
				b = b, I_sim = I_sim, alpha = alpha)
			cll_diff <- cll_prop$cll - cll$cll
			log_proposal_ratio <- 0
			for(i in 1:n_pars_f){
				log_proposal_ratio <- log_proposal_ratio + log(a_prop[i]) - log(a[i])
			}
			mh_ar <- cll_diff + log_prior_diff + log_proposal_ratio
			#       print("a1")
			#       print(cll_prop)
			if(is.na(mh_ar)==FALSE){
			if(log(runif(1)) < mh_ar){
				a = a_prop
				if(promotion_time$distribution == 'gamma_mixture'){
					w = w_prop
					p = p_prop
				}
				cll <- cll_prop
				a_all_accept_rate <- a_all_accept_rate + 1
				lw <- lw_prop
			}
			}
			a_mcmc[iter,] <- a
			
		}		


	       #------------------------------------------------
		#       Metropolis-Hastings proposal for b = (b0,b1,b2): multivariate normal proposal distribution
		#------------------------------------------------
		        b_prop <- b + b_prop_sd*rnorm(nCov)
		        cll_prop <- complete_log_likelihood_general(y = y, X = X, Censoring_status = Censoring_status, 
		        	g = g, lambda = lambda, log_f = lw$log_f,  log_F = lw$log_F,
		        	b= b_prop, I_sim = I_sim, alpha = alpha)
		        cll_diff <- cll_prop$cll - cll$cll
	#               print("5")
	#               print(cll_diff)
		        
		        log_prior_diff <- -alpha*0.5 * mahalanobis(b_prop, mu_b, Sigma)^2 + alpha*0.5 * mahalanobis(b, mu_b, Sigma)^2
		        log_proposal_ratio <- 0
		        mh_ar <- cll_diff + log_prior_diff + log_proposal_ratio
		#       print("b")
		#       print(cll_prop)
		        if(is.na(mh_ar)==FALSE){
		        if(log(runif(1)) < mh_ar){
		                b = b_prop
		                cll <- cll_prop
		                b_accept_rate <- b_accept_rate + 1
		        }
		        }
		        b_mcmc[iter,] <- b
		
			mh_iter <- mh_iter + 1
		}else{

			theta_previous <- as.numeric(c(g, lambda, a, b)	)


			gradient_vec <- derivative(f, var = theta_previous, 
					params = list(family = promotion_time$distribution, 
					K = promotion_time$K,
					I_sim = I_sim)
				)    
			#derv_complete_log_likelihood_fotis(y = y, X = X, Censoring_status = Censoring_status,
			#	theta_previous[1],theta_previous[2],theta_previous[3],
			#	theta_previous[4],theta_previous[-(1:4)],I_sim)

			# ftiaxe kai auto!
			lpg <- 0 #log_prior_gradient(g, lambda, a1, a2, b , mu_g, s2_g,  a_l, b_l, a_1, b_1, a_2, b_2, mu_b, Sigma)

			gradient_vec <- alpha*gradient_vec + alpha*lpg

			
			mean_prop <- theta_previous + tau_mala*gradient_vec
			sd_prop <- sqrt(2*tau_mala)
			theta_prop <- rnorm(nPars, mean_prop, sd_prop)
			rejectedMove = FALSE
			if(sum(is.na(theta_prop))){rejectedMove = TRUE}else{
			if(min(theta_prop[2:(2+n_pars_f)]) < 0){rejectedMove = TRUE}}
			if(rejectedMove == FALSE){
				log_prop_density <- sum(dnorm(theta_prop, mean_prop, sd = sd_prop, log = TRUE))
				gradient_vec_prop <- derivative(f, var = theta_prop, 
					params = list(family = promotion_time$distribution, 
					K = promotion_time$K,
					I_sim = I_sim)
				)    
				#derv_complete_log_likelihood_fotis(y = y, X = X, 
				#	Censoring_status = Censoring_status, 
				#	theta_prop[1],theta_prop[2],theta_prop[3], theta_prop[4],
				#	theta_prop[-(1:4)],I_sim)
			# ftiaxto
				lpg <- 0 #log_prior_gradient(theta_prop[1],theta_prop[2],theta_prop[3],
					#theta_prop[4], theta_prop[-(1:4)] , 
					#mu_g, s2_g,  a_l, b_l, a_1, b_1, a_2, b_2, mu_b, Sigma)

				gradient_vec_prop <- alpha*gradient_vec_prop + alpha*lpg
				
				mean_prop_rev <- theta_prop + tau_mala * gradient_vec_prop
				log_prop_density_reverse <- sum(dnorm(theta_previous, mean_prop_rev, sd = sd_prop, log = TRUE))
				b_prop = theta_prop[-(1:(2+n_pars_f ))]

				a_prop <- theta_prop[(3:(2+n_pars_f))]
				for(i in 1:n_pars_f){
					a_prop[i] = min(c(a_prop[i],c_max))
				}
				if(promotion_time$distribution == 'exponential'){
					lw_prop <- log_weibull(y, a1 = a_prop[1], a2 = 1,  c_under = 1e-9)
					log_prior_diff <- alpha*diff(log_inv_gamma_kernel(c(a[1], a_prop[1]), 
						promotion_time$prior_parameters[1,1], promotion_time$prior_parameters[1,2]))					
				}
				if(promotion_time$distribution == 'weibull'){
					lw_prop <- log_weibull(y, a1 = a_prop[1], a2 = a_prop[2],  c_under = 1e-9)
					log_prior_diff <- 0
					for(i in 1:n_pars_f){
					log_prior_diff <- log_prior_diff + alpha*diff(log_inv_gamma_kernel(c(a[i], a_prop[i]), 
						promotion_time$prior_parameters[i,1], promotion_time$prior_parameters[i,2]))				
					}
				}

				if(promotion_time$distribution == 'gamma'){
					lw_prop <- log_gamma(y, a1 = a_prop[1], a2 = a_prop[2],  c_under = 1e-9)
					log_prior_diff <- 0
					for(i in 1:n_pars_f){
					log_prior_diff <- log_prior_diff + alpha*diff(log_inv_gamma_kernel(c(a[i], a_prop[i]), 
						promotion_time$prior_parameters[i,1], promotion_time$prior_parameters[i,2]))
					}					
				}

				if(promotion_time$distribution == 'gompertz'){
					lw_prop <- log_gompertz(y, a1 = a_prop[1], a2 = a_prop[2],  c_under = 1e-9)
					log_prior_diff <- 0
					for(i in 1:n_pars_f){
					log_prior_diff <- log_prior_diff + alpha*diff(log_inv_gamma_kernel(c(a[i], a_prop[i]), 
						promotion_time$prior_parameters[i,1], promotion_time$prior_parameters[i,2]))
					}					
				}


				if(promotion_time$distribution == 'logLogistic'){
					lw_prop <- log_logLogistic(y, a1 = a_prop[1], a2 = a_prop[2],  c_under = 1e-9)
					log_prior_diff <- 0
					for(i in 1:n_pars_f){
					log_prior_diff <- log_prior_diff + alpha*diff(log_inv_gamma_kernel(c(a[i], a_prop[i]), 
						promotion_time$prior_parameters[i,1], promotion_time$prior_parameters[i,2]))
					}					
				}
				if(promotion_time$distribution == 'lomax'){
					lw_prop <- log_lomax(y, a1 = a_prop[1], a2 = a_prop[2],  c_under = 1e-9)
					log_prior_diff <- 0
					for(i in 1:n_pars_f){
					log_prior_diff <- log_prior_diff + alpha*diff(log_inv_gamma_kernel(c(a[i], a_prop[i]), 
						promotion_time$prior_parameters[i,1], promotion_time$prior_parameters[i,2]))				
					}					
				}
				if(promotion_time$distribution == 'dagum'){
					lw_prop <- log_dagum(y, a1 = a_prop[1], a2 = a_prop[2],  a3 = a_prop[3], c_under = 1e-9)
					log_prior_diff <- 0				
					for(i in 1:n_pars_f){
					log_prior_diff <- log_prior_diff +  alpha*diff(log_inv_gamma_kernel(c(a[i], a_prop[i]), 
						promotion_time$prior_parameters[i,1], promotion_time$prior_parameters[i,2]))	
					}			
				}
				if(promotion_time$distribution == 'gamma_mixture'){
					w_prop <- a_prop[-(1:(2*K))]
					w_prop2 <- c(w_prop, 1)
					p_prop <- w_prop2/sum(w_prop2)
					a1_prop = a_prop[a1_ind]
					a2_prop = a_prop[a2_ind]		
					lw_prop <- log_gamma_mixture(y, a1 = a1_prop, a2 = a2_prop, p = p_prop, c_under = 1e-9)
					log_prior_diff <- 0				
					for(i in 1:(2*K)){
						k <- which(a1_ind == i)
						if(length(k) == 0){
							k <- which(a2_ind == i)
						}
						j = (i+1)%%2 + 1
	#					k <- floor((1:(2*K))/K + 0.5)[i]
	#					j = i%%2 + 1
						log_prior_diff <- log_prior_diff + alpha*diff(log_inv_gamma_kernel(c(a[i], a_prop[i]), 
							promotion_time$prior_parameters[j,1,k], promotion_time$prior_parameters[j,2,k]))
					}
					log_prior_diff <-  log_prior_diff + 
						alpha*log_dirichlet_pdf(a_vec, p_prop) - alpha*log_dirichlet_pdf(a_vec, p)

				}

				cll_prop <- complete_log_likelihood_general(y = y, X = X, 
						Censoring_status = Censoring_status, 
						g = theta_prop[1], lambda = theta_prop[2], 
						log_f = lw_prop$log_f, log_F = lw_prop$log_F, 
						b = b_prop, I_sim = I_sim, alpha = alpha)				
				#complete_log_likelihood(y = y, X = X, Censoring_status = Censoring_status, 
				#	g = theta_prop[1], lambda = theta_prop[2], 
				#	a1 = theta_prop[3], a2 = theta_prop[4], b = theta_prop[-(1:4)], I_sim = I_sim, alpha = alpha)
				cll_diff <- cll_prop$cll - cll$cll
				log_proposal_ratio <- log_prop_density_reverse - log_prop_density
				log_prior_diff <- log_prior_diff + 
					diff(dnorm(c(g, theta_prop[1]), mean = mu_g, sd = sqrt(s2_g), log = TRUE))
				log_prior_diff <- log_prior_diff + diff(log_inv_gamma_kernel(c(lambda, theta_prop[2]), a_l, b_l))
				log_prior_diff <- log_prior_diff -0.5 * mahalanobis(b_prop, mu_b, Sigma)^2 + 0.5 * mahalanobis(b, mu_b, Sigma)^2
				mh_ar <- cll_diff + alpha*log_prior_diff + log_proposal_ratio

				if(is.na(mh_ar)==FALSE){
					if(log(runif(1)) < mh_ar){
						g = theta_prop[1]
						lambda = theta_prop[2]
						a = a_prop
						b = b_prop
						if(promotion_time$distribution == 'gamma_mixture'){
							w = w_prop
							p = p_prop
						}
						cll <- cll_prop
						lw <- lw_prop
						mala_accept_rate <- mala_accept_rate + 1
					}
				}

			}
			g_mcmc[iter] <- g
			lambda_mcmc[iter] <- lambda
			a_mcmc[iter,] <- a
			b_mcmc[iter,] <- b
			mala_iter <- mala_iter + 1
		}

   	
	#--------------------------------------------
	#	Gibbs step for latent variables	
	#--------------------------------------------
		logS <- cll$logS 
		logP0 <- cll$logP0
#		print(cll)
		susceptible_prob <- 1/(1 + (exp(logS - logP0) - 1)^{-alpha}) # 
		na_ind1 <- which(is.na(susceptible_prob))
		na_ind2 <- which(logS - logP0 <= 0)
		na_ind <- union(na_ind1, na_ind2)
		
		if(length(na_ind) > 0){
	#		print(length(na_ind))
	#		print(logS[na_ind]-logP0[na_ind])
			susceptible_prob[na_ind] <- 1/(1 + (1e-6)^{-alpha})
		}
		# latent binary variables (I = 1 => susceptible, I = 0 => cured)
		I_sim <- rep(1, n)
		# for i in D1: I_i = 1
		I_sim[D0] <- rbinom(n = length(D0), size = 1, prob = susceptible_prob[D0])
		I_sim_values_D0[iter, ] <- I_sim[D0]
		cll <- complete_log_likelihood_general(y = y, X = X, Censoring_status = Censoring_status, 
				g = g, lambda = lambda, 
				log_f = lw$log_f,  log_F = lw$log_F,
				b = b, I_sim = I_sim, alpha = alpha)
#		print("6")
#		print(cll$cll)

		log_prior_density <- alpha*log_gamma_mix_kernel(g,mu_g, s2_g)  + 
				alpha*log_inv_gamma_kernel(lambda, a_l, b_l) -
				alpha*0.5 * mahalanobis(b, mu_b, Sigma)^2

		if(promotion_time$distribution == 'gamma_mixture'){
			a1 <- a[a1_ind]
			a2 <- a[a2_ind]			
			for(k in 1:K){
				log_prior_density <- log_prior_density + sum(alpha*log_inv_gamma_kernel(c(a1[k], a2[k]), 
							promotion_time$prior_parameters[,1,k], promotion_time$prior_parameters[,2,k]))
			}
			a_vec <- rep(promotion_time$dirichlet_concentration_parameter, K)
			log_prior_density <- log_prior_density + log_dirichlet_pdf(a_vec, p)
			
		}else{
			for(i in 1:n_pars_f){
				log_prior_density <- log_prior_density + alpha*log_inv_gamma_kernel(a[i], 
						promotion_time$prior_parameters[i,1], promotion_time$prior_parameters[i,2])
			}
		}
				
		cllValues[iter] <- cll$cll
		lpd[iter] <- log_prior_density
		if(iter %% 100 == 0 & verbose){
			cat(paste0("* iteration = ", iter, ".  Current ergodic means (after discarding 30% of iterations):"), "\n")			
			burn <- floor(0.3*iter)
			ergMean <- c(mean(g_mcmc[(burn+1):iter]), mean(lambda_mcmc[(burn+1):iter]), 
				colMeans(as.matrix(a_mcmc[(burn+1):iter,])), 
				colMeans(as.matrix(b_mcmc[(burn+1):iter,])))
				names(ergMean) <- c("gamma", "lambda", paste0('a',1:n_pars_f), paste0('b',bRange))
			print(round(ergMean,2))
			cat(paste0("    g_accept_rate = ", round(100*g_accept_rate/mh_iter, 2), "%.   "), "\n")
                        cat(paste0("    lambda_accept_rate = ", round(100*lambda_accept_rate/mh_iter, 2), "%.   "), "\n")
                        for(i in 1:n_pars_f){
	                        cat(paste0("    a",i,"_accept_rate = ", round(100*a_accept_rate[i]/mh_single, 2), "%.   "), "\n")                                          
                        
                        }

                        cat(paste0("    a_all_accept_rate = ", round(100*a_all_accept_rate/mh_all, 2), "%.   "), "\n")                                          
                        
                        cat(paste0("    b_accept_rate = ", round(100*b_accept_rate/mh_iter, 2), "%.   "), "\n")                                                                                               
			cat(paste0("    mala_accept_rate = ", round(100*mala_accept_rate/mala_iter, 2), "%.   "), "\n")
			cat(paste0("\n"))
			if(plot){
				par(mfrow = c(2,3))		
				plot(g_mcmc[1:iter], type = "l", xlab = "MCMC iteration", ylab = bquote(gamma))
				plot(lambda_mcmc[1:iter], type = "l", xlab = "MCMC iteration", ylab = bquote(lambda))
				if(promotion_time$distribution == 'gamma_mixture'){
				matplot(as.matrix(a_mcmc[1:iter,1:(2*K)]), type = "l", xlab = "MCMC iteration", ylab = 'promotion time', 
					main = paste0(promotion_time$distribution, ' distr. parameters'))				
				lText <- paste0('a',1:(n_pars_f - K + 1))
				legend("topright", lText, col = 1:(n_pars_f-K+1), lty = 1:n_pars_f)
				props <- t(apply(cbind(a_mcmc[1:iter,(2*K+1):n_pars_f],1), 1, function(x)x/sum(x)))
				matplot(as.matrix(props), type = "l", xlab = "MCMC iteration", ylab = 'promotion time', 
					main = paste0(promotion_time$distribution, ' distr. parameters'))				
				lText <- paste0('p',1:K)
				legend("topright", lText, col = 1:K, lty = 1:K)

				}else{
				matplot(as.matrix(a_mcmc[1:iter,]), type = "l", xlab = "MCMC iteration", ylab = 'promotion time', 
					main = paste0(promotion_time$distribution, ' distr. parameters'))
					lText <- paste0('a',1:n_pars_f)
					legend("topright", lText, col = 1:n_pars_f, lty = 1:n_pars_f)

				}
				matplot(as.matrix(b_mcmc[1:iter,]), type = "l", xlab = "MCMC iteration", ylab = 'regression coefficients')														
				lText <- paste0('b',bRange)
				legend("topright", lText, col = 1:nCov, lty = 1:nCov)
				plot(cllValues[1:iter], type = "l", xlab = "MCMC iteration", ylab = "complete log-likelihood")
			}
		}

	}
		result <- vector("list", length = 5)
		colnames(b_mcmc) <- paste0('b',bRange, '_mcmc')
		colnames(a_mcmc) <- paste0('a',1:n_pars_f, '_mcmc')		
		mcmc_sample <- as.mcmc(cbind(g_mcmc, lambda_mcmc, a_mcmc, b_mcmc))
		latent_status_D0 <- I_sim_values_D0
		colnames(latent_status_D0) <- D0
		complete_log_likelihood <- cllValues
		result[[1]] <- mcmc_sample
		result[[2]] <- complete_log_likelihood
		result[[3]] <- c(g_accept_rate/mh_iter, lambda_accept_rate/mh_iter, 
				a_accept_rate/mh_iter, b_accept_rate/mh_iter, mala_accept_rate/mala_iter)
		result[[4]] <- latent_status_D0
		result[[5]] <- lpd
		names(result[[3]]) <- c("gamma", "lambda", paste0('a',1:n_pars_f, '_mcmc'), "betas","mala")
		names(result) <- c("mcmc_sample", "complete_log_likelihood", 
				"acceptance_rates", "latent_status_censored", "log_prior_density")
	}else{
		result <- vector("list", length = 5)
		complete_log_likelihood <- cllValues

		result[[2]] <- complete_log_likelihood
		result[[5]] <- lpd
		names(result) <- c("mcmc_sample", "complete_log_likelihood", 
				"acceptance_rates", "latent_status_censored", "log_prior_density")

	}
	return(result)	
}





cure_rate_MC3 <- function( y, X, Censoring_status, nChains = 12, 
				mcmc_cycles = 15000, 
				alpha = NULL, 
				nCores = 8, 
				sweep = 5,
				mu_g = 1, s2_g = 1, a_l = 2.1, b_l = 1.1, 
				mu_b = rep(0,dim(X)[2]), Sigma = 100*diag(dim(X)[2]),
				g_prop_sd = 0.045, 
				lambda_prop_scale = 0.03, 
				b_prop_sd = rep(0.022, dim(X)[2]), 
				initialValues = NULL, 
				plot = TRUE,
				adjust_scales = FALSE,
				verbose = FALSE, tau_mala = 0.000015, mala = 0.15, 
				promotion_time = list(distribution = 'weibull', 
						prior_parameters = matrix(rep(c(2.1, 1.1), 2), byrow = TRUE, 2, 2),
						prop_scale = c(0.1, 0.2)
					),
					single_MH_in_f = 0.2
					){
	X <- as.matrix(X)

#	specify default priors per distribution family
	if(is.null(promotion_time$prior_parameters)){
		if(promotion_time$distribution == 'exponential'){
			promotion_time$prior_parameters = matrix(rep(c(2.1, 1.1), 1), byrow = TRUE, 1, 2)
			if(is.null(promotion_time$prop_scale)){
				promotion_time$prop_scale = rep(0.1, dim(promotion_time$prior_parameters)[1])
			}
		}

		if(any(promotion_time$distribution %in% c('weibull', 'gamma', 'lomax', 'gompertz', 'logLogistic'))){
			promotion_time$prior_parameters = matrix(rep(c(2.1, 1.1), 2), byrow = TRUE, 2, 2)
			if(is.null(promotion_time$prop_scale)){
				promotion_time$prop_scale = rep(0.1, dim(promotion_time$prior_parameters)[1])
			}
		}

		if(promotion_time$distribution == 'dagum'){
			promotion_time$prior_parameters = matrix(rep(c(2.1, 1.1), 3), byrow = TRUE, 3, 2)
			if(is.null(promotion_time$prop_scale)){
				promotion_time$prop_scale = rep(0.1, dim(promotion_time$prior_parameters)[1])
			}
		}

		if(promotion_time$distribution == 'gamma_mixture'){
			if(is.null(promotion_time$K)){
				K = 2
				promotion_time$K = K
				cat('WARNING: the number of mixture components (K) is not specified. I will set it to K = 2','\n')
			}else{
				K = promotion_time$K
			}
			promotion_time$prior_parameters = array(data = NA, dim = c(2,2,K))
			for(k in 1:K){
				promotion_time$prior_parameters[,,k] = matrix(rep(c(2.1, 1.1), 2), byrow = TRUE, 2, 2)			
			}
			if(is.null(promotion_time$prop_scale)){
				promotion_time$prop_scale = rep(0.1, K*2 + K - 1)
			}
			if(is.null(promotion_time$dirichlet_concentration_parameter)){
				promotion_time$dirichlet_concentration_parameter = 1
			}
		}


	}
#	

	if(plot & verbose){
		oldpar <- par(no.readonly = TRUE)
		on.exit(par(oldpar)) 
		dev.off()
		dev.new(width=9, height=9, unit="in")						
	}
	if(min(y) < 0){stop("Negative values are not allowed in the response variable.")}				
	if(nChains < 2){
		print("Only one chains is considered. Redirecting to the `cure_rate_metropolis_hastings_cpp`...")
		
		run <- cure_rate_mcmc(y, X, Censoring_status,  m = mcmc_cycles, alpha = 1,
						mu_g = mu_g, s2_g = s2_g, a_l = a_l, b_l = b_l, 
						promotion_time = promotion_time,
						mu_b = mu_b, Sigma = Sigma,
						g_prop_sd = g_prop_sd, 
						lambda_prop_scale = lambda_prop_scale, 
						b_prop_sd = b_prop_sd, 
								initialValues = initialValues, 
								plot = plot,
								verbose = verbose,
								tau_mala = tau_mala, mala = mala, single_MH_in_f = single_MH_in_f
							)
		return(run)
		
	}
	if(is.null(alpha)){
		alpha <- 1/1.001^{(1:nChains)^3 - 1}
		if(verbose){
		cat("Default (inverse) temperatures: ", "\n")
		print(alpha)
		}
	}else{
		if(length(alpha)!= nChains){
			stop("number of temperatures does not match number of chains.")
		}
		if(alpha[1] != 1){
			stop("alpha[1] should be equal to 1.")
		}
	}
	nCov <- dim(X)[2]
	n_pars_f <- dim(promotion_time$prior_parameters)[1]
	if(promotion_time$distribution == 'gamma_mixture'){
		K = promotion_time$K
		if(dim(promotion_time$prior_parameters)[3] != K){stop("incosistent number of mixture components")}		
		n_pars_f = K * n_pars_f + K - 1# note that we include all mixing weights
	}
	
	nPars <- nCov + 2 + n_pars_f
	D0 <- which(Censoring_status == 0)
	nCensored <- length(D0)
	parallel_mcmc_samples <- array(data = NA, dim = c(1, nPars + nCensored, nChains))
	target_mcmc <- array(data = NA, dim = c(mcmc_cycles, nPars + nCensored))

	if(length(unique(X[,1])) == 1){
		b_names <- paste0('b',1:nCov - 1,'_mcmc')
	}else{
		b_names <- paste0('b',1:nCov,'_mcmc')	
	}
	a_names <- paste0('a',1:n_pars_f, '_mcmc')

	dimnames(parallel_mcmc_samples) = list(paste0("current_value"),
		c("g_mcmc", "lambda_mcmc", a_names, b_names,paste0("status_",D0)),
						paste0("chain_",1:nChains))
		dimnames(parallel_mcmc_samples)[[3]][1]	<- "chain_1_(target_chain)"
		dimnames(target_mcmc) = list(paste0("cycle_",1:mcmc_cycles),
			c("g_mcmc", "lambda_mcmc", a_names, b_names,paste0("status_",D0)))

	cllValues <- numeric(mcmc_cycles)


	g_prop_sd <- rep(g_prop_sd, nChains)
	lambda_prop_scale <- rep(lambda_prop_scale, nChains)
	b_prop_sd <- array(data = b_prop_sd, dim = c(nChains, nCov))
	a_prop_scale <- matrix(	promotion_time$prop_scale, n_pars_f, nChains, byrow = TRUE )
	promotion_time_all <- vector('list', length = nChains)
	for(i in 1:nChains){
		promotion_time_all[[i]] <- promotion_time
	}

	if(adjust_scales == TRUE){
		cat("Making a warm-up run in order to adjust the proposal scale per chain...","\n")
	      chain_iter <- 1

#	      g_prop_sd[chain_iter] = g_prop_sd 
#		lambda_prop_scale[chain_iter] = lambda_prop_scale
#		a1_prop_scale[chain_iter] = a1_prop_scale
#		a2_prop_scale[chain_iter] = a2_prop_scale
#		b_prop_sd[chain_iter,] = b_prop_sd

	      criterion <- TRUE
		adjFactor <- 0.75
		iter <- 0
	 
	      while(criterion & (iter < 20)){

		  run <- cure_rate_mcmc(  y, X, Censoring_status,  
		                        m = 1000, 
				mu_g = mu_g, s2_g = s2_g, a_l = a_l, b_l = b_l, 
				promotion_time = promotion_time_all[[chain_iter]],
				mu_b = mu_b, Sigma = Sigma,
		                        alpha = alpha[chain_iter],
		                        g_prop_sd = g_prop_sd[chain_iter], 
		                        lambda_prop_scale = lambda_prop_scale[chain_iter], 
		                        b_prop_sd = b_prop_sd[chain_iter,], 
		                        plot = FALSE, verbose = FALSE, tau_mala = tau_mala, 
		                        mala = mala, single_MH_in_f = single_MH_in_f)

			if(mala>0){
				indUp <- which((run$acceptance_rates > 0.6) == TRUE)
				indDown <- which((run$acceptance_rates < 0.4) == TRUE)
			}else{
				indUp <- which((run$acceptance_rates[-(2 + n_pars_f + 2)] > 0.3) == TRUE)
				indDown <- which((run$acceptance_rates[-(2 + n_pars_f + 2)] < 0.15) == TRUE)			
			}
		        if(length(indUp) > 0){
		                criterion = TRUE; 
		                if(1 %in% indUp){g_prop_sd[chain_iter] <- g_prop_sd[chain_iter]/adjFactor}
		                if(2 %in% indUp){lambda_prop_scale[chain_iter] <- lambda_prop_scale[chain_iter]/adjFactor}      
				for(i in 1:n_pars_f){
			                if(2 + i %in% indUp){a_prop_scale[i, chain_iter] <- a_prop_scale[i,chain_iter]/adjFactor}              				
				}
		                if(2 + n_pars_f + 1 %in% indUp){b_prop_sd[chain_iter,] <- b_prop_sd[chain_iter,]/adjFactor}
		                if(mala>0){
      		                	if(6 %in% indUp){tau_mala <- tau_mala/adjFactor}
      		                }
		                criterion1 = TRUE                               
		                
		        }else{criterion1 = FALSE}
		        if(length(indDown) > 0){
		                criterion = TRUE
		                if(1 %in% indDown){g_prop_sd[chain_iter] <- g_prop_sd[chain_iter] * adjFactor}
		                if(2 %in% indDown){lambda_prop_scale[chain_iter] <- lambda_prop_scale [chain_iter]* adjFactor}    
				for(i in 1:n_pars_f){
					if(2 + i %in% indDown){a_prop_scale[i,chain_iter] <- a_prop_scale[i,chain_iter] * adjFactor}
				}
		                if(2 + n_pars_f + 1 %in% indDown){b_prop_sd[chain_iter,] <- b_prop_sd[chain_iter,] * adjFactor}
		                if(mala>0){
		      			if(6 %in% indDown){tau_mala <- tau_mala*adjFactor}
		      		}
		                criterion2 = TRUE
		        }else{criterion2 = FALSE}
		        if((criterion1 == FALSE) & (criterion2 == FALSE)){criterion = FALSE}
		        iter <- iter + 1
		        #print(iter)

		}
		cat(paste0("    proposal scale for chain: ",chain_iter),"\n")
		print(c(g_prop_sd[chain_iter], lambda_prop_scale[chain_iter], a_prop_scale[,chain_iter], b_prop_sd[chain_iter,], tau_mala), digits = 4)

	for(chain_iter in 2:nChains){
	# initial metropolis-hastings proposal parameters (they will be adjusted)
	      g_prop_sd[chain_iter] = g_prop_sd[chain_iter - 1]
		lambda_prop_scale[chain_iter] = lambda_prop_scale[chain_iter-1]
		a_prop_scale[,chain_iter] = a_prop_scale[,chain_iter-1]
		b_prop_sd[chain_iter,] = b_prop_sd[chain_iter-1,]
		promotion_time_all[[chain_iter]]$prop_scale <- a_prop_scale[,chain_iter]
	      criterion <- TRUE
		adjFactor <- 0.75
		iter <- 0
	 
	      while(criterion & (iter < 5)){

		  run <- cure_rate_mcmc(  y, X, Censoring_status,  
		                        m = 500, 
				mu_g = mu_g, s2_g = s2_g, a_l = a_l, b_l = b_l, 
				promotion_time = promotion_time_all[[chain_iter]],
				mu_b = mu_b, Sigma = Sigma,
				alpha = alpha[chain_iter],
		                        g_prop_sd = g_prop_sd[chain_iter], 
		                        lambda_prop_scale = lambda_prop_scale[chain_iter], 
		                        b_prop_sd = b_prop_sd[chain_iter,], 
		                        plot = FALSE, verbose = FALSE, tau_mala = tau_mala, 
		                        mala = mala, single_MH_in_f = single_MH_in_f)

			if(mala>0){
				indUp <- which((run$acceptance_rates > 0.6) == TRUE)
				indDown <- which((run$acceptance_rates < 0.4) == TRUE)
			}else{
				indUp <- which((run$acceptance_rates[-(2 + n_pars_f + 2)] > 0.3) == TRUE)
				indDown <- which((run$acceptance_rates[-(2 + n_pars_f + 2)] < 0.15) == TRUE)			
			}

		        if(length(indUp) > 0){
		                criterion = TRUE; 
		                if(1 %in% indUp){g_prop_sd[chain_iter] <- g_prop_sd[chain_iter]/adjFactor}
		                if(2 %in% indUp){lambda_prop_scale[chain_iter] <- lambda_prop_scale[chain_iter]/adjFactor}      
				for(i in 1:n_pars_f){
					if(2 + i %in% indUp){a_prop_scale[i, chain_iter] <- a_prop_scale[i,chain_iter]/adjFactor}              		
				}		                
		                if(2 + n_pars_f + 1 %in% indUp){b_prop_sd[chain_iter,] <- b_prop_sd[chain_iter,]/adjFactor}
		                if(mala>0){
      		                	if(6 %in% indUp){tau_mala <- tau_mala/adjFactor}
      		                }
		                
		                criterion1 = TRUE                               
		                
		        }else{criterion1 = FALSE}
		        if(length(indDown) > 0){
		                criterion = TRUE
		                if(1 %in% indDown){g_prop_sd[chain_iter] <- g_prop_sd[chain_iter] * adjFactor}
		                if(2 %in% indDown){lambda_prop_scale[chain_iter] <- lambda_prop_scale [chain_iter]* adjFactor}    
		                for(i in 1:n_pars_f){
			                if(2 + i %in% indDown){a_prop_scale[i,chain_iter] <- a_prop_scale[i,chain_iter] * adjFactor}            
		                }
		                if(2 + n_pars_f + 1 %in% indDown){b_prop_sd[chain_iter,] <- b_prop_sd[chain_iter,] * adjFactor}                            
		                if(mala>0){
		      			if(6 %in% indDown){tau_mala <- tau_mala*adjFactor}
		      		}
		                
		                criterion2 = TRUE
		        }else{criterion2 = FALSE}
		        if((criterion1 == FALSE) & (criterion2 == FALSE)){criterion = FALSE}
		        iter <- iter + 1
		        #print(iter)

		}
		cat(paste0("    proposal scale for chain: ", chain_iter),"\n")
		print(c(g_prop_sd[chain_iter], lambda_prop_scale[chain_iter], a_prop_scale[,chain_iter], b_prop_sd[chain_iter,], tau_mala), digits = 4)
	}
	cat("   done Metropolis-Hastings proposal scale adjustements!","\n")
	}
	
	if(verbose){
		cat(paste0("Running MC^3 with ", nChains, " heated chains...", "\n"))
	}
	swap_accept_per_chain <- numeric(nChains - 1)
	n_attempts_per_chain <- numeric(nChains - 1)
	all_cll_values <- array(NA, dim = c(nChains, mcmc_cycles))
		swap_accept <- 0

		# initialization
		cycle <- 1
		registerDoParallel(nCores)
		parLoop <- foreach(chain_iter = 1:nChains #, .export = ls(envir = globalenv())
		) %dopar% {
#			if(chain_iter == 1){perform_mala = mala}else{
#				if(alpha[chain_iter] == 1){
#					perform_mala = mala
#				}else{perform_mala = 0}
#			}
			mcmc_chain <- cure_rate_mcmc(  y, X, Censoring_status,  
				                m = sweep, 
						mu_g = mu_g, s2_g = s2_g, a_l = a_l, b_l = b_l, 
						promotion_time = promotion_time_all[[chain_iter]],
						mu_b = mu_b, Sigma = Sigma,
				                alpha = alpha[chain_iter],
				                g_prop_sd = g_prop_sd[chain_iter], 
				                lambda_prop_scale = lambda_prop_scale[chain_iter], 
				                b_prop_sd = b_prop_sd[chain_iter,], plot = FALSE, 
				                tau_mala = tau_mala, mala = mala, single_MH_in_f = single_MH_in_f)
						
		}
		stopImplicitCluster()
		for(chain_iter in 1:nChains){
			parallel_mcmc_samples[1, ,chain_iter] <- c(parLoop[[chain_iter]]$mcmc_sample[sweep,], 
										parLoop[[chain_iter]]$latent_status_censored[sweep,])
			target_mcmc[cycle,] <- parallel_mcmc_samples[1,,1]
			all_cll_values[chain_iter,cycle] <- parLoop[[chain_iter]]$complete_log_likelihood[sweep]
		}
		cllValues[1] <- parLoop[[1]]$complete_log_likelihood[sweep]

		for(cycle in 2:mcmc_cycles){

			#attempt chain-swap
#			myPair <- sample(nChains,2,replace = FALSE)
#			j1 <- myPair[1]
#			j2 <- myPair[2]

			j1 <- sample(nChains - 1 ,1,replace = FALSE)
			j2 <- j1 + 1
			n_attempts_per_chain[j1] <- n_attempts_per_chain[j1] + 1

			log_posterior1 <- parLoop[[j1]]$complete_log_likelihood[sweep] + parLoop[[j1]]$log_prior_density[sweep]
			log_posterior2 <- parLoop[[j2]]$complete_log_likelihood[sweep] + parLoop[[j2]]$log_prior_density[sweep]
			denom <- log_posterior1 + log_posterior2
#			print(denom)
			state_j1 <- as.list(parLoop[[j1]]$mcmc_sample[sweep,])
			state_j1[["I_sim_D0"]] <- parLoop[[j1]]$latent_status_censored[sweep,]

			state_j2 <- as.list(parLoop[[j2]]$mcmc_sample[sweep,])
			state_j2[["I_sim_D0"]] <- parLoop[[j2]]$latent_status_censored[sweep,]
			nom1 <- cure_rate_mcmc(  y, X, Censoring_status,  
								m = 1, # only the log-posterior is computed at the initial values
								mu_g = mu_g, s2_g = s2_g, a_l = a_l, b_l = b_l, 
						promotion_time = promotion_time_all[[j1]],
						mu_b = mu_b, Sigma = Sigma,
								alpha = alpha[j1], plot = FALSE, 
								initialValues = state_j2, tau_mala = tau_mala, mala = mala)
			nom2 <- cure_rate_mcmc(  y, X, Censoring_status,  
								m = 1, # only the log-posterior is computed at the initial values
								alpha = alpha[j2], plot = FALSE, 
								mu_g = mu_g, s2_g = s2_g, a_l = a_l, b_l = b_l, 
						promotion_time = promotion_time_all[[j2]],
						mu_b = mu_b, Sigma = Sigma,
								initialValues = state_j1, tau_mala = tau_mala, mala = mala)
			nom <- nom1$complete_log_likelihood + nom1$log_prior_density +
				 nom2$complete_log_likelihood + nom2$log_prior_density	
	#		print(nom)
	#		print(nom-denom)
			mmm <- log(runif(1)) < nom - denom
			if (is.na(mmm) == FALSE){
		 	if(mmm){
		 		# chain switching accepted
		 		tmp <- parallel_mcmc_samples[1,,j2]
		 		parallel_mcmc_samples[1,,j2] <- parallel_mcmc_samples[1,,j1]
		 		parallel_mcmc_samples[1,,j1] <- tmp
				swap_accept <- 	swap_accept + 1 
				swap_accept_per_chain[j1] <- swap_accept_per_chain[j1] + 1
		 	}
		 	}

		 	# run mcmc for each chain now for sweep iterations, initialized by the previous values.
			registerDoParallel(nCores)
			parLoop <- foreach(chain_iter = 1:nChains#, .export = ls(envir = globalenv())
			) %dopar% {
				initialValues <-  as.list(parallel_mcmc_samples[1, 1:nPars,chain_iter])
				initialValues[["I_sim_D0"]] <- parallel_mcmc_samples[1, -(1:nPars),chain_iter]
#				if(chain_iter == 1){perform_mala = mala}else{
#					if(alpha[chain_iter] == 1){
#						perform_mala = mala
#					}else{perform_mala = 0}
#				}
				
				mcmc_chain <- cure_rate_mcmc(  y, X, Censoring_status,  
						        m = sweep, 
						        alpha = alpha[chain_iter],
							mu_g = mu_g, s2_g = s2_g, a_l = a_l, b_l = b_l, 
							promotion_time = promotion_time_all[[chain_iter]],
							mu_b = mu_b, Sigma = Sigma,
						        g_prop_sd = g_prop_sd[chain_iter], 
						        lambda_prop_scale = lambda_prop_scale[chain_iter], 
						        b_prop_sd = b_prop_sd[chain_iter,], plot = FALSE,
						        initialValues = initialValues, 
						        tau_mala = tau_mala, mala = mala, single_MH_in_f = single_MH_in_f)
							
			}
			stopImplicitCluster()
 			for(chain_iter in 1:nChains){
				parallel_mcmc_samples[1, ,chain_iter] <- c(parLoop[[chain_iter]]$mcmc_sample[sweep,], 
											parLoop[[chain_iter]]$latent_status_censored[sweep,])
				target_mcmc[cycle,] <- 	parallel_mcmc_samples[1, ,1]
				all_cll_values[chain_iter,cycle] <- parLoop[[chain_iter]]$complete_log_likelihood[sweep]														
			}
			cllValues[cycle] <- parLoop[[1]]$complete_log_likelihood[sweep]
 			if(cycle %% 500 == 0 && verbose == TRUE){
 				cat(paste0("** MCMC cycle: ", cycle, " completed."),"\n")
				cat(paste0("   Chain-swap rate (overall): ", round(100*swap_accept/cycle,2),"%."),"\n")
				cat(paste0("                 (per chain): ", paste0(round(100*swap_accept_per_chain/(n_attempts_per_chain),2),"%",collapse = ", ")),"\n")
				cat("   Current ergodic means and median (after discarding 30% of iterations):","\n")
		                burn <- floor(0.3*cycle)
		                ergMean <- colMeans(target_mcmc[(burn+1):cycle,1:nPars])
		                        names(ergMean) <- c("gamma", "lambda", a_names, b_names)
		                ergMed <- apply(target_mcmc[(burn+1):cycle,1:nPars], 2, median)
		                        names(ergMed) <- c("gamma", "lambda", a_names, b_names)
		                tmpMat <- rbind(ergMean, ergMed)
		                rownames(tmpMat) <- c("mean", "median")
		                print(round(tmpMat,2))

					if(plot){
					par(mfrow = c(2,2))		
					plot(target_mcmc[1:cycle, 1], 
						type = "l", xlab = "MCMC iteration", ylab = bquote(gamma))
					plot(target_mcmc[1:cycle, 2], 
						type = "l", xlab = "MCMC iteration", ylab = bquote(lambda))

					if(promotion_time$distribution == 'gamma_mixture'){
					K = promotion_time$K
					#matplot(as.matrix(target_mcmc[1:cycle,3:(3+2*K-1)]), type = "l", 
					#xlab = "MCMC iteration", ylab = 'promotion time', 
					#	main = paste0(promotion_time$distribution, ' distr. parameters'))				
					#lText <- paste0('a',1:(n_pars_f + 1 - K))
					#legend("topright", lText, col = 1:(n_pars_f + 1-K), lty = 1:n_pars_f)
					props <- t(apply(cbind(target_mcmc[1:cycle,(2 + 2*K+1):(2+n_pars_f)], 1), 
						1, function(x)x/sum(x)))
					matplot(as.matrix(props), type = "l", xlab = "MCMC iteration", ylab = 'promotion time', 
						main = paste0(promotion_time$distribution, ' mixing proportions'))				
					lText <- paste0('p',1:K)
					legend("topright", lText, col = 1:K, lty = 1:K)
					}else{
					matplot(as.matrix(target_mcmc[1:cycle,3:(2+n_pars_f)]), type = "l", 
						xlab = "MCMC iteration", ylab = '', 
						main = paste0(promotion_time$distribution, ' distr. parameters'))
						lText <- paste0('a',1:n_pars_f)
						legend("topright", lText, col = 1:n_pars_f, lty = 1:n_pars_f)
					}


						
					matplot(as.matrix(target_mcmc[1:cycle, (2+n_pars_f + 1):nPars]), col = 1:nCov, lty = 1:nCov,
						type = "l", xlab = "MCMC iteration", ylab = 'regression coefficients')														
					lText <- b_names
					legend("topright", lText, col = 1:nCov, lty = 1:nCov)
					#plot(cllValues[1:cycle], type = "l", xlab = "MCMC iteration", ylab = "complete log-likelihood",
					#ylim = quantile(cllValues[1:cycle],c(0.01,0.999)))
					#matplot(t(all_cll_values[,501:cycle]), type='l', col  = heat.colors(nChains) )
					#	points(all_cll_values[1,501:cycle], type = 'b', 
					#		col = heat.colors(nChains)[1])					
				}
				
			}


		}
		result <- vector("list", length = 9)
		if(promotion_time$distribution == 'gamma_mixture'){
		# get mixing proportions
			K = promotion_time$K
			props <- t(apply(cbind(target_mcmc[1:cycle,(2 + 2*K+1):(2+n_pars_f)], 1), 
				1, function(x)x/sum(x)))
			colnames(props) <- paste0('w',1:K)
			target_mcmc[,(2 + 2*K+1):(2+n_pars_f)] <- props[,-K]
			colnames(target_mcmc)[(2 + 2*K+1):(2+n_pars_f)] <- paste0('w',1:(K-1))
		}		
		result[[1]] <- as.mcmc(target_mcmc[,1:nPars])
		result[[2]] <- target_mcmc[,-(1:nPars)]
		result[[3]] <- cllValues
		result[[4]] <- swap_accept_per_chain/(n_attempts_per_chain + 0.001)
		result[[5]] <- all_cll_values
		result[[6]] <- vector('list', length = 10)
		result[[6]][[1]] = y
		result[[6]][[2]] = X
		result[[6]][[3]] = Censoring_status
		result[[6]][[4]] = mu_g
		result[[6]][[5]] = s2_g
		result[[6]][[6]] = a_l
		result[[6]][[7]] = b_l
		result[[6]][[8]] = mu_b
		result[[6]][[9]] = Sigma
		result[[6]][[10]] = promotion_time
		names(result[[6]]) <- c('y', 'X', 'Censoring_status', 'mu_g', 's2_g', 'a_l', 'b_l', 'mu_b', 'Sigma', 'promotion_time')
#####################################################################################
# logP computation
	ct = exp(exp(-1))
	log_S_p <- function(g, lambda, log_F, b, x){
		theta <- exp(x %*% b)
		return(-log(1 + g * theta * ct^{g*theta} * exp(log_F)^lambda)/g)

	}

	log_f_p <- function(g, lambda, log_f, log_F, b, logS){
		# logS = log_S_p(tau = tau, g = g, lambda = lambda, a1 = a1, a2 = a2, b0 = b0, b1 = b1, b2 = b2)
		log_theta <- X %*% b
		return(
		        (1 + g) * logS + log(lambda) + log_theta +
		        g*exp(log_theta)*log(ct) + 
		        (lambda - 1)*log_F + #### NOTE: this causes the log -> -inf, so now it regulated.
		        log_f
		)
	 }

	log_inv_gamma_kernel <- function(x, a, b){
		if(min(c(x,a,b)) < 0){
		stop("input should be positive")
		}
		return(-b/x - (a+1) * log(x))
	}

	log_prior_gamma <- function(g, mu_g, s2_g){
		# mu_g is a_g
		# s2_g is b_g
		return(mu_g * log(s2_g) - 2*lgamma(mu_g) + (mu_g - 1)*log(abs(g)) - s2_g*abs(g))
	}


	logL <- logP <- numeric(mcmc_cycles)
	if(length(unique(X[,1])) == 1){
		bIndices <- paste0('b',1:nCov - 1,'_mcmc')	
	}else{
		bIndices <- paste0('b',1:nCov,'_mcmc')
	}
	
	for(iter in 1:mcmc_cycles){
		if(promotion_time$distribution == 'exponential'){
			lw <- log_weibull(y, a1 = result[[1]][iter,'a1_mcmc'], a2 = 1,  c_under = 1e-9)
		}
		if(promotion_time$distribution == 'weibull'){
			lw <- log_weibull(y, a1 = result[[1]][iter,'a1_mcmc'], a2 = result[[1]][iter,'a2_mcmc'],  c_under = 1e-9)
		}
		if(promotion_time$distribution == 'gamma'){
			lw <- log_gamma(y, a1 = result[[1]][iter,'a1_mcmc'], a2 = result[[1]][iter,'a2_mcmc'],  c_under = 1e-9)
		}
		if(promotion_time$distribution == 'gompertz'){
			lw <- log_gompertz(y, a1 = result[[1]][iter,'a1_mcmc'], a2 = result[[1]][iter,'a2_mcmc'],  c_under = 1e-9)
		}
		if(promotion_time$distribution == 'logLogistic'){
			lw <- log_logLogistic(y, a1 = result[[1]][iter,'a1_mcmc'], a2 = result[[1]][iter,'a2_mcmc'],  c_under = 1e-9)
		}
		if(promotion_time$distribution == 'lomax'){
			lw <- log_lomax(y, a1 = result[[1]][iter,'a1_mcmc'], a2 = result[[1]][iter,'a2_mcmc'],  c_under = 1e-9)
		}
		if(promotion_time$distribution == 'dagum'){
			lw <- log_dagum(y, a1 = result[[1]][iter,'a1_mcmc'], a2 = result[[1]][iter,'a2_mcmc'], 
				a3 = result[[1]][iter,'a3_mcmc'], c_under = 1e-9)
		}
		if(promotion_time$distribution == 'gamma_mixture'){
			K = promotion_time$K
			a = result[[1]][iter,3:(3+n_pars_f - 1)]
			w <- c(a[-(1:(2*K))], 1)
			p <- w/sum(w)
			a1_ind <- seq(1, 2*K, by = 2)
			a2_ind <- seq(2, 2*K, by = 2)		
			a1 = a[a1_ind]
			a2 = a[a2_ind]		
			lw <- log_gamma_mixture(y, a1 = a1, a2 = a2, p = p, c_under = 1e-9)
		}

		logS <- log_S_p(g = result[[1]][iter,'g_mcmc'], 
			lambda = result[[1]][iter,'lambda_mcmc'], 
			log_F = lw$log_F,
			b = result[[1]][iter, bIndices], 
			x = X)

		logf <- log_f_p(
			g = result[[1]][iter,'g_mcmc'], 
			lambda = result[[1]][iter,'lambda_mcmc'], 
			log_f = lw$log_f, log_F = lw$log_F,
			b = result[[1]][iter, bIndices], 
			logS = logS
			)
		logL[iter] <- sum(Censoring_status * logf) + sum((1-Censoring_status)*logS)
		g = result[[1]][iter,'g_mcmc']
		lambda = result[[1]][iter,'lambda_mcmc']
		a <- result[[1]][iter, 3:(3+n_pars_f - 1)]
		
		b = result[[1]][iter, bIndices]
		log_prior_density <- log_prior_gamma(g, mu_g, s2_g) + 
			        log_inv_gamma_kernel(lambda, a_l, b_l) #+
#			        log_inv_gamma_kernel(a1, a_1, b_1) +
#			        log_inv_gamma_kernel(a2, a_2, b_2) -
#			        0.5 * mahalanobis(b, mu_b, Sigma)^2
#		print(logL[iter])
#		print(log_prior_density)
#		print('')
		logP[iter] <- logL[iter] + log_prior_density
		if(iter %% 500 == 0 && verbose == TRUE){
		cat(paste0('          Computing log-posterior at iteration: ', iter),'\r')
		}
	}	
	cat('\n')
	result[[7]] <- logP
	n <- dim(X)[1]
	n_parameters <- dim(X)[2] + 2 + n_pars_f
	BIC <- -2 * max(logL, na.rm = TRUE) + n_parameters * log(n)
	map_index <- which.max(logP)
	map_estimate <- result[[1]][map_index, ]
#        names(map_estimate) <- unlist(lapply(strsplit(names(map_estimate), split = '_'), function(x)x[1]))
#        names(map_estimate)[1] <- 'gamma'
 #       names(map_estimate)[(n_pars - nCov + 1):n_pars] <- paste0("beta_", colnames(x$input_data_and_model_prior$X))
#	colnames(result[[1]]) <- names(map_estimate)
	result[[8]] <- map_estimate
	result[[9]] <- BIC
#####################################################################################
	 	names(result) <- c("mcmc_sample", "latent_status_censored", "complete_log_likelihood","swap_accept_rate",'all_cll_values', 'input_data_and_model_prior', 'log_posterior', 'map_estimate', 'BIC')
	 	class(result) <- c('list', 'bayesCureModel')
	 	return(result)
}



#' @export
summary.bayesCureModel <- function(object, burn = NULL, gamma_mix = TRUE, K_gamma = 3, K_max = 3, fdr = 0.1, covariate_levels = NULL, yRange = NULL, alpha = 0.1, ...){

	if(is.null(burn)){
		burn = floor(dim(object$mcmc_sample)[1]/3)
		cat(paste0('By default, I will discard the first one third of the mcmc sample as burn-in period.\n Alternatively, you may set the "burn" parameter to another value.'),'\n')
	}else{
		if(burn > dim(object$mcmc_sample)[1] - 1){stop('burn in period not valid.')}
		if(burn < 0){stop('burn in period not valid.')}		
	}
	y = object$input_data_and_model_prior$y
	X = object$input_data_and_model_prior$X
	Censoring_status = object$input_data_and_model_prior$Censoring_status
	promotion_time = object$input_data_and_model_prior$promotion_time
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


	log_inv_gamma_kernel <- function(x, a, b){
		if(min(c(x,a,b)) < 0){
		stop("input should be positive")
		}
		return(-b/x - (a+1) * log(x))
	}

	log_prior_gamma <- function(g, mu_g, s2_g){
		# mu_g is a_g
		# s2_g is b_g
		return(mu_g * log(s2_g) - 2*lgamma(mu_g) + (mu_g - 1)*log(abs(g)) - s2_g*abs(g))
	}
	
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

	
	a <- numeric(n_pars_f)

	log_S_p <- function(g, lambda, log_F, b, x){
		theta <- exp(x %*% b)
		return(-log(1 + g * theta * ct^{g*theta} * exp(log_F)^lambda)/g)

	}

	log_p0 <- function(g, b, x){
		theta <- exp(x %*% b)
		return(-(log(1 + g*theta*ct^(g*theta)))/g)
	}


	c_under <- 10^{-9}
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



	m <- dim(retained_mcmc)[1]
	ind <- 1:m
	
	logL <- logP <- numeric(m)
	tz <- 0
	if(length(unique(X[,1])) == 1){
		bIndices <- paste0('b',1:nCov - 1,'_mcmc')	
	}else{
		bIndices <- paste0('b',1:nCov,'_mcmc')
	}
	
	


	
	
	n <- dim(X)[1]
	n_parameters <- dim(x)[2] + 2 + n_pars_f
	BIC <- object$BIC
	logP <- object$log_posterior[-(1:burn)]
	map_estimate <- object$map_estimate
	# computation of DIC2, celeux 2006
#	D_theta_hat <- max(logL, na.rm = TRUE)
#	p_D <- -2*mean(logL, na.rm = TRUE) + 2*D_theta_hat
#	dic2 <- -2 * D_theta_hat + 2 * p_D
	hdis <- NULL
	if(length(logP) < 100){K_max = 1}
	trans_logP <- log(-min(logP)+logP+abs(mean(logP))+0.0001)
#	set.seed(1)
	hh <- Mclust(trans_logP, G = 1:K_max, verbose = FALSE)
	ind <- which(hh$classification == hh$G)
#	if(map_index %in% ind){
#		ind <- ind
#	}else{
#		ind <- sort(union(ind, map_index))
#	}
	main_mode_index <- ind

	if(is.null(alpha)==FALSE){
#	if(main_mode == FALSE){
		hdis <- fpf <- plot(object, alpha = alpha, plot_graphs = FALSE, burn = burn, index_of_main_mode = NULL)
#	}else{
#		hdis <- fpf <- plot(object, alpha = alpha, plot_graphs = FALSE, burn = burn, index_of_main_mode = main_mode_index)
#	}
	}


	dic_main <- NULL
#	if(main_mode){
		latent_cured_status <- 1 - colMeans(object$latent_status_censored[burn + ind,])
#	}else{
#		latent_cured_status <- 1 - colMeans(object$latent_status_censored[-(1:burn),])
#		cat('WARNING: main_mode is set to FALSE. If minor modes are present, the estimate of latent cured status may be conservative.', '\n')
#	}
	p_cured_given_tau <- p_cured_given_tau_values <- NULL 
	if(is.null(covariate_levels) == FALSE){
	if(is.matrix(covariate_levels) == FALSE){
		covariate_levels = as.matrix(t(covariate_levels))
	}
	nLevels <- dim(covariate_levels)[1]
	if(length(covariate_levels[1,])!=dim(X)[2]){
		stop("the length of distinct covariate_levels should be equal to the number of columns in the design matrix (X).")
	}
	if(is.null(yRange)){
		tau_values <- seq(min(y), max(y), length = 100)
	}else{
		tau_values = seq(min(yRange), max(yRange), length = 100)
	}
	p_cured_given_tau <- array(NA, dim = c(length(tau_values), nLevels))
	p_cured_given_tau_values <- array(data = 0, dim = c(length(tau_values), dim(retained_mcmc)[1], nLevels))



	
	
	for(iter in 1:dim(retained_mcmc)[1]){
		g = retained_mcmc[iter,'g_mcmc']
		lambda = retained_mcmc[iter,'lambda_mcmc']
		b = retained_mcmc[iter, bIndices]

		if(promotion_time$distribution == 'exponential'){
			lw <- log_weibull(y = tau_values, a1 = retained_mcmc[iter,'a1_mcmc'], 
			a2 = 1,  c_under = 1e-9)
		}

		if(promotion_time$distribution == 'weibull'){
			lw <- log_weibull(y = tau_values, a1 = retained_mcmc[iter,'a1_mcmc'], 
				a2 = retained_mcmc[iter,'a2_mcmc'],  c_under = 1e-9)
		}

		if(promotion_time$distribution == 'gamma'){
			lw <- log_gamma(y = tau_values, a1 = retained_mcmc[iter,'a1_mcmc'], 
				a2 = retained_mcmc[iter,'a2_mcmc'],  c_under = 1e-9)
		}

		if(promotion_time$distribution == 'gompertz'){
			lw <- log_gompertz(y = tau_values, a1 = retained_mcmc[iter,'a1_mcmc'], 
				a2 = retained_mcmc[iter,'a2_mcmc'],  c_under = 1e-9)
		}


		if(promotion_time$distribution == 'logLogistic'){
			lw <- log_logLogistic(y = tau_values, a1 = retained_mcmc[iter,'a1_mcmc'], 
				a2 = retained_mcmc[iter,'a2_mcmc'],  c_under = 1e-9)
		}

		if(promotion_time$distribution == 'lomax'){
			lw <- log_lomax(y = tau_values, a1 = retained_mcmc[iter,'a1_mcmc'], 
				a2 = retained_mcmc[iter,'a2_mcmc'],  c_under = 1e-9)
		}
		if(promotion_time$distribution == 'dagum'){
			lw <- log_dagum(y = tau_values, a1 = retained_mcmc[iter,'a1_mcmc'], a2 = retained_mcmc[iter,'a2_mcmc'], 
				a3 = retained_mcmc[iter,'a3_mcmc'], c_under = 1e-9)
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
			lw <- log_gamma_mixture(y = tau_values, a1 = a1, a2 = a2, p = p, c_under = 1e-9)
		}



						
		log_F = lw$log_F
		i <- 0
		for(tau in tau_values){
			i <- i + 1
			for(j in 1:nLevels){
			p_cured_given_tau_values[i, iter,j] <- exp(log_p0(g = g, 
							b = b,
							x = covariate_levels[j,]) - 
						log_S_p(g = g, 
							lambda = lambda,
							log_F = log_F[i], 
							b = b,
							x = covariate_levels[j,]
							)
						)
			}
			#p_cured_given_tau[i] <- p_cured_given_tau[i] +  p_cured_given_tau_values[i, iter]           				
		}
		if(iter == map_index){
			for(j in 1:nLevels){
				p_cured_given_tau[,j] <- p_cured_given_tau_values[,iter,j]
			}
		}
	}
	}
#	p_cured_given_tau <- p_cured_given_tau/length(ind)

	

#############################FDR control
	posterior_probs <- 1 - latent_cured_status
        p <- matrix(1 - posterior_probs, ncol = 1)
        perm <- order(p,decreasing = TRUE)
        orderedP <- p[perm,]
	K <- dim(p)[1]
	myList <-  1 - orderedP[1]
	k <- 1 
	criterion <- myList
	while ((criterion < fdr) & (k < length(orderedP))){
		k <- k + 1 
		myList <- myList + 1 - orderedP[k]
		criterion <- myList/k
	}
	if(k > 1){
		ind <- perm[1:(k-1)]
	}else{
		ind <- c()
	}
	cured_at_given_FDR = rep('susceptible', table(object$input_data_and_model_prior$Censoring_status)[1])
	cured_at_given_FDR[ind] = 'cured'
	names(cured_at_given_FDR) <- which(object$input_data_and_model_prior$Censoring_status == 0)
#############################################################################
	
#	hdis <- NULL
	results <- vector('list', length = 6)
#	results[[1]] <- logP
#	results[[2]] <- BIC
#	results[[3]] <- dic2	
#	results[[4]] <- dic_main


	results[[1]] <- map_estimate
	results[[2]] <- hdis
	results[[3]] <- latent_cured_status
	results[[4]] <- cured_at_given_FDR
	results[[6]] <- main_mode_index
	if(is.null(covariate_levels)){
cat('                           MCMC summary','\n')
	myDF <- data.frame(MAP_estimate = round(map_estimate,2), HPD_interval = character(n_parameters))
	for(i in 1:n_parameters){
		nIntervals <- dim(hdis[[i]])[1]
		myInt <- paste0('(', paste0(round(hdis[[i]][1,], 2), collapse=', '), ")")
		if(nIntervals > 1){
			for(j in 2:nIntervals){
				myInt <- paste0(myInt, "U(" ,paste0(round(hdis[[i]][j,], 2), collapse=', '), ")")
			}
		}
		myDF$HPD_interval[i] <- myInt
	}
	print(myDF)
	cat('\n')

cat(paste0('Among ', length(latent_cured_status) ,' cencored observations, I found ', sum(cured_at_given_FDR == 'cured'), ' cured subjects (FDR = ', fdr, ').'))
	cat('\n')
	}
	if(is.null(covariate_levels)==FALSE){
		p_cured_output <- vector('list', length = 6)
		p_cured_output[[1]] <- 	p_cured_given_tau_values
		p_cured_output[[2]] <- p_cured_given_tau
		p_cured_output[[3]] <- tau_values
		p_cured_output[[4]] <- covariate_levels
		p_cured_output[[5]] <- main_mode_index	
		if(is.null(colnames(object$input_data_and_model_prior$X))){
			colnames(object$input_data_and_model_prior$X) <- paste0('x_', 1:nCov)
		}
		p_cured_output[[6]] <- colnames(object$input_data_and_model_prior$X)
		names(p_cured_output) <- c('mcmc', 'map', 'tau_values', 'covariate_levels', 'index_of_main_mode', 'Xnames')
		results[[5]] <- p_cured_output
	}
	names(results) <- c('map_estimate', 'highest_density_indervals', 'latent_cured_status', 'cured_at_given_FDR', 'p_cured_output', 'index_of_main_mode')
	return(results)

}



#' @export

plot.bayesCureModel <- function(x, burn = NULL, alpha = 0.05, gamma_mix = TRUE, K_gamma = 5, plot_graphs = TRUE, bw = 'nrd0', what = NULL, p_cured_output = NULL, index_of_main_mode = NULL, ...){
	retained_mcmc = x$mcmc_sample
	map_estimate = x$map_estimate
#	if(plot_graphs){
#	oldpar <- par(no.readonly = TRUE)
#	on.exit(par(oldpar)) 
#	}

	if(is.null(burn)){
		burn = floor(dim(x$mcmc_sample)[1]/3)
	}else{
		if(burn > dim(x$mcmc_sample)[1] - 1){stop('burn in period not valid.')}
		if(burn < 0){stop('burn in period not valid.')}		
	}
	if(burn > 0){
	retained_mcmc = retained_mcmc[-(1:burn),]	
	}
	nPars <- dim(retained_mcmc)[2]
	n_pars_f <- dim(x$input_data_and_model_prior$promotion_time$prior_parameters)[1]
	if(x$input_data_and_model_prior$promotion_time$distribution == 'gamma_mixture'){
		K = x$input_data_and_model_prior$promotion_time$K	
		n_pars_f = K * n_pars_f + K - 1
	}
	hdis <- vector('list', length = nPars)
	m <- dim(retained_mcmc)[1]
	if(is.null(index_of_main_mode)==FALSE){
		ind <- index_of_main_mode
		K_gamma = 1
	}else{
		ind <- 1:m
	}

	myXlim <- matrix(NA, nPars, 2)
	for(i in 1:nPars){
		myXlim[i,] <- quantile(retained_mcmc[ind,i],probs = c(0.001,0.999))# c(-5,2)
	}

	
	varnames <- numeric(nPars)
	varnames[1:2] <- as.expression(c(
		bquote(gamma), 
		bquote(lambda)))
	b_ind <- as.numeric(unlist(lapply(strsplit(unlist(lapply(strsplit(colnames(retained_mcmc)[-(1:(2+n_pars_f))], 
		split = '_'), function(x)x[1])), split = 'b'), function(x)x[2])))

	for(i in 3:(2+n_pars_f)){
		varnames[i] <- as.expression(bquote(alpha[.(i-2)]))
	}


	for(i in (3+n_pars_f):nPars){
		varnames[i] <- as.expression(bquote(beta[.(b_ind[i-(2+n_pars_f)])]))
	}


	if(x$input_data_and_model_prior$promotion_time$distribution == 'gamma_mixture'){
	# get mixing proportions
		wRange <- (2 + 2*K+1):(2+n_pars_f)
		i <- 0
		for (j in wRange){
			i <- i + 1
			varnames[j] <- as.expression(bquote(w[.(i)]))
		}
	}		


	hdi_alpha = alpha
	if(is.null(what)){
		what <- 1:nPars
	}
	if(what[1] == 'cured_prob'){
		plot(
			p_cured_output$tau_values, 
			p_cured_output$map[,1], 
			xlab = "t", 
			ylab = bquote(hat("P")[theta]*"(cured|T">="t)"), 
			type = 'n', 
			...
		)
		xti <- par('xaxp')
		rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "gray90")
		for(j in (0:10)/10){
		abline (h = j, lwd = 2, col = 'white')
		}
		xxseq <- seq(xti[1], xti[2], by = xti[3])
		for(j in xxseq){
			abline (v = j, lwd = 2, col = 'white')
		}
		nLevels <- dim(p_cured_output$covariate_levels)[1]
		for(j in 1:nLevels){
			points(
				p_cured_output$tau_values, 
				p_cured_output$map[,j],
				type = 'l', 
				col = j + 1,
				...
			)
			
			points(p_cured_output$tau_values, p_cured_output$map[,j], type = 'l', lwd = 2, col = j+1)
			if(is.null(alpha) == FALSE){
			p_cured_given_tau_low <- p_cured_given_tau_up <- numeric(length(p_cured_output$tau_values))
			for(k in 1:length(p_cured_output$tau_values)){
				p_cured_given_tau_low[k] <- hdi(p_cured_output$mcmc[k,p_cured_output$index_of_main_mode,j], credMass = 1 - alpha)[1]
				p_cured_given_tau_up[k] <- hdi(p_cured_output$mcmc[k,p_cured_output$index_of_main_mode,j], credMass = 1 - alpha)[2]        
			}
			Ti <- length(p_cured_output$tau_values)
			rgb.val <- col2rgb(j+1)/255
			polygon(x = c(p_cured_output$tau_values, p_cured_output$tau_values[Ti:1]), 
				y = c(p_cured_given_tau_low, p_cured_given_tau_up[Ti:1]), 
				col= rgb(rgb.val[1], rgb.val[2], rgb.val[3],0.1), 
#				col = j+1,
#				density = 10,
#				angle = 90*(j+1),
				border = NA
				)
			}
		}
		lText <- paste0(apply(round(p_cured_output$covariate_levels,2), 1, function(y)paste0(y, collapse=', ')))
		legend('bottomright', col = 2:(nLevels+1), lty = 1, 
			title = paste0('covariate levels\n',paste0(p_cured_output$Xnames, collapse=', ')), lText)
#		p_cured_given_tau_low <- p_cured_given_tau_up <- numeric(length(p_cured_output$tau_values))
#		for(j in 1:nLevels){
#		for(k in 1:length(p_cured_output$tau_values)){
#			p_cured_given_tau_low[k] <- hdi(p_cured_output$mcmc[k,p_cured_output$index_of_main_mode,j], credMass = 1 - alpha)[1]
#			p_cured_given_tau_up[k] <- hdi(p_cured_output$mcmc[k,p_cured_output$index_of_main_mode,j], credMass = 1 - alpha)[2]        
#		}
#		polygon(x = c(p_cured_output$tau_values, p_cured_output$tau_values[Ti:1]), 
#			y = c(p_cured_given_tau_low, p_cured_given_tau_up[Ti:1]), 
#			lty = 2)
#		}

	
	}else{
	for(i in what){
#		pdf(file = paste0("../img/recidivism_new_data_parameter_",i,".pdf"), width = 12, height = 3)
		if(i == 1){
			shouldIask = FALSE
		}else{shouldIask = TRUE; if(length(what) == 1){shouldIask = FALSE} }
		if(plot_graphs & length(what) > 1){
		par(ask = shouldIask)	
		}
		xxx <- retained_mcmc[ind,i]

		myD <- density(xxx, bw = bw)
#		if(i < 5){
#			myD <- density(x, bw = 'bcv')			
#		}
		if(i == 1){
		if(gamma_mix){
		fit <- Mclust(xxx,G=1:K_gamma, modelNames = "V", verbose = FALSE)
		k <- fit$G
		if( k > 1){
			multMode = TRUE
			mu <- fit$parameters$mean
			w <- fit$parameters$pro
			s2 <- fit$parameters$variance$sigmasq
			dd <- range(xxx) + 0.01*c(-1,1)
			xvals <- seq(dd[1],dd[2], length = 512)
			yvals <- w[1]*dnorm(xvals,mean = mu[1], sd = sqrt(s2[1]))
			for(j in 2:k){
				yvals <- yvals + w[j]*dnorm(xvals,mean = mu[j], sd = sqrt(s2[j]))
			}
			myD <- vector("list", length = 2)
			names(myD) <- c('x','y')
			myD$x <- xvals
			myD$y <- yvals
			class(myD) <- 'density'
		}
		
		}}
		if(is.null(index_of_main_mode)==FALSE){
		allow_split = FALSE
		}else{
		allow_split = TRUE
		}
		#if(i %in% c(2,3,4)){allow_split = FALSE}
		hdi_95 <- hdi(myD,allowSplit=allow_split, credMass = 1 - hdi_alpha)				


		if(is.null(dim(hdi_95))){hdi_length = diff(hdi_95)}else{
		hdi_length = sum(apply(hdi_95,1,diff))}
		hdis[[i]] <- hdi_95
		if(is.null(dim(hdi_95))){hdi_95 = matrix(hdi_95,1,2)}
		if(plot_graphs){
		plot(myD, xlab = varnames[i], xlim = myXlim[i,], ...)
		}
#		if(i == 1){
#					if(plot_graphs){
#			legend('topleft', c('estimate (map)', paste0(100*(1-hdi_alpha), '% HDI')), 
#				col = c('red','papayawhip'),lty = 1,lwd = c(1,10))
#				}
#		}
		#hist(mcmcmc16$mcmc_sample[ind,i], xlab = varnames[i], main = '');abline(v = truePars[i], col = 1, lwd = 2)

		for(j in 1:dim(hdi_95)[1]){
			#abline(v = hdi_95[j,], lty = 2, col = 1+j)
			x_dens1 <- which.min(abs(hdi_95[j,1] - myD$x))
			x_dens2 <- which.min(abs(hdi_95[j,2] - myD$x))
			y_dens <- myD$y[c(x_dens1,x_dens2)]
					if(plot_graphs){
			polygon(c(myD$x[c(x_dens1:x_dens2,x_dens2:x_dens1)]),
				c(rep(0,x_dens2 - x_dens1+1),myD$y[c(x_dens2:x_dens1)]),col='papayawhip')
				}
			#points(c(hdi_95[j,1],hdi_95[j,1]),c(0,y_dens[1]), type = 'l')	
			#points(c(hdi_95[j,2],hdi_95[j,2]),c(0,y_dens[2]), type = 'l')			
		}
		hdi_length = sum(apply(hdi_95,1,diff))
		if(is.null(map_estimate) == FALSE){
				if(plot_graphs){
			abline(v = c(map_estimate[i]), col = 2, lwd = 2, lty = 1)
				}
		}
#		abline(v = truePars[i], col = 'green', lwd = 2, lty = 2)
#		dev.off()
		
	}
	if(plot_graphs){
	par(ask=FALSE)
	}
	names(hdis) <- colnames(x$mcmc_sample)
	if(plot_graphs == FALSE){
	return(hdis)	
	}
	}

}

#' @export
print.bayesCureModel <- function(x, ...){
       if( 'bayesCureModel' %in% class(x) ){
                cat("\n")
                cat(paste0("* Run information:"),"\n")
                cat(paste0("      Fitted model: `", x$input_data_and_model_prior$promotion_time$distribution,"'\n"))
                cat(paste0("      BIC: ", round(x$BIC, 3),"\n"))
                cat(paste0("      MCMC cycles: ", dim(x$mcmc_sample)[1],"\n"))
                cat(paste0("      Number of parallel heated chains: ", 1+length(x$swap_accept_rate),"\n"))
                ss <- summary(x$swap_accept_rate)
                cat(paste0("      Swap rates of adjacent chains: ","\n"))                                                
                print(ss[c(1,3,6)], digits = 2)
                cat("\n")
                cat(paste0("* Maximum A Posteriori (MAP) estimate of parameters"),"\n")
                map <- x$map_estimate
                n_pars <- length(x$map_estimate)
                nCov <- dim(x$input_data_and_model_prior$X)[2]              
                txt <- colnames(x$input_data_and_model_prior$X)
		names(map)[(n_pars - nCov + 1):n_pars] <- paste0(names(x$map_estimate)[(n_pars - nCov + 1):n_pars], ' [',txt,']')
                myD <- t(t(map))
                colnames(myD) <- 'MAP estimate'
		print(myD)
                cat('\n')

        }else{
                cat(paste("    The input is not in class `bayesCureModel`"),'\n')
        }
}


