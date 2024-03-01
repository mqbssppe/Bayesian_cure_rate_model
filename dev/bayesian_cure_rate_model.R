library("Rcpp")
library("RcppArmadillo")
sourceCpp("cll_reparameterization_fix.cpp")
library('coda')
library('doParallel')
library('foreach')
library('HDInterval')
library('mclust')

#myData must contain the columns "Covariate1", "Covariate2" and "Y"
#We also need  D1, D0, and I_sim
#dercomp takes on {1,2,3,4,5,6,7}->{g,lambda,a1,a2,b0,b1,b2} depending on which derivative must be computed 

derv_complete_log_likelihood_fotis <- function(y, X, Censoring_status, g, lambda, a1, a2, b, I_sim){
	X <- as.matrix(X)
	nCov <- dim(X)[2] 
	n <- dim(X)[1]
	nPars <- 4 + nCov
	D1 <- which(Censoring_status == 1)	# fixed
	D0 <- which(Censoring_status == 0)	# fixed
	
	Dp0b <- array(data = NA, dim = c(n,nCov))
	Dfpb <- array(data = NA, dim = c(n,nCov))
	DSpb <- array(data = NA, dim = c(n,nCov))
	parg1 <- matrix(NA, length(D1), nCov)
	parg2 <- matrix(NA, length(D0), nCov)	
	parg3 <- matrix(NA, length(D0), nCov)
	com1 <- com2 <- com3 <- numeric(nCov)
	
	gradientVec <- numeric(nPars)
	ct = exp(exp(-1))
	x <- X
	tau <- y


	theta1 <- exp(x %*% b)
	p00 <- (1 + g*theta1*ct^(g*theta1))^(-1/g)
	FW <- 1 - exp(-(a1*tau)^a2)
	fW <- (a2/tau)*exp(-(a1*tau)^a2)*(a1*tau)^a2
	Sp <- (1 + g *theta1* ct^(g*theta1) * FW^lambda)^(-1/g)
	fp <- lambda*ct^(g*theta1)*theta1*(1 + g *theta1* ct^(g*theta1) * FW^lambda)^(-1/g-1)*FW^(lambda-1)*fW
	fp0 <- (1 + g * theta1*ct^(g*theta1) * FW^lambda)^(-1/g-1)

	DFa1 <- a2*tau*((a1*tau)^(a2-1))*exp(-(a1*tau)^a2) 
	DFa2 <- log(a1*tau)*((a1*tau)^a2)*exp(-(a1*tau)^a2) 
	Dfa1 <- ((a1*tau)^(a2-1))*a2^2*exp(-(a1*tau)^a2)*(1-(a1*tau)^a2) 
	Dfa2 <- (1/tau)*exp(-(a1*tau)^a2)*(a1*tau)^a2*(1+a2*log(a1*tau)-a2*(a1*tau)^a2*log(a1*tau))

	Dp0g<-   p00*(1/g)*(log(1 + g*theta1*ct^(g*theta1))/g-(theta1*ct^(g*theta1))*(1+g*exp(-1)*theta1)/(1 + g*theta1*ct^(theta1*g)))
	for(j in 1:nCov){
		Dp0b[,j] <- -p00*ct^(g*theta1)*((1+g*exp(-1)*theta1)/(1 + g*theta1*ct^(theta1*g)))*theta1*x[,j]
	}


	Dlogfp0lam <- (-1/g-1)*(g*theta1*ct^(g*theta1)*FW^lambda*log(FW))/(1 + g *theta1* ct^(g*theta1) * FW^lambda)
	Dlogfp0gam <-(1/g^2)*log((1 + g *theta1* ct^(g*theta1) * FW^lambda))+(-1/g-1)*theta1* ct^(g*theta1) * FW^lambda*(1+g*exp(-1)*theta1)/(1 + g*theta1*ct^(theta1*g)*FW^lambda)
	Dlogfp0a1 <- (-1/g-1)*(lambda*g*theta1*ct^(g*theta1)*FW^(lambda-1)*DFa1)/(1+g*theta1*ct^(g*theta1)*FW^(lambda))
	Dlogfp0a2 <- (-1/g-1)*(lambda*g*theta1*ct^(g*theta1)*FW^(lambda-1)*DFa2)/(1+g*theta1*ct^(g*theta1)*FW^(lambda))
	Dlogfp0theta <- (-1/g-1)*(g*ct^(g*theta1)*FW^(lambda))*(1+g*exp(-1)*theta1)/(1+g*theta1*ct^(g*theta1)*FW^(lambda))

	Dfplam <- fp*(lambda^(-1)+log(FW)+Dlogfp0lam)
	Dfpg <- fp*(theta1*exp(-1)+Dlogfp0gam)
	Dfpa1 <- fp*((lambda-1)*DFa1/FW+Dfa1/fW+Dlogfp0a1)
	Dfpa2 <- fp*((lambda-1)*DFa2/FW+Dfa2/fW+Dlogfp0a2)
	for(j in 1:nCov){
		Dfpb[,j] <- fp*(theta1^(-1)+g*exp(-1)+Dlogfp0theta)*theta1*x[,j]
	}

	DSpg <- Sp*
	(1/g^2*log((1 + g *theta1* ct^(g*theta1) * FW^lambda))-((theta1* ct^(g*theta1) * FW^lambda)/(g))*(1+g*theta1*exp(-1))/((1 + g *theta1* ct^(g*theta1) * FW^lambda)))
	DSplam <- -fp0*theta1*ct^(g*theta1)*FW^lambda*log(FW)
	DSpa1 <- -fp0*theta1*ct^(g*theta1)*lambda*FW^(lambda-1)*(DFa1)
	DSpa2 <- -fp0*theta1*ct^(g*theta1)*lambda*FW^(lambda-1)*(DFa2)
	for(j in 1:nCov){
		DSpb[,j] <- -fp0*ct^(g*theta1)*FW^(lambda)*(1+g*theta1*exp(-1))*theta1*x[,j]
	}


	parg1c<-fp[D1]
	parg2c<-p00[D0]
	parg3c<-Sp[D0]

	parg11<-Dfpg[D1]
	parg21<-Dp0g[D0]
	parg31<-DSpg[D0]

	parg12<-Dfplam[D1]
	parg22<-0
	parg32<-DSplam[D0]  

	parg13<-Dfpa1[D1]
	parg23<-0
	parg33<-DSpa1[D0]  

	parg14<-Dfpa2[D1]
	parg24<-0
	parg34<-DSpa2[D0]    

	for(j in 1:nCov){
		parg1[, j] <- Dfpb[D1, j]
		parg2[, j] <- Dp0b[D0, j]
		parg3[, j] <- DSpb[D0, j]   
	}

#	parg15<-Dfpb0[D1]
#	parg25<-Dp0b0[D0]
#	parg35<-DSpb0[D0]   
	  
#	parg16<-Dfpb1[D1]
#	parg26<-Dp0b1[D0]
#	parg36<-DSpb1[D0] 

#	parg17<-Dfpb2[D1]
#	parg27<-Dp0b2[D0]
#	parg37<-DSpb2[D0]    


	com11 <- sum(parg11/parg1c)
	com21 <- sum(parg21/parg2c)
	com31 <- sum(I_sim[D0]*((parg31*parg2c-parg3c*parg21)/(parg2c*(parg3c-parg2c))))

	com12 <- sum(parg12/parg1c)
	com22 <- sum(parg22/parg2c)
	com32 <- sum(I_sim[D0]*((parg32*parg2c-parg3c*parg22)/(parg2c*(parg3c-parg2c))))

	com13 <- sum(parg13/parg1c)
	com23 <- sum(parg23/parg2c)
	com33 <- sum(I_sim[D0]*((parg33*parg2c-parg3c*parg23)/(parg2c*(parg3c-parg2c))))

	com14 <- sum(parg14/parg1c)
	com24 <- sum(parg24/parg2c)
	com34 <- sum(I_sim[D0]*((parg34*parg2c-parg3c*parg24)/(parg2c*(parg3c-parg2c))))

	for(j in 1:nCov){
		com1[j] <- sum(parg1[,j]/parg1c)
		com2[j] <- sum(parg2[,j]/parg2c)
		com3[j] <- sum(I_sim[D0]*((parg3[,j]*parg2c-parg3c*parg2[,j])/(parg2c*(parg3c-parg2c))))
	}

#	com15 <- sum(parg15/parg1c)
#	com25 <- sum(parg25/parg2c)
#	com35 <- sum(I_sim[D0]*((parg35*parg2c-parg3c*parg25)/(parg2c*(parg3c-parg2c))))

#	com16 <- sum(parg16/parg1c)
#	com26 <- sum(parg26/parg2c)
#	com36 <- sum(I_sim[D0]*((parg36*parg2c-parg3c*parg26)/(parg2c*(parg3c-parg2c))))

#	com17 <- sum(parg17/parg1c)
#	com27 <- sum(parg27/parg2c)
#	com37 <- sum(I_sim[D0]*((parg37*parg2c-parg3c*parg27)/(parg2c*(parg3c-parg2c))))


	gradientVec[1:4] <- c(
		com11 + com21 + com31, 
		com12 + com22 + com32,
		com13 + com23 + com33,
		com14 + com24 + com34)

	gradientVec[-(1:4)] <- com1 + com2 + com3


#	gradientVec[5] <- com15 + com25 + com35
#	gradientVec[6] <- com16 + com26 + com36
#	gradientVec[7] <- com17 + com27 + com37


		
	return(gradientVec)
}


#################################################
#########################test####################


log_prior_gradient <- function(g, lambda, a1, a2, b ,mu_g, s2_g,  a_l, b_l, a_1, b_1, a_2, b_2, mu_b, Sigma){
	
	log_gradient <- numeric(4+length(b))
	log_gradient[1] <- if(g < 0){
		(mu_g - 1)/g + s2_g
	}else{
		(mu_g - 1)/g - s2_g
	}
	log_gradient[2] <- -(a_l+1)/lambda + b_l/(lambda^2)
	log_gradient[3] <- -(a_1+1)/a1 + b_1/(a1^2)
	log_gradient[4] <- -(a_2+1)/a2 + b_2/(a2^2)		
	log_gradient[5:(4+length(b))] <- (mu_b - b)/diag(Sigma)
	return(log_gradient)
}


cure_rate_metropolis_hastings_cpp <- function( y, X, Censoring_status,  m, alpha = 1,
				mu_g = 0.2, s2_g = 0.1, a_l = 2.001, b_l = 1, 
				a_1 = 2.001, b_1 = 1, a_2 = 2.001, b_2 = 1, 
				mu_b = NULL, Sigma = NULL,
				g_prop_sd = 0.045, 
				lambda_prop_scale = 0.03, 
				a1_prop_scale = 0.2, 
				a2_prop_scale = 0.03, 
				b_prop_sd = NULL, 
						initialValues = NULL, 
						plot = TRUE,
						verbose = TRUE,
						tau_mala = 0.000015, mala = 0.33
					){
#	nCov <- dim(myData)[2] - 1
	X <- as.matrix(X)
	nCov <- dim(X)[2]
	nPars <- nCov + 4
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
	#	a1 ~ IG(a_1, b_1)
	#	a2 ~ IG(a_2, b_2)
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


	g_mcmc <- lambda_mcmc <- a1_mcmc <- a2_mcmc <- numeric(m)
	b_mcmc <- array(data = NA, dim = c(m,nCov))
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
		a1 = rgamma(1, shape = 1, rate = 1)
		a2 =  rgamma(1, shape = 1, rate = 1)
	}else{
		g = initialValues[['g_mcmc']]
		b <- numeric(nCov)
		for(j in 1:nCov){
			b[j] <- initialValues[[4+j]]
		}
		lambda = initialValues[['lambda_mcmc']]								
		a1 = initialValues[['a1_mcmc']]
		a2 = initialValues[['a2_mcmc']]						
	}
	


	g_mcmc[1] <- g
	b_mcmc[1,] <- b
	a1_mcmc[1] <- a1
	a2_mcmc[1] <- a2
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


	cll <- complete_log_likelihood(y = y, X = X, Censoring_status = Censoring_status, g = g, lambda = lambda, a1 = a1, a2 = a2, b = b, I_sim = I_sim, alpha = alpha)
	cllValues[1] <- cll$cll
	##########################################################edw stamatisa################################################3
#	print("0")
#	print(cllValues[1])
#	print(c(g, lambda, a1, a2, b0, b1, b2))

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
				alpha*log_inv_gamma_kernel(lambda, a_l, b_l) +
				alpha*log_inv_gamma_kernel(a1, a_1, b_1) +
				alpha*log_inv_gamma_kernel(a2, a_2, b_2) -
				alpha*0.5 * mahalanobis(b, mu_b, Sigma)^2

	lpd[1] <- log_prior_density

	mala_accept_rate <- 0
       g_accept_rate <- lambda_accept_rate <- a1_accept_rate <- a2_accept_rate <- b_accept_rate <- 0



	if(plot){
	dev.new(width=18, height=8, unit="in")
	}
	######################################
	# MCMC sampler
	######################################
	if(m > 1){
	mh_iter <- 0
	mala_iter <- 0
	for(iter in 2:m){

		joe <- runif(1)
		if(joe < 1 - mala){
	   #------------------------------------------------
		#       Metropolis-Hastings proposal for gamma: normal proposal distribution
		#------------------------------------------------
		        g_prop <- rnorm(1, mean = g, sd = g_prop_sd)
		        cll_prop <- complete_log_likelihood(y = y, X = X, Censoring_status = Censoring_status, 
		        		g = g_prop, lambda = lambda, a1 = a1, a2 = a2, b = b, I_sim = I_sim, alpha = alpha)
		        cll_diff <- cll_prop$cll - cll$cll
	#               print("1")
	#               print(cll_diff)
		        log_prior_diff <- alpha*diff(log_gamma_mix_kernel(c(g, g_prop), mu_g, s2_g))
		        log_proposal_ratio <- 0
		        mh_ar <- cll_diff + log_prior_diff + log_proposal_ratio
		#       print("gamma")
		#       print(cll_prop)
		        if(is.na(mh_ar)==F){
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
		        cll_prop <- complete_log_likelihood(y = y, X = X, Censoring_status = Censoring_status, 
		        	g = g, lambda = lambda_prop, a1 = a1, a2 = a2, 
		        	b = b, I_sim = I_sim, alpha = alpha)
		        cll_diff <- cll_prop$cll - cll$cll
	#               print("2")
	#               print(cll_diff)         
		        log_prior_diff <- alpha*diff(log_inv_gamma_kernel(c(lambda, lambda_prop), a_l, b_l))
		        log_proposal_ratio <- log(lambda_prop) - log(lambda)
		        mh_ar <- cll_diff + log_prior_diff + log_proposal_ratio
		#       print("lambda")
		#       print(cll_prop)

		        if(is.na(mh_ar)==F){
		        if(log(runif(1)) < mh_ar){
		                lambda = lambda_prop
		                cll <- cll_prop
		                lambda_accept_rate <- lambda_accept_rate + 1
		        }
		        }
		        lambda_mcmc[iter] <- lambda
		
	     #------------------------------------------------
		#       Metropolis-Hastings proposal for a1: lognormal proposal distribution
		#------------------------------------------------
		        a1_prop = rlnorm(1, meanlog = log(a1), sdlog = a1_prop_scale)
		        cll_prop <- complete_log_likelihood(y = y, X = X, Censoring_status = Censoring_status, 
		        	g = g, lambda = lambda, a1 = a1_prop, a2 = a2, 
		        	b = b, I_sim = I_sim, alpha = alpha)
		        cll_diff <- cll_prop$cll - cll$cll
	#               print("3")
	#               print(cll_diff)
		        
		        log_prior_diff <- alpha*diff(log_inv_gamma_kernel(c(a1, a1_prop), a_1, b_1))
		        log_proposal_ratio <- log(a1_prop) - log(a1)
		        mh_ar <- cll_diff + log_prior_diff + log_proposal_ratio
		#       print("a1")
		#       print(cll_prop)
		        if(is.na(mh_ar)==F){
		        if(log(runif(1)) < mh_ar){
		                a1 = a1_prop
		                cll <- cll_prop
		                a1_accept_rate <- a1_accept_rate + 1
		        }
		        }
		        a1_mcmc[iter] <- a1

		#------------------------------------------------
		#       Metropolis-Hastings proposal for a2: lognormal proposal distribution
		#------------------------------------------------
		        a2_prop = rlnorm(1, meanlog = log(a2), sdlog = a2_prop_scale)
		        cll_prop <- complete_log_likelihood(y = y, X = X, Censoring_status = Censoring_status, 
		        	g = g, lambda = lambda, a1 = a1, a2 = a2_prop, 
		        	b = b, I_sim = I_sim, alpha = alpha)
		        cll_diff <- cll_prop$cll - cll$cll
	#               print("4")
	#               print(cll_diff)
		        
		        log_prior_diff <- alpha*diff(log_inv_gamma_kernel(c(a2, a2_prop), a_2, b_2))
		        log_proposal_ratio <- log(a2_prop) - log(a2)
		        mh_ar <- cll_diff + log_prior_diff + log_proposal_ratio
		#       print("a2")
		#       print(cll_prop)
		        if(is.na(mh_ar)==F){
		        if(log(runif(1)) < mh_ar){
		                a2 = a2_prop
		                cll <- cll_prop
		                a2_accept_rate <- a2_accept_rate + 1
		        }
		        }
		        a2_mcmc[iter] <- a2

	       #------------------------------------------------
		#       Metropolis-Hastings proposal for b = (b0,b1,b2): multivariate normal proposal distribution
		#------------------------------------------------
		        b_prop <- b + b_prop_sd*rnorm(nCov)
		        cll_prop <- complete_log_likelihood(y = y, X = X, Censoring_status = Censoring_status, 
		        	g = g, lambda = lambda, a1 = a1, a2 = a2, 
		        	b= b_prop, I_sim = I_sim, alpha = alpha)
		        cll_diff <- cll_prop$cll - cll$cll
	#               print("5")
	#               print(cll_diff)
		        
		        log_prior_diff <- -alpha*0.5 * mahalanobis(b_prop, mu_b, Sigma)^2 + alpha*0.5 * mahalanobis(b, mu_b, Sigma)^2
		        log_proposal_ratio <- 0
		        mh_ar <- cll_diff + log_prior_diff + log_proposal_ratio
		#       print("b")
		#       print(cll_prop)
		        if(is.na(mh_ar)==F){
		        if(log(runif(1)) < mh_ar){
		                b = b_prop
		                cll <- cll_prop
		                b_accept_rate <- b_accept_rate + 1
		        }
		        }
		        b_mcmc[iter,] <- b
		
			mh_iter <- mh_iter + 1
		}else{
			theta_previous <- c(g, lambda, a1, a2, b)		


			gradient_vec <- derv_complete_log_likelihood_fotis(y = y, X = X, Censoring_status = Censoring_status,
				theta_previous[1],theta_previous[2],theta_previous[3],
				theta_previous[4],theta_previous[-(1:4)],I_sim)

			lpg <- log_prior_gradient(g, lambda, a1, a2, b , mu_g, s2_g,  a_l, b_l, a_1, b_1, a_2, b_2, mu_b, Sigma)

			gradient_vec <- alpha*gradient_vec + alpha*lpg

			
			mean_prop <- theta_previous + tau_mala*gradient_vec
			sd_prop <- sqrt(2*tau_mala)
			theta_prop <- rnorm(nPars, mean_prop, sd_prop)
			rejectedMove = FALSE
			if(sum(is.na(theta_prop))){rejectedMove = TRUE}else{
			if(min(theta_prop[2:4]) < 0){rejectedMove = TRUE}}
			if(rejectedMove == FALSE){
				log_prop_density <- sum(dnorm(theta_prop, mean_prop, sd = sd_prop, log = TRUE))
				gradient_vec_prop <- derv_complete_log_likelihood_fotis(y = y, X = X, 
					Censoring_status = Censoring_status, 
					theta_prop[1],theta_prop[2],theta_prop[3], theta_prop[4],
					theta_prop[-(1:4)],I_sim)
				
				lpg <- log_prior_gradient(theta_prop[1],theta_prop[2],theta_prop[3],
					theta_prop[4], theta_prop[-(1:4)] , 
					mu_g, s2_g,  a_l, b_l, a_1, b_1, a_2, b_2, mu_b, Sigma)

				gradient_vec_prop <- alpha*gradient_vec_prop + alpha*lpg
				
				mean_prop_rev <- theta_prop + tau_mala * gradient_vec_prop
				log_prop_density_reverse <- sum(dnorm(theta_previous, mean_prop_rev, sd = sd_prop, log = TRUE))
				cll_prop <- complete_log_likelihood(y = y, X = X, Censoring_status = Censoring_status, 
					g = theta_prop[1], lambda = theta_prop[2], 
					a1 = theta_prop[3], a2 = theta_prop[4], b = theta_prop[-(1:4)], I_sim = I_sim, alpha = alpha)
				cll_diff <- cll_prop$cll - cll$cll
				log_proposal_ratio <- log_prop_density_reverse - log_prop_density
				log_prior_diff <- diff(dnorm(c(g, theta_prop[1]), mean = mu_g, sd = sqrt(s2_g), log = TRUE))
				log_prior_diff <- log_prior_diff + diff(log_inv_gamma_kernel(c(lambda, theta_prop[2]), a_l, b_l))
				log_prior_diff <- log_prior_diff + diff(log_inv_gamma_kernel(c(a1, theta_prop[3]), a_1, b_1))
				log_prior_diff <- log_prior_diff + diff(log_inv_gamma_kernel(c(a2, theta_prop[4]), a_2, b_2))
				b_prop <- theta_prop[-(1:4)]
				log_prior_diff <- log_prior_diff -0.5 * mahalanobis(b_prop, mu_b, Sigma)^2 + 0.5 * mahalanobis(b, mu_b, Sigma)^2
				mh_ar <- cll_diff + alpha*log_prior_diff + log_proposal_ratio

				if(is.na(mh_ar)==F){
					if(log(runif(1)) < mh_ar){
						g = theta_prop[1]
						lambda = theta_prop[2]
						a1 = theta_prop[3]
						a2 = theta_prop[4]
						b = b_prop
						cll <- cll_prop
						mala_accept_rate <- mala_accept_rate + 1
					}
				}

			}
			g_mcmc[iter] <- g
			lambda_mcmc[iter] <- lambda
			a1_mcmc[iter] <- a1
			a2_mcmc[iter] <- a2
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
		cll <- complete_log_likelihood(y = y, X = X, Censoring_status = Censoring_status, 
				g = g, lambda = lambda, a1 = a1, a2 = a2, b = b, I_sim = I_sim, alpha = alpha)
#		print("6")
#		print(cll$cll)

		log_prior_density <- alpha*log_gamma_mix_kernel(g,mu_g, s2_g)  + 
					alpha*log_inv_gamma_kernel(lambda, a_l, b_l) +
					alpha*log_inv_gamma_kernel(a1, a_1, b_1) +
					alpha*log_inv_gamma_kernel(a2, a_2, b_2) -
					alpha*0.5 * mahalanobis(b, mu_b, Sigma)^2
		cllValues[iter] <- cll$cll
		lpd[iter] <- log_prior_density
		if(iter %% 10000 == 0 & verbose){
			cat(paste0("* iteration = ", iter, ".  Current ergodic means (after discarding 30% of iterations):"), "\n")			
			burn <- floor(0.3*iter)
			ergMean <- c(mean(g_mcmc[(burn+1):iter]), mean(lambda_mcmc[(burn+1):iter]), 
				mean(a1_mcmc[(burn+1):iter]), mean(a2_mcmc[(burn+1):iter]), 
				colMeans(as.matrix(b_mcmc[(burn+1):iter,])))
				names(ergMean) <- c("gamma", "lambda", "a1", "a2", paste0('b',bRange))
			print(round(ergMean,2))
			cat(paste0("    g_accept_rate = ", round(100*g_accept_rate/mh_iter, 2), "%.   "), "\n")
                        cat(paste0("    lambda_accept_rate = ", round(100*lambda_accept_rate/mh_iter, 2), "%.   "), "\n")
                        cat(paste0("    a1_accept_rate = ", round(100*a1_accept_rate/mh_iter, 2), "%.   "), "\n")                                          
                        cat(paste0("    a2_accept_rate = ", round(100*a2_accept_rate/mh_iter, 2), "%.   "), "\n")
                        cat(paste0("    b_accept_rate = ", round(100*b_accept_rate/mh_iter, 2), "%.   "), "\n")                                                                                               
			cat(paste0("    mala_accept_rate = ", round(100*mala_accept_rate/mala_iter, 2), "%.   "), "\n")
			cat(paste0("\n"))
			if(plot){
				par(mfrow = c(2,3))		
				plot(g_mcmc[1:iter], type = "l", xlab = "MCMC iteration", ylab = bquote(gamma))
				plot(lambda_mcmc[1:iter], type = "l", xlab = "MCMC iteration", ylab = bquote(lambda))
				plot(a1_mcmc[1:iter], type = "l", xlab = "MCMC iteration", ylab = bquote(alpha[1]))						
				plot(a2_mcmc[1:iter], type = "l", xlab = "MCMC iteration", ylab = bquote(alpha[2]))
				matplot(as.matrix(b_mcmc[1:iter,]), type = "l", xlab = "MCMC iteration", ylab = 'regression coefficients')														
				lText <- paste0('b',bRange)
				legend("topright", lText, col = 1:nCov, lty = 1:nCov)
				plot(cllValues[1:iter], type = "l", xlab = "MCMC iteration", ylab = "complete log-likelihood")
			}
		}

	}
		result <- vector("list", length = 5)
		colnames(b_mcmc) <- paste0('b',bRange, '_mcmc')
		mcmc_sample <- as.mcmc(cbind(g_mcmc, lambda_mcmc, a1_mcmc, a2_mcmc, b_mcmc))
		latent_status_D0 <- I_sim_values_D0
		colnames(latent_status_D0) <- D0
		complete_log_likelihood <- cllValues
		result[[1]] <- mcmc_sample
		result[[2]] <- complete_log_likelihood
		result[[3]] <- c(g_accept_rate/mh_iter, lambda_accept_rate/mh_iter, 
				a1_accept_rate/mh_iter, a2_accept_rate/mh_iter, b_accept_rate/mh_iter, mala_accept_rate/mala_iter)
		result[[4]] <- latent_status_D0
		result[[5]] <- lpd
		names(result[[3]]) <- c("gamma", "lambda", "a1", "a2", "betas","mala")
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






cure_rate_MC3 <- function( y, X, Censoring_status, nChains = 16, 
				mcmc_cycles = 15000, 
				alpha = NULL, 
				nCores = 8, 
				sweep = 10,
				mu_g = 0.2, s2_g = 0.1, a_l = 2.001, b_l = 1, 
				a_1 = 2.001, b_1 = 1, a_2 = 2.001, b_2 = 1, 
				mu_b = rep(0,dim(X)[2]), Sigma = 100*diag(dim(X)[2]),
				g_prop_sd = 0.045, 
				lambda_prop_scale = 0.03, 
				a1_prop_scale = 0.2, 
				a2_prop_scale = 0.03, 
				b_prop_sd = rep(0.022, dim(X)[2]), 
				initialValues = NULL, 
				plot = TRUE,
				adjust_scales = TRUE,
				adjust_alpha = FALSE,
				verbose = TRUE, tau_mala = 0.000015, mala = 0.33
					){
	X <- as.matrix(X)
	if(min(y) < 0){stop("Negative values are not allowed in the response variable.")}				
	if(nChains < 2){stop("at least two chains should be considered.")}
	if(is.null(alpha)){
		alpha <- 1/1.001^{(1:nChains)^2.5 - 1}	# temperature per chain (the target chain always corresponds to chain 1)
		cat("Default (inverse) temperatures: ", "\n")
		print(alpha)
	}else{
		if(length(alpha)!= nChains){
			stop("number of temperatures does not match number of chains.")
		}
		if(alpha[1] != 1){
			stop("alpha[1] should be equal to 1.")
		}
	}
	nCov <- dim(X)[2]
	nPars <- 4 + nCov
	D0 <- which(Censoring_status == 0)
	nCensored <- length(D0)
	parallel_mcmc_samples <- array(data = NA, dim = c(1, nPars + nCensored, nChains))
	target_mcmc <- array(data = NA, dim = c(mcmc_cycles, nPars + nCensored))

	if(length(unique(X[,1])) == 1){
		b_names <- paste0('b',1:nCov - 1,'_mcmc')
	}else{
		b_names <- paste0('b',1:nCov,'_mcmc')	
	}

	dimnames(parallel_mcmc_samples) = list(paste0("current_value"),
		c("g_mcmc", "lambda_mcmc", "a1_mcmc", "a2_mcmc", b_names,paste0("status_",D0)),
						paste0("chain_",1:nChains))
		dimnames(parallel_mcmc_samples)[[3]][1]	<- "chain_1_(target_chain)"
		dimnames(target_mcmc) = list(paste0("cycle_",1:mcmc_cycles),
			c("g_mcmc", "lambda_mcmc", "a1_mcmc", "a2_mcmc", b_names,paste0("status_",D0)))

	cllValues <- numeric(mcmc_cycles)


	g_prop_sd <- rep(g_prop_sd, nChains)
	lambda_prop_scale <- rep(lambda_prop_scale, nChains)
	a1_prop_scale <- rep(a1_prop_scale, nChains)
	a2_prop_scale <- rep(a2_prop_scale, nChains)
	b_prop_sd <- array(data = b_prop_sd, dim = c(nChains, nCov))

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

		  run <- cure_rate_metropolis_hastings_cpp(  y, X, Censoring_status,  
		                        m = 1000, 
				mu_g = mu_g, s2_g = s2_g, a_l = a_l, b_l = b_l, 
				a_1 = a_1, b_1 = b_1, a_2 = a_2, b_2 = b_2, 
				mu_b = mu_b, Sigma = Sigma,
		                        alpha = alpha[chain_iter],
		                        g_prop_sd = g_prop_sd[chain_iter], 
		                        lambda_prop_scale = lambda_prop_scale[chain_iter], 
		                        a1_prop_scale = a1_prop_scale[chain_iter], 
		                        a2_prop_scale = a2_prop_scale[chain_iter], 
		                        b_prop_sd = b_prop_sd[chain_iter,], 
		                        plot = FALSE, verbose = FALSE, tau_mala = tau_mala, mala = mala)

			if(mala>0){
				indUp <- which((run$acceptance_rates > 0.6) == TRUE)
				indDown <- which((run$acceptance_rates < 0.4) == TRUE)
			}else{
				indUp <- which((run$acceptance_rates[-6] > 0.3) == TRUE)
				indDown <- which((run$acceptance_rates[-6] < 0.15) == TRUE)			
			}
		        if(length(indUp) > 0){
		                criterion = TRUE; 
		                if(1 %in% indUp){g_prop_sd[chain_iter] <- g_prop_sd[chain_iter]/adjFactor}
		                if(2 %in% indUp){lambda_prop_scale[chain_iter] <- lambda_prop_scale[chain_iter]/adjFactor}      
		                if(3 %in% indUp){a1_prop_scale[chain_iter] <- a1_prop_scale[chain_iter]/adjFactor}              
		                if(4 %in% indUp){a2_prop_scale[chain_iter] <- a2_prop_scale[chain_iter]/adjFactor}                      
		                if(5 %in% indUp){b_prop_sd[chain_iter,] <- b_prop_sd[chain_iter,]/adjFactor}
		                if(mala>0){
      		                	if(6 %in% indUp){tau_mala <- tau_mala/adjFactor}
      		                }
		                criterion1 = TRUE                               
		                
		        }else{criterion1 = FALSE}
		        if(length(indDown) > 0){
		                criterion = TRUE
		                if(1 %in% indDown){g_prop_sd[chain_iter] <- g_prop_sd[chain_iter] * adjFactor}
		                if(2 %in% indDown){lambda_prop_scale[chain_iter] <- lambda_prop_scale [chain_iter]* adjFactor}    
		                if(3 %in% indDown){a1_prop_scale[chain_iter] <- a1_prop_scale[chain_iter] * adjFactor}            
		                if(4 %in% indDown){a2_prop_scale[chain_iter] <- a2_prop_scale[chain_iter] * adjFactor}                    
		                if(5 %in% indDown){b_prop_sd[chain_iter,] <- b_prop_sd[chain_iter,] * adjFactor}
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
		print(c(g_prop_sd[chain_iter], lambda_prop_scale[chain_iter], a1_prop_scale[chain_iter], a2_prop_scale[chain_iter], b_prop_sd[chain_iter,], tau_mala), digits = 4)

	for(chain_iter in 2:nChains){
	# initial metropolis-hastings proposal parameters (they will be adjusted)
	      g_prop_sd[chain_iter] = g_prop_sd[chain_iter - 1]
		lambda_prop_scale[chain_iter] = lambda_prop_scale[chain_iter-1]
		a1_prop_scale[chain_iter] = a1_prop_scale[chain_iter-1]
		a2_prop_scale[chain_iter] = a2_prop_scale[chain_iter-1]
		b_prop_sd[chain_iter,] = b_prop_sd[chain_iter-1,]

	      criterion <- TRUE
		adjFactor <- 0.75
		iter <- 0
	 
	      while(criterion & (iter < 5)){

		  run <- cure_rate_metropolis_hastings_cpp(  y, X, Censoring_status,  
		                        m = 500, 
				mu_g = mu_g, s2_g = s2_g, a_l = a_l, b_l = b_l, 
				a_1 = a_1, b_1 = b_1, a_2 = a_2, b_2 = b_2, 
				mu_b = mu_b, Sigma = Sigma,
				alpha = alpha[chain_iter],
		                        g_prop_sd = g_prop_sd[chain_iter], 
		                        lambda_prop_scale = lambda_prop_scale[chain_iter], 
		                        a1_prop_scale = a1_prop_scale[chain_iter], 
		                        a2_prop_scale = a2_prop_scale[chain_iter], 
		                        b_prop_sd = b_prop_sd[chain_iter,], 
		                        plot = FALSE, verbose = FALSE, tau_mala = tau_mala, mala = mala)

			if(mala>0){
				indUp <- which((run$acceptance_rates > 0.6) == TRUE)
				indDown <- which((run$acceptance_rates < 0.4) == TRUE)
			}else{
				indUp <- which((run$acceptance_rates[-6] > 0.3) == TRUE)
				indDown <- which((run$acceptance_rates[-6] < 0.15) == TRUE)			
			}

		        if(length(indUp) > 0){
		                criterion = TRUE; 
		                if(1 %in% indUp){g_prop_sd[chain_iter] <- g_prop_sd[chain_iter]/adjFactor}
		                if(2 %in% indUp){lambda_prop_scale[chain_iter] <- lambda_prop_scale[chain_iter]/adjFactor}      
		                if(3 %in% indUp){a1_prop_scale[chain_iter] <- a1_prop_scale[chain_iter]/adjFactor}              
		                if(4 %in% indUp){a2_prop_scale[chain_iter] <- a2_prop_scale[chain_iter]/adjFactor}                      
		                if(5 %in% indUp){b_prop_sd[chain_iter,] <- b_prop_sd[chain_iter,]/adjFactor}
		                if(mala>0){
      		                	if(6 %in% indUp){tau_mala <- tau_mala/adjFactor}
      		                }
		                
		                criterion1 = TRUE                               
		                
		        }else{criterion1 = FALSE}
		        if(length(indDown) > 0){
		                criterion = TRUE
		                if(1 %in% indDown){g_prop_sd[chain_iter] <- g_prop_sd[chain_iter] * adjFactor}
		                if(2 %in% indDown){lambda_prop_scale[chain_iter] <- lambda_prop_scale [chain_iter]* adjFactor}    
		                if(3 %in% indDown){a1_prop_scale[chain_iter] <- a1_prop_scale[chain_iter] * adjFactor}            
		                if(4 %in% indDown){a2_prop_scale[chain_iter] <- a2_prop_scale[chain_iter] * adjFactor}                    
		                if(5 %in% indDown){b_prop_sd[chain_iter,] <- b_prop_sd[chain_iter,] * adjFactor}                            
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
		print(c(g_prop_sd[chain_iter], lambda_prop_scale[chain_iter], a1_prop_scale[chain_iter], a2_prop_scale[chain_iter], b_prop_sd[chain_iter,], tau_mala), digits = 4)
	}
	cat("   done Metropolis-Hastings proposal scale adjustements!","\n")
	}
	if(adjust_alpha){
		alpha0 <- 1.05
		alpha <- 1/alpha0^{1:nChains - 1}

		cat(paste0("   Initial temperatures: "),"\n")
		print(alpha)

		mcmc_cyc <- 100
		if(mcmc_cyc > 4){mcmc_cyc <- 2*mcmc_cyc}
		if(mcmc_cyc > 8){mcmc_cyc <- 2*mcmc_cyc}		
		
	
	
		cat(paste0("Adjusting temperatures..."),"\n")
		run <- cure_rate_MC3_cpp(  y, X, Censoring_status, nChains = nChains, 
				mcmc_cycles = mcmc_cyc, 
				alpha = alpha, 
				nCores = nCores, 
				sweep = sweep,
				mu_g = mu_g, s2_g = s2_g, a_l = a_l, b_l = b_l, 
				a_1 = a_1, b_1 = b_1, a_2 = a_2, b_2 = b_2, 
				mu_b = mu_b, Sigma = Sigma,
		                lambda_prop_scale = lambda_prop_scale[chain_iter], 
		                a1_prop_scale = a1_prop_scale[chain_iter], 
		                a2_prop_scale = a2_prop_scale[chain_iter], 
		                b_prop_sd = b_prop_sd[chain_iter,],
		                initialValues = NULL, 
				plot = FALSE,
				adjust_scales = FALSE,
				adjust_alpha = FALSE,
				verbose = FALSE, tau_mala = tau_mala, mala = mala
			)
		
		criterion <- criterion_up <- criterion_down <- FALSE
		if ( min(run[["swap_accept_rate"]]) < 0.35 ){
			criterion <- TRUE
			criterion_down <- TRUE
		}
		if ( min(run[["swap_accept_rate"]]) > 0.75 ){
			criterion <- TRUE
			criterion_up <- TRUE
		}
		cat(paste0("Swap accept rates:"), "\n")
		print(run[["swap_accept_rate"]])
		iter <- 0
		while( criterion & (iter < 30)){
		
			if(criterion_down){
				alpha0 <- max(c(alpha0 - 0.01, 1.001))
			}
			if(criterion_up){
				alpha0 <- alpha0 + 0.01
			}
			alpha <- 1/alpha0^{1:nChains - 1}
			cat(paste0("   Adjusted temperatures: "),"\n")
			print(alpha)

			run <- cure_rate_MC3_cpp(  y, X, Censoring_status, nChains = nChains, 
				mcmc_cycles = mcmc_cyc, 
				alpha = alpha, 
				nCores = nCores, 
				sweep = sweep,
				mu_g = mu_g, s2_g = s2_g, a_l = a_l, b_l = b_l, 
				a_1 = a_1, b_1 = b_1, a_2 = a_2, b_2 = b_2, 
				mu_b = mu_b, Sigma = Sigma,
		                lambda_prop_scale = lambda_prop_scale[chain_iter], 
		                a1_prop_scale = a1_prop_scale[chain_iter], 
		                a2_prop_scale = a2_prop_scale[chain_iter], 
		                b_prop_sd = b_prop_sd[chain_iter,],
		                initialValues = NULL, 
				plot = FALSE,
				adjust_scales = FALSE,
				adjust_alpha = FALSE,
				verbose = FALSE, tau_mala = tau_mala, mala = mala
			)
			print(run[["swap_accept_rate"]])
			criterion <- criterion_up <- criterion_down <- FALSE
			if ( min(run[["swap_accept_rate"]]) < 0.25 ){
				criterion <- TRUE
				criterion_down <- TRUE
			}
			if ( min(run[["swap_accept_rate"]]) > 0.65 ){
				criterion <- TRUE
				criterion_up <- TRUE
			}
			iter <- iter + 1
			
		}
		cat("   done temperature adjustements!","\n")
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
		parLoop <- foreach(chain_iter = 1:nChains, .export = ls(envir = globalenv())) %dopar% {
#			if(chain_iter == 1){perform_mala = mala}else{
#				if(alpha[chain_iter] == 1){
#					perform_mala = mala
#				}else{perform_mala = 0}
#			}
			mcmc_chain <- cure_rate_metropolis_hastings_cpp(  y, X, Censoring_status,  
				                m = sweep, 
						mu_g = mu_g, s2_g = s2_g, a_l = a_l, b_l = b_l, 
						a_1 = a_1, b_1 = b_1, a_2 = a_2, b_2 = b_2, 
						mu_b = mu_b, Sigma = Sigma,
				                alpha = alpha[chain_iter],
				                g_prop_sd = g_prop_sd[chain_iter], 
				                lambda_prop_scale = lambda_prop_scale[chain_iter], 
				                a1_prop_scale = a1_prop_scale[chain_iter], 
				                a2_prop_scale = a2_prop_scale[chain_iter], 
				                b_prop_sd = b_prop_sd[chain_iter,], plot = FALSE, 
				                tau_mala = tau_mala, mala = mala)
						
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

			log_posterior1 <- parLoop[[j1]]$complete_log_likelihood[sweep] #+ parLoop[[j1]]$log_prior_density[sweep]
			log_posterior2 <- parLoop[[j2]]$complete_log_likelihood[sweep] #+ parLoop[[j2]]$log_prior_density[sweep]
			denom <- log_posterior1 + log_posterior2
#			print(denom)
			state_j1 <- as.list(parLoop[[j1]]$mcmc_sample[sweep,])
			state_j1[["I_sim_D0"]] <- parLoop[[j1]]$latent_status_censored[sweep,]

			state_j2 <- as.list(parLoop[[j2]]$mcmc_sample[sweep,])
			state_j2[["I_sim_D0"]] <- parLoop[[j2]]$latent_status_censored[sweep,]
			nom1 <- cure_rate_metropolis_hastings_cpp(  y, X, Censoring_status,  
								m = 1, # only the log-posterior is computed at the initial values
								mu_g = mu_g, s2_g = s2_g, a_l = a_l, b_l = b_l, 
						a_1 = a_1, b_1 = b_1, a_2 = a_2, b_2 = b_2, 
						mu_b = mu_b, Sigma = Sigma,
								alpha = alpha[j1], plot = FALSE, 
								initialValues = state_j2, tau_mala = tau_mala, mala = mala)
			nom2 <- cure_rate_metropolis_hastings_cpp(  y, X, Censoring_status,  
								m = 1, # only the log-posterior is computed at the initial values
								alpha = alpha[j2], plot = FALSE, 
								mu_g = mu_g, s2_g = s2_g, a_l = a_l, b_l = b_l, 
						a_1 = a_1, b_1 = b_1, a_2 = a_2, b_2 = b_2, 
						mu_b = mu_b, Sigma = Sigma,
								initialValues = state_j1, tau_mala = tau_mala, mala = mala)
			nom <- nom1$complete_log_likelihood + #nom1$log_prior_density 
				+ nom2$complete_log_likelihood 
				#+ nom2$log_prior_density	
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
			parLoop <- foreach(chain_iter = 1:nChains, .export = ls(envir = globalenv())) %dopar% {
				initialValues <-  as.list(parallel_mcmc_samples[1, 1:nPars,chain_iter])
				initialValues[["I_sim_D0"]] <- parallel_mcmc_samples[1, -(1:nPars),chain_iter]
#				if(chain_iter == 1){perform_mala = mala}else{
#					if(alpha[chain_iter] == 1){
#						perform_mala = mala
#					}else{perform_mala = 0}
#				}
				
				mcmc_chain <- cure_rate_metropolis_hastings_cpp(  y, X, Censoring_status,  
						        m = sweep, 
						        alpha = alpha[chain_iter],
							mu_g = mu_g, s2_g = s2_g, a_l = a_l, b_l = b_l, 
							a_1 = a_1, b_1 = b_1, a_2 = a_2, b_2 = b_2, 
							mu_b = mu_b, Sigma = Sigma,
						        g_prop_sd = g_prop_sd[chain_iter], 
						        lambda_prop_scale = lambda_prop_scale[chain_iter], 
						        a1_prop_scale = a1_prop_scale[chain_iter], 
						        a2_prop_scale = a2_prop_scale[chain_iter], 
						        b_prop_sd = b_prop_sd[chain_iter,], plot = FALSE,
						        initialValues = initialValues, 
						        tau_mala = tau_mala, mala = mala)
							
			}
			stopImplicitCluster()
 			for(chain_iter in 1:nChains){
				parallel_mcmc_samples[1, ,chain_iter] <- c(parLoop[[chain_iter]]$mcmc_sample[sweep,], 
											parLoop[[chain_iter]]$latent_status_censored[sweep,])
				target_mcmc[cycle,] <- 	parallel_mcmc_samples[1, ,1]
				all_cll_values[chain_iter,cycle] <- parLoop[[chain_iter]]$complete_log_likelihood[sweep]														
			}
			cllValues[cycle] <- parLoop[[1]]$complete_log_likelihood[sweep]
 			if(cycle %% 1000 == 0 && verbose == TRUE){
 				cat(paste0("** MCMC cycle: ", cycle, " completed."),"\n")
				cat(paste0("   Chain-swap rate (overall): ", round(100*swap_accept/cycle,2),"%."),"\n")
				cat(paste0("                 (per chain): ", paste0(round(100*swap_accept_per_chain/(n_attempts_per_chain),2),"%",collapse = ", ")),"\n")
				cat("   Current ergodic means and median (after discarding 30% of iterations):","\n")
		                burn <- floor(0.3*cycle)
		                ergMean <- colMeans(target_mcmc[(burn+1):cycle,1:nPars])
		                        names(ergMean) <- c("gamma", "lambda", "a1", "a2", b_names)
		                ergMed <- apply(target_mcmc[(burn+1):cycle,1:nPars], 2, median)
		                        names(ergMed) <- c("gamma", "lambda", "a1", "a2", b_names)
		                tmpMat <- rbind(ergMean, ergMed)
		                rownames(tmpMat) <- c("mean", "median")
		                print(round(tmpMat,2))

					if(plot){
					par(mfrow = c(2,3))		
					plot(target_mcmc[1:cycle, 1], 
						type = "l", xlab = "MCMC iteration", ylab = bquote(gamma))
					plot(target_mcmc[1:cycle, 2], 
						type = "l", xlab = "MCMC iteration", ylab = bquote(lambda))
					plot(target_mcmc[1:cycle, 3], 
						type = "l", xlab = "MCMC iteration", ylab = bquote(alpha[1]),
						ylim = quantile(target_mcmc[1:cycle, 3],c(0.001,0.99)))						
					plot(target_mcmc[1:cycle, 4], 
						type = "l", xlab = "MCMC iteration", ylab = bquote(alpha[2]))
					matplot(as.matrix(target_mcmc[1:cycle, 5:nPars]), col = 1:nCov, lty = 1:nCov,
						type = "l", xlab = "MCMC iteration", ylab = 'regression coefficients')														
					lText <- b_names
					legend("topright", lText, col = 1:nCov, lty = 1:nCov)
					#plot(cllValues[1:cycle], type = "l", xlab = "MCMC iteration", ylab = "complete log-likelihood",
					#ylim = quantile(cllValues[1:cycle],c(0.01,0.999)))
					matplot(t(all_cll_values[,501:cycle]), type='l', col  = heat.colors(nChains) )
						points(all_cll_values[1,501:cycle], type = 'b', 
							col = heat.colors(nChains)[1])					
				}
				
			}


		}
		result <- vector("list", length = 5)
		result[[1]] <- as.mcmc(target_mcmc[,1:nPars])
		result[[2]] <- target_mcmc[,-(1:nPars)]
		result[[3]] <- cllValues
		result[[4]] <- swap_accept_per_chain/(n_attempts_per_chain + 0.001)
		result[[5]] <- all_cll_values
	 	names(result) <- c("mcmc_sample", "latent_status_censored", "complete_log_likelihood","swap_accept_rate",'all_cll_values')
	 	return(result)
}


compute_map_and_hdis <- function(y, X, Censoring_status, retained_mcmc, prior_parameters = NULL, gamma_mix = TRUE, K_gamma = 2){

	if(is.null(prior_parameters)){
		mu_g = 1
		s2_g = 1
		a_l = 2.001
		b_l = 1 
		a_1 = 2.001
		b_1 = 1
		a_2 = 2.001
		b_2 = 1
		mu_b = rep(0,dim(X)[2])
		Sigma = 100*diag(dim(X)[2])
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
	
	ct = exp(exp(-1))
	X <- as.matrix(X)
	x <- X
	nCov <- dim(x)[2]
	log_S_p <- function(tau, g, lambda, a1, a2, b, x){
		theta <- exp(x %*% b)
		return(-log(1 + g * theta * ct^{g*theta} * (1 - exp(-(a1*tau)^a2))^lambda)/g)
	}

	c_under <- 10^{-9}
	log_f_p <- function(tau, g, lambda, a1, a2, b, logS){
		# logS = log_S_p(tau = tau, g = g, lambda = lambda, a1 = a1, a2 = a2, b0 = b0, b1 = b1, b2 = b2)
		one_minus_exp <- apply(cbind((1 - exp(-(a1*tau)^a2)), c_under), 1, max)
		log_weibull_dens <- log(a1) + log(a2)  -(a1*tau)^a2 + (a2 - 1)*log(a1*tau)
		log_theta <- x %*% b
		return(
		        (1 + g) * logS + log(lambda) + log_theta +
		        g*exp(log_theta)*log(ct) + 
		        (lambda - 1)*log(one_minus_exp) + #### NOTE: this causes the log -> -inf, so now it regulated.
		        log_weibull_dens
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
	
	for(iter in ind){


		logS <- log_S_p(tau = y, 
			g = retained_mcmc[iter,'g_mcmc'], 
			lambda = retained_mcmc[iter,'lambda_mcmc'], 
			a1 = retained_mcmc[iter,'a1_mcmc'], 
			a2 = retained_mcmc[iter,'a2_mcmc'], 
			b = retained_mcmc[iter, bIndices], 
			x = x)

		logf <- log_f_p(tau = y, 
			g = retained_mcmc[iter,'g_mcmc'], 
			lambda = retained_mcmc[iter,'lambda_mcmc'], 
			a1 = retained_mcmc[iter,'a1_mcmc'], 
			a2 = retained_mcmc[iter,'a2_mcmc'], 
			b = retained_mcmc[iter, bIndices], 
			logS = logS
			)
		tz <- tz + 1
		logL[tz] <- sum(Censoring_status * logf) + sum((1-Censoring_status)*logS)
		g = retained_mcmc[iter,'g_mcmc']
		lambda = retained_mcmc[iter,'lambda_mcmc']
		a1 = retained_mcmc[iter,'a1_mcmc']
		a2 = retained_mcmc[iter,'a2_mcmc'] 
		b = retained_mcmc[iter, bIndices]
		log_prior_density <- log_prior_gamma(g, mu_g, s2_g) + 
			        log_inv_gamma_kernel(lambda, a_l, b_l) +
			        log_inv_gamma_kernel(a1, a_1, b_1) +
			        log_inv_gamma_kernel(a2, a_2, b_2) -
			        0.5 * mahalanobis(b, mu_b, Sigma)^2

		logP[tz] <- logL[tz] + log_prior_density
		if(iter %% 100 == 0){
		cat(paste0('iteration: ', iter),'\r')
		}
		
	}
	cat('\n')
	n <- dim(X)[1]
	n_parameters <- dim(x)[2] + 4
	BIC <- -2 * max(logL) + n_parameters * log(n)
	map_estimate <- retained_mcmc[which.max(logP), ]
	hdis <- fpf <- plot.bayesCureModel(retained_mcmc = retained_mcmc, alpha = 0.05, plot = FALSE)
	results <- vector('list', length = 4)
	results[[1]] <- logP
	results[[2]] <- BIC
	results[[3]] <- map_estimate
	results[[4]] <- hdis
	names(results) <- c('log_posterior', 'bic','map_estimate', 'highest_density_indervals')
	return(results)

}

plot.bayesCureModel <- function(retained_mcmc, map_estimate = NULL, alpha = 0.05, gamma_mix = TRUE, K_gamma = 2, plot = TRUE, bw = 'ucv'){
	nPars <- dim(retained_mcmc)[2]
	myXlim <- matrix(NA, nPars, 2)
	for(i in 1:nPars){
		myXlim[i,] <- quantile(retained_mcmc[,i],probs = c(0.001,0.999))# c(-5,2)
	}
	hdis <- vector('list', length = nPars)
	m <- dim(retained_mcmc)[1]
	ind <- 1:m
	
	varnames <- numeric(nPars)
	varnames[1:4] <- as.expression(c(
		bquote(gamma), 
		bquote(lambda), 
		bquote(alpha[1]), 
		bquote(alpha[2])))
	b_ind <- as.numeric(unlist(lapply(strsplit(unlist(lapply(strsplit(colnames(retained_mcmc)[-(1:4)], 
		split = '_'), function(x)x[1])), split = 'b'), function(x)x[2])))
	for(i in 5:nPars){
		varnames[i] <- as.expression(bquote(beta[.(b_ind[i-4])]))
	}

	hdi_alpha = alpha
	for(i in 1:nPars){
#		pdf(file = paste0("../img/recidivism_new_data_parameter_",i,".pdf"), width = 12, height = 3)
		if(i == 1){
			shouldIask = FALSE
		}else{shouldIask = TRUE}
		if(plot){
		par(mar = c(4,4,0.1,0.5), ask = shouldIask)	
		}
		x <- retained_mcmc[ind,i]

		myD <- density(x, bw = bw)
#		if(i < 5){
#			myD <- density(x, bw = 'bcv')			
#		}
		if(i == 1){
		if(gamma_mix){
		fit <- Mclust(x,G=1:K_gamma, modelNames = "V")
		k <- fit$G
		if( k > 1){
			multMode = TRUE
			mu <- fit$parameters$mean
			w <- fit$parameters$pro
			s2 <- fit$parameters$variance$sigmasq
			dd <- range(x) + 0.01*c(-1,1)
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
		allow_split = TRUE
		if(i %in% c(2,3,4)){allow_split = FALSE}
		hdi_95 <- hdi(myD,allowSplit=allow_split, credMass = 1 - hdi_alpha)				


		if(is.null(dim(hdi_95))){hdi_length = diff(hdi_95)}else{
		hdi_length = sum(apply(hdi_95,1,diff))}
		hdis[[i]] <- hdi_95
		if(is.null(dim(hdi_95))){hdi_95 = matrix(hdi_95,1,2)}
		if(plot){
		plot(myD, xlab = varnames[i], main = '', xlim = myXlim[i,])
		}
		if(i == 1){
					if(plot){
			legend('topleft', c('estimate (map)', paste0(100*(1-hdi_alpha), '% HDI')), 
				col = c('red','papayawhip'),lty = 1,lwd = c(1,10))
				}
		}
		#hist(mcmcmc16$mcmc_sample[ind,i], xlab = varnames[i], main = '');abline(v = truePars[i], col = 1, lwd = 2)

		for(j in 1:dim(hdi_95)[1]){
			#abline(v = hdi_95[j,], lty = 2, col = 1+j)
			x_dens1 <- which.min(abs(hdi_95[j,1] - myD$x))
			x_dens2 <- which.min(abs(hdi_95[j,2] - myD$x))
			y_dens <- myD$y[c(x_dens1,x_dens2)]
					if(plot){
			polygon(c(myD$x[c(x_dens1:x_dens2,x_dens2:x_dens1)]),
				c(rep(0,x_dens2 - x_dens1+1),myD$y[c(x_dens2:x_dens1)]),col='papayawhip')
				}
			#points(c(hdi_95[j,1],hdi_95[j,1]),c(0,y_dens[1]), type = 'l')	
			#points(c(hdi_95[j,2],hdi_95[j,2]),c(0,y_dens[2]), type = 'l')			
		}
		hdi_length = sum(apply(hdi_95,1,diff))
		if(is.null(map_estimate) == FALSE){
				if(plot){
			abline(v = c(map_estimate[i]), col = 2, lwd = 2, lty = 1)
				}
		}
#		abline(v = truePars[i], col = 'green', lwd = 2, lty = 2)
#		dev.off()
		
	}
			if(plot){
	par(ask=FALSE)
	}
	names(hdis) <- c('g', 'lambda', 'a1', 'a2',paste0('b',b_ind))
	return(hdis)

}





