library('calculus')


f <- function(g, l, a1, a2, b0, b1, b2)complete_log_likelihood(y = y, X = X1, Censoring_status=stat, g, l,a1, a2 , c(b0, b1, b2), I_sim = sim[,5], alpha = 1)$cll

derivative(f, var = c(g = 1, l = 1, a1 = 1, a2 = 1,b0 = 0, b1 = 1, b2 = 1))

derv_complete_log_likelihood_fotis(y = y, X = X1, Censoring_status = stat, g = 1, lambda = 1, a1 = 1, a2 = 1, b = c(0,1,1), I_sim = sim[,5])



# with the generic functions


F_parameters <- list(
	K = 2,
	p = c(0.5,0.5),
	gamma_parameters = matrix(c(1,1,6,1), 2, 2, byrow =TRUE)
)

#
sim <- sim_fotis_model(n = 1000, lambda1 = 6, F_parameters = F_parameters, ab = 0.05, ranunL=0, ranunU=1, seed = 1, discrcov = 1)


y <- sim[,1]
X1 <- cbind(1, sim[,3:4])
stat <- sim[,2]

I_sim = sim[,5]
alpha = 1
f_gamma <- function(g, lambda, a1, a2, b0, b1, b2){
	lw <- log_gamma(y = y, a1 = a1, a2 = a2, c_under = 1e-9)
	b = c(b0, b1, b2)
	return(complete_log_likelihood_general(y = y, X = X1, Censoring_status = stat, g = g, lambda = lambda, 
		log_f = lw$log_f, log_F = lw$log_F, b = b, I_sim = I_sim, alpha = alpha)$cll
		)
	}
derivative(f_gamma, var = c(g = 1, lambda = 1, a1 = 1, a2 = 1,b0 = 0, b1 = 1, b2 = 1))	
	
f_exponential <- function(g, lambda, a1, b0, b1, b2){
	lw <- log_weibull(y, a1 = a1, a2 = 1,  c_under = 1e-9)
	b = c(b0, b1, b2)
	return(complete_log_likelihood_general(y = y, X = X1, Censoring_status = stat, g = g, lambda = lambda, 
		log_f = lw$log_f, log_F = lw$log_F, b = b, I_sim = I_sim, alpha = alpha)$cll
		)
	}
derivative(f_exponential, var = c(g = 1, lambda = 1, a1 = 1, b0 = 0, b1 = 1, b2 = 1))	

Κ <- 2
w <- a[-(1:(2*K))]
w2 <- c(w, 1)
p <- w2/sum(w2)

a1_ind <- seq(1, 2*K, by = 2)
a2_ind <- seq(2, 2*K, by = 2)		
a1 = a[a1_ind]
a2 = a[a2_ind]		

f_gamma_mixture <- function(g, lambda, a11, a12, a21, a22, p1, b0, b1, b2){
#			a11, a12,..., a1K denotes gamma-shape parameter per component (1, 2, ... ,K)
#			a21, a22,..., a2K denotes gamma-rate parameter per component (1, 2, ... ,K)
	a1 = c(a11, a12)
	a2 = c(a21, a22)
	p = c(p1, 1 - p1)	
	lw <- log_gamma_mixture(y, a1 = a1, a2 = a2, p = p, c_under = 1e-9)
	b = c(b0, b1, b2)
	return(complete_log_likelihood_general(y = y, X = X1, Censoring_status = stat, g = g, lambda = lambda, 
		log_f = lw$log_f, log_F = lw$log_F, b = b, I_sim = I_sim, alpha = alpha)$cll
		)
	}
derivative(f_gamma_mixture, var = c(g = 1, 
				lambda = 6, 
				a11 = F_parameters[['gamma_parameters']][1,1], # shape component 1
				a12 = F_parameters[['gamma_parameters']][2,1], # shape component 2
				a21 = F_parameters[['gamma_parameters']][1,2], # rate component 1
				a22 = F_parameters[['gamma_parameters']][2,2], # rate component 2
				p1=0.5, 
				b0 = 1.5, b1 = 1.5, b2 = -0.8))    



f_gamma_mixture2 <- function(g, lambda, a1, a2, p1, b0, b1, b2){
#			a11, a12,..., a1K denotes gamma-shape parameter per component (1, 2, ... ,K)
#			a21, a22,..., a2K denotes gamma-rate parameter per component (1, 2, ... ,K)
	p = c(p1, 1-p1)
	lw <- log_gamma_mixture(y, a1 = a1, a2 = a2, p = p, c_under = 1e-9)
	b = c(b0, b1, b2)
	return(complete_log_likelihood_general(y = y, X = X1, Censoring_status = stat, g = g, lambda = lambda, 
		log_f = lw$log_f, log_F = lw$log_F, b = b, I_sim = I_sim, alpha = alpha)$cll
		)
	}



f <- function(par_vec, family, K, I_sim){
#			par_vec should be gamma, lambda, a1, a2, p, b
#			a11, a12,..., a1K denotes gamma-shape parameter per component (1, 2, ... ,K)
#			a21, a22,..., a2K denotes gamma-rate parameter per component (1, 2, ... ,K)
	g = par_vec[1]
	lambda = par_vec[2]	
	if(family == 'weibull'){
		lw <- log_weibull(y, a1 = par_vec[3], a2 = par_vec[4],  c_under = 1e-9)	
		b = par_vec[-(1:4)]
	}

	if(family == 'exponential'){
		lw <- log_weibull(y, a1 = par_vec[3], a2 = 1,  c_under = 1e-9)	
		b = par_vec[-(1:3)]
	}


	if(family == 'gamma'){
		lw <- log_gamma(y = y, a1 = par_vec[3], a2 = par_vec[4], c_under = 1e-9)	
		b = par_vec[-(1:4)]
	}
	
	if(family == 'gamma_mixture'){
		p = c(par_vec[(2+K*2 + 1):(2+K*2 + K-1)], 1-sum(par_vec[(2+K*2 + 1):(2+K*2 + K-1)]))
		lw <- log_gamma_mixture(y, a1 = par_vec[3:(2+K)], a2 = par_vec[(K+3):(2*K+2)], p = p, c_under = 1e-9)
		b = par_vec[-(1:(3*K+1))]
	}
	
	return(complete_log_likelihood_general(y = y, X = X, Censoring_status = stat, g = g, lambda = lambda, 
		log_f = lw$log_f, log_F = lw$log_F, b = b, I_sim = I_sim, alpha = alpha)$cll
		)
}

my_vec0 <- c(1, 6, 1, 1.5, 1.5, -0.8)
f(my_vec0, family = 'exponential', I_sim = sim[,5])	
derivative(f, var = my_vec0, params = list(family = 'exponential', I_sim = sim[,5]))    


my_vec <- c(1, 6, 1, 2, 1.5, 1.5, -0.8)
f(my_vec, family = 'gamma', I_sim = sim[,5])	
derivative(f, var = my_vec, params = list(family = 'gamma', I_sim = sim[,5]))    

f(my_vec, family = 'weibull', I_sim = sim[,5])	
derivative(f, var = my_vec, params = list(family = 'weibull', I_sim = sim[,5]))    


my_vec <- c(1, 6, F_parameters[['gamma_parameters']][ ,1], F_parameters[['gamma_parameters']][ ,2], 0.5, 1.5, 1.5, -0.8)
f(my_vec, family = 'gamma_mixture', K = 2, I_sim = sim[,5])	
derivative(f, var = my_vec, params = list(family = 'gamma_mixture', K = 2, I_sim = sim[,5]))    

my_vec2 <- c(1, 6, c(F_parameters[['gamma_parameters']][ ,1],1), c(F_parameters[['gamma_parameters']][ ,2],3), c(0.5, 0.1), 1.5, 1.5, -0.8)
f(my_vec2, family = 'gamma_mixture', K = 3, I_sim = sim[,5])	
derivative(f, var = my_vec2, params = list(family = 'gamma_mixture', K = 3, I_sim = sim[,5]))    





