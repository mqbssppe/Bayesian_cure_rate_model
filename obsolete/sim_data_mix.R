
library(pracma)



sim_fotis_model <- function(n, gamma1 = 1, lambda1 = 1.5, F_parameters, truebetas = c(1.5,1.5,-0.8), ab = 0.85, ranunL=0, ranunU=1, seed = NULL, discrcov, miny = 1e-6){
	if(is.null(seed) == FALSE){
		set.seed(seed)	}

  #??he first entries of truePars should be (in order): gamma, lambda, a1, a2, bvector (vector of regression coefficients, including a constant; 
  #requires at least two covariates, wherein the last is discrete. The discrete covariate should be unifom on {0,1,...,discrcov})
  
#  if(length(truePars)< 7){stop("At least 7 arguments should be provided in truePars, since we need at least two covariates.")}
  if(discrcov< 1){stop("The discrete covariate should be unifom on {0,1,...,discrcov}")}
  
#  gamma1 = truePars[1]
#	lambda1 = truePars[2]
#	a1 = truePars[3]
#	a2 = truePars[4]
#	truebetas=truePars[5:length(truePars)]


	nComponents <- F_parameters[['K']]
	mixing_proportions <- F_parameters[['p']]
	a_parameters <- F_parameters[['gamma_parameters']]


	my_objective_function <- function (x, a){
		k <- 1
		my_cdf <- -a +  mixing_proportions[k] * pgamma(x, shape = a_parameters[k,1], rate = a_parameters[k,2])
		if(nComponents > 1){
			for(k in 2:nComponents){
				my_cdf <- my_cdf + mixing_proportions[k] * pgamma(x, shape = a_parameters[k,1], rate = a_parameters[k,2])
			}
		}
		return(my_cdf)
		
	} 
	invert_cdf <- function(a_quant)return(uniroot(my_objective_function, c(0, 100), tol = 0.00001, a = a_quant)$root)
	



#	if(a1 < 0){stop("a1 < 0")}
#	if(a2 < 0){stop("a2 < 0")}	
	if(lambda1 < 0){stop("lambda < 0")}	
	if(ab < 0){stop("ab < 0")}		
	
	# keep true values for plotting purposes
	gTrue = gamma1
	lambdaTrue = lambda1 
#	a1True = a1
#	a2True =  a2
	betastrue=truebetas
	
 #The magic constant 
	cc=exp(exp(-1))

 #Simulation of censoring times
	if(ab > 0){
		lcenso=ab
		ctimes=rexp(n,rate=lcenso)
	}else{
#		ctimes = rep(1,n)
		ctimes = rchisq(n, df = 1)
	}

	
  #Simulate the design matrix: simulation of covariates where the first column has all its entries equal to one. 
	#All the covariates are continuous uniform on (ranunL,ranunU), except the last one whici is discrete uniform


	covardata=matrix(ncol=length(truebetas),nrow=n)
	
	
	for(j in 1:n){
	  covardata[j,1]=1}

	
	
	nnn=length(truebetas)-1
	for(i in 2:nnn){
	  for(j in 1:n){
	  covardata[j,i]=runif(1,ranunL,ranunU)}
	}
	
	
	  for(j in 1:n){
	    covardata[j,length(truebetas)]=sample(0:discrcov, 1, replace = TRUE)}

	
#Save the values of the regression model for each subject 
	
regrmodel=c()
for(i in 1:n){
  regrmodel[i]=exp(sum(truebetas*covardata[i,]))
  
}
	

#Computation of cure rate for each subject 
	pcure=c()
	for(i in 1:n){
	  pcure[i]=(1+gamma1*regrmodel[i]*cc^(regrmodel[i]*gamma1))^(-1/gamma1)
	  }
	


	#Simulation of cure status  for each subject 

	curestat=c()
	for(i in 1:n)
	{curestat=c(curestat,rbinom(1, 1, max(c(pcure[i],10^{-10}))))}

	
	#simulation of survival time of susceptibles, for Weibull lifetimes for Wj, given above
	

	omoio=runif(n,0,1)
	omionew=c()
	for(i in 1:n){
	  omionew[i]=(1-omoio[i])*(1-(1+gamma1*regrmodel[i]*(cc^(regrmodel[i]*gamma1)))^(-1/gamma1))+((1+gamma1*regrmodel[i]*(cc^(regrmodel[i]*gamma1)))^(-1/gamma1))
	}
	
	#Inverse of Population Survival for gamma
	Fpinvg=c()
	
#	if (promodis == 1){
#	  for(i in 1:n){
#	    Fpinvg[i]=	((-log(1-((omionew[i]^(-gamma1)-1)/(gamma1*regrmodel[i]*cc^(regrmodel[i]*gamma1)))^(1/lambda1)))^(1/a2))/a1
#	  }} 
#	else{
#	for(i in 1:n){
#	  Fpinvg[i]=qgamma(,shape = a1,rate =a2)
#}}
	
	for(i in 1:n){
		Fpinvg[i] =  invert_cdf(((omionew[i]^(-gamma1)-1)/(gamma1*regrmodel[i]*cc^(regrmodel[i]*gamma1)))^(1/lambda1))
		Fpinvg[i] = max(c(Fpinvg[i], miny))
	}
	
	


	nocov=ncol(covardata)-1
	nocov2=2+nocov
	data <- matrix( ncol=2+ncol(covardata) + 1,nrow = n)
#	data[,2+ncol(covardata) + 1] <- NA
	data[,3:nocov2] = covardata[,2:ncol(covardata)]
	
	
	for(i in 1:n){
	  if(curestat[i] == 1){
	    data[i,1]=ctimes[i]
	    data[i,2]=0
	    data[i, 2+ncol(covardata) + 1] = ctimes[i]
	  }else{
	    w=Fpinvg[i]
	    data[i,1]=min(w,ctimes[i])
	    data[i, 2+ncol(covardata) + 1] = w
	    if(min(w,ctimes[i])==ctimes[i]){
	      data[i,2]=0
	    }else{
	      data[i,2]=1}}}
	data[,ncol(data)-1] <- curestat

	mean(data[,ncol(data)-1])
	mean(1-data[,2])
	mean(data[,1])

	#Theoretical values, under the above set of values 
	Prob0=c()
	for(i in 1:n){
	Prob0[i]<-(1+gamma1*regrmodel[i]*(cc^(regrmodel[i]*gamma1)))^(-1/gamma1)}
	
	# Y_i = min(T_i, C_i), observed data
	# Censoring_status: 0 = censored, 1 = full time
	# Covariates: continuous and discrete uniform
	# cured_status: 1 = cured, 0 = susceptible (1 - I sto paper).
	names11=c()
	for(i in 1:nocov){
	  names11=c(names11,paste("Covariate",i,sep=""))}
		names11=c("Y","Censoring_status", names11, "cured_status", 'uncensored_time')
		colnames(data) <- names11
	
		return(data)
}



