
library(pracma)
sim_fotis_model <- function(n, truePars, ab = 0, ranunL=0, ranunU=1, seed = NULL){
	if(is.null(seed) == FALSE){
		set.seed(seed)
	}
#	truePars should contain: gamma, lambda, a1, a2, b0, b1, b2
	if(length(truePars) != 7){stop("7 arguments should be provided in truePars.")}
	gamma1 = truePars[1]
	lambda1 = truePars[2]

	a1 = truePars[3]
	a2 = truePars[4]

	b0 = truePars[5]
	b1 = truePars[6]
	b2 = truePars[7]
	if(a1 < 0){stop("a1 < 0")}
	if(a2 < 0){stop("a2 < 0")}	
	if(lambda1 < 0){stop("lambda < 0")}	
	if(ab < 0){stop("ab < 0")}		
	
	# keep true values for plotting purposes
	gTrue = gamma1
	b0True = b0
	b1True = b1
	b2True = b2
	lambdaTrue = lambda1 
	a1True = a1
	a2True =  a2


	#the magic constant 
	cc=exp(exp(-1))

	#simulation of censoring times
	if(ab > 0){
		lcenso=ab
		ctimes=rexp(n,rate=lcenso)
	}else{
#		ctimes = rep(1,n)
		ctimes = rchisq(n, df = 1)
	}

	#simulation of covariates
	probc1=0.5
#	ranunL=0
#	ranunU=1
	cov1=rbinom(n,1,probc1)
#	cov1=runif(n,ranunL,ranunU)
	cov2=runif(n,ranunL,ranunU)

	#the CDF of promotion time, i.e. the common distribution of Wj
	FW<-function(t1){
	  valFW<-1-exp(-(a1*t1)^a2)
	return(valFW)}

	#the survival population function
	Sp<-function(t1,x1,x2){
	  ccocm<-(1+gamma1*exp(b0+b1*x1+b2*x2)*cc^(gamma1*exp(b0+b1*x1+b2*x2))*FW(t1)^lambda1)^(-1/gamma1)
	  return(ccocm)}


	#computation of cure rate for each subject 
	pcure=c()
	for(i in 1:n)
	{pcure[i]=(1+gamma1*exp(b0+b1*cov1[i]+b2*cov2[i])*cc^(exp(b0+b1*cov1[i]+b2*cov2[i])*gamma1))^(-1/gamma1)}



	#simulation of cure status

	curestat=c()
	for(i in 1:n)
	{curestat=c(curestat,rbinom(1, 1, max(c(pcure[i],10^{-10}))))}


	#simulation of survival time of susceptibles, for Weibull lifetimes for Wj, given above
	#Inverse of Population Survival and CD Function

	Spinv<-function(u1,x1,x2){
	  ccocm<-((-log(1-((u1^(-gamma1)-1)/(gamma1*exp(b0+b1*x1+b2*x2)*cc^(exp(b0+b1*x1+b2*x2)*gamma1)))^(1/lambda1)))^(1/a2))/a1
	 return(ccocm)}

	Fpinv<-function(u1,x1,x2){
	  cdcd<-Spinv((1-u1)*(1-(1+gamma1*exp(b0+b1*cov1[i]+b2*cov2[i])*(cc^(exp(b0+b1*cov1[i]+b2*cov2[i])*gamma1)))^(-1/gamma1))+((1+gamma1*exp(b0+b1*cov1[i]+b2*cov2[i])*(cc^(exp(b0+b1*cov1[i]+b2*cov2[i])*gamma1)))^(-1/gamma1))
		    ,x1,x2)
	  return(cdcd) }



	data<-matrix(ncol=5,nrow=n)

	for(i in 1:n){
	  data[i,3]=cov1[i]
	  data[i,4]=cov2[i]
	  if(curestat[i] == 1){
	    data[i,1]=ctimes[i]
	    data[i,2]=0
	  }else{
	    w=Fpinv(runif(1,0,1),cov1[i],cov2[i])
	    data[i,1]=min(w,ctimes[i])
	    if(min(w,ctimes[i])==ctimes[i]){
	      data[i,2]=0
	    }else{
	      data[i,2]=1}}}
	data[,5] <- curestat

	mean(data[,5])
	mean(1-data[,2])
	mean(data[,1])

	#Theoretical values, under the above set of values 

	Prob0<-function(x1,x2){
	  cscs<-(1+gamma1*exp(b0+b1*x1+b2*x2)*(cc^(exp(b0+b1*x1+b2*x2)*gamma1)))^(-1/gamma1)
	  return(cscs)
	}

	avercurerate<-function(r1){
	  xsxs<-(1/(ranunU-ranunL))*(Prob0(0,r1)*(1-probc1)+Prob0(1,r1)*(probc1))
	  return(xsxs)}

	#the average theoretical and observed cure rate

	integrate(avercurerate,ranunL,ranunU)
	mean(data[,5])#must be around 5%


	#########

#	if(ab > 0){
#	averSP<-function(t1,r1){
#	  xsxs<-((1/(ranunU-ranunL))*(Sp(t1,0,r1)*(1-probc1)+Sp(t1,1,r1)*(probc1)))*lcenso*exp(-lcenso*t1)
#	  return(xsxs)}
#
#	integral2(averSP,0,50,ranunL,ranunU, reltol = 1e-10)
#	mean(1-data[,2]) #must be around 15%
#	}
	#########

	mean(data[,1]) #must be around 0.20

	# Y_i = min(T_i, C_i), observed data
	# Censoring_status: 0 = censored, 1 = full time
	# Covariate 1: binary
	# Covariate 2: cont
	# cured_status: 1 = cured, 0 = susceptible (1 - I sto paper).
	colnames(data) <- c("Y", "Censoring_status", "Covariate1", "Covariate2", "cured_status")
	

	return(data)
}


