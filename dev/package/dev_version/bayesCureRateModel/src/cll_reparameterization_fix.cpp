#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;


// c_under denotes a small positive threshold used for underflows in the cdf computation. 
// the default value should be c_under = pow(10.0,-9.0);

// [[Rcpp::export]]
List log_weibull(arma::vec y, double a1, double a2, double c_under){
	int i, n = y.n_elem;
	arma::vec log_f(n), log_F(n);
	double ct4;
	ct4 = a2*log(a1) + log(a2);

	for(i = 0; i < n; i++){
		log_F(i) = log(std::max(1.0 - exp(-pow(a1*y(i),a2)), c_under));
		log_f(i) = ct4 - pow(a1*y(i),a2) + (a2 - 1.0)*log(y(i));
	}
	
	List result;
	result["log_f"] = log_f;
	result["log_F"] = log_F;
	return(result);
}


// [[Rcpp::export]]
List complete_log_likelihood_general(arma::vec y, arma::mat X, IntegerVector Censoring_status, double g, double lambda, arma::vec log_f, arma::vec log_F, arma::vec b, IntegerVector I_sim, double alpha){
//	log_f: denotes the log_pdf of the promotion times, evaluated for each observation (i = 1,...n) for a given set of parameters
//	log_F: denotes the log_cdf of the promotion times, evaluated for each observation (i = 1,...n) for a given set of parameters
	int i, j, k, n = X.n_rows, nstiles = X.n_cols; 
	arma::vec exp_term(n), log_p0(n), log_S_p(n), log_theta(n);
	double cll, ct, ct3, sum_log_f_p, log_lambda, g_log_ct, sum_1_minus_I_log_p0_D0, LastTerm, log_ct, yy, log_g, y2, c_under;
	log_ct = exp(-1.0);
	ct = exp(log_ct);
	sum_log_f_p = 0.0;
	log_lambda = log(lambda);
	c_under = pow(10.0,-9.0);

	sum_1_minus_I_log_p0_D0 = 0.0;
	LastTerm = 0.0;
	cll = 0.0;
//	Rcout <<  X(0,0) << " " << X(0,1) << " " << X(0,2)<< "\n";				
	log_theta = X * b; 
	if(g < 0){	
//	if(1 < 2){	
		for(i = 0; i < n; i++){
//			log_theta(i) = b0 + b1*myData(i,2) + b2*myData(i,3);
			exp_term(i) = exp(log_theta(i));
			g_log_ct = g * exp_term(i);
			ct3 =  g_log_ct * pow(ct, g_log_ct);
			g_log_ct = g_log_ct*log_ct;
			log_p0(i) = -log(1.0 + ct3)/g; 
//			Rcout <<  log_p0(i) << " ";		
//			one_minus_exp(i) = std::max(1.0 - exp(-pow(a1*y(i),a2)), c_under);
			log_S_p(i) = - log(1.0 + ct3 * pow(exp(log_F(i)),lambda))/g;
	//		Rcout <<  log_S_p(i) << " ";				
			//Rcout <<  one_minus_exp(i) << " ";				
			if(Censoring_status(i) == 1){
				//D1 only here
//				log_weibull_dens(i) = ct4 - pow(a1*y(i),a2) + (a2 - 1.0)*log(y(i));			
	//			computing sum(log_f_p[D1])
				sum_log_f_p = sum_log_f_p +  (1.0 + g) * log_S_p(i) + log_lambda + log_theta(i) + g_log_ct + (lambda - 1.0)*log_F(i) + log_f(i);
	//			Rcout <<  sum_log_f_p << "\n";			
			}else{
			//D0 only here
	//			computing sum((1-I_sim[D0])*logP0[D0])
				sum_1_minus_I_log_p0_D0 = sum_1_minus_I_log_p0_D0 + (1.0 - I_sim(i))*log_p0(i);
				// LastTerm
				if(I_sim(i) == 1){
					LastTerm = LastTerm + log(std::max((exp(log_S_p(i)) - exp(log_p0(i))), c_under));	
				}
			}
		}
	}else{
		log_g = log(g);
//		logsumexps
		for(i = 0; i < n; i++){
//			log_theta(i) = b0 + b1*myData(i,2) + b2*myData(i,3);
			exp_term(i) = exp(log_theta(i));
			g_log_ct = g * exp_term(i);
			yy = log_theta(i) + log_g + g_log_ct * log_ct;
			log_p0(i) = -(yy+log(1.0 + exp(-yy)))/g; 
//			one_minus_exp(i) = std::max(1.0 - exp(-pow(a1*y(i),a2)), c_under);
			y2 = yy + lambda*log_F(i);
///////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////
//			Rcout <<  log_p0(i) << " ";		
			log_S_p(i) = - (y2 + log(1.0+exp(-y2)))/g;
			g_log_ct = g_log_ct*log_ct;
	//		Rcout <<  log_S_p(i) << " ";				
			//Rcout <<  one_minus_exp(i) << " ";				
			if(Censoring_status(i) == 1){
				//D1 only here
//				log_weibull_dens(i) = ct4 - pow(a1*y(i),a2) + (a2 - 1.0)*log(y(i));			
	//			computing sum(log_f_p[D1])
				sum_log_f_p = sum_log_f_p +  (1.0 + g) * log_S_p(i) + log_lambda + log_theta(i) + g_log_ct + (lambda - 1.0)*log_F(i) + log_f(i);
	//			Rcout <<  sum_log_f_p << "\n";			
			}else{
			//D0 only here
	//			computing sum((1-I_sim[D0])*logP0[D0])
				sum_1_minus_I_log_p0_D0 = sum_1_minus_I_log_p0_D0 + (1.0 - I_sim(i))*log_p0(i);
				// LastTerm
				if(I_sim(i) == 1){
					LastTerm = LastTerm + log(std::max((exp(log_S_p(i)) - exp(log_p0(i))), c_under));	
				}
			}
		}


			
	}
	cll = alpha * sum_log_f_p + alpha*sum_1_minus_I_log_p0_D0 + alpha*LastTerm;
	List result;
	result["cll"] = cll;
	result["logS"] = log_S_p;
	result["logP0"] = log_p0;		
	return(result);
}




