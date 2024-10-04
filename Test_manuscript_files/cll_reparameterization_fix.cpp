#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export]]
List complete_log_likelihood(arma::mat myData, double g, double lambda, double a1, double a2, double b0, double b1, double b2, IntegerVector I_sim, double alpha){
	int i, j, k, n = myData.n_rows;; 
	arma::vec exp_term(n), log_p0(n), log_S_p(n), one_minus_exp(n),log_weibull_dens(n), log_theta(n);
	double cll, ct, ct3, ct4, c_under, sum_log_f_p, log_lambda, g_log_ct, sum_1_minus_I_log_p0_D0, LastTerm, log_ct, y, log_g, y2;	 	
	log_ct = exp(-1.0);
	ct = exp(log_ct);
	c_under = pow(10.0,-9.0);
	ct4 = a2*log(a1) + log(a2);
	sum_log_f_p = 0.0;
	log_lambda = log(lambda);
	
	sum_1_minus_I_log_p0_D0 = 0.0;
	LastTerm = 0.0;
	cll = 0.0;
	if(g < 0){	
//	if(1 < 2){	
		for(i = 0; i < n; i++){
			log_theta(i) = b0 + b1*myData(i,2) + b2*myData(i,3);
			exp_term(i) = exp(log_theta(i));
			g_log_ct = g * exp_term(i);
			ct3 =  g_log_ct * pow(ct, g_log_ct);
			g_log_ct = g_log_ct*log_ct;
			log_p0(i) = -log(1.0 + ct3)/g; 
//			Rcout <<  log_p0(i) << " ";		
			one_minus_exp(i) = std::max(1.0 - exp(-pow(a1*myData(i,0),a2)), c_under);
			log_S_p(i) = - log(1.0 + ct3 * pow(one_minus_exp(i),lambda))/g;
	//		Rcout <<  log_S_p(i) << " ";				
			//Rcout <<  one_minus_exp(i) << " ";				
			if(myData(i,1) == 1){
				//D1 only here
				log_weibull_dens(i) = ct4 - pow(a1*myData(i,0),a2) + (a2 - 1.0)*log(myData(i,0));			
	//			computing sum(log_f_p[D1])
				sum_log_f_p = sum_log_f_p +  (1.0 + g) * log_S_p(i) + log_lambda + log_theta(i) + g_log_ct + (lambda - 1.0)*log(one_minus_exp(i)) + log_weibull_dens(i);
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
			log_theta(i) = b0 + b1*myData(i,2) + b2*myData(i,3);
			exp_term(i) = exp(log_theta(i));
			g_log_ct = g * exp_term(i);
			y = log_theta(i) + log_g + g_log_ct * log_ct;
			log_p0(i) = -(y+log(1.0 + exp(-y)))/g; 
			one_minus_exp(i) = std::max(1.0 - exp(-pow(a1*myData(i,0),a2)), c_under);
			y2 = y + lambda*log(one_minus_exp(i));
///////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////
//			Rcout <<  log_p0(i) << " ";		
			log_S_p(i) = - (y2 + log(1.0+exp(-y2)))/g;
			g_log_ct = g_log_ct*log_ct;
	//		Rcout <<  log_S_p(i) << " ";				
			//Rcout <<  one_minus_exp(i) << " ";				
			if(myData(i,1) == 1){
				//D1 only here
				log_weibull_dens(i) = ct4 - pow(a1*myData(i,0),a2) + (a2 - 1.0)*log(myData(i,0));			
	//			computing sum(log_f_p[D1])
				sum_log_f_p = sum_log_f_p +  (1.0 + g) * log_S_p(i) + log_lambda + log_theta(i) + g_log_ct + (lambda - 1.0)*log(one_minus_exp(i)) + log_weibull_dens(i);
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
