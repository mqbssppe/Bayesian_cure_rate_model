# Bayesian_cure_rate_model
Bayesian inference and cure rate modeling

Code for the paper **Papastamoulis and Milienos (2023)**. *Bayesian inference and cure rate
modeling for event history data*. [arXiv:2310.06926](https://arxiv.org/abs/2310.06926)

## Illustrative examples

The file [example.R](https://github.com/mqbssppe/Bayesian_cure_rate_model/blob/main/example.R) generates synthetic data and then applies the proposed methodology. 

The file [Recidivism_Iowa_5000.txt](https://github.com/mqbssppe/Bayesian_cure_rate_model/blob/main/example.R) contains the recidivism dataset used in our paper. 

## Required R packages
- [Rcpp](https://CRAN.R-project.org/package=Rcpp)
- [RcppArmadillo](https://CRAN.R-project.org/package=RcppArmadillo)
- [coda](https://CRAN.R-project.org/package=coda)
- [doParallel](https://CRAN.R-project.org/package=doParallel)
- [foreach](https://CRAN.R-project.org/package=foreach)
- [pracma](https://CRAN.R-project.org/package=pracma)

## Updates	

The package `bayesCureRateModel` (version 1.1) is now available on [CRAN](https://CRAN.R-project.org/package=bayesCureRateModel ) 

Developer versions will be available here 

* [version 1.1](https://github.com/mqbssppe/Bayesian_cure_rate_model/tree/main/dev/package/version_1.1): fixed a bug on Windows 


