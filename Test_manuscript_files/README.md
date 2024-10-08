# Bayesian_cure_rate_model
Bayesian inference and cure rate modeling

Code for the paper **Papastamoulis and Milienos (2024)**. [*Bayesian inference and cure rate
modeling for event history data*](https://doi.org/10.1007/s11749-024-00942-w). TEST. 

## Illustrative examples

The file [example.R](https://github.com/mqbssppe/Bayesian_cure_rate_model/blob/main/Test_manuscript_files/example.R) generates synthetic data and then applies the proposed methodology. 

The file [Recidivism_Iowa_5000.txt](https://github.com/mqbssppe/Bayesian_cure_rate_model/blob/main/Test_manuscript_files/Recidivism_Iowa_5000.txt) contains the recidivism dataset used in our paper. 

## Required R packages
- [Rcpp](https://CRAN.R-project.org/package=Rcpp)
- [RcppArmadillo](https://CRAN.R-project.org/package=RcppArmadillo)
- [coda](https://CRAN.R-project.org/package=coda)
- [doParallel](https://CRAN.R-project.org/package=doParallel)
- [foreach](https://CRAN.R-project.org/package=foreach)
- [pracma](https://CRAN.R-project.org/package=pracma)

## Note 
These files are based on a preliminary version of our code. For a more user-friendly implementation and various extensions of the method, see the R package `bayesCureRateModel`. 


