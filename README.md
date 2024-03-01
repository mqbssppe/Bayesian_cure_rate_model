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

## See dev version 

The [https://github.com/mqbssppe/Bayesian_cure_rate_model/tree/main/dev](developer version) of the package is now available. Main features
* allow general number of covariates (with/out constant terms)
* dedicated functions for plotting and summarizing the output.

## Session info

```
> sessionInfo()         
R version 4.2.1 (2022-06-23)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 22.04.2 LTS

Matrix products: default
BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.10.0
LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.10.0

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
 [3] LC_TIME=el_GR.UTF-8        LC_COLLATE=en_US.UTF-8    
 [5] LC_MONETARY=el_GR.UTF-8    LC_MESSAGES=en_US.UTF-8   
 [7] LC_PAPER=el_GR.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=el_GR.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] parallel  stats     graphics  grDevices utils     datasets  methods  
[8] base     

other attached packages:
[1] pracma_2.3.8             doParallel_1.0.17        iterators_1.0.14        
[4] foreach_1.5.2            coda_0.19-4              RcppArmadillo_0.10.5.0.0
[7] Rcpp_1.0.8.3            

loaded via a namespace (and not attached):
[1] compiler_4.2.1   codetools_0.2-18 grid_4.2.1       lattice_0.20-44 
```
