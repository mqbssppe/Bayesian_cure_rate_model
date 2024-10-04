


R> library(RcppArmadillo)
R> RcppArmadillo.package.skeleton(name = "bayesCureRateModel", code_files = "bayesian_cure_rate_model.R", example_code=FALSE)

T> cp fabMix.cpp fabMix/src/fabMix.cpp

meta edit description kai namespace

R> library(Rcpp)
R> compileAttributes("bayesCureRateModel")

meta edit to man

meta R CMD buil check etc

