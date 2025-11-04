# bayesCureRateModel: Bayesian Cure Rate Modeling for Time-to-Event Data

A fully Bayesian approach in order to estimate a general family of cure rate models under the presence of covariates, see [*Papastamoulis and Milienos, 2024a*](https://doi.org/10.1007/s11749-024-00942-w) and [Papastamoulis and Milienos, 2024b](https://arxiv.org/abs/2409.10221). 

The promotion time can be modelled 
* parametrically using typical distributional assumptions for time to event data (including the Weibull, Exponential, Gompertz, log-Logistic distributions), or 
* semiparametrically using finite mixtures of distributions.
 
In both cases, user-defined families of distributions are allowed under some specific requirements. Posterior inference is carried out by constructing a Metropolis-coupled Markov chain Monte Carlo (MCMC) sampler, which combines Gibbs sampling for the latent cure indicators and Metropolis-Hastings steps with Langevin diffusion dynamics for parameter updates. The main MCMC algorithm is embedded within a parallel tempering scheme by considering heated versions of the target posterior distribution.

The R package `bayesCureRateModel` package is available on [CRAN](https://CRAN.R-project.org/package=bayesCureRateModel). The latest version is 1.5 (1/11/2025). 

## References

**Papastamoulis P and Milienos FS (2024a)**. [*Bayesian inference and cure rate
modeling for event history data*](https://doi.org/10.1007/s11749-024-00942-w). TEST. 

**Papastamoulis P and Milienos FS (2024b)**. [*bayesCureRateModel: Bayesian Cure Rate Modeling for Time to Event Data in R*](https://arxiv.org/abs/2409.10221). arXiv pre-print. 






