\name{compute_fdr_tpr}
\alias{compute_fdr_tpr}
\title{
Compute the achieved FDR and TPR
}
\description{
This help function computes the achieved False Discovery Rate (FDR) and True Positive Rate (TPR). Useful for simulation studies where the ground truth classification of subjects in susceptibles and cured items is known. 
}
\usage{
compute_fdr_tpr(true_latent_status, posterior_probs, 
	myCut = c(0.01, 0.05, 0.1, 0.15))
}
\arguments{
  \item{true_latent_status}{
	a binary vector containing the true latent status: 1 should correspond to the positive instances ("cured") and 0 to the negative ("susceptibles"). 
}
  \item{posterior_probs}{
	a numeric vector with entries between 0 and 1 containing the scores (posterior probabilities) of being positive ("cured") for each item.
}
  \item{myCut}{
	Vector containing the desired nominal FDR levels.
}
}
\value{
This function will return for every nominal FDR level the following quantities:
\item{achieved_fdr}{the achieved false discovery rate.}
\item{tpr}{the true positive rate.}
\item{nominal_fdr}{the nominal FDR level.}
}
\author{
Panagiotis Papastamoulis
}
\examples{
set.seed(1)
v1 <- sample(0:1, size = 100, replace=TRUE, prob=c(0.8,0.2) )
v2 <- runif(100)
compute_fdr_tpr(true_latent_status = v1, posterior_probs = v2)
}

