\name{print.bayesCureModel}
\alias{print.bayesCureModel}
\title{
Print method
}
\description{
This function prints a summary of objects returned by the \code{cure_rate_MC3} function.
}
\usage{
\method{print}{bayesCureModel}(x, ...)
}
\arguments{
  \item{x}{
An object of class \code{bayesCureModel}, which is returned by the \code{cure_rate_MC3} function.
}
  \item{...}{
        ignored.
}
}
\value{No return value, called for side effects.}
\details{
The function prints some basic information for a \code{cure_rate_MC3}, such as the MAP estimate of model parameters and the value of Bayesian information criterion.
}
\author{
Panagiotis Papastamoulis
}

