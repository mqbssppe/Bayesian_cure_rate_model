\name{Surv}
\alias{Surv}
\title{Create a Survival Object}
\description{
  Create a survival object for use in survival analysis, as imported from the \pkg{survival} package.
}
\usage{
Surv(time, time2, event, type = c("right", "left", "interval", 
	"counting", "interval2"), origin = 0)
}
\arguments{
  \item{time}{The follow-up time. For counting process data, this is the start time.}
  \item{time2}{The end time for counting process or interval-censored data. Optional for right-censored data.}
  \item{event}{The event indicator, normally 0=alive/censored, 1=dead/event occurred. Can also be a factor or logical.}
  \item{type}{Type of censoring. Options are \code{"right"}, \code{"left"}, \code{"interval"}, \code{"counting"} or \code{"interval2"}.}
  \item{origin}{Starting point for time scale. Default is 0. Only used for type \code{"counting"}.}
}
\details{
  The \code{Surv} function creates an object of class \code{"Surv"}, which is used to represent survival data. Depending on the arguments, the object can represent different types of censoring.

  \itemize{
    \item Right-censoring: one time and event indicator.
    \item Left-censoring: similar to right-censoring but event=1 for censored.
    \item Interval-censoring: requires both \code{time} and \code{time2}.
    \item Counting process: both \code{time} and \code{time2} used to specify start and stop times.
  }

  The resulting object is used as a response in survival regression models and estimation functions. 
}
\value{
  An object of class \code{"Surv"} which is used as a response in survival models.
}
\note{
The implementation in the \pkg{bayesCureRateModel} package only supports right-censored data. The binary censoring indicators are interpreted as a time-to-event (1) or as a censoring time (0).
}
\seealso{
  \code{\link{cure_rate_MC3}}
}
\references{
Therneau T (2026). A Package for Survival Analysis in R. R package version 3.8-5, https://CRAN.R-project.org/package=survival. 
}

\examples{
# Right-censored survival data
Surv(5, 1)
Surv(c(5, 10), c(1, 0))
}
\keyword{survival}

