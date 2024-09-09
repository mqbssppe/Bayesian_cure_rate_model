    \name{marriage_dataset}
     \docType{data}
     \alias{marriage_dataset}
     \title{Marriage data}
     \description{The variable of interest (\code{time})  corresponds to the duration (in years) of first marriage for 1500 individuals. The available covariates are: 
     \describe{
     \item{\code{age} }{ age of respondent (in years) at the time of first marriage. The values are standardized (sample mean and variance equal to 0 and 1, respectively). }
     \item{\code{kids} }{ factor: whether there were kids during the first marriage (\code{"yes"}) or not (\code{"no"}).}
     \item{\code{race} }{ race of respondent with levels: \code{"black"}, \code{"hispanic"} and \code{"other"}.}
     }
     Among the 1500 observations, there are 1018 censoring times (\code{censoring = 0}) and 482 divorces (\code{censoring = 1}). Source: National Longitudinal Survey of Youth 1997 (NLSY97).}
     \usage{data(marriage_dataset)}
     \format{Time-to-event data.}
     \keyword{datasets}
     \references{
Bureau of Labor Statistics, U.S. Department of Labor. National Longitudinal Survey of Youth 1997 cohort, 1997-2022 (rounds 1-20). Produced and distributed by the Center for Human Resource Research (CHRR), The Ohio State University. Columbus, OH: 2023.
}

