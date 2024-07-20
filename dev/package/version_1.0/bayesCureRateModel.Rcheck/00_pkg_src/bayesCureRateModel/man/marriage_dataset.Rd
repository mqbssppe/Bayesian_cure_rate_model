    \name{marriage_dataset}
     \docType{data}
     \alias{marriage_dataset}
     \title{Marriage data}
     \description{The variable of interest (\code{time})  corresponds to the duration (in years) of first marriage for 1500 individuals. The available covariates are: 
     \describe{
     \item{\code{age} }{ age of respondent (in years) at the time of fist marriage.}
     \item{\code{kids} }{ factor: whether there were kids during the first marriage (1) or not (0).}
     \item{\code{race} }{ race of respondent decoded as: black (1), hispanic (2) and non-black/non-hispanic (4).}
     }
     Among the 1500 observations, there are 1018 censoring times (\code{censoring = 0}) and 482 divorces (\code{censoring = 1}). Source: National Longitudinal Survey of Youth 1997 (NLSY97).}
     \usage{data(marriage_dataset)}
     \format{Time-to-event data.}
     \keyword{datasets}
     \references{
Bureau of Labor Statistics, U.S. Department of Labor. National Longitudinal Survey of Youth 1997 cohort, 1997-2022 (rounds 1-20). Produced and distributed by the Center for Human Resource Research (CHRR), The Ohio State University. Columbus, OH: 2023.
}

