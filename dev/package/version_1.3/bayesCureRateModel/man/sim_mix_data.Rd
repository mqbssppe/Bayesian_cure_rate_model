    \name{sim_mix_data}
     \docType{data}
     \alias{sim_mix_data}
     \title{Simulated dataset}
     \description{A synthetic dataset generated from a bimodal promotion time distribution. The available covariates are: 
     \describe{
     \item{\code{x1} }{ continuous. }
     \item{\code{x2} }{ factor with three levels.}
     }
     Among the 500 observations, there are 123 censoring times (\code{censoring = 0}) and 377 "events" (\code{censoring = 1}). The true status (cured or susceptible) is contained in the column \code{true_status} and contains 59 cured and 441 susceptible subjects.}
     \usage{data(sim_mix_data)}
     \format{Time-to-event data.}
     \keyword{datasets}

