nChains <- 12
adjust = TRUE
mcmcIterations <- 10000
alpha <- 1/1.00001^{(1:nChains)^2.5 - 1}


# run an exponential model
promotion_time = list(distribution = 'exponential', 
                prior_parameters = matrix(rep(c(2.1, 1.1), 1), byrow = TRUE, 1, 2),
                prop_scale = c(0.1)
        )


run_exp <- cure_rate_MC3_general( 
        y = y, 
        X = X1,
        Censoring_status = stat,                
        nChains = nChains,              # number of heated chains (denoted as `C` in the paper)
        mcmc_cycles = mcmcIterations,   # number of MCMC cycles (denoted as `N` in the paper)
        nCores = 8,                     # number of cores to use for parallel processing
        sweep = 10,                     # number of iterations per MCMC cycle (denoted  as `m_1` in the paper)
        alpha = alpha,                  # a decreasing sequence of temperatures per chain 
#                                               (starting from 1 - the target chain is the first one). 
#
#                               PRIOR PARAMETERS        
#######################################################################################################
        mu_g = 1, s2_g = 1,     # this denotes a_g and b_g in the paper: gamma ~ f_g(a_g, b_g
        a_l = 2.1, b_l = 1.1,   # mean = 1, var = 10    lambda ~ IG(a_l, b_l)
        promotion_time = promotion_time,
########################################################################################################                
        plot = TRUE,            # opens plot device for plotting on the run the MCMC output
        adjust_scales = adjust,   # uses an initial phase in order to adjust the MH scale per parameter
        tau_mala = 0.0001,      # initial value of the MALA scale parameter - it will be adjusted as well.
        mala = 0,
        single_MH_in_f = 0.5
)
burn = floor(mcmcIterations/10)
getSummary_exp <- compute_map_and_hdis_general(y = y, X = X1, Censoring_status = stat, 
        retained_mcmc = run_exp$mcmc_sample[-(1:burn), ], prior_parameters = NULL, promotion_time = promotion_time)





# run a gamma model (including const. term)
promotion_time = list(distribution = 'gamma', 
                prior_parameters = matrix(rep(c(2.1, 1.1), 2), byrow = TRUE, 2, 2),
                prop_scale = c(0.1, 0.2)
        )

run_gamma <- cure_rate_MC3_general( 
        y = y, 
        X = X1,
        Censoring_status = stat,                
        nChains = nChains,              # number of heated chains (denoted as `C` in the paper)
        mcmc_cycles = mcmcIterations,   # number of MCMC cycles (denoted as `N` in the paper)
        nCores = 8,                     # number of cores to use for parallel processing
        sweep = 10,                     # number of iterations per MCMC cycle (denoted  as `m_1` in the paper)
        alpha = alpha,                  # a decreasing sequence of temperatures per chain 
#                                               (starting from 1 - the target chain is the first one). 
#
#                               PRIOR PARAMETERS        
#######################################################################################################
        mu_g = 1, s2_g = 1,     # this denotes a_g and b_g in the paper: gamma ~ f_g(a_g, b_g
        a_l = 2.1, b_l = 1.1,   # mean = 1, var = 10    lambda ~ IG(a_l, b_l)
        promotion_time = promotion_time,
########################################################################################################                
        plot = TRUE,            # opens plot device for plotting on the run the MCMC output
        adjust_scales = adjust,   # uses an initial phase in order to adjust the MH scale per parameter
        tau_mala = 0.0001,      # initial value of the MALA scale parameter - it will be adjusted as well.
        mala = 0,
        single_MH_in_f = 0.5
)
getSummary_gamma <- compute_map_and_hdis_general(y = y, X = X1, Censoring_status = stat, 
        retained_mcmc = run_gamma$mcmc_sample[-(1:burn), ], prior_parameters = NULL, promotion_time = promotion_time)


# run a weibull model (including const. term)

promotion_time = list(distribution = 'weibull', 
                prior_parameters = matrix(rep(c(2.1, 1.1), 2), byrow = TRUE, 2, 2),
                prop_scale = c(0.1, 0.2)
        )


run_weibull <- cure_rate_MC3_general( 
        y = y, 
        X = X1,
        Censoring_status = stat,                
        nChains = nChains,              # number of heated chains (denoted as `C` in the paper)
        mcmc_cycles = mcmcIterations,   # number of MCMC cycles (denoted as `N` in the paper)
        nCores = 8,                     # number of cores to use for parallel processing
        sweep = 10,                     # number of iterations per MCMC cycle (denoted  as `m_1` in the paper)
        alpha = alpha,                  # a decreasing sequence of temperatures per chain 
#                                               (starting from 1 - the target chain is the first one). 
#
#                               PRIOR PARAMETERS        
#######################################################################################################
        mu_g = 1, s2_g = 1,     # this denotes a_g and b_g in the paper: gamma ~ f_g(a_g, b_g
        a_l = 2.1, b_l = 1.1,   # mean = 1, var = 10    lambda ~ IG(a_l, b_l)
        promotion_time = promotion_time,
########################################################################################################                
        plot = TRUE,            # opens plot device for plotting on the run the MCMC output
        adjust_scales = adjust,   # uses an initial phase in order to adjust the MH scale per parameter
        tau_mala = 0.0001,      # initial value of the MALA scale parameter - it will be adjusted as well.
        mala = 0, single_MH_in_f = 0.5
)


getSummary_weibull <- compute_map_and_hdis(y = y, X = X1, Censoring_status = stat, 
        retained_mcmc = run_weibull$mcmc_sample[-(1:burn), ], prior_parameters = NULL)

#getSummary_weibull0 <- compute_map_and_hdis(y = y, X = X1, Censoring_status = stat, 
 #       retained_mcmc = run_weibull_0$mcmc_sample[-(1:burn), ], prior_parameters = NULL)


# run a gamma mixture model


K = 2
prior_parameters = array(data = NA, dim = c(2,2,K))
prior_parameters[,,1] = matrix(rep(c(2.1, 1.1), 2), byrow = TRUE, 2, 2)
prior_parameters[,,2] = matrix(rep(c(2.1, 1.1), 2), byrow = TRUE, 2, 2)

promotion_time = list(distribution = 'gamma_mixture',
                K = 2, 
                prior_parameters = prior_parameters,
                prop_scale = c(0.1, 0.1, 0.1, 0.1, 0.1 ),# a11, a21, ..., a1K, a2K, w_1, ...w_K
                dirichlet_concentration_parameter = 1
        )

run_gamma_mix2 <- cure_rate_MC3_general( 
        y = y, 
        X = X1,
        Censoring_status = stat,                
        nChains = nChains,              # number of heated chains (denoted as `C` in the paper)
        mcmc_cycles = mcmcIterations,   # number of MCMC cycles (denoted as `N` in the paper)
        nCores = 8,                     # number of cores to use for parallel processing
        sweep = 10,                     # number of iterations per MCMC cycle (denoted  as `m_1` in the paper)
        alpha = alpha,                  # a decreasing sequence of temperatures per chain 
#                                               (starting from 1 - the target chain is the first one). 
#
#                               PRIOR PARAMETERS        
#######################################################################################################
        mu_g = 1, s2_g = 1,     # this denotes a_g and b_g in the paper: gamma ~ f_g(a_g, b_g
        a_l = 2.1, b_l = 1.1,   # mean = 1, var = 10    lambda ~ IG(a_l, b_l)
        promotion_time = promotion_time,
########################################################################################################                
        plot = TRUE,            # opens plot device for plotting on the run the MCMC output
        adjust_scales = adjust,  # uses an initial phase in order to adjust the MH scale per parameter
        tau_mala = 0.0001,      # initial value of the MALA scale parameter - it will be adjusted as well.
        mala = 0, 
        single_MH_in_f = 0.5
)

getSummary_gamma_mix2 <- compute_map_and_hdis_general(y = y, X = X1, Censoring_status = stat, 
        retained_mcmc = run_gamma_mix2$mcmc_sample[-(1:burn), ], prior_parameters = NULL, promotion_time = promotion_time)


# run a gamma mixture model

K = 3
prior_parameters = array(data = NA, dim = c(2,2,K))
prior_parameters[,,1] = matrix(rep(c(2.1, 1.1), 2), byrow = TRUE, 2, 2)
prior_parameters[,,2] = matrix(rep(c(2.1, 1.1), 2), byrow = TRUE, 2, 2)
prior_parameters[,,3] = matrix(rep(c(2.1, 1.1), 2), byrow = TRUE, 2, 2)



promotion_time = list(distribution = 'gamma_mixture',
                K = 3, 
                prior_parameters = prior_parameters,
                prop_scale = rep(0.1, 8),# a11, a21, ..., a1K, a2K, w_1, ...w_K
                dirichlet_concentration_parameter = 1
        )




run_gamma_mix3 <- cure_rate_MC3_general( 
        y = y, 
        X = X1,
        Censoring_status = stat,                
        nChains = nChains,              # number of heated chains (denoted as `C` in the paper)
        mcmc_cycles = mcmcIterations,     # number of MCMC cycles (denoted as `N` in the paper)
        nCores = 8,                     # number of cores to use for parallel processing
        sweep = 10,                     # number of iterations per MCMC cycle (denoted  as `m_1` in the paper)
        alpha = alpha,                  # a decreasing sequence of temperatures per chain 
#                                               (starting from 1 - the target chain is the first one). 
#
#                               PRIOR PARAMETERS        
#######################################################################################################
        mu_g = 1, s2_g = 1,     # this denotes a_g and b_g in the paper: gamma ~ f_g(a_g, b_g
        a_l = 2.1, b_l = 1.1,   # mean = 1, var = 10    lambda ~ IG(a_l, b_l)
        promotion_time = promotion_time,
########################################################################################################                
        plot = TRUE,            # opens plot device for plotting on the run the MCMC output
        adjust_scales = FALSE,  # uses an initial phase in order to adjust the MH scale per parameter
        tau_mala = 0.0001,      # initial value of the MALA scale parameter - it will be adjusted as well.
        mala = 0, 
        single_MH_in_f = 0.5
)

getSummary_gamma_mix3 <- compute_map_and_hdis_general(y = y, X = X1, Censoring_status = stat, 
        retained_mcmc = run_gamma_mix3$mcmc_sample[-(1:burn), ], prior_parameters = NULL, promotion_time = promotion_time)


