######################.
##### Regression #####
######################.

if (cfg$which_sim=="regression") {

  one_simulation <- function() {

    # Generate dataset
    dat_obj <- generate_data(n=L$n, distr_A=L$distr_A,
                             theta_true=L$theta_true, sigma=L$sigma)
    dat <- dat_obj$dat
    theta_0 <- dat_obj$theta_0

    # Estimate regression function
    theta_n <- est_curve(dat, L$reg_type)$theta_n

    # Get estimated values
    sim_res <- list()
    for (val in round(seq(0,1,0.02),2)) {
      sim_res[[paste0("est_",format(val,nsmall=2))]] <- theta_n(val)
      sim_res[[paste0("theta_",format(val,nsmall=2))]] <- theta_0(val)
    }

    return(sim_res)

  }

}



###################.
##### Density #####
###################.

if (cfg$which_sim=="density") {

  one_simulation <- function() {

    # # Generate dataset
    # dat_obj <- generate_data(n=L$n, distr_A=L$distr_A,
    #                          theta_true=L$theta_true, sigma=L$sigma)
    # dat <- dat_obj$dat
    # theta_0 <- dat_obj$theta_0
    #
    # # Estimate regression function
    # theta_n <- est_curve(dat, L$reg_type)$theta_n
    #
    # # Get estimated values
    # sim_res <- list()
    # for (val in round(seq(0,1,0.02),2)) {
    #   sim_res[[paste0("est_",format(val,nsmall=2))]] <- theta_n(val)
    #   sim_res[[paste0("theta_",format(val,nsmall=2))]] <- theta_0(val)
    # }
    #
    # return(sim_res)

  }

}
