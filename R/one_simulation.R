#' Run a single simulation
#'
#' @return A list
one_simulation <- function() {

  # Generate dataset
  dat_obj <- generate_data(n=L$n, distr_A=L$distr_A,
                           theta_true=L$theta_true, sigma=L$sigma)
  dat <- dat_obj$dat
  theta_0 <- dat_obj$theta_0

  # Estimate regression function
  reg <- est_curve(dat, L$reg_type)

  # Get estimated values
  sim_res <- list()
  for (val in round(seq(0,1,0.1),1)) {
    sim_res[[paste0("est_",format(val,nsmall=1))]] <- reg(val)
    sim_res[[paste0("theta_",format(val,nsmall=1))]] <- theta_0(val)
  }

  return(sim_res)

}
