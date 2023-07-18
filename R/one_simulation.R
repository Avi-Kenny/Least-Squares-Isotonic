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



########################.
##### Dissertation #####
########################.

if (cfg$which_sim=="dissertation") {

  one_simulation <- function() {

    # Generate dataset
    batch({
      dat_obj <- generate_data(n=L$n, distr_A=L$distr_A,
                               theta_true=L$theta_true, sigma=L$sigma)
    })
    dat <- dat_obj$dat
    theta_0 <- dat_obj$theta_0

    # Estimate regression function
    theta_n <- est_curve(dat, "Iso GCM")$theta_n
    theta_s <- est_curve(dat, "Iso CLS")$theta_n

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

    # Generate data from an exponential distribution
    x <- rexp(n=L$n, rate=1)

    # True density function
    theta_0 <- function(x) { dexp(x=x, rate=1) }

    # Estimate density
    {
      # Phi_n is just the identity function
      Phi_n <- function(x) { x }

      # Gamma_n is the (inverted) ecdf
      F_n <- ecdf(x)
      Gamma_n <- function(x) { -1 * F_n(x) }

      # Psi_n equals Gamma_n
      Psi_n <- Gamma_n
    }

    # Create CUSUM diagram
    cusum <- dplyr::arrange(data.frame(x=c(0,x), y=c(0,Gamma_n(x))), x)

    if (L$type=="GCM") {

      GCM <- fdrtool::gcmlcm(x=cusum$x, y=cusum$y, type="gcm")
      dGCM <- approxfun(
        x = GCM$x.knots[-length(GCM$x.knots)],
        y = GCM$slope.knots,
        method = "constant",
        rule = 2,
        f = 0
      )
      theta_n <- function(x) { -1 * dGCM(Phi_n(x)) }

    } else if (L$type=="CLS") {

      CLS <- simest::cvx.lse.reg(t=cusum$x, z=cusum$y)
      pred_x <- seq(0, max(x), length.out=1000)
      pred_y <- predict(CLS, newdata=pred_x)
      dCLS <- Vectorize(function(x) {
        index <- which.min(abs(x-CLS$x.values))
        return(CLS$deriv[index])
      })
      theta_n <- function(x) { -1 * dCLS(Phi_n(x)) }

    } else {
      stop("Invalid type")
    }

    # Get estimated values
    sim_res <- list()
    for (val in round(seq(0,1,0.02),2)) {
      sim_res[[paste0("est_",format(val,nsmall=2))]] <- theta_n(val)
      sim_res[[paste0("theta_",format(val,nsmall=2))]] <- theta_0(val)
    }

    return(sim_res)

  }

}
