#' Generate data
#'
#' @param n Sample size
#' @param distr_A Marginal distribution of A; one of c("Unif(0,1)")
#' @param theta_true True regression function; one of c("identity", "square")
#' @param sigma Standard deviation of Y (!!!!! allow for heteroskedasticity)
#' @return A dataframe of X and Y variables
generate_data <- function(n, distr_A, theta_true, sigma) {

  # Generate A values
  if (distr_A=="Unif(0,1)") {
    a <- runif(n)
  } else if (distr_A=="N(0.5,0.09)") {
    a <- rtruncnorm(n, a=0, b=1, mean=0.5, sd=0.3)
  } else {
    stop("distr_A not valid")
  }

  # Generate Y values
  if (theta_true=="identity") {
    theta_0 <- function(x) { x }
  } else if (theta_true=="square") {
    theta_0 <- function(x) { x^2 }
  } else if (theta_true=="flat") {
    theta_0 <- function(x) { 0 }
  } else {
    stop("theta_true not valid")
  }
  y <- theta_0(a) + rnorm(n, sd=sigma)

  return(list(
    dat = data.frame(a=a, y=y),
    theta_0 = theta_0
  ))

}
