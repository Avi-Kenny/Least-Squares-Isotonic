#' Estimate the regression function
#'
#' @param dat Dataset returned by generate_dataset
#' @param type Type of estimator. One of c("Linear", "Quadratic", "Iso GCM",
#'     "Iso GCM2", "Iso LS")
#' @return A regression estimator function
est_curve <- function(dat, type) {

  # Parametric (linear)
  if (type=="Linear") {
    model <- lm(y~a, dat=dat)
    cf <- as.numeric(coefficients(model))
    reg <- function(x) { cf[1] + cf[2]*x }
  }

  # Parametric (quadratic)
  if (type=="Quadratic") {
    model <- lm(y~a+I(a^2), dat=dat)
    cf <- as.numeric(coefficients(model))
    reg <- function(x) { cf[1] + cf[2]*x + cf[3]*x^2 }
  }

  # Standard isotonic regression estimator
  if (type %in% c("Iso GCM", "Iso LS")) {

    # Estimate primitive
    Gamma_n <- Vectorize(function(x) { mean(dat$y * as.integer(dat$a<=x)) })

    # Estimate empirical CDF
    Phi_n <- ecdf(dat$a)

    # Create CUSUM diagram
    cusum <- arrange(data.frame(x=c(0,Phi_n(dat$a)), y=c(0,Gamma_n(dat$a))), x)

    # Take the GCM and its derivative
    if (type=="Iso GCM") {
      GCM <- gcmlcm(x=cusum$x, y=cusum$y, type="gcm")
      dGCM <- approxfun(
        x = GCM$x.knots[-length(GCM$x.knots)],
        y = GCM$slope.knots,
        method = "constant",
        rule = 2,
        f = 0
      )
      reg <- function(x) { dGCM(Phi_n(x)) }
    }

    # Take the convex least squares line and its derivative
    if (type=="Iso LS") {
      LS <- cvx.lse.reg(t=cusum$x, z=cusum$y)
      pred_x <- round(seq(0,1,0.001),3)
      pred_y <- predict(LS, newdata=pred_x)
      dLS <- Vectorize(function(x) {
        width <- 0.01
        x1 <- x - width/2; x2 <- x + width/2;
        if (x1<0) { x2 <- x2 - x1; x1 <- 0; }
        if (x2>1) { x1 <- x1 - x2 + 1; x2 <- 1; }
        x1 <- round(x1,3); x2 <- round(x2,3);
        ind1 <- which(pred_x==x1); ind2 <- which(pred_x==x2);
        y1 <- pred_y[ind1]; y2 <- pred_y[ind2];
        return((y2-y1)/width)
      })
      reg <- function(x) { dLS(Phi_n(x)) }
    }

  }

  # Isotonic regression estimator, using Iso package
  # Note: this estimator jumps at the midpoint
  if (type=="Iso GCM2") {
    dat <- arrange(dat,a)
    reg <- Vectorize(function(x) {
      index <- which.min(abs(x-dat$a))
      pred <- Iso::pava(y=dat$y)
      return(pred[index])
    })
  }

  return(reg)

}
