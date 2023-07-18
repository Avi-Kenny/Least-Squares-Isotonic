#' Estimate the regression function
#'
#' @param dat Dataset returned by generate_dataset
#' @param type Type of estimator. One of c("Linear", "Quadratic", "Iso GCM",
#'     "Iso GCM2", "Iso CLS")
#' @param return_Theta_n Boolean; return primitive estimator
#' @param return_cusum Boolean; return CUSUM diagram
#' @return A regression estimator function (along with a primitive estimator,
#'     when applicable)
est_curve <- function(dat, type, return_Theta_n=F, return_cusum=F) {

  # Set placeholders
  Theta_n <- NA
  cusum <- NA

  # Parametric (linear)
  if (type=="Linear") {
    model <- lm(y~a, dat=dat)
    cf <- as.numeric(coefficients(model))
    theta_n <- function(x) { cf[1] + cf[2]*x }
  }

  # Parametric (quadratic)
  if (type=="Quadratic") {
    model <- lm(y~a+I(a^2), dat=dat)
    cf <- as.numeric(coefficients(model))
    theta_n <- function(x) { cf[1] + cf[2]*x + cf[3]*x^2 }
  }

  if (type %in% c("Iso GCM", "Iso GCM2", "Iso CLS")) {

    # Estimate primitive
    Gamma_n <- Vectorize(function(x) { mean(dat$y * as.integer(dat$a<=x)) })

    # Estimate empirical CDF
    Phi_n <- ecdf(dat$a)

    # Create CUSUM diagram
    cusum <- arrange(data.frame(x=c(0,Phi_n(dat$a)), y=c(0,Gamma_n(dat$a))), x)

    # Compute the GCM
    if (type %in% c("Iso GCM", "Iso GCM2")) {
      GCM <- gcmlcm(x=cusum$x, y=cusum$y, type="gcm")
      Theta_n <- approxfun(
        x = GCM$x.knots,
        y = GCM$y.knots,
        method = "linear",
        rule = 2
      )
    }

    # Take the GCM and its derivative
    if (type=="Iso GCM") {
      dGCM <- approxfun(
        x = GCM$x.knots[-length(GCM$x.knots)],
        y = GCM$slope.knots,
        method = "constant",
        rule = 2,
        f = 0
      )
      theta_n <- function(x) { dGCM(Phi_n(x)) }
    }

    # Isotonic regression estimator, using Iso package
    # Note: this estimator jumps at the midpoint between two points
    if (type=="Iso GCM2") {
      dat <- arrange(dat,a)
      theta_n <- Vectorize(function(x) {
        index <- which.min(abs(x-dat$a))
        pred <- Iso::pava(y=dat$y)
        return(pred[index])
      })
    }

    # Take the convex least squares line and its derivative
    if (type=="Iso CLS") {
      CLS <- cvx.lse.reg(t=cusum$x, z=cusum$y)
      pred_x <- round(seq(0,1,0.001),3)
      pred_y <- predict(CLS, newdata=pred_x)
      dCLS <- Vectorize(function(x) {
        width <- 0.01
        x1 <- x - width/2; x2 <- x + width/2;
        if (x1<0) { x2 <- x2 - x1; x1 <- 0; }
        if (x2>1) { x1 <- x1 - x2 + 1; x2 <- 1; }
        x1 <- round(x1,3); x2 <- round(x2,3);
        ind1 <- which(pred_x==x1); ind2 <- which(pred_x==x2);
        y1 <- pred_y[ind1]; y2 <- pred_y[ind2];
        return((y2-y1)/width)
      })
      theta_n <- function(x) { dCLS(Phi_n(x)) }
      Theta_n <- Vectorize(function(x) {
        ind <- which.min(abs(x-pred_x))
        return(pred_y[ind])
      })
    }

  }

  # Temporary: look at difference between two estimators
  if (type=="difference") {

    # !!!!! Check this code to ensure consistency with above

    # Gamma_n <- Vectorize(function(x) { mean(dat$y * as.integer(dat$a<=x)) })
    # Phi_n <- ecdf(dat$a)
    # cusum <- arrange(data.frame(x=c(0,Phi_n(dat$a)), y=c(0,Gamma_n(dat$a))), x)
    # CLS <- cvx.lse.reg(t=cusum$x, z=cusum$y)
    # pred_x <- round(seq(0,1,0.001),3)
    # pred_y <- predict(CLS, newdata=pred_x)
    # dCLS <- Vectorize(function(x) {
    #   width <- 0.01
    #   x1 <- x - width/2; x2 <- x + width/2;
    #   if (x1<0) { x2 <- x2 - x1; x1 <- 0; }
    #   if (x2>1) { x1 <- x1 - x2 + 1; x2 <- 1; }
    #   x1 <- round(x1,3); x2 <- round(x2,3);
    #   ind1 <- which(pred_x==x1); ind2 <- which(pred_x==x2);
    #   y1 <- pred_y[ind1]; y2 <- pred_y[ind2];
    #   return((y2-y1)/width)
    # })
    # theta_n1 <- function(x) { dCLS(Phi_n(x)) }
    #
    # # Iso GCM2
    # dat <- arrange(dat,a)
    # theta_n2 <- Vectorize(function(x) {
    #   index <- which.min(abs(x-dat$a))
    #   pred <- Iso::pava(y=dat$y)
    #   return(pred[index])
    # })
    #
    # theta_n <- function (x) { theta_n1(x) - theta_n2(x) }

  }

  res <- list(theta_n=theta_n)
  if (return_Theta_n) { res$Theta_n <- Theta_n }
  if (return_cusum) { res$cusum <- cusum }
  return(res)

}
