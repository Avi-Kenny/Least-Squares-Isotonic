######################.
###### Debugging #####
######################.

if (F) {

  L <- list(n=40, distr_A="Unif(0,1)", theta_true="identity", sigma=0.2)
  dat <- generate_data(n=L$n, distr_A=L$distr_A,
                           theta_true=L$theta_true, sigma=L$sigma)$dat

  Gamma_n <- Vectorize(function(x) { mean(dat$y * as.integer(dat$a<=x)) })
  Phi_n <- ecdf(dat$a)
  cusum <- arrange(data.frame(x=c(0,Phi_n(dat$a)), y=c(0,Gamma_n(dat$a))), x)

  # Compute the GCM
  GCM <- gcmlcm(x=cusum$x, y=cusum$y, type="gcm")
  GCM_fn <- approxfun(
    x = GCM$x.knots,
    y = GCM$y.knots,
    method = "linear",
    rule = 2
  )

  # Compute CLS
  CLS <- cvx.lse.reg(t=cusum$x, z=cusum$y)
  pred_x <- round(seq(0,1,0.001),3)
  pred_y <- predict(CLS, newdata=pred_x)
  CLS_fn <- Vectorize(function(x) {
    index <- which.min(abs(x-pred_x))
    return(pred_y[index])
  })

  # ggplot(data.frame(x=dat$a, y=dat$y), aes(x=x, y=y)) + geom_point() + labs(title="orig data")
  grid <- seq(0,1,0.001)
  plot_df <- data.frame(
    x = rep(grid,2),
    y = c(GCM_fn(grid), CLS_fn(grid)),
    which = rep(c("GCM", "CLS"), each=length(grid))
  )
  ggplot(plot_df, aes(x=x, y=y, group=which, color=which)) +
    geom_line() +
    geom_point(
      aes(x=x, y=y),
      data = data.frame(x=cusum$x, y=cusum$y),
      alpha = 0.5,
      inherit.aes = F
    ) +
    labs(title="cusum")
  print(CLS_fn(0.5)-GCM_fn(0.5))

}

######################################################.
###### Mini simulation to investigate edge issue #####
######################################################.

library(SimEngine)

# Function to generate data
generate_data <- function(n, edge=0) {

  A <- runif(n)
  if (edge!=0) { A <- A * rbinom(n, size=1, prob=edge) }
  W <- runif(n)
  Y <- 2*A + 3*W + rnorm(n)

  return(data.frame(A=A, W=W, Y=Y))

}

# Define other functions
In <- as.integer
mu_0 <- function(a,w) { 2*a + 3*w }
theta_0 <- function(a) {
  mean(sapply(seq(0,1,0.001), function(w) { mu_0(a,w) }))
}
mu_n <- mu_0
g_n <- function(a,w) { rep(1,length(a)) }

# Simulation
sim <- new_sim()
sim %<>% set_levels(n=c(50), edge=c(0,0.8))
sim %<>% set_config(num_sim=100)
sim %<>% set_script(function() {

  # print("rep") # !!!!!

  # Generate dataset
  d <- generate_data(n=L$n, edge=L$edge)

  # Compute ecdf transformation
  F_n <- ecdf(d$A)

  # Compute one-step estimator
  Gamma_n <- function(a) {
    piece_1 <- mean( In(d$A<=a)*((d$Y-mu_n(d$A,d$W))/(g_n(d$A,d$W))) )
    piece_2 <- 0
    for (i in c(1:L$n)) {
      for (j in c(1:L$n)) {
        piece_2 <- piece_2 + In(d$A[i]<=a)*mu_n(d$A[i],d$W[j])
      }
    }
    piece_2 <- piece_2/(L$n)^2
    return(piece_1+piece_2)
  }

  # Perform GCM procedure
  gcm_x_vals <- sapply(sort(unique(d$A)), F_n)
  indices_to_keep <- !base::duplicated(gcm_x_vals)
  gcm_x_vals <- gcm_x_vals[indices_to_keep]
  gcm_y_vals <- sapply(sort(unique(d$A))[indices_to_keep], Gamma_n)
  if (!any(gcm_x_vals==0)) {
    gcm_x_vals <- c(0, gcm_x_vals)
    gcm_y_vals <- c(0, gcm_y_vals)
  }
  gcm <- fdrtool::gcmlcm(x=gcm_x_vals, y=gcm_y_vals, type="gcm")
  dGCM <- approxfun(
    x = gcm$x.knots[-length(gcm$x.knots)],
    y = gcm$slope.knots,
    method = "constant",
    rule = 2,
    f = 0
  )
  theta_n <- function(a) { dGCM(F_n(a)) }

  return (list(
    "theta_n_0.0" = theta_n(0),
    "theta_0_0.0" = theta_0(0),
    "theta_n_0.2" = theta_n(0.2),
    "theta_0_0.2" = theta_0(0.2),
    "theta_n_0.5" = theta_n(0.5),
    "theta_0_0.5" = theta_0(0.5),
    "F_n_0.0" = F_n(0),
    "F_n_0.1" = F_n(0.1),
    "F_n_0.3" = F_n(0.3)
  ))

})
sim %<>% run()

sim %>% summarize(
  list(stat="mean", x="theta_n_0.0"),
  list(stat="mean", x="theta_0_0.0"),
  list(stat="mean", x="theta_n_0.2"),
  list(stat="mean", x="theta_0_0.2"),
  list(stat="mean", x="theta_n_0.5"),
  list(stat="mean", x="theta_0_0.5"),
  list(stat="mean", x="F_n_0.0"),
  list(stat="mean", x="F_n_0.1"),
  list(stat="mean", x="F_n_0.3")
)


# ggplot(df, aes(x=x1, y=y1, color=factor(grp2))) + geom_line()





chk(18)
if (p$convex_type=="GCM") {
} else if (p$convex_type=="CLS") {
}
chk(19)

# Construct Grenander-based r_Mn estimator (and truncate to lie within [0,1])
