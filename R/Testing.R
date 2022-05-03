
# Test different regression methods
{

  # Generate dataset
  L <- list(
    n = 50,
    distr_A = "N(0.5,0.04)", # "Unif(0,1)"
    theta_true = "identity"
  )
  dat_obj <- generate_data(n=L$n, distr_A=L$distr_A, theta_true=L$theta_true)
  dat <- dat_obj$dat
  theta_0 <- dat$theta_0

  # Estimate regression
  reg_Linear <- est_curve(dat, type="Linear")
  reg_IsoGCM <- est_curve(dat, type="Iso GCM")
  reg_IsoGCM2 <- est_curve(dat, type="Iso GCM2")
  reg_IsoLS <- est_curve(dat, type="Iso LS")

  # Plot results
  grid <- seq(0,1,0.001)
  df_plot <- data.frame(
    a = rep(grid, 4),
    y = c(reg_Linear(grid), reg_IsoGCM(grid), reg_IsoGCM2(grid), reg_IsoLS(grid)),
    which = rep(c("Linear", "Iso GCM", "Iso GCM2", "Iso LS"), each=length(grid))
  )
  ggplot(df_plot, aes(x=a, y=y, color=which)) +
    geom_point(aes(x=a, y=y), data=dat, inherit.aes=F, alpha=0.5) +
    geom_line() +
    labs(color="Estimator")

}
