# Title: "Alternative least squares isotonic regression estimator"
# Author: Avi Kenny

##################.
##### CONFIG #####
##################.

# To run multiple sims/analyses concurrently, ONLY change Slurm commands

# Set global config
cfg <- list(
  main_task = "run", # run update
  which_sim = "dissertation (density)", # "regression" "density" "dissertation (regression)" "dissertation (density)"
  level_set_which = "level_set_density_2", # level_set_regression_1 level_set_regression_2 level_set_density_1 level_set_density_2
  num_sim = 10000,
  pkgs = c("dplyr", "Iso", "fdrtool", "simest", "tidyr", "truncnorm",
           "modeest"),
  pkgs_nocluster = c("ggplot2", "ChernoffDist"),
  parallel = "none",
  n_cores = 500,
  stop_at_error = FALSE
)

# Set cluster config
cluster_config <- list(
  # js = "ge",
  # dir = paste0("/home/users/avikenny/Desktop/", Sys.getenv("project"))
  js = "slurm",
  dir = paste0("/home/akenny/", Sys.getenv("project"))
)



#################.
##### SETUP #####
#################.

# Set local vs. cluster variables
if (Sys.getenv("USERDOMAIN")=="AVI-KENNY-T460") {
  # Local
  setwd(paste0("C:/Users/avike/OneDrive/Desktop/Avi/Biostats + Research/Resear",
               "ch/Marco Carone/Project - LSI/Least-Squares-Isotonic/R"))
  load_pkgs_local <- TRUE
} else {
  # Cluster
  setwd(paste0(cluster_config$dir, "/R"))
  if (cfg$main_task %in% c("run", "update")) {
    load_pkgs_local <- FALSE
  } else {
    load_pkgs_local <- TRUE
  }
}

# Load packages (if running locally)
if (load_pkgs_local) {
  for (pkg in c(cfg$pkgs,cfg$pkgs_nocluster)) {
    suppressMessages({ do.call("library", list(pkg)) })
  }
}

# Load SimEngine + functions
{
  library(SimEngine)
  source("one_simulation.R", local=TRUE)
  source("generate_data.R", local=TRUE)
  source("est_curve.R", local=TRUE)
}



##########################################################.
##### MAIN: Set level sets for different simulations #####
##########################################################.

if (Sys.getenv("sim_run") %in% c("first", "")) {

  # Regression: main
  level_set_regression_1 <- list(
    n = c(100,200,400,800,1600,3200),
    distr_A = c("Unif(0,1)", "N(0.5,0.09)"), # "Unif(0,1)", "N(0.5,0.09)"
    theta_true = c("identity"), # "identity", "square", "constant"
    reg_type = c("Iso CLS", "Iso GCM"), # "Linear", "Iso GCM", "Iso GCM2", "Iso CLS"
    sigma = 0.2
  )

  # Regression: dissertation
  level_set_regression_2 <- level_set_regression_1
  level_set_regression_2$reg_type <- NULL

  # # Regression: difference between estimators
  # level_set_regression_3 <- level_set_regression_1
  # level_set_regression_3$reg_type <- "difference"

  # Regression: main
  level_set_density_1 <- list(
    n = c(100,200,400,800,1600,3200),
    type = c("GCM", "CLS")
  )

  # Regression: dissertation
  level_set_density_2 <- level_set_density_1
  level_set_density_2$type <- NULL

  level_set <- get(cfg$level_set_which)

}



##########################################.
##### MAIN: Setup and run simulation #####
##########################################.

# Set global constants
C <- list() # placeholder

if (cfg$main_task=="run") {

  run_on_cluster(

    first = {

      # Simulation setup
      sim <- new_sim()
      sim %<>% set_config(
        num_sim = cfg$num_sim,
        parallel = cfg$parallel,
        n_cores = cfg$n_cores,
        stop_at_error = cfg$stop_at_error,
        # batch_levels = c("n", "distr_A", "theta_true", "sigma"),
        # return_batch_id = T,
        # seed = 123, # !!!!!
        packages = cfg$pkgs
      )
      sim <- do.call(set_levels, c(list(sim), level_set))
      if (!is.null(cfg$keep)) { sim %<>% set_levels(.keep=cfg$keep) }

      # Simulation script
      sim %<>% set_script(one_simulation)

    },

    main = {
      sim %<>% run()
    },

    last = {

      sim %>% SimEngine::summarize() %>% print()

    },

    cluster_config = cluster_config

  )

} else if (cfg$main_task=="update") {

  update_sim_on_cluster(

    first = {
      sim <- readRDS(paste0(cluster_config$dir,"/sim.rds"))
      sim <- do.call(set_levels, c(list(sim), level_set))
    },

    main = {
      sim %<>% update_sim()
    },

    last = {},

    cluster_config = cluster_config

  )

} else {

  source(cfg$main_task, local=TRUE)

}



####################################.
##### VIZ: Bias, Variance, MSE #####
####################################.

if (F) {

  # Set this manually
  w <- list(
    sim_type = "density",
    print_or_save = "save", # "print" "save"
    theta_true = "identity", # "identity" "constant"
    scaled = T,
    zoomed = F
  )

  # Summarize results
  summ_bias <- summ_var <- summ_sd <- summ_mse <- list()
  for (i in c(1:51)) {
    m <- format(round(i/50-0.02,2), nsmall=2)
    summ_bias[[length(summ_bias)+1]] <- list(
      stat = "bias",
      name = paste0("bias_n__",m),
      estimate = paste0("theta_n_",m),
      truth = paste0("theta_0_",m)
    )
    summ_bias[[length(summ_bias)+1]] <- list(
      stat = "bias",
      name = paste0("bias_s__",m),
      estimate = paste0("theta_s_",m),
      truth = paste0("theta_0_",m)
    )
    summ_var[[length(summ_var)+1]] <- list(
      stat = "var",
      name = paste0("var_n__",m),
      x = paste0("theta_n_",m)
    )
    summ_var[[length(summ_var)+1]] <- list(
      stat = "var",
      name = paste0("var_s__",m),
      x = paste0("theta_s_",m)
    )
    summ_sd[[length(summ_sd)+1]] <- list(
      stat = "sd",
      name = paste0("sd_n__",m),
      x = paste0("theta_n_",m)
    )
    summ_sd[[length(summ_sd)+1]] <- list(
      stat = "sd",
      name = paste0("sd_s__",m),
      x = paste0("theta_s_",m)
    )
  }
  summ <- do.call(SimEngine::summarize,
                  c(list(sim), summ_bias, summ_var, summ_sd))
  if (w$sim_type=="regression") {
    summ %<>% filter(distr_A=="Unif(0,1)")
    # summ %<>% rename("Estimator"=reg_type)
  }
  summ$n_reps <- NULL

  if (w$sim_type=="regression") {
    cols <- c("level_id","n","distr_A","theta_true","sigma")
  } else if (w$sim_type=="density") {
    cols <- c("level_id","n")
    # cols <- c("level_id","n","type")
  }
  p_data <- pivot_longer(
    data = summ,
    cols = -cols,
    names_to = c("stat","point"),
    names_sep = "__"
  )
  p_data %<>% mutate(point=as.numeric(point))

  # Add `Estimator` and `stat` columns
  p_data %<>% mutate(
    Estimator = case_when(
      stat %in% c("bias_n", "var_n", "sd_n") ~ "GCM",
      stat %in% c("bias_s", "var_s", "sd_s") ~ "CLS",
      TRUE ~ "Error"
    ),
    stat = substr(stat,1,1)
  )

  # Add scaled bias/var
  if (w$sim_type=="regression") {

    p_data %<>% mutate(
      value2 = case_when(
        theta_true=="identity" & stat %in% c("b", "s") ~ n^(1/3)*value,
        theta_true=="identity" & stat=="v" ~ n^(2/3)*value,
        theta_true=="constant" & stat %in% c("b", "s") ~ n^(1/2)*value,
        theta_true=="constant" & stat=="v" ~ n*value
      )
    )

  } else {

    p_data %<>% mutate(
      value2 = case_when(
        stat %in% c("b", "s") ~ n^(1/3)*value,
        stat=="v" ~ n^(2/3)*value,
      )
    )

  }

  if (w$sim_type=="regression") {

    # Manipulations specific to `w` plot config object
    w$factor <- ifelse(!w$scaled, "", ifelse(
      w$theta_true=="constant", "n^{1/2}", "n^{1/3}")
    )
    p_data %<>% filter(theta_true==w$theta_true)

  } else {

    w$factor <- ifelse(!w$scaled, "", "n^{1/3}")

  }

  if (w$scaled) { p_data$value <- p_data$value2 }
  p_data$value2 <- NULL

  # This step is necessary because of outliers
  if (w$sim_type=="density") { p_data %<>% filter(point!=0) }

  # Y axis limits
  # Note: not using zooming for now
  if (w$zoomed==T) {
    ylim_b <- filter(p_data, stat=="b" & Estimator=="CLS")$value %>%
      abs() %>% max() %>% (function(x) { c(-2*x, 2*x) })
    ylim_v <- filter(p_data, stat=="v" & Estimator=="CLS")$value %>%
      max() %>% (function(x) { c(0, 2*x) })
    ylim_s <- filter(p_data, stat=="s" & Estimator=="CLS")$value %>%
      max() %>% (function(x) { c(0, 2*x) })
  } else {
    ylim_b <- filter(p_data, stat=="b")$value %>%
      abs() %>% max() %>% (function(x) { c(-1*x, x) })
    ylim_v <- filter(p_data, stat=="v")$value %>%
      max() %>% (function(x) { c(0, x) })
    ylim_s <- filter(p_data, stat=="s")$value %>%
      max() %>% (function(x) { c(0, x) })
  }

  n_levels <- paste("n =", 100*c(1,2,4,8,16,32))

  # Bias plot
  # Export: 10" x 6"
  # Note: change "value" to "value2" for scaled plots
  plot_b <- ggplot(
    filter(p_data, stat=="b"),
    aes(x=point, y=value, color=Estimator, group=Estimator)
  ) +
    geom_line() +
    # facet_grid(rows=dplyr::vars(distr_A), cols=dplyr::vars(n)) +
    facet_wrap(~factor(paste("n =",n), levels=n_levels)) +
    theme(legend.position="bottom") +
    ylim(ylim_b) +
    labs(title = unname(latex2exp::TeX(paste0("$Bias(", w$factor,"\\theta_n)$, ",
                                              w$theta_true, " function"))),
         x = "X",
         y = "Bias")

  # Variance plot
  # Export: 10" x 6"
  plot_v <- ggplot(
    filter(p_data, stat=="v"),
    aes(x=point, y=value, color=Estimator, group=Estimator)
  ) +
    geom_line() +
    # facet_grid(rows=dplyr::vars(distr_A), cols=dplyr::vars(n)) +
    facet_wrap(~factor(paste("n =",n), levels=n_levels)) +
    theme(legend.position="bottom") +
    ylim(ylim_v) +
    labs(title = unname(latex2exp::TeX(paste0("$Var(", w$factor,"\\theta_n)$, ",
                                              w$theta_true, " function"))),
         x = "X",
         y = "Variance")

  # Standard deviation plot
  # Export: 10" x 6"
  plot_s <- ggplot(
    filter(p_data, stat=="s"),
    aes(x=point, y=value, color=Estimator, group=Estimator)
  ) +
    geom_line() +
    # facet_grid(rows=dplyr::vars(distr_A), cols=dplyr::vars(n)) +
    facet_wrap(~factor(paste("n =",n), levels=n_levels)) +
    theme(legend.position="bottom") +
    ylim(ylim_s) +
    labs(title = unname(latex2exp::TeX(paste0("$SD(", w$factor,"\\theta_n)$, ",
                                              w$theta_true, " function"))),
         x = "X",
         y = "Standard deviation")

  if (w$print_or_save=="print") {
    for (p in c("b", "v", "s")) { print(get(paste0("plot_",p))) }
  } else if (w$print_or_save=="save") {
    for (p in c("b", "v", "s")) {
      sc <- ifelse(w$scaled, "scaled", "unscaled")
      wh <- ifelse(p=="b", "Bias", ifelse(p=="v", "Variance", "SD"))
      zm <- ifelse(w$zoomed, "zoomed", "unzoomed")
      ggsave(
        filename = paste0(wh, " (", w$theta_true, ",", sc, ",", zm, ").pdf"),
        plot = get(paste0("plot_",p)), device="pdf", width=10, height=6
      )
    }
  }

}



##############################.
##### VIZ: Distributions #####
##############################.

if (F) {

  # Set this manually
  sim_type <- "density"

  # Code applicable to all plots
  p_data <- sim$results
  if (sim_type=="regression") {
    p_data %<>% filter(theta_true=="identity" & distr_A=="Unif(0,1)")
  }
  n_levels <- paste("n =", 100*c(1,2,4,8,16,32))

  # Create plotting dataset 1
  p_data1 <- p_data %>%
    subset(select=c(n, theta_n_0.30, theta_s_0.30, theta_0_0.30)) %>%
    mutate(
      GCM = n^(1/3)*(theta_n_0.30-theta_0_0.30),
      CLS = n^(1/3)*(theta_s_0.30-theta_0_0.30)
    ) %>%
    subset(select=c(n, GCM, CLS)) %>%
    pivot_longer(
      cols = c(GCM,CLS),
      names_to = c("Estimator")
    )

  # Plot of n^(1/3) (theta_n-theta_0)
  if (sim_type=="regression") {
    sigma <- sim$results$sigma[1]
    tau_0 <- (4*sigma^2)^(1/3)
    grid <- seq(-1,1,0.01)
  } else if (sim_type=="density") {
    point <- 0.3
    # tau_0 <- (4*exp(-2*point))^(1/3)
    dns <- function(x) { dbeta(x=x, shape1=1, shape2=5) }
    dens_0.3 <- dns(0.3)
    deriv_0.3 <- (dns(0.301)-dns(0.299))/(0.301-0.299)
    tau_0 <- (-4*dens_0.3*deriv_0.3)^(1/3)
    grid <- seq(-6,6,0.01)
  }
  ggplot(p_data1, aes(x=value, color=Estimator)) +
    geom_density() +
    geom_area(
      aes(x=x, y=y),
      data = data.frame(
        x = grid,
        y = dChern(grid/tau_0)/tau_0
      ),
      inherit.aes = F,
      alpha = 0.2
    ) +
    facet_wrap(~factor(paste("n =",n), levels=n_levels)) +
    labs(
      title = unname(latex2exp::TeX("$n^{1/3}(\\theta_n-\\theta_0)$")),
      x = "X",
      y = "Estimated density"
    ) +
    theme(legend.position="bottom")

  # Create plotting dataset 2
  p_data2 <- p_data %>%
    subset(select=c(n, Gamma_n_0.30, Gamma_s_0.30)) %>%
    mutate(diff = n^(1/2)*(Gamma_n_0.30-Gamma_s_0.30)) %>%
    subset(select=c(n, diff))

  # Plot of n^(1/2) (Gamma_n-Gamma_n*)
  ggplot(p_data2, aes(x=diff)) +
    geom_density(fill="forestgreen", alpha=0.5, color="white") +
    facet_wrap(~factor(paste("n =",n), levels=n_levels)) +
    labs(
      title = unname(latex2exp::TeX(
        "$n^{1/2}(\\Gamma_n(0.3)-\\Gamma_n^*(0.3))$"
      )),
      x = NULL,
      y = "Estimated density"
    )

  # !!!!! Marco new plot request (GCM vs Gamma_n)
  if (F) {

    p_data6 <- p_data %>%
      subset(select=c(n, Gamma_n_0.30, GCM_n_0.30)) %>%
      mutate(diff = n^(1/2)*(Gamma_n_0.30-GCM_n_0.30)) %>%
      subset(select=c(n, diff))

    # Plot of n^(1/2) (GCM_n-Gamma_n)
    ggplot(p_data6, aes(x=diff)) +
      geom_density(fill="forestgreen", alpha=0.5, color="white") +
      facet_wrap(~factor(paste("n =",n), levels=n_levels)) +
      labs(
        title = unname(latex2exp::TeX(
          "$n^{1/2}(\\Gamma_n(0.3)-\\GCM_n(0.3))$"
        )),
        x = NULL,
        y = "Estimated density"
      )


  }

  # Create plotting dataset 3
  p_data3 <- p_data %>%
    subset(select=c(n, theta_n_0.30, theta_s_0.30)) %>%
    mutate(diff = n^(1/3)*(theta_n_0.30-theta_s_0.30)) %>%
    subset(select=c(n, diff))

  # Plot of n^(1/3) (theta_n-theta_s)
  ggplot(p_data3, aes(x=diff)) +
    geom_density(fill="forestgreen", alpha=0.5, color="white") +
    facet_wrap(~factor(paste("n =",n), levels=n_levels)) +
    labs(
      title = unname(latex2exp::TeX("$n^{1/3}(\\theta_n-\\theta_n^*)$")),
      x = NULL,
      y = "Estimated density"
    )

  # Create plotting dataset 4
  p_data4 <- sim$results
  if (sim_type=="regression") {
    p_data4 %<>% filter(theta_true=="identity" & distr_A=="Unif(0,1)")
  }
  p_data4 %<>%
    filter(rep_id<=1000) %>%
    subset(select=c(n, theta_n_0.30, theta_s_0.30)) %>%
    rename(
      "GCM" = theta_n_0.30,
      "CLS" = theta_s_0.30
    )

  # Plot of n^(1/3) (theta_n-theta_0)
  if (sim_type=="regression") {
    pt_center <- 0.3
  } else {
    pt_center <- dbeta(x=0.3, shape1=1, shape2=5)
  }
  ggplot(p_data4, aes(x=GCM, y=CLS)) +
    geom_abline(slope=1, intercept=0, color="orange") +
    geom_point(alpha=0.1, pch=16) +
    geom_point(
      data = data.frame(CLS=pt_center, GCM=pt_center),
      color = "maroon",
      pch = 8,
      size = 2
    ) +
    facet_wrap(~factor(paste("n =",n), levels=n_levels)) +
    labs(
      title = unname(latex2exp::TeX(
        "$\\theta_n$ (GCM) vs. $\\theta_n^* (CLS)$"
      )),
      x = "GCM",
      y = "CLS"
    )

  # Create plotting dataset 5
  # NOT INCLUDED IN DISSERTATION
  p_data5 <- p_data %>%
    subset(select=c(n, GCM_n_0.30, GCM_s_0.30)) %>%
    mutate(diff = n^(1/2)*(GCM_n_0.30-GCM_s_0.30)) %>%
    subset(select=c(n, diff))

  # Plot of n^(1/2) (GCM(Gamma_n)-GCM(Gamma_n*))
  # NOT INCLUDED IN DISSERTATION
  ggplot(p_data5, aes(x=diff)) +
    geom_density(fill="forestgreen", alpha=0.5, color="white") +
    facet_wrap(~factor(paste("n =",n), levels=n_levels)) +
    labs(
      title = unname(latex2exp::TeX(
        "$n^{1/2}(GCM(\\Gamma_n)(0.3)-GCM(\\Gamma_n)^*(0.3))$"
      )),
      x = NULL,
      y = "Estimated density"
    )

}



#####################################.
##### VIZ: Concept illlustraion #####
#####################################.

if (F) {

  # Generate dataset
  set.seed(180)
  dat_obj <- generate_data(n=20, distr_A="Unif(0,1)",
                           theta_true="identity", sigma=1) # !!!!! identity
  dat <- dat_obj$dat
  theta_0 <- dat_obj$theta_0

  # Estimate regression function
  res_gcm <- est_curve(dat, "Iso GCM2", return_Gamma_n=T, return_cusum=T,
                       return_GCM=T)
  res_ls <- est_curve(dat, "Iso CLS", return_CLS=T)
  theta_gcm <- res_gcm$theta_n
  theta_ls <- res_ls$theta_n
  Gamma_gcm <- res_gcm$GCM
  Gamma_ls <- res_ls$CLS
  cusum <- res_gcm$cusum

  # Create plot: theta_n
  which_plot <- 1 # c(1,2,3)
  grid <- seq(0,1,0.01)
  df_plot1 <- data.frame(
    x = rep(grid, 3),
    y = c(theta_gcm(grid), theta_ls(grid), theta_0(grid)),
    which = rep(c("theta_n (GCM)", "theta_n (CLS)", "theta_0"), each=length(grid))
  )
  if (which_plot==1) {
    df_plot1 %<>% filter(which=="theta_0")
  } else if (which_plot==2) {
    df_plot1 %<>% filter(which!="theta_n (CLS)")
  }
  ggplot(df_plot1, aes(x=x, y=y, color=which)) +
    geom_line() +
    geom_point(
      aes(x=x, y=y),
      data.frame(x=dat$a, y=dat$y),
      inherit.aes=F,
      alpha=0.5
    ) +
    scale_color_manual(values=c("theta_n (CLS)"="salmon",
                                "theta_n (GCM)"="deepskyblue",
                                "theta_0"="darkolivegreen4")) +
    labs(color="Function", x="X", y="Y") +
    theme(legend.position="bottom")

  # Create plot: Gamma_n
  which_plot <- 2 # c(1,2)
  df_plot2 <- data.frame(
    x = rep(grid, 2),
    y = c(Gamma_gcm(grid), Gamma_ls(grid)),
    which = rep(c("GCM", "CLS"), each=length(grid))
  )
  if (which_plot==1) { df_plot2 %<>% filter(which=="GCM") }
  ggplot(df_plot2, aes(x=x, y=y, color=which)) +
    geom_line() +
    geom_point(
      aes(x=x, y=y),
      data.frame(x=cusum$x, y=cusum$y),
      inherit.aes=F,
      alpha=0.5
    ) +
    scale_color_manual(values=c("CLS"="salmon",
                                "GCM"="deepskyblue")) +
    labs(color="Function", x="F_n(X)", y="Gamma_n(X)") + # title="Primitive"
    theme(legend.position="bottom")

}


