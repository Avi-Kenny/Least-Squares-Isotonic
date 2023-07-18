# Title: "Alternative least squares isotonic regression estimator"
# Author: Avi Kenny

##################.
##### CONFIG #####
##################.

# To run multiple sims/analyses concurrently, ONLY change Slurm commands

# Set global config
cfg <- list(
  main_task = "run", # run update
  which_sim = "dissertation", # regression density dissertation
  level_set_which = "level_set_regression_2", # level_set_regression_1 level_set_density_1
  num_sim = 3,
  pkgs = c("dplyr", "Iso", "fdrtool", "simest", "tidyr", "truncnorm",
           "modeest"),
  pkgs_nocluster = c("ggplot2"),
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
    n = c(100,200,400,800,1600), # ,3200
    distr_A = c("Unif(0,1)"), # "Unif(0,1)", "N(0.5,0.09)"
    theta_true = c("identity"), # "identity", "square", "constant"
    reg_type = c("Iso CLS", "Iso GCM"), # "Linear", "Iso GCM", "Iso GCM2", "Iso CLS"
    sigma = 0.2
  )

  # Regression: dissertation
  level_set_regression_2 <- level_set_regression_1
  level_set_regression_2$n <- c(100,200) # !!!!! TEMP
  level_set_regression_2$reg_type <- NULL

  # # Regression: difference between estimators
  # level_set_regression_3 <- level_set_regression_1
  # level_set_regression_3$reg_type <- "difference"

  # Regression: main
  level_set_density_1 <- list(
    n = c(100,200,400,800), # 1600,3200
    type = c("GCM", "CLS")
    # distr_A = c("Unif(0,1)", "N(0.5,0.09)"),
    # theta_true = c("identity", "constant"), # "square"
    # reg_type = c("Iso CLS", "Iso GCM2"), # "Linear", "Iso GCM", "Iso GCM2", "Iso CLS"
    # sigma = 0.2
  )

  level_set <- get(cfg$level_set_which)

}



##########################################.
##### MAIN: Setup and run simulation #####
##########################################.

# Set global constants
C <- list() # placeholde

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
        batch_levels = c("n", "distr_A", "theta_true", "sigma"),
        return_batch_id = T,
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
    which_sim = cfg$which_sim,
    print_or_save = "print", # "print" "save"
    theta_true = "identity", # "identity" "constant"
    scaled = T,
    zoomed = T
  )

  # Summarize results
  summ_bias <- summ_var <- summ_sd <- summ_mse <- list()
  for (i in c(1:51)) {
    m <- format(round(i/50-0.02,2), nsmall=2)
    summ_bias[[i]] <- list(
      stat = "bias",
      name = paste0("bias_",m),
      estimate = paste0("est_",m),
      truth = paste0("theta_",m)
    )
    summ_var[[i]] <- list(
      stat = "var",
      name = paste0("var_",m),
      x = paste0("est_",m)
    )
    summ_mse[[i]] <- list(
      stat = "mse",
      name = paste0("mse_",m),
      estimate = paste0("est_",m),
      truth = paste0("theta_",m)
    )
  }
  summ <- do.call(SimEngine::summarize,
                  c(list(sim), summ_bias, summ_var, summ_mse))
  if (w$which_sim=="regression") { summ %<>% rename("Estimator"=reg_type) }
  summ$n_reps <- NULL

  if (w$which_sim=="regression") {
    cols <- c("level_id","n","distr_A","theta_true","Estimator","sigma")
  } else if (w$which_sim=="density") {
    cols <- c("level_id","n","type")
  }
  p_data <- pivot_longer(
    data = summ,
    cols = -cols,
    names_to = c("stat","point"),
    names_sep = "_"
  )
  p_data %<>% mutate(point=as.numeric(point))

  if (w$which_sim=="regression") {

    # Add scaled bias/var/MSE
    p_data %<>% mutate(
      value2 = case_when(
        theta_true=="identity" & stat %in% c("bias") ~ n^(1/3)*value,
        theta_true=="identity" & stat %in% c("var","mse") ~ n^(2/3)*value,
        theta_true=="constant" & stat %in% c("bias") ~ n^(1/2)*value,
        theta_true=="constant" & stat %in% c("var","mse") ~ n^(1/1)*value
      )
    )

    # Manipulations specific to `w` plot config object
    w$factor <- ifelse(!w$scaled, "", ifelse(
      w$theta_true=="constant", "n^{1/2}", "n^{1/3}")
    )
    p_data %<>% filter(theta_true==w$theta_true)

    if (w$scaled) { p_data$value <- p_data$value2 }

    # Y axis limits
    if (w$zoomed==T) {
      ylim_b <- filter(p_data, stat=="bias" & Estimator=="Rearrangement")$value %>%
        abs() %>% max() %>% (function(x) { c(-1.2*x, 1.2*x) })
      # ylim_b <- filter(p_data, stat=="bias" & Estimator=="Iso CLS")$value %>%
      #   abs() %>% max() %>% (function(x) { c(-2*x, 2*x) })
      ylim_v <- filter(p_data, stat=="var" & Estimator=="Iso CLS")$value %>%
        max() %>% (function(x) { c(0, 2*x) })
      ylim_m <- filter(p_data, stat=="mse" & Estimator=="Iso CLS")$value %>%
        max() %>% (function(x) { c(0, 2*x) })
    } else {
      ylim_b <- filter(p_data, stat=="bias")$value %>%
        abs() %>% max() %>% (function(x) { c(-1*x, x) })
      ylim_v <- filter(p_data, stat=="var")$value %>%
        max() %>% (function(x) { c(0, x) })
      ylim_m <- filter(p_data, stat=="mse")$value %>%
        max() %>% (function(x) { c(0, x) })
    }

  }

  # !!!!!
  ggplot(
    filter(p_data, stat=="bias"),
    aes(x=point, y=value, color=type, group=type)
  ) +
    geom_line() +
    facet_grid(cols=dplyr::vars(n)) +
    theme(legend.position="bottom") +
    labs(
      # title = unname(latex2exp::TeX(paste0("$Bias(", w$factor,"\\theta_n)$, ",
                                              # w$theta_true, " function"))),
         # x = "A",
         y = "Bias")


  # Bias plot
  # Export: 10" x 6"
  # Note: change "value" to "value2" for scaled plots
  plot_b <- ggplot(
    filter(p_data, stat=="bias"),
    aes(x=point, y=value, color=Estimator, group=Estimator)
  ) +
    geom_line() +
    facet_grid(rows=dplyr::vars(distr_A), cols=dplyr::vars(n)) +
    theme(legend.position="bottom") +
    ylim(ylim_b) +
    labs(title = unname(latex2exp::TeX(paste0("$Bias(", w$factor,"\\theta_n)$, ",
                                              w$theta_true, " function"))),
         x = "A",
         y = "Bias")

  # Variance plot
  # Export: 10" x 6"
  # Note: change "value" to "value2" for scaled plots
  plot_v <- ggplot(
    filter(p_data, stat=="var"),
    aes(x=point, y=value, color=Estimator, group=Estimator)
  ) +
    geom_line() +
    facet_grid(rows=dplyr::vars(distr_A), cols=dplyr::vars(n)) +
    theme(legend.position="bottom") +
    ylim(ylim_v) +
    labs(title = unname(latex2exp::TeX(paste0("$Var(", w$factor,"\\theta_n)$, ",
                                              w$theta_true, " function"))),
         x = "A",
         y = "Variance")

  # MSE plot
  # Export: 10" x 6"
  # Note: change "value" to "value2" for scaled plots
  plot_m <- ggplot(
    filter(p_data, stat=="mse"),
    aes(x=point, y=value, color=Estimator, group=Estimator)
  ) +
    geom_line() +
    facet_grid(rows=dplyr::vars(distr_A), cols=dplyr::vars(n)) +
    theme(legend.position="bottom") +
    ylim(ylim_m) +
    labs(title = unname(latex2exp::TeX(paste0("$MSE(", w$factor,"\\theta_n)$, ",
                                              w$theta_true, " function"))),
         x = "A",
         y = "MSE")

  if (w$print_or_save=="print") {
    for (p in c("b", "v", "m")) { print(get(paste0("plot_",p))) }
  } else if (w$print_or_save=="save") {
    for (p in c("b", "v", "m")) {
      sc <- ifelse(w$scaled, "scaled", "unscaled")
      wh <- ifelse(p=="b", "Bias", ifelse(p=="v", "Variance", "MSE"))
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

  # Manipulate data
  p_data <- pivot_longer(
    data = filter(sim$results, theta_true=="identity" & distr_A=="Unif(0,1)"),
    cols = -c(sim_uid,level_id,sim_id,n,distr_A,theta_true,reg_type,sigma,
              runtime),
    names_to = c("stat","point"),
    names_sep = "_"
  )
  p_data %<>% subset(select=-c(sim_uid,level_id,sim_id,sigma,
                               runtime,theta_true))
  p_data %<>% mutate(
    point = as.numeric(point),
    # value = n^(1/3)*(value-point)
    value = n^(1/3)*value # "difference"
  )
  p_data %<>% filter(stat=="est" & point %in% c(0.1,0.5))

  # # Create reference distribution
  # sigma <- sim$results$sigma[1]
  # sd <- (4*sigma^2)^(1/3) * 0.52

  # # Plot of n^(1/3) (theta_n-theta_0)
  # ggplot(p_data, aes(x=value, color=reg_type)) +
  #   labs(
  #     title = unname(latex2exp::TeX("$n^{1/3}(\\theta_n-\\theta_0)$")),
  #     color = "Estimator",
  #     x = NULL,
  #     y = NULL
  #   ) +
  #   geom_area(
  #     aes(x=x, y=y),
  #     data = data.frame(x=seq(-1,1,0.1), y=dnorm(seq(-1,1,0.1), sd=sd)),
  #     inherit.aes = F,
  #     alpha = 0.2
  #   ) +
  #   geom_density() +
  #   facet_grid(rows=dplyr::vars(point), cols=dplyr::vars(n)) +
  #   xlim(-1,1)

  # Plot of n^(1/3) (theta1_n-theta2_n)
  ggplot(p_data, aes(x=value, color=reg_type)) +
    labs(
      title = unname(
        latex2exp::TeX("$n^{1/3}(\\theta_n^{CLS}-\\theta_n^{GCM})$")
      ),
      color = "Estimator",
      x = NULL,
      y = NULL) +
    geom_density() +
    facet_grid(rows=dplyr::vars(point), cols=dplyr::vars(n))
    # xlim(-1,1)

}



#####################################.
##### VIZ: Concept illlustraion #####
#####################################.

if (F) {

  # Generate dataset
  dat_obj <- generate_data(n=25, distr_A="Unif(0,1)",
                           theta_true="constant", sigma=0.5)
  dat <- dat_obj$dat
  theta_0 <- dat_obj$theta_0

  # Estimate regression function
  res_gcm <- est_curve(dat, "Iso GCM2", return_Theta_n=T, return_cusum=T)
  res_ls <- est_curve(dat, "Iso CLS", return_Theta_n=T, return_cusum=T)
  theta_gcm <- res_gcm$theta_n
  theta_ls <- res_ls$theta_n
  Theta_gcm <- res_gcm$Theta_n
  Theta_ls <- res_ls$Theta_n
  cusum <- res_gcm$cusum

  # Create plot: theta_n
  grid <- seq(0,1,0.01)
  df_plot1 <- data.frame(
    x = rep(grid, 3),
    y = c(theta_gcm(grid), theta_ls(grid), theta_0(grid)),
    which = rep(c("theta_n (GCM)", "theta_n (CLS)", "theta_0"), each=length(grid))
  )
  ggplot(df_plot1, aes(x=x, y=y, color=which)) +
    geom_line() +
    geom_point(
      aes(x=x, y=y),
      data.frame(x=dat$a, y=dat$y),
      inherit.aes=F,
      alpha=0.5
    ) +
    scale_color_manual(values=c("theta_n (CLS)"="deepskyblue",
                                "theta_n (GCM)"="salmon",
                                "theta_0"="darkolivegreen4")) +
    labs(color="Function", x="X", y="Y") +
    theme(legend.position="bottom")

  # Create plot: Theta_n
  df_plot2 <- data.frame(
    x = rep(grid, 2),
    y = c(Theta_gcm(grid), Theta_ls(grid)),
    which = rep(c("GCM", "CLS"), each=length(grid))
  )
  ggplot(df_plot2, aes(x=x, y=y, color=which)) +
    geom_line() +
    geom_point(
      aes(x=x, y=y),
      data.frame(x=cusum$x, y=cusum$y),
      inherit.aes=F,
      alpha=0.5
    ) +
    scale_color_manual(values=c("CLS"="deepskyblue",
                                "GCM"="salmon")) +
    labs(color="Function", x="i/n", y="Gamma_n(i)") + # title="Primitive"
    theme(legend.position="bottom")

}












