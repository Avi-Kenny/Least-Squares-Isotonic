# Title: "Alternative least squares isotonic regression estimator"
# Author: Avi Kenny

##################.
##### CONFIG #####
##################.

# To run multiple sims/analyses concurrently, ONLY change Slurm commands

# Set global config
cfg <- list(
  main_task = "run", # run update analysis.R
  which_sim = "estimation", # "estimation" "edge" "testing" "Cox"
  level_set_which = "level_set_estimation_1", # level_set_estimation_1 level_set_testing_1 level_set_Cox_1
  # keep = c(1:3,7:9,16:18,22:24),
  num_sim = 1000,
  pkgs = c("dplyr", "Iso", "fdrtool", "simest", "tidyr", "truncnorm",
           "modeest"),
  pkgs_nocluster = c("ggplot2"),
  parallel = "none",
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



##########################################.
##### MAIN: Setup and run simulation #####
##########################################.

if (cfg$main_task=="run") {

  run_on_cluster(

    first = {

      # Simulation setup
      sim <- new_sim()
      sim %<>% set_config(
        num_sim = cfg$num_sim,
        parallel = cfg$parallel,
        stop_at_error = cfg$stop_at_error,
        packages = cfg$pkgs
      )
      sim %<>% set_levels(
        n = c(100,200,400,800,1600,3200),
        distr_A = c("Unif(0,1)", "N(0.5,0.09)"),
        theta_true = c("identity", "constant"), # "constant", "identity", "square"
        reg_type = c("Iso LS", "Iso GCM2"),
        # reg_type = c("Linear", "Iso GCM", "Iso GCM2", "Iso LS", "difference"),
        sigma = 0.2
      )
      if (!is.null(cfg$keep)) { sim %<>% set_levels(.keep=cfg$keep) }

      # Simulation script
      sim %<>% set_script(one_simulation)

    },

    main = {
      sim %<>% run()
    },

    last = {

      sim %>% summarize() %>% print()

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
    print_or_save = "save", # "print" "save"
    theta_true = "constant", # "identity" "constant"
    scaled = T,
    zoomed = T
  )

  # Summarize results
  summ_bias <- summ_var <- summ_sd <- summ_mse <- list()
  for (i in c(1:51)) {
    m <- format(round(i/50-0.02,2), nsmall=2)
    summ_bias[[i]] <- list(
      name = paste0("bias_",m),
      estimate = paste0("est_",m),
      truth = paste0("theta_",m)
    )
    summ_var[[i]] <- list(
      name = paste0("var_",m),
      x = paste0("est_",m)
    )
    summ_mse[[i]] <- list(
      name = paste0("mse_",m),
      estimate = paste0("est_",m),
      truth = paste0("theta_",m)
    )
  }
  summ <- summarize(sim, bias=summ_bias, var=summ_var, mse=summ_mse)
  summ %<>% rename("Estimator"=reg_type)

  p_data <- pivot_longer(
    data = summ,
    cols = -c(level_id,n,distr_A,theta_true,Estimator,sigma),
    names_to = c("stat","point"),
    names_sep = "_"
  )
  p_data %<>% mutate(point=as.numeric(point))

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
  if (w$scaled) {
    p_data$value <- p_data$value2
  }

  # Y axis limits
  if (w$zoomed==T) {
    ylim_b <- filter(p_data, stat=="bias" & Estimator=="Iso LS")$value %>%
      abs() %>% max() %>% (function(x) { c(-2*x, 2*x) })
    ylim_v <- filter(p_data, stat=="var" & Estimator=="Iso LS")$value %>%
      max() %>% (function(x) { c(0, 2*x) })
    ylim_m <- filter(p_data, stat=="mse" & Estimator=="Iso LS")$value %>%
      max() %>% (function(x) { c(0, 2*x) })
  } else {
    ylim_b <- filter(p_data, stat=="bias")$value %>%
      abs() %>% max() %>% (function(x) { c(-1*x, x) })
    ylim_v <- filter(p_data, stat=="var")$value %>%
      max() %>% (function(x) { c(0, x) })
    ylim_m <- filter(p_data, stat=="mse")$value %>%
      max() %>% (function(x) { c(0, x) })
  }

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
    data = sim$results,
    cols = -c(sim_uid,level_id,sim_id,n,distr_A,theta_true,reg_type,sigma,
              runtime),
    names_to = c("stat","point"),
    names_sep = "_"
  )
  p_data %<>% subset(select=-c(sim_uid,level_id,sim_id,sigma,
                               runtime,theta_true))
  p_data %<>% filter(stat=="est" & point %in% c(0.1,0.5))
  p_data %<>% mutate(
    point = as.numeric(point),
    value = n^(1/3)*(value-point)
    # value = n^(1/3)*value # "difference"
  )

  # Create reference distribution
  sigma <- sim$results$sigma[1]
  sd <- (4*sigma^2)^(1/3) * 0.52

  # Plot of n^(1/3) (theta_n-theta_0)
  ggplot(p_data, aes(x=value, color=reg_type)) +
    labs(
      title = unname(latex2exp::TeX("$n^{1/3}(\\theta_n-\\theta_0)$")),
      color = "Estimator",
      x = NULL,
      y = NULL
    ) +
    geom_area(
      aes(x=x, y=y),
      data = data.frame(x=seq(-1,1,0.1), y=dnorm(seq(-1,1,0.1), sd=sd)),
      inherit.aes = F,
      alpha = 0.2
    ) +
    geom_density() +
    facet_grid(rows=dplyr::vars(point), cols=dplyr::vars(n)) +
    xlim(-1,1)

  # Plot of n^(1/3) (theta1_n-theta2_n)
  ggplot(p_data, aes(x=value, color=reg_type)) +
    labs(
      title = unname(
        latex2exp::TeX("$n^{1/3}(\\theta_n^{LS}-\\theta_n^{GCM})$")
      ),
      color = "Estimator",
      x = NULL,
      y = NULL
    ) +
    geom_density() +
    facet_grid(rows=dplyr::vars(point), cols=dplyr::vars(n))
    # xlim(-1,1)

}
