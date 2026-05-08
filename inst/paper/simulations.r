# ============================================================
# Monte Carlo study for the Contracted-Beta (CB) distribution
# ============================================================

library(CB)
library(dplyr)
library(tidyr)
library(patchwork)

# ------------------------------------------------------------
# One Monte Carlo replicate
# ------------------------------------------------------------

mc_one_cb <- function(n, alpha_true, beta_true) {
  
  x <- rCB(n = n, alpha = alpha_true, beta = beta_true)
  
  # ----------------------------------------------------------
  # Moment estimator
  # ----------------------------------------------------------
  
  fit_mom <- try(fitCB_mom(x), silent = TRUE)
  
  if (!inherits(fit_mom, "try-error") &&
      isTRUE(fit_mom$success) &&
      !is.null(fit_mom$par) &&
      all(is.finite(fit_mom$par))) {
    
    alpha_mom <- unname(fit_mom$par["alpha"])
    beta_mom  <- unname(fit_mom$par["beta"])
    success_mom <- TRUE
    
  } else {
    
    alpha_mom <- NA_real_
    beta_mom  <- NA_real_
    success_mom <- FALSE
  }
  
  # ----------------------------------------------------------
  # Maximum likelihood estimator
  # ----------------------------------------------------------
  
  fit_mle <- try(
    fitCB_mle(
      x,
      methods = "L-BFGS-B",
      multistart = FALSE
    ),
    silent = TRUE
  )
  
  if (inherits(fit_mle, "try-error") ||
      is.null(fit_mle$par) ||
      any(!is.finite(fit_mle$par)) ||
      fit_mle$convergence != 0) {
    
    return(data.frame(
      alpha_true = alpha_true,
      beta_true = beta_true,
      n = n,
      alpha_mle = NA_real_,
      beta_mle = NA_real_,
      alpha_mom = alpha_mom,
      beta_mom = beta_mom,
      se_alpha = NA_real_,
      se_beta = NA_real_,
      cp_alpha = NA_real_,
      cp_beta = NA_real_,
      success_mle = FALSE,
      success_mom = success_mom
    ))
  }
  
  alpha_mle <- unname(fit_mle$par["alpha"])
  beta_mle  <- unname(fit_mle$par["beta"])
  
  se_alpha <- unname(fit_mle$se["alpha"])
  se_beta  <- unname(fit_mle$se["beta"])
  
  cp_alpha <- NA_real_
  cp_beta  <- NA_real_
  
  if (is.finite(se_alpha) && se_alpha > 0) {
    li_alpha <- alpha_mle - 1.96 * se_alpha
    ls_alpha <- alpha_mle + 1.96 * se_alpha
    cp_alpha <- as.numeric(alpha_true >= li_alpha && alpha_true <= ls_alpha)
  }
  
  if (is.finite(se_beta) && se_beta > 0) {
    li_beta <- beta_mle - 1.96 * se_beta
    ls_beta <- beta_mle + 1.96 * se_beta
    cp_beta <- as.numeric(beta_true >= li_beta && beta_true <= ls_beta)
  }
  
  data.frame(
    alpha_true = alpha_true,
    beta_true = beta_true,
    n = n,
    alpha_mle = alpha_mle,
    beta_mle = beta_mle,
    alpha_mom = alpha_mom,
    beta_mom = beta_mom,
    se_alpha = se_alpha,
    se_beta = se_beta,
    cp_alpha = cp_alpha,
    cp_beta = cp_beta,
    success_mle = TRUE,
    success_mom = success_mom
  )
}

# ------------------------------------------------------------
# Full Monte Carlo study
# ------------------------------------------------------------

mc_CB_study <- function(R = 1000,
                        n_values = c(50, 100, 200, 300, 500),
                        alpha_values = c(0.5, 2, 4),
                        beta_values = c(1, 3, 6),
                        seed = 123,
                        verbose = TRUE) {
  
  set.seed(seed)
  
  scenarios <- expand.grid(
    alpha_true = alpha_values,
    beta_true = beta_values,
    n = n_values
  )
  
  results <- vector("list", nrow(scenarios) * R)
  id <- 1L
  
  for (i in seq_len(nrow(scenarios))) {
    
    alpha_true <- scenarios$alpha_true[i]
    beta_true  <- scenarios$beta_true[i]
    n          <- scenarios$n[i]
    
    if (isTRUE(verbose)) {
      message(
        "Scenario ", i, "/", nrow(scenarios),
        ": alpha = ", alpha_true,
        ", beta = ", beta_true,
        ", n = ", n
      )
    }
    
    for (r in seq_len(R)) {
      
      aux <- mc_one_cb(
        n = n,
        alpha_true = alpha_true,
        beta_true = beta_true
      )
      
      aux$replicate <- r
      results[[id]] <- aux
      id <- id + 1L
    }
  }
  
  bind_rows(results) %>%
    select(
      replicate,
      alpha_true,
      beta_true,
      n,
      alpha_mle,
      beta_mle,
      alpha_mom,
      beta_mom,
      se_alpha,
      se_beta,
      cp_alpha,
      cp_beta,
      success_mle,
      success_mom
    )
}

# ------------------------------------------------------------
# Summary table:
# Bias and RMSE for MLE and MoM; CP only for MLE
# ------------------------------------------------------------

summarize_CB_mc <- function(mc_results) {
  
  mc_results %>%
    group_by(alpha_true, beta_true, n) %>%
    summarise(
      Bias_alpha_MLE = mean(alpha_mle - alpha_true, na.rm = TRUE),
      Bias_beta_MLE  = mean(beta_mle  - beta_true,  na.rm = TRUE),
      RMSE_alpha_MLE = sqrt(mean((alpha_mle - alpha_true)^2, na.rm = TRUE)),
      RMSE_beta_MLE  = sqrt(mean((beta_mle  - beta_true)^2,  na.rm = TRUE)),
      
      Bias_alpha_MoM = mean(alpha_mom - alpha_true, na.rm = TRUE),
      Bias_beta_MoM  = mean(beta_mom  - beta_true,  na.rm = TRUE),
      RMSE_alpha_MoM = sqrt(mean((alpha_mom - alpha_true)^2, na.rm = TRUE)),
      RMSE_beta_MoM  = sqrt(mean((beta_mom  - beta_true)^2,  na.rm = TRUE)),
      
      CP_alpha_MLE = mean(cp_alpha, na.rm = TRUE),
      CP_beta_MLE  = mean(cp_beta,  na.rm = TRUE),
      
      SuccessRate_MLE = mean(success_mle, na.rm = TRUE),
      SuccessRate_MoM = mean(success_mom, na.rm = TRUE),
      .groups = "drop"
    )
}

# ------------------------------------------------------------
# Function to assemble 1 x 3 panels
# ------------------------------------------------------------

build_simulation_panel <- function(mc_results,
                                   alpha_val,
                                   beta_vals,
                                   scenario_labels) {
  
  plots <- vector("list", length(beta_vals))
  
  for (i in seq_along(beta_vals)) {
    
    beta_val <- beta_vals[i]
    scen_lab <- scenario_labels[i]
    
    plots[[i]] <- plot_box_cb_free_y(
      mc_results = mc_results,
      alpha_val = alpha_val,
      beta_val = beta_val
    ) +
      ggtitle(
        bquote(
          .(scen_lab) * ": " ~
            alpha == .(alpha_val) * "," ~
            beta == .(beta_val)
        )
      )
  }
  
  wrap_plots(plots, nrow = 1) &
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        size = 12
      )
    )
}

# ------------------------------------------------------------
# Run the Monte Carlo study
# ------------------------------------------------------------

mc_res <- mc_CB_study(
  R = 1000,
  n_values = c(50, 100, 200, 300, 500),
  alpha_values = c(0.5, 2, 4),
  beta_values = c(1, 3, 6),
  seed = 123,
  verbose = TRUE
)

# ------------------------------------------------------------
# Monte Carlo summaries
# ------------------------------------------------------------

mc_summary_full <- summarize_CB_mc(mc_res)

mc_summary_alpha05 <- mc_summary_full %>%
  filter(alpha_true == 0.5) %>%
  mutate(across(where(is.numeric), ~ round(.x, 5)))

mc_summary_alpha2 <- mc_summary_full %>%
  filter(alpha_true == 2) %>%
  mutate(across(where(is.numeric), ~ round(.x, 5)))

mc_summary_alpha4 <- mc_summary_full %>%
  filter(alpha_true == 4) %>%
  mutate(across(where(is.numeric), ~ round(.x, 5)))

as.data.frame(mc_summary_alpha05)

as.data.frame(mc_summary_alpha2)

as.data.frame(mc_summary_alpha4)


