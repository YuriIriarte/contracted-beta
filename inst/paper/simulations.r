# ============================================================
# Monte Carlo study for the canonical CB distribution
# ============================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

# ------------------------------------------------------------
# One Monte Carlo replicate
# ------------------------------------------------------------

mc_one_cb <- function(n, alpha_true, beta_true) {
  
  x <- rCB(n = n, alpha = alpha_true, beta = beta_true)
  
  fit <- try(
    fitCB_mle(
      x,
      methods = "L-BFGS-B",
      multistart = FALSE
    ),
    silent = TRUE
  )
  
  if (inherits(fit, "try-error") ||
      is.null(fit$par) ||
      any(!is.finite(fit$par)) ||
      fit$convergence != 0) {
    
    return(data.frame(
      alpha_true = alpha_true,
      beta_true = beta_true,
      n = n,
      alpha_hat = NA_real_,
      beta_hat = NA_real_,
      se_alpha = NA_real_,
      se_beta = NA_real_,
      cp_alpha = NA_real_,
      cp_beta = NA_real_,
      success = FALSE
    ))
  }
  
  alpha_hat <- unname(fit$par["alpha"])
  beta_hat  <- unname(fit$par["beta"])
  
  se_alpha <- unname(fit$se["alpha"])
  se_beta  <- unname(fit$se["beta"])
  
  cp_alpha <- NA_real_
  cp_beta  <- NA_real_
  
  if (is.finite(se_alpha) && se_alpha > 0) {
    li_alpha <- alpha_hat - 1.96 * se_alpha
    ls_alpha <- alpha_hat + 1.96 * se_alpha
    cp_alpha <- as.numeric(alpha_true >= li_alpha && alpha_true <= ls_alpha)
  }
  
  if (is.finite(se_beta) && se_beta > 0) {
    li_beta <- beta_hat - 1.96 * se_beta
    ls_beta <- beta_hat + 1.96 * se_beta
    cp_beta <- as.numeric(beta_true >= li_beta && beta_true <= ls_beta)
  }
  
  data.frame(
    alpha_true = alpha_true,
    beta_true = beta_true,
    n = n,
    alpha_hat = alpha_hat,
    beta_hat = beta_hat,
    se_alpha = se_alpha,
    se_beta = se_beta,
    cp_alpha = cp_alpha,
    cp_beta = cp_beta,
    success = TRUE
  )
}

# ------------------------------------------------------------
# Full Monte Carlo study
# ------------------------------------------------------------

mc_CB_study <- function(R = 1000,
                        n_values = c(50, 100, 200, 300, 500),
                        alpha_values = c(1, 4, 8),
                        beta_values = c(1, 4, 8),
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
    
    if (verbose) {
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
      alpha_hat,
      beta_hat,
      se_alpha,
      se_beta,
      cp_alpha,
      cp_beta,
      success
    )
}

# ------------------------------------------------------------
# Summary table: mean, bias, RMSE, coverage probability
# ------------------------------------------------------------

summarize_CB_mc <- function(mc_results) {
  
  mc_results %>%
    group_by(alpha_true, beta_true, n) %>%
    summarise(
      mean_alpha = mean(alpha_hat, na.rm = TRUE),
      mean_beta  = mean(beta_hat, na.rm = TRUE),
      
      Bias_alpha = mean(alpha_hat - alpha_true, na.rm = TRUE),
      Bias_beta  = mean(beta_hat - beta_true, na.rm = TRUE),
      
      RMSE_alpha = sqrt(mean((alpha_hat - alpha_true)^2, na.rm = TRUE)),
      RMSE_beta  = sqrt(mean((beta_hat - beta_true)^2, na.rm = TRUE)),
      
      CP_alpha = mean(cp_alpha, na.rm = TRUE),
      CP_beta  = mean(cp_beta, na.rm = TRUE),
      
      SuccessRate = mean(success, na.rm = TRUE),
      .groups = "drop"
    )
}

# ------------------------------------------------------------
# Boxplots of the parameter estimates
# ------------------------------------------------------------

plot_box_cb_free_y <- function(mc_results, alpha_val, beta_val) {
  
  df <- mc_results %>%
    filter(success == TRUE,
           alpha_true == alpha_val,
           beta_true == beta_val) %>%
    select(alpha_true, beta_true, n, alpha_hat, beta_hat) %>%
    pivot_longer(
      cols = c(alpha_hat, beta_hat),
      names_to = "Parameter",
      values_to = "Estimate"
    ) %>%
    mutate(
      Parameter = recode(
        Parameter,
        alpha_hat = "hat(alpha)",
        beta_hat  = "hat(beta)"
      ),
      Panel = paste0(Parameter),
      n = factor(n)
    )
  
  ggplot(df, aes(x = n, y = Estimate)) +
    geom_boxplot(
      width = 0.60,
      fill = "white",
      color = "black",
      outlier.size = 1.1,
      outlier.alpha = 0.55
    ) +
    facet_wrap(
      ~ Panel,
      scales = "free_y",
      nrow = 1,
      labeller = label_parsed
    ) +
    labs(
      x = "Sample size (n)",
      y = "Estimated values"
    ) +
    theme_bw() +
    theme(
      strip.background = element_rect(fill = "grey85", color = NA),
      strip.text = element_text(face = "bold", size = 14),
      axis.title.x = element_text(size = 14, face = "bold"),
      axis.title.y = element_text(size = 14, face = "bold"),
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12),
      panel.grid.major = element_line(color = "grey82"),
      panel.grid.minor = element_blank(),
      panel.spacing = unit(1.2, "lines")
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
  seed = 123
)

# ------------------------------------------------------------
# Monte Carlo summaries
# ------------------------------------------------------------

mc_summary_full <- summarize_CB_mc(mc_res)

mc_summary_full %>%
  filter(alpha_true == 0.5) %>%
  mutate(across(where(is.numeric), ~ round(.x, 5)))

mc_summary_full %>%
  filter(alpha_true == 2) %>%
  mutate(across(where(is.numeric), ~ round(.x, 5)))

mc_summary_full %>%
  filter(alpha_true == 4) %>%
  mutate(across(where(is.numeric), ~ round(.x, 5)))

# ------------------------------------------------------------
# Boxplots
# ------------------------------------------------------------

p1 <- plot_box_cb_free_y(mc_res, alpha_val = 0.5, beta_val = 1) +
  ggtitle(expression(paste("Scenario A: ", alpha == 0.5, ", ", beta == 1)))

p2 <- plot_box_cb_free_y(mc_res, alpha_val = 0.5, beta_val = 3) +
  ggtitle(expression(paste("Scenario B: ", alpha == 0.5, ", ", beta == 3)))

p3 <- plot_box_cb_free_y(mc_res, alpha_val = 0.5, beta_val = 6) +
  ggtitle(expression(paste("Scenario C: ", alpha == 0.5, ", ", beta == 6)))

library(patchwork)

fig <- (p1 | p2 | p3) &
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 12
    )
  )

fig

p4 <- plot_box_cb_free_y(mc_res, alpha_val = 2, beta_val = 1) +
  ggtitle(expression(paste("Scenario D: ", alpha == 2, ", ", beta == 1)))

p5 <- plot_box_cb_free_y(mc_res, alpha_val = 2, beta_val = 3) +
  ggtitle(expression(paste("Scenario E: ", alpha == 2, ", ", beta == 3)))

p6 <- plot_box_cb_free_y(mc_res, alpha_val = 2, beta_val = 6) +
  ggtitle(expression(paste("Scenario F: ", alpha == 2, ", ", beta == 6)))

fig <- (p4 | p5 | p6) &
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 12
    )
  )

fig

p7 <- plot_box_cb_free_y(mc_res, alpha_val = 4, beta_val = 1) +
  ggtitle(expression(paste("Scenario G: ", alpha == 4, ", ", beta == 1)))

p8 <- plot_box_cb_free_y(mc_res, alpha_val = 4, beta_val = 3) +
  ggtitle(expression(paste("Scenario H: ", alpha == 4, ", ", beta == 3)))

p9 <- plot_box_cb_free_y(mc_res, alpha_val = 4, beta_val = 6) +
  ggtitle(expression(paste("Scenario I: ", alpha == 4, ", ", beta == 6)))

fig <- (p7 | p8 | p9) &
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 12
    )
  )

fig




