# ------------------------------------------------------------
# Required packages
# ------------------------------------------------------------
# install.packages("ggplot2")
# install.packages("dplyr")
# install.packages("tidyr")
# install.packages("CB")   
# install.packages(betareg)

library(ggplot2)
library(dplyr)
library(tidyr)
library(CB)
library(betareg)

# ------------------------------------------------------------
# Food data
# ------------------------------------------------------------

data(FoodExpenditure)

foo <- FoodExpenditure$food

inc <- FoodExpenditure$income

x <- foo/inc

# ------------------------------------------------------------
# Fit competing models
# ------------------------------------------------------------

cmp <- compareCBmodels(x)

summary(cmp)

# Comparison of fitted models for unit data
# ---------------------------------------------
#
# Call:
# compareCBmodels(x = x)
#
# Sample size: 38 
# 
# Model comparison:
#  Model npar  logLik      AIC     AICc      BIC convergence
#     CB    2 36.3303 -68.6606 -68.3178 -65.3855           0
#      B    2 35.3464 -66.6929 -66.3500 -63.4177           0
#      K    2 33.4891 -62.9782 -62.6353 -59.7030           0
#
# Parameter estimates:
#  Model parameter estimate      se
#     CB     alpha   8.0923  2.5750
#     CB      beta  11.2016  3.7639
#      B    shape1   6.0716  1.3586
#      B    shape2  14.8221  3.3988
#      K         a   2.9546  0.3692
#      K         b  26.9654 10.8271
#
# Best model according to AIC:   CB 
# Best model according to AICc:  CB 
# Best model according to BIC:   CB 

goffit <- gofCB_boot(x = x, 
                     B = 1000, 
                     method = "L-BFGS-B", 
                     multistart = FALSE, 
                     seed = 2026)

summary(goffit)

# Bootstrap goodness-of-fit for the CB distribution
# ------------------------------------------------
#
# Call:
# gofCB_boot(x = x, B = 1000, seed = 2026, method = "L-BFGS-B", 
#     multistart = FALSE)
#
# Sample size: 38 
# Bootstrap replicates: 1000 
# Valid bootstrap replicates: 1000 
#
# Parameter estimates:
#  Parameter Estimate
#      alpha   8.0923
#       beta  11.2016
#
# Fit statistics:
#   logLik      AIC      BIC
#  36.3303 -68.6606 -65.3855
#
# Goodness-of-fit statistics:
#              Test Statistic p_value
#  Anderson-Darling    0.3077  0.6464
#  Cramer-von Mises    0.0403  0.7163

# ------------------------------------------------------------
# Fitted histogram
# ------------------------------------------------------------

dens_grid <- seq(0.01, 0.75, length.out = 600)

dens_df <- data.frame(
  x = dens_grid,
  B = dbeta(dens_grid, shape1 = 6.0716, shape2 = 14.8221),
  K = CB:::dKum(dens_grid, a = 2.9546, b = 26.9654),
  CB = dCB(dens_grid, alpha = 8.0923, beta = 11.2016)
) %>%
  pivot_longer(
    cols = c(B, K, CB),
    names_to = "Distribution",
    values_to = "Density"
  ) %>%
  mutate(
    Distribution = factor(
      Distribution,
      levels = c("CB", "B", "K")
    )
  )

ggplot(data.frame(x = x), aes(x = x)) +
  geom_histogram(
    aes(y = after_stat(density)),
    breaks = seq(0, 0.7, by = 0.05),
    fill = "white",
    color = "black",
    alpha = 0.65
  ) +
  geom_line(
    data = dens_df,
    aes(x = x, y = Density,
        linetype = Distribution,
        linewidth = Distribution),
    color = "black"
  ) +
  scale_linetype_manual(
    values = c(
      CB = "solid",
      B = "dashed",
      K = "dotdash"
    )
  ) +
  scale_linewidth_manual(
    values = c(
      CB = 1.1,
      B = 0.9,
      K = 0.9
    )
  ) +
  labs(
    x = "Food expenditure share",
    y = "Density function",
    linetype = "Distribution",
    linewidth = "Distribution"
  ) +
  guides(
    linewidth = "none",
    linetype = guide_legend(
      override.aes = list(linewidth = c(1.1, 0.9, 0.9))
    )
  ) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid.major = element_line(color = "grey85", linewidth = 0.35),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7),
    
    axis.title = element_text(size = 15, face = "bold"),
    axis.text = element_text(size = 12, color = "black"),
    
    legend.position = c(0.74, 0.98),
    legend.justification = c(0, 1),
    legend.background = element_rect(
      fill = scales::alpha("white", 0.85),
      color = "grey70",
      linewidth = 0.3
    ),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    legend.key.width = unit(1.6, "cm"),
    legend.key.height = unit(0.45, "cm")
  )

# ------------------------------------------------------------
# cumulatives distributions
# ------------------------------------------------------------

alpha_cb <- 8.0923
beta_cb <- 11.2016

grid_x <- seq(0.01, 0.75, length.out = 600)

cb_df <- data.frame(
  x = grid_x,
  CDF = pCB(grid_x, alpha = alpha_cb, beta = beta_cb)
)

ggplot() +
  stat_ecdf(
    data = data.frame(x = x),
    aes(x = x, linetype = "Empirical"),
    geom = "step",
    linewidth = 1,
    color = "black"
  ) +
  geom_line(
    data = cb_df,
    aes(x = x, y = CDF, linetype = "CB"),
    linewidth = 1,
    color = "black"
  ) +
  scale_linetype_manual(
    values = c(
      Empirical = "solid",
      CB = "longdash"
    )
  ) +
  labs(
    x = "Food expenditure share",
    y = "Cumulative function",
    linetype = "Distribution"
  ) +
  coord_cartesian(xlim = c(0.01, 0.75), ylim = c(0, 1)) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid.major = element_line(color = "grey85", linewidth = 0.35),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7),
    axis.title = element_text(size = 15, face = "bold"),
    axis.text = element_text(size = 12, color = "black"),
    legend.position = c(0.02, 0.98),
    legend.justification = c(0, 1),
    legend.background = element_rect(
      fill = scales::alpha("white", 0.85),
      color = "grey70",
      linewidth = 0.3
    ),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    legend.key.width = unit(1.6, "cm"),
    legend.key.height = unit(0.45, "cm")
  )


