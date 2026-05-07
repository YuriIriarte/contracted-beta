# ------------------------------------------------------------
# Required packages
# ------------------------------------------------------------
# install.packages("ggplot2")
# install.packages("dplyr")
# install.packages("tidyr")
# install.packages("CB")   
# install.packages(datasets)

library(ggplot2)
library(dplyr)
library(tidyr)
library(CB)
library(datasets)

# ------------------------------------------------------------
# Measurements on Petroleum Rock Samples (shape variable)
# ------------------------------------------------------------

x <- rock[,3]

# ------------------------------------------------------------
# Fit competing models
# ------------------------------------------------------------

cmp <- compareCBmodels(x)

summary(cmp)

# Call:
# compareCBmodels(x = x)
#
# Sample size: 48 
# 
# Model comparison:
#  Model npar  logLik       AIC      AICc       BIC convergence
#     CB    2 57.2306 -110.4611 -110.1945 -106.7187           0
#      B    2 55.6002 -107.2004 -106.9338 -103.4580           0
#      K    2 52.4915 -100.9831 -100.7164  -97.2407           0
#
# Parameter estimates:
#  Model parameter estimate      se
#     CB     alpha   8.5960  2.4256
#     CB      beta  18.7245  5.5210
#      B    shape1   5.9418  1.1814
#      B    shape2  21.2057  4.3469
#      K         a   2.7187  0.2935
#      K         b  44.6604 17.5757
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
# Sample size: 48 
# Bootstrap replicates: 1000 
# Valid bootstrap replicates: 1000 
#
# Parameter estimates:
#  Parameter Estimate
#      alpha   8.5959
#       beta  18.7245
#
# Fit statistics:
#   logLik       AIC       BIC
#  57.2306 -110.4611 -106.7187
#
# Goodness-of-fit statistics:
#              Test Statistic p_value
#  Anderson-Darling    0.4989  0.2557
#  Cramer-von Mises    0.0837  0.2428

# ------------------------------------------------------------
# Fitted histogram
# ------------------------------------------------------------

dens_grid <- seq(0.01, 0.55, length.out = 600)

dens_df <- data.frame(
  x = dens_grid,
  B = dbeta(dens_grid, shape1 = 5.9418, shape2 = 21.2057),
  K = CB:::dKum(dens_grid, a = 2.7187, b = 44.6604),
  CB = dCB(dens_grid, alpha = 8.5960, beta = 18.7245)
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
    breaks = seq(0, 0.7, by = 0.065),
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
    x = "Shape",
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

alpha_cb <- 8.5960
beta_cb <- 18.7245

grid_x <- seq(0.01, 0.55, length.out = 600)

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
    x = "Shape",
    y = "Cumulative function",
    linetype = "Distribution"
  ) +
  coord_cartesian(xlim = c(0.01, 0.6), ylim = c(0, 1)) +
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


