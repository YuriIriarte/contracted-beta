#' Compare Fitted Models for Unit Data
#'
#' Fits and compares the Beta, Kumaraswamy, and Contracted-Beta
#' (CB) distributions for data supported on the unit interval.
#'
#' The comparison is based on the maximized log-likelihood, Akaike information
#' criterion (AIC), corrected Akaike information criterion (AICc), and Bayesian
#' information criterion (BIC). Parameter estimates and approximate standard
#' errors are also returned for each fitted model.
#'
#' @param x Numeric vector of observations in \eqn{(0,1)}.
#' @param suppress_warnings Logical. If \code{TRUE}, warnings generated during
#' model fitting are suppressed.
#'
#' @return An object of class \code{"compareCBmodels"}.
#'
#' @examples
#' set.seed(123)
#' x <- rCB(n = 200, alpha = 2, beta = 5)
#' cmp <- compareCBmodels(x)
#' summary(cmp)
#'
#' @export
compareCBmodels <- function(x, suppress_warnings = TRUE) {
  
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  
  if (length(x) < 2L) {
    stop("x must have length >= 2.", call. = FALSE)
  }
  
  if (any(x <= 0 | x >= 1)) {
    stop("All observations must lie in (0,1).", call. = FALSE)
  }
  
  fit_try <- function(expr) {
    if (isTRUE(suppress_warnings)) {
      suppressWarnings(try(expr, silent = TRUE))
    } else {
      try(expr, silent = TRUE)
    }
  }
  
  n <- length(x)
  fits <- list()
  rows <- list()
  id <- 1L
  
  add_row <- function(model, fit) {
    data.frame(
      Model = model,
      npar = fit$k,
      logLik = fit$logLik,
      AIC = fit$AIC,
      AICc = if (!is.null(fit$AICc)) fit$AICc else {
        fit$AIC + (2 * fit$k * (fit$k + 1)) / max(n - fit$k - 1, 1)
      },
      BIC = fit$BIC,
      convergence = fit$convergence,
      success = fit$success,
      stringsAsFactors = FALSE
    )
  }
  
  fit_cb <- fit_try(fitCB_mle(x))
  
  if (!inherits(fit_cb, "try-error")) {
    fits$cb <- fit_cb
    rows[[id]] <- add_row("CB", fit_cb)
    id <- id + 1L
  }
  
  fit_beta <- fit_try(fitBeta_mle(x))
  
  if (!inherits(fit_beta, "try-error")) {
    fits$beta <- fit_beta
    rows[[id]] <- add_row("B", fit_beta)
    id <- id + 1L
  }
  
  fit_kum <- fit_try(fitKum_mle(x))
  
  if (!inherits(fit_kum, "try-error")) {
    fits$kumaraswamy <- fit_kum
    rows[[id]] <- add_row("K", fit_kum)
    id <- id + 1L
  }
  
  if (length(rows) == 0L) {
    stop("No model could be fitted.", call. = FALSE)
  }
  
  comparison <- do.call(rbind, rows)
  rownames(comparison) <- NULL
  
  comparison <- comparison[order(comparison$AIC), ]
  rownames(comparison) <- NULL
  
  estimates <- list()
  
  add_estimates <- function(model, fit) {
    data.frame(
      Model = model,
      parameter = names(fit$par),
      estimate = as.numeric(fit$par),
      se = as.numeric(fit$se),
      stringsAsFactors = FALSE
    )
  }
  
  if (!is.null(fits$cb)) {
    estimates[[length(estimates) + 1L]] <- add_estimates("CB", fits$cb)
  }
  
  if (!is.null(fits$beta)) {
    estimates[[length(estimates) + 1L]] <- add_estimates("B", fits$beta)
  }
  
  if (!is.null(fits$kumaraswamy)) {
    estimates[[length(estimates) + 1L]] <- add_estimates("K", fits$kumaraswamy)
  }
  
  estimates <- do.call(rbind, estimates)
  rownames(estimates) <- NULL
  
  out <- list(
    comparison = comparison,
    estimates = estimates,
    fits = fits,
    n = n,
    best_AIC = comparison$Model[which.min(comparison$AIC)],
    best_AICc = comparison$Model[which.min(comparison$AICc)],
    best_BIC = comparison$Model[which.min(comparison$BIC)],
    suppress_warnings = suppress_warnings,
    call = match.call()
  )
  
  class(out) <- c("compareCBmodels", "list")
  out
}