#' @export
summary.compareCBmodels <- function(object, digits = 4, ...) {
  
  out <- list(
    call = object$call,
    n = object$n,
    comparison = object$comparison,
    estimates = object$estimates,
    best_AIC = object$best_AIC,
    best_AICc = object$best_AICc,
    best_BIC = object$best_BIC
  )
  
  class(out) <- c("summary.compareCBmodels", "list")
  out
}

#' @export
print.summary.compareCBmodels <- function(x, digits = 4, ...) {
  
  cat("\nComparison of fitted models for unit data\n")
  cat(strrep("-", 45), "\n\n", sep = "")
  
  if (!is.null(x$call)) {
    cat("Call:\n")
    print(x$call)
    cat("\n")
  }
  
  cat("Sample size:", x$n, "\n\n")
  
  # ----------------------------------------------------------
  # Comparison table
  # ----------------------------------------------------------
  
  cat("Model comparison:\n")
  
  comp <- x$comparison
  
  if ("success" %in% names(comp)) {
    comp$success <- NULL
  }
  
  comp$logLik <- round(comp$logLik, digits)
  comp$AIC    <- round(comp$AIC, digits)
  
  if ("AICc" %in% names(comp)) {
    comp$AICc <- round(comp$AICc, digits)
  }
  
  comp$BIC    <- round(comp$BIC, digits)
  
  print(comp, row.names = FALSE)
  
  # ----------------------------------------------------------
  # Estimates
  # ----------------------------------------------------------
  
  cat("\nParameter estimates:\n")
  
  est <- x$estimates
  est$estimate <- round(est$estimate, digits)
  est$se       <- round(est$se, digits)
  
  print(est, row.names = FALSE)
  
  # ----------------------------------------------------------
  # Best models
  # ----------------------------------------------------------
  
  cat("\nBest model according to AIC:  ", x$best_AIC, "\n")
  
  if (!is.null(x$best_AICc)) {
    cat("Best model according to AICc: ", x$best_AICc, "\n")
  }
  
  cat("Best model according to BIC:  ", x$best_BIC, "\n")
  
  invisible(x)
}