#' @export
summary.fitCB_mom <- function(object, digits = 4, ...) {
  
  par_tab <- data.frame(
    Parameter = names(object$par),
    Estimate = as.numeric(object$par),
    row.names = NULL
  )
  
  mom_tab <- data.frame(
    Moment = names(object$moments),
    Value = as.numeric(object$moments),
    row.names = NULL
  )
  
  out <- list(
    call = object$call,
    coefficients = par_tab,
    moments = mom_tab,
    success = object$success
  )
  
  class(out) <- c("summary.fitCB_mom", "list")
  out
}


#' @export
print.summary.fitCB_mom <- function(x, digits = 4, ...) {
  
  cat("\n")
  cat("Moment estimation for the CB distribution\n")
  cat(strrep("-", 41), "\n\n", sep = "")
  
  cat("Call:\n")
  print(x$call)
  
  cat("\nParameter estimates:\n")
  coef_tab <- x$coefficients
  coef_tab$Estimate <- round(coef_tab$Estimate, digits)
  print(coef_tab, row.names = FALSE)
  
  cat("\nMoment-based quantities:\n")
  mom_tab <- x$moments
  mom_tab$Value <- round(mom_tab$Value, digits)
  print(mom_tab, row.names = FALSE)
  
  cat("\nSuccessful fit:", x$success, "\n")
  
  invisible(x)
}