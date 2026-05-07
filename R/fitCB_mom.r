#' Moment Estimation for the Contracted-Beta Distribution
#'
#' Estimates the parameters of the Contracted-Beta (CB)
#' distribution by matching the first two theoretical raw moments with
#' their empirical counterparts.
#'
#' The estimator is obtained in closed form by exploiting the relationship
#' between the moments of the Contracted-Beta distribution and those of
#' the Beta baseline distribution.
#'
#' @param x Numeric vector of observations in \eqn{(0,1)}.
#'
#' @return An object of class \code{"fitCB_mom"}, which is a list containing:
#' \itemize{
#'   \item \code{par}: estimated parameters \eqn{(\alpha, \beta)}.
#'   \item \code{moments}: empirical and transformed moment quantities used in the estimation.
#'   \item \code{success}: logical indicator of whether the estimates are valid.
#'   \item \code{call}: matched function call.
#' }
#'
#' @details
#' Let \eqn{T} follow a Contracted-Beta distribution obtained from a
#' Beta\eqn{(\alpha,\beta)} baseline. The raw moments of \eqn{T} can be
#' expressed as scaled versions of the corresponding Beta moments.
#'
#' The estimation proceeds by:
#' \enumerate{
#'   \item Matching the first two empirical moments of \eqn{T}.
#'   \item Recovering the corresponding Beta moments via known scaling factors.
#'   \item Solving the resulting system to obtain \eqn{\alpha} and \eqn{\beta}.
#' }
#'
#' This yields a fast and fully explicit estimator, suitable for initialization
#' in maximum likelihood estimation procedures.
#'
#' @examples
#' set.seed(123)
#' x <- rCB(n = 200, alpha = 2, beta = 5)
#'
#' fit_mom <- fitCB_mom(x)
#' summary(fit_mom)
#'
#' @export
fitCB_mom <- function(x) {
  
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  
  if (length(x) < 2L) {
    stop("x must have length >= 2.", call. = FALSE)
  }
  if (any(x <= 0 | x >= 1)) {
    stop("x must lie strictly in (0,1).", call. = FALSE)
  }
  
  m1 <- mean(x)
  m2 <- mean(x^2)
  
  c1 <- log(2)
  c2 <- 1 / 2
  
  m1_beta <- m1 / c1
  m2_beta <- m2 / c2
  
  v_beta <- m2_beta - m1_beta^2
  
  phi <- m1_beta * (1 - m1_beta) / v_beta - 1
  
  alpha_hat <- m1_beta * phi
  beta_hat  <- (1 - m1_beta) * phi
  
  success <- is.finite(alpha_hat) && is.finite(beta_hat) &&
    alpha_hat > 0 && beta_hat > 0 &&
    m1_beta > 0 && m1_beta < 1 &&
    v_beta > 0
  
  out <- list(
    par = c(alpha = alpha_hat, beta = beta_hat),
    success = success,
    moments = c(
      mean_T = m1,
      second_T = m2,
      mean_beta = m1_beta,
      second_beta = m2_beta,
      var_beta = v_beta
    ),
    call = match.call()
  )
  
  class(out) <- c("fitCB_mom", "list")
  out
}