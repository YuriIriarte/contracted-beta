#' Density of the Contracted-Beta Distribution
#'
#' Computes the probability density function of the 
#' Contracted-Beta (CB) distribution.
#'
#' The CB distribution is defined through the stochastic representation
#' \deqn{T = X / U,}
#' where \eqn{X \sim Beta(\alpha,\beta)}, \eqn{U \sim Uniform(1,2)},
#' and \eqn{X} and \eqn{U} are independent.
#'
#' @param x Numeric vector of quantiles.
#' @param alpha Positive shape parameter.
#' @param beta Positive shape parameter.
#' @param log Logical. If \code{TRUE}, the log-density is returned.
#'
#' @return A numeric vector with the density or log-density values.
#'
#' @details
#' The density is evaluated using a numerically stable implementation based on
#' the log-density. This is useful when the density involves small differences
#' between Beta distribution functions.
#'
#' @examples
#' dCB(0.4, alpha = 2, beta = 5)
#' dCB(c(0.2, 0.5, 0.8), alpha = 2, beta = 5)
#' dCB(0.4, alpha = 2, beta = 5, log = TRUE)
#'
#' curve(dCB(x, alpha = 2, beta = 5), from = 0, to = 1)
#'
#' @export
dCB <- function(x, alpha, beta, log = FALSE) {
  
  logdens <- log_dCB_core(
    x = x,
    alpha = alpha,
    beta = beta
  )
  
  if (log) {
    return(logdens)
  }
  
  dens <- exp(logdens)
  dens[!is.finite(dens)] <- 0
  
  dens
}