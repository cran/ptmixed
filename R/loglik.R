#' Loglikelihood of Poisson-Tweedie generalized linear mixed model with random intercept
#'
#' Evaluates the loglikelihood of a Poisson-Tweedie generalized linear mixed model 
#' with random intercept, using the adaptive Gauss-Hermite quadrature rule.
#' @param beta Vector of regression coefficients
#' @param D Deviance
#' @param a Power parameter (must be < 1)
#' @param Sigma A matrix with the variance of the random intercept
#' @param y Response vector (discrete)
#' @param X Design matrix for the fixed effects
#' @param Z Design matrix for the random effects
#' @param t Time variable
#' @param id Id indicator (it should be numeric)
#' @param offset Offset term to be added to the linear predictor
#' @param GHk Number of quadrature points (default is 10)
#' @param tol Tolerance value for the evaluation of the probability mass function of the Poisson-Tweedie distribution
#' @return The loglikelihood value obtained using a Gauss-Hermite quadrature 
#' approximation with \code{GHk} quadrature points.
#' @export
#' @author Mirko Signorelli
#' @seealso \code{\link{ptmixed}} and the examples therein

loglik.pt.1re <- function(beta, D, a, Sigma, y, X, Z, t, id, offset = NULL, GHk = 10, tol = 1e-323) {
  requireNamespace('tweeDEseq')
  requireNamespace('mvtnorm')
  with.offset = !is.null(offset)
  RE.size = 1
  GHs <- GHpoints.pt.1re(y, id, X, Z, beta, D, a, Sigma, offset, RE.size = RE.size,
                         GHk, tol = tol) # added offset argument
  if (with.offset) Delta.ij <- c(X %*% beta + c(offset))
  else Delta.ij <- c(X %*% beta)
  # calculate log Lik
  log.p.b <- t(sapply(GHs$scaled.b, mvtnorm::dmvnorm, mean = rep(0, RE.size), 
                      sigma = Sigma, log = TRUE))
  aux = dim(GHs$scaled.b[[1]])[1]
  n = length(unique(id))
  mus = matrix(NA, n*t, aux) 
  k = 1
  for (i in 1:n) {
    d.ij = Delta.ij[k:(k+t-1)]
    Z.ij = as.matrix(Z[k:(k+t-1),])
    b.values = GHs$scaled.b[[i]]
    for (j in 1:aux) {
      mus[k:(k+t-1), j] = d.ij + Z.ij %*% b.values[j,]
    }
    k = k+t
  }
  mus = exp(mus)
  #
  ff <- function (yy) {
    log.p.y.b <- rowsum(dpt.yvec.mumatr(yvec = yy, mu.matr = mus, D, a, 
                                        log = T, tol = tol), id)
    #    log(mapply(dPT, x = yy, mu = mus, D = D, a = a, tol = 1e-500))
    c((GHs$det.Bs * exp(log.p.y.b + log.p.b)) %*% GHs$wGH)
  }
  p.y <- ff(y)
  p.y[which(p.y < 1e-323)] = 1e-323 # otherwise it's too small for log(p.y)
  #
  sum(log(p.y), na.rm = TRUE) 
}