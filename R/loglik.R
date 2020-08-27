#' Loglikelihood of Poisson-Tweedie generalized linear mixed model with random intercept
#'
#' Evaluates the loglikelihood of a Poisson-Tweedie generalized linear mixed model 
#' with random intercept, using the adaptive Gauss-Hermite quadrature rule.
#' @param beta Vector of regression coefficients
#' @param D Dispersion parameter (must be > 1)
#' @param a Power parameter (must be < 1)
#' @param Sigma A matrix with the variance of the random intercept
#' @param y Response vector (discrete)
#' @param X Design matrix for the fixed effects
#' @param Z Design matrix for the random effects
#' @param id Id indicator (it should be numeric)
#' @param offset Offset term to be added to the linear predictor
#' @param GHk Number of quadrature points (default is 5)
#' @param tol Tolerance value for the evaluation of the probability mass function of the Poisson-Tweedie distribution
#' @param GHs Quadrature points at which to evaluate the loglikelihood. If \code{NULL} (default value), 
#' the GH quadrature points are computed internally
#' @return The loglikelihood value obtained using a Gauss-Hermite quadrature
#' approximation with \code{GHk} quadrature points.
#' @export
#' @author Mirko Signorelli
#' @references Signorelli, M., Spitali, P., Tsonaka, R. (2020). Poisson-Tweedie 
#' mixed-effects model: a flexible approach for the analysis of longitudinal RNA-seq
#' data. Statistical Modelling. URL: https://doi.org/10.1177/1471082X20936017
#' @seealso \code{\link{ptmixed}} and the examples therein

loglik.pt.1re = function(beta, D, a, Sigma, y, X, Z, id, offset = NULL, 
                          GHk = 5, tol = 1e-323, GHs = NULL) {
  requireNamespace('tweeDEseq')
  requireNamespace('mvtnorm')
  with.offset = !is.null(offset)
  RE.size = 1
  numid = as.numeric(as.factor(id))
  if (is.null(GHs)) {
    GHs = GHpoints.pt.1re(y = y, id = id, X = X, Z = Z, 
                       beta = beta, D = D, a = a, Sigma = Sigma, 
                       offset = offset, RE.size = RE.size,
                       GHk = GHk, tol = tol)
  }
  if (with.offset) Delta.ij = c(X %*% beta + c(offset))
  else Delta.ij = c(X %*% beta)
  # calculate log Lik
  log.p.b = sapply(GHs$scaled.b, mvtnorm::dmvnorm, mean = rep(0, RE.size), 
                      sigma = Sigma, log = TRUE)
  if (GHk > 1) log.p.b = t(log.p.b)
  else if (GHk == 1) log.p.b = as.matrix(log.p.b)
  aux = dim(GHs$scaled.b[[1]])[1]
  n = length(unique(id))
  mus = matrix(NA, length(y), aux) 
  for (i in 1:n) {
    pos = (numid == i)
    d.ij = Delta.ij[pos]
    Z.ij = as.matrix(Z[pos,])
    b.values = GHs$scaled.b[[i]]
    for (j in 1:aux) {
      mus[pos, j] = d.ij + Z.ij %*% b.values[j,]
    }
  }
  mus = exp(mus)
  #
  ff = function (yy) {
    log.p.y.b = rowsum(dpt.yvec.mumatr(yvec = yy, mu.matr = mus, D, a, 
                                        log = T, tol = tol), id)
    c((GHs$det.Bs * exp(log.p.y.b + log.p.b)) %*% GHs$wGH)
  }
  p.y = ff(y)
  p.y[which(p.y < 1e-323)] = 1e-323 # otherwise it's too small for log(p.y)
  #
  sum(log(p.y), na.rm = TRUE) 
}
