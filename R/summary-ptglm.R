#' Summarizing Poisson-Tweedie and negative binomial GLM estimation results
#' 
#' Provides parameter estimates, standard errors and univariate Wald test
#' for the Poisson-Tweedie and the negative binomial generalized linear models
#' (fitted through \code{ptglm} and \code{nbglm} respectively)
#' 
#' @param object an object of class \code{ptglm} 
#' (obtained from \code{ptglm} or \code{nbglm}).
#' @param ... Further arguments passed to or from other methods.
#' @return A list with the following elements: \code{logl}, \code{coefficients}, 
#' \code{D}, \code{a}
#' @export
#' @author Mirko Signorelli
#' @seealso \code{\link{ptglm}}, \code{\link{nbglm}} and the examples therein
 
summary.ptglm = function(object, ...) {
  ncov = length(object$mle) - 2
  beta = object$mle[1:ncov]
  D = object$mle[ncov+1]
  a = object$mle[ncov+2]
  coef.table = cbind(beta)
  colnames(coef.table) = 'Estimate'
  # compute asymptotic SEs
  requireNamespace('matrixcalc')
  check.pdef1 = matrixcalc::is.positive.definite(object$fisher.info, tol = 1e-10)
  last.eig = tail(eigen(object$fisher.info)$values, 1)
  if (check.pdef1 & last.eig > 0.01) {
    inv.hess = solve(object$fisher.info)
    se.beta = sqrt(diag(inv.hess)[1:ncov])
  }
  else {
    red.hess = object$fisher.info[1:ncov, 1:ncov]
    check.pdef2 = matrixcalc::is.positive.definite(red.hess, tol = 1e-10)
    last.eig.red = tail(eigen(red.hess)$values, 1)
    if (check.pdef2 & last.eig.red > 0.01) {
      warning('Full Fisher information matrix not positive definite. Standard errors are estimated
              using the hessian of the profile likelihood of beta | (D, a, sigma2)')
      inv.hess = solve(red.hess)
      se.beta = sqrt(diag(inv.hess))
    }
    else stop('Fisher information matrix is not positive definite. Standard errors cannot be computed.')
  }
  z.score = beta/se.beta
  p = 2*pnorm(abs(z.score), lower.tail = F)
  coef.table = cbind(beta, se.beta, z.score, p)
  colnames(coef.table) = c('Estimate', 'Std. error', 'z', 'p.value')
  cat(paste('Loglikelihood:', round(object$logl, 3), '\n'))
  cat('Parameter estimates:\n')
  print(coef.table)
  cat('\n')
  cat(paste('Dispersion =', round(D,2), '\n'))
  cat(paste('Power =', round(a,2), '\n'))
  out = list('logl' = round(object$logl,2), 'coefficients' = coef.table,
             'D'= D, 'a' = a)
}
