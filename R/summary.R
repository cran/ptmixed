#' Summary results for Poisson-Tweedie generalized linear mixed model
#'
#' Provides parameter estimates, standard errors and univariate Wald test
#' for the Poisson-Tweedie generalized linear mixed model
#' 
#' @param object an object of class \code{ptmm} (obtained from \code{pt.mixed}).
#' @param wald logical. If \code{TRUE}, standard errors and univariate Wald test are computed. Default is \code{TRUE}.
#' @param ... Further arguments passed to or from other methods.
#' @return A list with the following elements: \code{logl}, \code{coefficients}, 
#' \code{D}, \code{a}, \code{sigma2}
#' @export
#' @author Mirko Signorelli
#' @seealso \code{\link{ptmixed}} and the examples therein
 
summary.ptmm = function(object, wald = T, ...) {
  if (wald) {
    if ('hessian' %in% ls(object) == F) stop('Hessian matrix not available in object. 
  Hessian is needed to return standard errors and Wald test. Use pt.model with option wald = F to obtain it.')
  }
  ncov = length(object$mle) - 3
  beta = object$mle[1:ncov]
  D = object$mle[ncov+1]; a = object$mle[ncov+2]; sigma2 = object$mle[ncov+3]
  coef.table = cbind(beta)
  colnames(coef.table) = 'Estimate'
  if (wald == T) {
    requireNamespace('matrixcalc')
    check.pdef1 = matrixcalc::is.positive.definite(object$hessian, tol = 1e-10)
    if (check.pdef1) {
      inv.hess = solve(object$hessian)
      se.beta = sqrt(diag(inv.hess)[1:ncov])
    }
    else {
      red.hess = object$hessian
      check.pdef2 = matrixcalc::is.positive.definite(red.hess, tol = 1e-10)
      if (check.pdef2) {
        warning('Full hessian not positive definite. Standard errors are estimated
                using the hessian of the profile likelihood of beta | (D, a, sigma2)')
        inv.hess = solve(red.hess)
        se.beta = sqrt(diag(red.hess))
      }
      if (!check.pdef2) stop('Hessian matrix is not positive definite. Standard errors cannot be computed.')
    }
    z.score = beta/se.beta
    p = (1-pnorm(abs(z.score)))
    coef.table = cbind(beta, se.beta, z.score, p)
    colnames(coef.table) = c('Estimate', 'Std. error', 'z', 'p.value')
  }
  cat(paste('Loglikelihood:', object$logl, '\n'))
  cat('Parameter estimates:\n')
  print(coef.table)
  cat('\n')
  cat(paste('Deviance =', round(D,2), '\n'))
  cat(paste('Power =', round(a,2), '\n'))
  cat(paste('Variance =', round(sigma2,2), '\n'))
  out = list('logl' = round(object$logl,2), 'coefficients' = coef.table,
             'D'= D, 'a' = a, 'sigma2' = sigma2)
}
