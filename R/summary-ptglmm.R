#' Summarizing Poisson-Tweedie mixed model estimation results
#' 
#' Provides parameter estimates, standard errors and univariate Wald test
#' for the Poisson-Tweedie generalized linear mixed model
#' 
#' @param object an object of class \code{ptglmm} (obtained from \code{ptmixed}).
#' @param wald logical. If \code{TRUE}, standard errors and univariate Wald test are computed. Default is \code{TRUE}.
#' @param ... Further arguments passed to or from other methods.
#' @return A list with the following elements: value of the loglikelihood at the MLE (\code{logl}),
#' table with maximum likelihood estimates of the regression coefficients, SEs and Wald tests \code{coefficients}, 
#' and maximum likelihood estimates of the other parameters (\code{D}, \code{a} and \code{sigma2})
#' @export
#' @author Mirko Signorelli
#' @seealso \code{\link{ptmixed}} and the examples therein
 
summary.ptglmm = function(object, wald = T, ...) {
  if (wald) {
    if ('fisher.info' %in% ls(object) == F) stop('Observed Fisher information matrix not available in object. 
  It is needed to return standard errors and Wald test. Set hessian = T in ptmixed() to obtain it.')
  }
  ncov = length(object$mle) - 3
  beta = object$mle[1:ncov]
  D = object$mle[ncov+1]; a = object$mle[ncov+2]; sigma2 = object$mle[ncov+3]
  coef.table = cbind(beta)
  colnames(coef.table) = 'Estimate'
  if (wald == T) {
    requireNamespace('matrixcalc')
    check.pdef1 = matrixcalc::is.positive.definite(object$fisher.info, tol = 1e-10)
    if (check.pdef1) {
      inv.hess = solve(object$fisher.info)
      se.beta = sqrt(diag(inv.hess)[1:ncov])
    }
    else {
      red.hess = object$fisher.info[1:ncov, 1:ncov]
      check.pdef2 = matrixcalc::is.positive.definite(red.hess, tol = 1e-10)
      if (check.pdef2) {
        warning('Full Fisher information matrix not positive definite. Standard errors are estimated
                using the hessian of the profile likelihood of beta | (D, a, sigma2)')
        inv.hess = solve(red.hess)
        se.beta = sqrt(diag(red.hess))
      }
      if (!check.pdef2) stop('Fisher information matrix is not positive definite. Standard errors cannot be computed.')
    }
    z.score = beta/se.beta
    p = (1-pnorm(abs(z.score)))
    coef.table = cbind(beta, se.beta, z.score, p)
    colnames(coef.table) = c('Estimate', 'Std. error', 'z', 'p.value')
  }
  cat(paste('Loglikelihood:', round(object$logl, 3), '\n'))
  cat('Parameter estimates:\n')
  print(coef.table)
  cat('\n')
  cat(paste('Dispersion =', round(D,2), '\n'))
  cat(paste('Power =', round(a,2), '\n'))
  cat(paste('Variance =', round(sigma2,2), '\n'))
  out = list('logl' = round(object$logl,2), 'coefficients' = coef.table,
             'D'= D, 'a' = a, 'sigma2' = sigma2)
}
