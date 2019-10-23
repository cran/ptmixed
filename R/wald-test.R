#' Wald test for regression coefficients
#'
#' Compute a multivariate Wald test for one of the following
#' models: Poisson-Tweedie GLMM, negative binomial GLMM, 
#' Poisson-Tweedie GLM, negative binomial GLM. The null
#' hypothesis has to be specified in the (matrix) form
#' $L b = k$
#' 
#' @param obj an object of class \code{ptglmm} (obtained from 
#' \code{ptmixed} or \code{nbmixed}) or \code{ptglm} (obtained from 
#' \code{ptglm} or \code{nbglm})
#' @param L a matrix used to define the hypothesis to test,
#' in the form $L b = k$
#' @param k a vector used to define the hypothesis to test,
#' in the form $L b = k$. Default is a null vector ($L b = 0$)
#' @return A data frame with the result of the test
#' @export
#' @author Mirko Signorelli
#' @examples
#' \donttest{
#' # generate data
#' set.seed(1234)
#' n = 50
#' group = rep(c(0,1), each = n/2)
#' age = rpois(n, lambda = 5)
#' beta = c(3, 0.05, 0.2, 0.03)
#' X = model.matrix(~group + age + age*group)
#' mu = exp(X %*% beta)
#' y = rep(NA, n) 
#' library(tweeDEseq)
#' for (i in 1:n) y[i] = rPT(1, mu = mu[i], D = 2, a = 0, max = 1000)
#' dataset = data.frame(y, group, age)
#' rm(list = setdiff(ls(), 'dataset'))
#' # estimate a negative binomial glm
#' fit1 = nbglm(formula = y ~ group + age + age*group, data = dataset)
#' # define L for beta2 = beta4 = 0
#' L = matrix(0, nrow = 2, ncol = 4)
#' L[1, 2] = L[2, 4] = 1              
#' # compute multivariate Wald test
#' wald.test(obj = fit1, L = L, k = NULL)
#' }

wald.test <- function (obj, L, k = NULL) {
  requireNamespace('aod')
  requireNamespace('matrixcalc')
  mle = obj$mle
  if (class(obj)[1] == 'ptglmm') p = length(mle) - 3
  else if (class(obj)[1] == 'ptglm') p = length(mle) - 2
  else stop('obj should be of class ptglmm or ptglm')
  beta.hat = mle[1:p]
  if (ncol(L) != p) stop('L is not conformable with b')
  if (is.null(k)) k = rep(0, nrow(L))
  if (nrow(L) != length(k)) stop('number of rows in L and K differs')
  check.pdef1 = matrixcalc::is.positive.definite(obj$fisher.info, tol = 1e-10)
  last.eig = tail(eigen(obj$fisher.info)$values, 1)
  if (check.pdef1 & last.eig > 0.01) {
    inv.hess = solve(obj$fisher.info)
    var.beta = inv.hess[1:p, 1:p]
  }
  else {
    red.hess = obj$fisher.info[1:p, 1:p]
    check.pdef2 = matrixcalc::is.positive.definite(red.hess, tol = 1e-10)
    last.eig.red = tail(eigen(red.hess)$values, 1)
    if (check.pdef2 & last.eig.red > 0.01) {
      warning('Full Fisher information matrix not positive definite. Standard errors are estimated
                using the hessian of the profile likelihood of beta | (D, a, sigma2)')
      inv.hess = solve(red.hess)
      var.beta = inv.hess[1:p, 1:p]
    }
    else stop('Fisher information matrix is not positive definite. Standard errors cannot be computed.')
  }
  test = aod::wald.test(b = beta.hat, Sigma = var.beta,
                        L = L, H0 = k)
  out = as.data.frame(t(test$result$chi2))
  return(out)
}
