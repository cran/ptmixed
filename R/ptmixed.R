#' Poisson-Tweedie generalized linear mixed model
#'
#' Estimates the Poisson-Tweedie generalized linear mixed model with random intercept.
#' Likelihood approximation for the model is based on the adaptive Gauss-Hermite quadrature rule.
#' 
#' @param fixef.formula A formula for the fixed effects part of the model. It should be in the form \code{y ~ x1 + x2}
#' @param id A variable to distinguish observations from the same subject.
#' @param offset An offset to be added to the linear predictor. Default is \code{NULL}.
#' @param data A data frame containing the variables declared in \code{fixef.formula}.
#' @param npoints Number of quadrature points employed in the adaptive quadrature. Default is 10.
#' @param hessian Logical value. If \code{TRUE}, the hessian matrix is evaluated at the MLE. Default is \code{TRUE}.
#' @param trace Logical value. If \code{TRUE}, additional information is printed during the optimization. Default is \code{TRUE}.
#' @param theta.start Numeric vector comprising initial parameter values for the
#' vector of regression coefficients, the deviance parameter, the power parameter 
#' and the variance of the random intercept (to be specified exactlyin this order!). 
#' Default is \code{NULL}: initial
#' parameter estimates are computed automatically by the function.
#' @param reltol Relative tolerance to be used in optim. Default to 1e-8
#' @param maxit Vector containing the maximum number of iterations used in optim by
#' the Nelder-Mead method and, if this fails, by the BFGS method
#' @return A list containing the following elements: maximum likelihood estimate (\code{mle});  value of the
#' loglikelihood at the mle (\code{logl}); \code{convergence} value (if 0, the optimization converged);
#' the hessian evaluated at the mle (\code{hessian}), if \code{hessian = T}.
#' @import stats
#' @importFrom utils tail
#' @export
#' @author Mirko Signorelli
#' @seealso \code{\link{summary.ptmm}}
#' @examples
#' # generate data
#' set.seed(123)
#' n = 6; t = 3
#' id = rep(1:n, each = t)
#' rand.int = rep(rnorm(n, sd = 0.7), each = t)
#' group = rep(c(0,1), each = n*t/2)
#' time = rep(0:(t-1), n)
#' offset = rnorm(n*t, sd = 0.3)
#' 
#' beta = c(3, 0.3, 0.1)
#' X = model.matrix(~group + time)
#' mu = exp(X %*% beta + rand.int + offset)
#' y = rep(NA, n*t)
#' library(tweeDEseq)
#' for (i in 1:(n*t)) y[i] = rPT(1, mu = mu[i], D = 2, a = -0.5, max = 1000)
#' 
#' data.long = data.frame(y, group, time, id, offset)
#' rm(list = setdiff(ls(), 'data.long'))
#' 
#' # estimate the model
#' # very quick (and inaccurate) example of execution 
#' # (see fit2 below for a much more accurate, but slower, estimation)
#' fit1 = ptmixed(fixef.formula = y ~ group + time, id = data.long$id,
#'               offset = data.long$offset, data = data.long, npoints = 3, 
#'               hessian = FALSE, trace = TRUE, reltol = 1e-3)
#' # print summary:
#' summary(fit1, wald = FALSE)
#' 
#' # more accurate estimation (increase npoints and maxit, reduce reltol,
#' # evaluate hessian matrix at mle) - this takes more time
#' \donttest{
#' fit2 = ptmixed(fixef.formula = y ~ group + time, id = data.long$id,
#'               offset = data.long$offset, data = data.long, npoints = 10, 
#'               hessian = TRUE, trace = TRUE, reltol = 1e-8)
#' # print and get summary:
#' x2 = summary(fit2, wald = TRUE)
#' ls(x2)
#' }

ptmixed = function(fixef.formula, id, offset = NULL,
                   data, npoints = 10, hessian = T, trace = T,
                   theta.start = NULL, reltol = 1e-8, maxit = c(1e4, 100)) {
  # preliminary checks:
  if (length(fixef.formula) != 3) stop('fixef.formula should be in the form y ~ x1 + x2 +...')
  t = dim(data)[1] / length(unique(id))
  if (t %% 1 !=0) stop('The dataset appears to be unbalanced. The code for
                       the unbalanced case is not yet implemented (check back soon!)')
  if (npoints == 1) stop('Use at least 2 quadrature points')
  if (npoints %%1 !=0) stop('npoints should be a natural number > 1')
  if (!is.null(maxit)) {
    if (length(maxit) == 1) maxit = c(maxit, 100)
  }
  # identify elements
  y = data[, all.vars(fixef.formula[[2]])]
  X = model.matrix(fixef.formula[-2], data = data)
  Z = as.matrix(model.matrix(~1, data = data))
  # fix id:
  id = as.numeric(as.factor(id))
  # optim control values:
  optim.control.nm = list(maxit = maxit[1], reltol = reltol)
  optim.control.bfgs = list(maxit = maxit[2], reltol = reltol)
  if (trace) {
    optim.control.nm = list(trace = 1, REPORT = 1, maxit = maxit[1], reltol = reltol)
    optim.control.bfgs = list(trace = 1, REPORT = 1, maxit = maxit[2], reltol = reltol)
  }
  # starting values:
  if (is.null(theta.start)) {
    fixef.glmmad = fixef.formula
    if (!is.null(offset)) {
      fixef.glmmad = update.formula(fixef.glmmad, ~ . + offset(offset))
      df.glmmad = cbind(data, 'id' = id, 'offset' = offset)
    }
    else df.glmmad = cbind(data, 'id' = id)
    theta.init = get.initial.theta(fixef.formula = fixef.glmmad, data = df.glmmad, 
                                   y = y, id = id)
    if (exp(tail(theta.init,1)) < 0.001) warning('Starting value of variance of random intercept is < 0.001.
      Consider using a GLM instead of a GLMM.')
  }
  else if (!is.null(theta.start)) {
    q = length(theta.start)
    p = dim(X)[2]
    if (q != p+3) stop('Wrong number of elements in theta.start')
    if (theta.start[p+1] <= 1) stop('Initial deviance value should be > 1')
    if (theta.start[p+2] >= 1) stop('Initial power value should be < 1')
    if (theta.start[p+3] <= 1e-6) stop('Initial variance value should be > 0')
    theta.init = c(theta.start[1:p], log(theta.start[p+1]-1),
                   log(1-theta.start[p+2]), log(theta.start[p+3]))
  }
  # rename correctly elements in theta.init
  names(theta.init) = c(colnames(X), 'D', 'a', 'sigma2')
  # negative loglikelihood:
  nlogl = function(theta, offset, y, X, Z, id, GHk) {
    p = length(theta)
    beta = theta[1:(p-3)]
    D = 1+exp(theta[p-2])
    a = 1-exp(theta[p-1])
    Sigma = as.matrix(exp(theta[p]))
    ll = loglik.pt.1re(beta, D, a, Sigma, y, X, Z, id, offset = offset, GHk)
    return(-ll)
  }
  # check that starting point is finite:
  try.init1 = try(nlogl(theta.init, offset, y, X, Z, id = id, GHk = npoints))
  if (trace) cat(paste('initial loglik value:', try.init1));cat('\n')
  if(is.infinite(try.init1)) {
    try.init2 = Inf
    while(is.infinite(try.init2)) {
      new.a.init = runif(1, -10, 1)
      p = length(theta.init)
      theta.init[p-1] = log(1 - new.a.init)
      #print(paste('proposed a.init:', new.a.init))
      try.init2 = try(nlogl(theta.init, offset, y, X, Z, id, GHk = npoints))
      print(paste('retry: initial loglik value:', try.init2))
    }
  }
  # optimization: try Nelder-Mead
  if (trace) cat('Beginning optimization with Nelder-Mead:')
  mle = try( optim(theta.init, nlogl, method = "Nelder-Mead", offset = offset,
                   y = y, X = X, Z = Z, id = id, GHk = npoints, hessian = F,
                   control = optim.control.nm) )
  # check convergence
  redo = redo2 = F
  check1 = exists('mle')
  if (!check1) redo = T
  # if it does not exist, redo
  if (check1) {
    check2 = inherits(mle, 'try-error')
    if (check2) redo = T
    # if it has error status, redo
    if (check2==F) {
      check3 = (mle$convergence == 0)
      if (!check3) redo = T
    }
  }
  if (trace & redo == F) print('Optimization was successful')
  # if it failed: try BFGS
  if (redo) {
    print('Optimization with Nelder-Mead was not successful. Trying BFGS...')
    if (trace) cat('Beginning optimization with BFGS:')
    mle = try( optim(theta.init, nlogl, method = "BFGS", offset = offset,
                     y = y, X = X, Z = Z, id = id, GHk = npoints, hessian = F,
                     control = optim.control.bfgs) )
    check1 = exists('mle')
    if (!check1) redo2 = T
    # if it does not exist, redo
    if (check1) {
      check2 = inherits(mle, 'try-error')
      if (check2) redo2 = T
      # if it has error status, redo
      if (check2==F) {
        check3 = (mle$convergence == 0)
        if (!check3) redo2 = T
      }
    }
    if (trace & redo2 == F) print('Second optimization was successful')
  }
  # if no convergence
  if (redo2) stop('Convergence not reached')
  # if convergence:
  if (!redo2 & hessian) {
    if (trace) cat('\n'); cat('Computing hessian...')
    requireNamespace('numDeriv')
    logl.hess = function(theta) {
      ll = nlogl(theta, offset = offset, y = y, X = X, Z = Z, 
                 id = id, GHk = npoints)
      return(ll)
    }
    hess = numDeriv::hessian(func = logl.hess, x = mle$par)
    if (trace) cat('... done')
  }
  # collect results:
  mle.est = mle$par
  ncov = length(mle.est) - 3
  mle.est[ncov + 1] = 1 + exp(mle.est[ncov + 1])
  mle.est[ncov + 2] = 1 - exp(mle.est[ncov + 2])
  mle.est[ncov + 3] = exp(mle.est[ncov + 3])
  names(mle.est)[1:3 + ncov] = c('D', 'a', 'sigma2')
  logl = - mle$value
  if (!redo2 & hessian == F) out = list('mle' = mle.est, 'logl' = logl,
                                        'convergence' = mle$convergence)
  if (!redo2 & hessian) out = list('mle' = mle.est, 'logl' = logl, 'hessian' = hess,
                                   'convergence' = mle$convergence)
  class(out)=c('ptmm', 'list')
  return(out)
}