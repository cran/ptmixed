#' Negative binomial generalized linear model
#'
#' Estimates a negative binomial generalized linear model.
#'
#' Maximum likelihood estimation of a negative binomial GLM
#' (the NB distribution is obtained as special case of the Poisson-Tweedie distribution when a = 0).
#' 
#' @param formula A formula for the fixed effects part of the model. It should be in the form \code{y ~ x1 + x2}
#' @param offset An offset to be added to the linear predictor. Default is \code{NULL}.
#' @param data A data frame containing the variables declared in \code{formula}.
#' @param maxit Vector containing the maximum number of iterations used in optim by
#' the BFGS method and, if this fails, by the Nelder-Mead method
#' @param trace Logical value. If \code{TRUE}, additional information is printed during the optimization. Default is \code{TRUE}.
#' @param theta.start Numeric vector comprising initial parameter values for the
#' vector of regression coefficients and the dispersion parameter
#' @return A list containing the following elements: function's call (\code{call}); 
#' maximum likelihood estimate (\code{mle});  value of the
#' loglikelihood at the mle (\code{logl}); \code{convergence} value (if 0, the optimization converged);
#' the observed Fisher information (\code{fisher.info}) and the starting values
#' used in the optimization (\code{theta.init})
#' @export
#' @author Mirko Signorelli
#' @references Signorelli, M., Spitali, P., Tsonaka, R. (2021). Poisson-Tweedie 
#' mixed-effects model: a flexible approach for the analysis of longitudinal RNA-seq
#' data. Statistical Modelling, 21 (6), 520-545. URL: https://doi.org/10.1177/1471082X20936017
#' @seealso \code{\link{ptmixed}} for the Poisson-Tweedie GLMM
#' @examples
#' data(df1, package = 'ptmixed')
#' 
#' # estimate the model
#' fit1 = nbglm(formula = y ~ group*time, data = df1)
#' 
#' # view model summary
#' summary(fit1)

nbglm = function(formula, offset = NULL, data, maxit = c(500, 1e5), 
                 trace = T, theta.start = NULL) {
  call = match.call()
  # preliminary checks:
  if (length(formula) != 3) stop('formula should be in the form y ~ x1 + x2 +...')
  # identify elements
  y = data[, all.vars(formula[[2]])]
  X = model.matrix(formula[-2], data = data)
  check0 = missing(offset)
  if (check0) offset = NULL
  if (!check0) {
    offset = data[ , deparse(substitute(offset))]
    data$offset = offset
  }
  # starting values:
  if (is.null(theta.start)) {
    if (!is.null(offset)) {
      poi = glm(formula = formula, offset = offset,
                family = poisson(), data = data)
    }
    if (is.null(offset)) poi = glm(formula = formula, family = poisson(), 
                                    data = data)
    beta.init = coef(poi)
    D.init = 3
    theta.start = c(beta.init, D.init)
  }
  p = length(theta.start)
  theta.init = c(theta.start[1:(p-1)], log(theta.start[p]) - 1)
  # define negative loglik
  nlogl = function(theta, y, X, offset) {
    requireNamespace('tweeDEseq')
    p = length(theta)
    beta = theta[1:(p-1)]
    D = exp(theta[p]) + 1
    mu = exp(X %*% beta)
    if (!is.null(offset)) mu = exp(X %*% beta + offset)
    dens = mapply(tweeDEseq::dPT, x = y, mu = mu, D = D, a = 0)
    dens[which(dens < 1e-323)] = 1e-323 # otherwise it's too small for log(dens)
    ll = sum(log(dens))
    return(-ll)
  }
  # check that starting point is finite:
  try.init1 = try(nlogl(theta = theta.init, y = y, X = X, offset = offset))
  if (trace) {
    cat('Starting values: \n')
    cat(round(theta.start, 3)); cat('\n')
    cat(paste('initial loglik value:', round(-try.init1, 2))); cat('\n')
  }
  if (is.infinite(try.init1)) {
    try.init2 = Inf
    while(is.infinite(try.init2)) {
      new.D.init = runif(1, 1.5, 5)
      p = length(theta.init)
      theta.init[p] = log(new.D.init - 1)
      try.init2 = try(nlogl(theta = theta.init, y = y, X = X, offset = offset))
      cat('\n')
      cat(paste('retry: initial loglik value:', round(-try.init2, 2)))
    }
  }
  
  # optimization
  if (maxit[1] > 0) {
    ctrl = list(maxit = maxit[1])
    if (trace) ctrl = list(trace = 1, REPORT = 1, maxit = maxit[1])
    mle = try( optim(theta.init, nlogl, method = "BFGS", 
                   y = y, X = X, offset = offset,
                   control = ctrl) )
  }
  else if (maxit[1] == 0) cat('Skipping BFGS optimization\n')
  
  do.nm = F
  temp1 = exists('mle')
  if (!temp1) do.nm = T
  if (temp1) {
    temp2 = inherits(mle, 'try-error')
    if (temp2) do.nm = T
    else if (!temp2) {
      temp3 = (mle$convergence == 0)
      if (!temp3) do.nm = T
    }
  }
  if (do.nm) {
    if (trace & maxit[1] !=0) cat('Optimization with BFGS failed. Trying with Nelder-Mead...')
    ctrl = list(maxit = maxit[2])
    if (trace) ctrl = list(trace = 1, REPORT = 1, maxit = maxit[2])
    mle = try( optim(theta.init, nlogl, method = "Nelder-Mead", 
                     y = y, X = X, offset = offset,
                     control = ctrl) )
  }
  converged = F
  temp1 = exists('mle')
  if (temp1) {
    temp2 = inherits(mle, 'try-error')
    if (!temp2) {
      temp3 = (mle$convergence == 0)
      if (temp3) converged = T
    }
  }
  if (!converged) stop('Maximization failed')
  if (converged) {
    if (trace) cat('Convergence reached. Computing hessian...\n')
    requireNamespace('numDeriv')
    nlogl.hess = function(theta) {
      nll = nlogl(theta, offset = offset, y = y, X = X)
      return(nll)
    }
    hess = numDeriv::hessian(func = nlogl.hess, x = mle$par)
    if (trace) cat('... done\n')
  }
  
  # collect results:
  mle.est = mle$par
  ncov = length(mle.est) - 1
  mle.est[ncov + 1] = 1 + exp(mle.est[ncov + 1])
  mle.est = c(mle.est, 0)
  names(mle.est)[1:2 + ncov] = c('D', 'a')
  logl = - mle$value
  out = list('call' = call, 'mle' = mle.est, 
             'logl' = logl, 'convergence' = mle$convergence,
             'fisher.info' = hess, 'theta.init' = theta.start)
  class(out)=c('ptglm', 'list')
  return(out)
}
