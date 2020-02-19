#' Negative binomial generalized linear mixed model
#'
#' Estimates the negative binomial generalized linear mixed model with random intercept
#' (here, the NB distribution is obtained as special case of the Poisson-Tweedie distribution when a = 0).
#' Likelihood approximation for the model is based on the adaptive Gauss-Hermite quadrature rule.
#' 
#' @param fixef.formula A formula for the fixed effects part of the model. It should be in the form \code{y ~ x1 + x2}
#' @param id A variable to distinguish observations from the same subject.
#' @param offset An offset to be added to the linear predictor. Default is \code{NULL}.
#' @param data A data frame containing the variables declared in \code{fixef.formula}.
#' @param npoints Number of quadrature points employed in the adaptive quadrature. Default is 5.
#' @param hessian Logical value. If \code{TRUE}, the hessian matrix is evaluated at the MLE to derive the observed Fisher information matrix. Default is \code{TRUE}.
#' @param trace Logical value. If \code{TRUE}, additional information is printed during the optimization. Default is \code{TRUE}.
#' @param theta.start Numeric vector comprising initial parameter values for the
#' vector of regression coefficients, the dispersion parameter (using the same parametrization of \code{ptmixed}) 
#' and the variance of the random intercept. 
#' Default is \code{NULL}: initial
#' parameter estimates are computed automatically by the function.
#' @param reltol Relative tolerance to be used in optim. Default to 1e-8
#' @param maxit Vector containing the maximum number of iterations used in optim by
#' the Nelder-Mead method and, if this fails, by the BFGS method
#' @param freq.updates Number of iterations after which the quadrature points are updated when the Nelder-Mead
#' algorithm is used for the optimization. Default value is 200. To update the quadrature points at every iteration 
#' (note that this may make the computation about 10x slower), set \code{freq.updates = 1} 
#' or \code{freq.updates = NA}. The function first tries to optimize the loglikelihood using the Nelder-Mead
#' algorithm, updating the quadrature points every \code{freq.updates} iterations. If this fails to converge,
#' a second attempt is made using the BFGS algorithm, for which the quadrature points are updated at every iteration.
#' @param min.var.init If the initial estimate of the variance of the random intercept is smaller than
#' this value, estimation is stopped and the user is advided to use the simpler Poisson-Tweedie GLM is used. Default is 1e-3.
#' @return A list containing the following elements: function's call (\code{call}); 
#' maximum likelihood estimate (\code{mle});  value of the
#' loglikelihood at the mle (\code{logl}); \code{convergence} value (if 0, the optimization converged);
#' the observed Fisher information (\code{fisher.info}), if \code{hessian = T}; the number of quadrature points 
#' used (\code{quad.points}) and the starting value used in the optimization (\code{theta.init}); 
#' relevant warnings (\code{warnings}).
#' @import stats
#' @importFrom utils tail
#' @export
#' @author Mirko Signorelli
#' @seealso \code{\link{summary.ptglmm}}, \code{\link{ranef}}
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
#' for (i in 1:(n*t)) y[i] = rPT(1, mu = mu[i], D = 2, a = 0, max = 1000)
#' 
#' data.long = data.frame(y, group, time, id, offset)
#' rm(list = setdiff(ls(), 'data.long'))
#' 
#' # 1) Quick example (5 quadrature points, hessian and SEs not computed)
#' # estimate the model
#' fit1 = nbmixed(fixef.formula = y ~ group + time, id = data.long$id,
#'               offset = data.long$offset, data = data.long, npoints = 5, 
#'               freq.updates = 200, hessian = FALSE, trace = TRUE)
#' # print summary:
#' summary(fit1, wald = FALSE)
#' 
#' \donttest{
#' # 2) Full computation, including hessian evaluation and using more quadrature points
#' # estimate the model
#' fit2 = nbmixed(fixef.formula = y ~ group + time, id = data.long$id,
#'               offset = data.long$offset, data = data.long, npoints = 10, 
#'               freq.updates = 200, hessian = TRUE, trace = TRUE)
#' # print and get summary:
#' results = summary(fit2, wald = TRUE)
#' ls(results)
#' # view table with estimates of regression coefficients, standard errors and Wald test:
#' results$coefficients
#' }

nbmixed = function(fixef.formula, id, offset = NULL,
                   data, npoints = 5, hessian = T, trace = T,
                   theta.start = NULL, reltol = 1e-8, maxit = c(1e4, 100),
                   freq.updates = 200, min.var.init = 1e-3) {
  call = match.call()
  warning.list = list()
  # preliminary checks:
  if (length(fixef.formula) != 3) stop('fixef.formula should be in the form y ~ x1 + x2 +...')
  t = dim(data)[1] / length(unique(id))
  if (t %% 1 !=0) stop('The dataset appears to be unbalanced. The code for
                       the unbalanced case is not implemented yet')
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
  # NB: updates at every iteration coded as NA (but it can also be provided = 1)
  if (freq.updates == 1) freq.updates = NA
  # optim control values:
  optim.control.nm = list(maxit = maxit[1], reltol = reltol)
  optim.control.bfgs = list(maxit = maxit[2], reltol = reltol)
  if (!is.na(freq.updates)) optim.control.nm = list(maxit = freq.updates, reltol = reltol)
  if (trace) {
    optim.control.nm = list(trace = 1, REPORT = 1, maxit = maxit[1], reltol = reltol)
    optim.control.bfgs = list(trace = 1, REPORT = 1, maxit = maxit[2], reltol = reltol)
    if (!is.na(freq.updates)) optim.control.nm = list(trace = 1, REPORT = 1, maxit = freq.updates, reltol = reltol)
  }
  # starting values:
  if (is.null(theta.start)) {
    fixef.glmmad = fixef.formula
    if (!is.null(offset)) {
      fixef.glmmad = update.formula(fixef.glmmad, ~ . + offset(offset))
      df.glmmad = cbind(data, 'id' = id, 'offset' = offset)
    }
    else df.glmmad = cbind(data, 'id' = id)
    temp = get.initial.theta(fixef.formula = fixef.glmmad, data = df.glmmad, 
                             y = y, id = id)
    q = length(temp$theta.init)
    theta.init = temp$theta.init[-(q-1)]
    warning.list = temp$warnings
  }
  else if (!is.null(theta.start)) {
    q = length(theta.start)
    p = dim(X)[2]
    if (q != p+2) stop('Wrong number of elements in theta.start')
    if (theta.start[p+1] <= 1) stop('Initial dispersion value should be > 1')
    theta.init = c(theta.start[1:p], log(theta.start[p+1]-1),
                   log(theta.start[p+2]))
  }
  # rename correctly elements in theta.init
  names(theta.init) = c(colnames(X), 'D', 'sigma2')
  # do a check on the initial starting value of the variance:
  include.ranef = T
  if (exp(tail(theta.init,1)) < min.var.init) {
    warning('Initial variance value < min.var.init. Estimating a GLM instead of a GLMM')
    warning.list = c(warning.list, 
    'Initial variance value < min.var.init. A GLM was estimated (hence sigma2 = 0)')
    include.ranef = F
  }
  # estimate NB glm
  if (!include.ranef) {
    mle = nbglm(fixef.formula, offset = offset, 
              data = data, maxit = c(500, 1e5), 
              trace = trace, theta.start = NULL)
  }
  # estimate NB glmm
  if (include.ranef) {
    # negative loglikelihood:
    nlogl = function(theta, offset, y, X, Z, id, GHk, GHs = NULL) {
      p = length(theta)
      beta = theta[1:(p-2)]
      D = 1+exp(theta[p-1])
      Sigma = as.matrix(exp(theta[p]))
      ll = loglik.pt.1re(beta, D, a = 0, Sigma, y, X, Z, id, 
                         offset = offset, GHk, GHs = GHs)
      return(-ll)
    }
    
    redo = redo2 = F
    # optimization: try Nelder-Mead
    ##stop = F
    ##niter = 0
    # starting values
    p = length(theta.init)
    theta.curr = theta.init
    beta.curr = theta.init[1:(p-2)]
    D.curr = 1+exp(theta.init[p-1])
    S.curr = as.matrix(exp(theta.init[p]))
    GHs = NULL
    if (trace) {
      cat('\n')
      cat(paste('Initial D =', round(D.curr, 2), 'S =', round(S.curr, 3)))
      cat('\n')
      cat('Beginning optimization with Nelder-Mead:'); cat('\n')
    }
    
    # optimization updating quadrature points in NM
    # every k iterations
    if (!is.na(freq.updates)) {
      stop = F
      niter = 0
      while (stop == F) {
        GH.up = try( GHpoints.pt.1re(y = y, id = id, X = X, Z = Z, 
                                   beta = beta.curr, D = D.curr, a = 0, Sigma = S.curr, 
                                   offset = offset, RE.size = 1, GHk = npoints, tol = 1e-323),
                    silent = T)
      
        if (exists('GH.up')) {
          if (!inherits(GH.up, 'try-error')) GHs = GH.up
          else warning('Quadrature points not updated')
        }
        mle = try( optim(theta.curr, nlogl, method = "Nelder-Mead", offset = offset,
                       y = y, X = X, Z = Z, id = id, GHk = npoints, hessian = F,
                       GHs = GHs, control = optim.control.nm) )
        # check convergence
        temp1 = exists('mle')
        if (!temp1) { # proceed with BFGS
          stop = T
          redo = T
        }
        else if (temp1) {
          temp2 = inherits(mle, 'try-error')
          if (temp2) { # proceed with BFGS
            stop = T
            redo = T
          }
          else if (!temp2) { # check if convergence reached
            niter = niter + mle$counts[1]
            temp3 = (mle$convergence == 0)
            if (temp3) stop = T
            else if (niter > maxit[1]) {
              stop = T
              redo = T
            }
            else {
              theta.curr = mle$par
              p = length(theta.curr)
              beta.curr = theta.curr[1:(p-2)]
              D.curr = 1+exp(theta.curr[p-1])
              S.curr = as.matrix(exp(theta.curr[p]))
              if (trace) print(paste('D =', round(D.curr, 2), 'S =', round(S.curr, 2)))
            }
          }
        }
      }
      cat('\n')
      cat(paste('Total number of iterations =', niter))
      cat('\n')
    }
    
    # optimization updating quadrature points at every
    # iteration
    if (is.na(freq.updates)) {
      mle = try( optim(theta.curr, nlogl, method = "Nelder-Mead", offset = offset,
                       y = y, X = X, Z = Z, id = id, GHk = npoints, hessian = F,
                       GHs = NULL, control = optim.control.nm) )
      # check convergence
      redo = T
      temp1 = exists('mle')
      if (temp1) {
        temp2 = inherits(mle, 'try-error')
        if (!temp2) {
          temp3 = (mle$convergence == 0)
          if (temp3) {
            redo = F
            theta.curr = mle$par
            p = length(theta.curr)
            beta.curr = theta.curr[1:(p-2)]
            D.curr = 1+exp(theta.curr[p-1])
            S.curr = as.matrix(exp(theta.curr[p]))
          }
        }
      }
    }
    
    # if it failed: try BFGS
    if (redo) {
      print('Optimization with Nelder-Mead was not successful. Trying BFGS...'); cat('\n')
      if (trace) cat('Beginning optimization with BFGS:'); cat('\n')
      mle = try( optim(theta.init, nlogl, method = "BFGS", offset = offset,
                       y = y, X = X, Z = Z, id = id, GHk = npoints, hessian = F,
                       control = optim.control.bfgs, GHs = NULL) )
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
      if (trace) cat('\n'); cat('Convergence reached. Computing hessian...'); cat('\n')
      requireNamespace('numDeriv')
      nlogl.hess = function(theta) {
        nll = nlogl(theta, offset = offset, y = y, X = X, Z = Z, 
                    id = id, GHk = npoints)
        return(nll)
      }
      hess = numDeriv::hessian(func = nlogl.hess, x = mle$par)
      if (trace) cat('... done')
    }
  }
  
  # collect results:
  if (include.ranef) {
    mle.est = mle$par
    ncov = length(mle.est) - 2
    mle.est[ncov + 1] = 1 + exp(mle$par[ncov + 1])
    mle.est[ncov + 2] = 0
    mle.est[ncov + 3] = exp(mle$par[ncov + 2])
    logl = - mle$value
  }
  if (!include.ranef) {
    mle.est = mle$mle
    mle.est = c(mle.est, 0)
    logl = mle$logl
  }
  ncov = length(mle.est)-3
  names(mle.est)[1:3 + ncov] = c('D', 'a', 'sigma2')
  # starting value:
  if (include.ranef) {
    # transform back initial starting value
    q = length(theta.init)
    theta.init[q-1] = 1+exp(theta.init[q-2])
    theta.init[q] = exp(theta.init[q])
  }
  if (!include.ranef) {
    theta.init = mle$theta.init
  }
  # output list:
  if (include.ranef) {
    if (!redo2 & hessian == F) out = list('call' = call,
                                          'mle' = mle.est, 'logl' = logl,
                                          'convergence' = mle$convergence, 'quad.points' = npoints,
                                          'theta.init' = theta.init, 'warnings' = warning.list)
    if (!redo2 & hessian) out = list('call' = call,
                                     'mle' = mle.est, 'logl' = logl, 'fisher.info' = hess,
                                     'convergence' = mle$convergence, 'quad.points' = npoints,
                                     'theta.init' = theta.init, 'warnings' = warning.list)
  }
  if (!include.ranef) {
    out = list('call' = call,
               'mle' = mle.est, 'logl' = logl, 'fisher.info' = mle$fisher.info,
               'convergence' = mle$convergence, 'quad.points' = NA,
               'theta.init' = mle$theta.init, 'warnings' = warning.list)
  }
  class(out)=c('ptglmm', 'list')
  return(out)
}
