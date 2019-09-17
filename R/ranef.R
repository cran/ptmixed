#' Compute random effects for Poisson-Tweedie mixed model
#'
#' Compute the empirical Bayes estimates of the random effects for 
#' the Poisson-Tweedie generalized linear mixed model
#' 
#' @param obj an object of class \code{ptglmm} (obtained from \code{ptmixed}).
#' @return A vector with the EB estimates of the random effects
#' @export
#' @author Mirko Signorelli
#' @seealso \code{\link{ptmixed}}
#' @examples
#' \donttest{
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
#' fit1 = ptmixed(fixef.formula = y ~ group + time, id = data.long$id,
#'               offset = data.long$offset, data = data.long, npoints = 5, 
#'               freq.updates = 200, hessian = FALSE, trace = TRUE)
#'               
#' # obtain random effect estimates
#' ranef(obj = fit1)
#' }

ranef <- function (obj) {
  if (class(obj)[1] != 'ptglmm') stop('obj should be of class ptglmm')
  # extract arguments from ptglmm object call and mle
  # get y, id, X, Z, offset:
  temp = as.list(obj$call)
  fixef.formula = temp$fixef.formula
  df = eval(temp$data)
  y = df[, all.vars(fixef.formula[[2]])]
  X = model.matrix(as.formula(fixef.formula[-2]), data = df)
  Z = as.matrix(model.matrix(~1, data = df))
  id = eval(temp$id)
  id = as.numeric(as.factor(id))
  offset = eval(temp$offset)
  # get beta, D, a, Sigma
  mle = obj$mle
  p = length(mle) - 3
  beta = mle[1:p]
  D = mle[p+1]
  a = mle[p+2]
  Sigma = mle[p+3]
  # other arguments
  RE.size = 1
  GHk = temp$npoints
  tol = 1e-323
  
  requireNamespace('tweeDEseq')
  #require(mvtnorm)
  with.offset = !is.null(offset)
  if (with.offset) Delta.ij <- c(X %*% beta + c(offset))
  else Delta.ij <- c(X %*% beta)
  n <- length(unique(id))
  GH <- gauher(GHk) # returns ascissae and weights for Gauss-Hermite quadrature
  b <- as.matrix(expand.grid(rep(list(GH$x), RE.size)), drop = F)
  wGH <- as.matrix(expand.grid(rep(list(GH$w), RE.size)), drop = F)
  wGH <- 2^(RE.size/2) * apply(wGH, 1, prod) * exp(rowSums(b * b))
  b <- sqrt(2) * b
  ###########
  fn <- function (b, y.i, delta.ij, Z.ij, tol) {
    log.p.b <- dnorm(b, 0, sd = sqrt(Sigma), log = T)
    p.yb = mapply(PT.logdens, x = y.i, mu = exp(delta.ij + Z.ij %*% b), 
                  D = D, a = a, tol = tol)
    log.p.yb <- sum(p.yb)
    - log.p.yb - log.p.b
  }
  gr <- function (b, y.i, delta.ij, Z.ij, tol) {
    cd(b, fn, y.i = y.i, delta.ij = delta.ij, Z.ij = Z.ij, tol = tol)
  }
  
  out = rep(NA, n)
  for (i in 1:n) {
    #print(i)
    id.i <- id == i
    opt <- try( optim(rep(0, RE.size), fn = fn, gr = gr,
                      y.i = y[id.i], delta.ij = Delta.ij[id.i], 
                      Z.ij = as.matrix(Z[id.i,]), tol = tol, method = "BFGS",
                      hessian = FALSE), silent = F )
    out[i] = opt$par
  }
  names(out) = unique(id)
  return(out)
}