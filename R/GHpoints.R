GHpoints.pt.1re <- function (y, id, X, Z, beta, D, a, Sigma, offset = NULL, RE.size = 1, GHk = 10, tol = 1e-323){
  requireNamespace('tweeDEseq')
  #require(mvtnorm)
  with.offset = !is.null(offset)
  if (with.offset) Delta.ij <- c(X %*% beta + c(offset)) #check this is ok
  else Delta.ij <- c(X %*% beta)
  n <- length(unique(id))
  GH <- gauher(GHk)
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
  scaled.b <- vector("list", n)
  inv.chol <- rep(NA, n)
  det.Bs <- numeric(n)
  for (i in 1:n) {
    #print(i)
    id.i <- id == i
    opt <- try( optim(rep(0, RE.size), fn, gr, y.i = y[id.i], delta.ij = Delta.ij[id.i], 
                      Z.ij = as.matrix(Z[id.i,]), tol = tol, method = "BFGS", #control = list(trace = 1),
                      hessian = TRUE), silent = F )
    chol = try(chol(opt$hessian))
    if (!inherits(chol, 'try-error')) {
      inv.chol[i] <- solve(chol(opt$hessian))
      scaled.b[[i]] <- t(opt$par + tcrossprod(inv.chol[i], b))
      det.Bs[i] <- det(as.matrix(inv.chol[i]))
    }
  }
  no.chol = which(is.na(inv.chol))
  if (length(no.chol >= 1)) {
    avg.chol = mean(inv.chol, na.rm = T)
    avg.det = mean(det.Bs, na.rm = T)
    for (j in no.chol) {
      scaled.b[[j]] <- t(opt$par + tcrossprod(avg.chol, b))
      det.Bs[j] = avg.det
    }
  }
  list(scaled.b = scaled.b, det.Bs = det.Bs, wGH = wGH, unscaled.b = b, 
       inv.chol = inv.chol, chol.diagnostic = length(no.chol))
}