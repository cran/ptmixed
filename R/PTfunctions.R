dpt.yvec.mumatr = function(yvec, mu.matr, D, a, log = T, tol = 1e-323) {
  requireNamespace('tweeDEseq')
  if (length(yvec) != dim(mu.matr)[1]) stop('length(yec) != dim(mumatr)[1]')
  nrows = dim(mu.matr)[1]
  ncols = dim(mu.matr)[2]
  out = matrix(NA, nrows, ncols)
  if (log == T) {
    for (j in 1:ncols) {
      out[,j] = mapply(PT.logdens, x=yvec, mu=mu.matr[,j], D = D, a = a, tol = tol)
    } 
  }
  else if (log == F) {
    for (j in 1:ncols) {
      out[,j] = mapply(tweeDEseq::dPT, yvec, mu.matr[,j], D = D, a = a, tol = tol)
    } 
  }
  return(out)
}

PT.logdens = function(x, mu, D, a, tol = 1e-323) {
  requireNamespace('tweeDEseq')
  out = tweeDEseq::dPT(x, mu, D, a, tol)
  out[out < 1e-323] = 1e-323
  return(log(out))
}