#' Simulate data from the Poisson-Tweedie generalized linear mixed model
#'
#' Simulates a dataset comprising t repeated measurements for n subjects
#' from a Poisson-Tweedie GLMM. Subjects are assumed to belong to two
#' different groups. The linear predictor comprises an intercept,
#' main effects of group and of time, and the interaction between time and
#' group; a random intercept; and, optionally, a normally-distributed 
#' offset term.
#' 
#' @param n number of subjects
#' @param t number of time points (0, 1, ..., t-1)
#' @param seed seed for random number generation
#' @param beta vector of regression coefficients, to be specified in
#' this order: intercept, group main effect, time main effect, group*time interaction
#' @param D dispersion parameter of the Poisson-Tweedie distribution (D > 1)
#' @param a power parameter of the Poisson-Tweedie distribution (a < 1)
#' @param sigma2 Variance of the subject-specific random intercept
#' @param offset Logical value. If \code{TRUE}, a normally-distributed offset is added to the linear predictor.
#' @return A list containing the following elements: a dataframe (\code{data})
#' containing the response y, the subject id, the group indicator and time; 
#' a vector with the true random intercept values (\code{true.randint}).
#' @export
#' @author Mirko Signorelli
#' @references Signorelli, M., Spitali, P., Tsonaka, R. (2021). Poisson-Tweedie 
#' mixed-effects model: a flexible approach for the analysis of longitudinal RNA-seq
#' data. Statistical Modelling, 21 (6), 520-545. URL: https://doi.org/10.1177/1471082X20936017
#' @examples
#' # simulate a simple, small dataset
#' example1 = simulate_ptglmm(n = 5, t = 2)
#' example1$data
#' # the function allows to set several different parameters
#' example2 = simulate_ptglmm(n = 20, t = 5, seed = 1,
#' beta = c(2.2, 1.2, 0.3, -0.5), D = 1.8, a = 0.5,
#' sigma2 = 0.7, offset = TRUE)
#' # view the distribution of the response variable:
#' pmf(example2$data$y)
#' # visualize the data with a trajectory plot:
#' make.spaghetti(x = time, y = y, id = id,
#' group = group, data = example2$data)

simulate_ptglmm = function(n = 20, t = 5, seed = 1, 
                           beta = c(3, 0, 0, 0.4), D = 1.5, a = -1, 
                           sigma2 = 0.8^2, offset = F) {
  requireNamespace('tweeDEseq')
  if (D < 1) stop('D should be >= 1')
  if (a > 1) stop('D should be <= 1')
  if (a == 1 & D > 1) stop('When a = 1 (Poisson distribution), D can only be = 1')
  if (sigma2 < 0) stop('sigma2 should be > 0')
  warn.text = 'is close to the boundary of the PT GLMM. This typically makes
  estimation of the PT GLMM challenging. Proceed with care.'
  if (D <= 1.1) warning(paste('D', warn.text))
  if (a >= 0.9) warning(paste('a', warn.text))
  if (sigma2 < 0.1) warning(paste('sigma2', warn.text))
  id = rep(1:n, each = t)
  n.ceil = ceiling(n/2)
  group = ifelse(id <= n.ceil, 0, 1)
  time = rep(0:(t-1), n)
  if (offset) off = rnorm(n*t, sd = 0.5)
  set.seed(seed)
  y = rep(NA, n*t)
  X = model.matrix(~ group + time + time*group)
  rand.int = rep(rnorm(n, sd = 0.8), each = t)
  eta = X %*% beta + rand.int
  if (offset) eta = eta + off
  for (i in 1:(n*t)) y[i] = tweeDEseq::rPT(1, mu = exp(eta[i]), 
                                D = D, a = a, max = 1000)
  out1 = data.frame(y, id, group, time)
  if (!offset) out1 = data.frame(y, id, group, time)
  if (offset) out1 = data.frame(y, id, group, time, 'offset' = off)
  out = list('data' = out1, 'true.randint' = unique(rand.int))
  return(out)
}