get.initial.theta = function(fixef.formula, data, y, id, time) {
  # NB: the offset should be added as offset(log.offset) in fixef.formula!
  requireNamespace('GLMMadaptive')
  # start with Poisson glmm:
  poi.glmm = try(GLMMadaptive::mixed_model(fixed = fixef.formula, random = ~ 1 | id, data = data, 
                         family = poisson(), nAGQ=21,
                         initial_values = list(betas = poisson())),
                 silent = T)
  # try NB glmm
  nb.glmm = try(
    GLMMadaptive::mixed_model(fixed = fixef.formula, random = ~ 1 | id, data = data, 
                family = negbin.family(), n_phis = 1, nAGQ=21,
                initial_values = list("betas" = poi.glmm$coefficients,
                                      "D" = poi.glmm$D)),
    silent = T)
  # initial beta and variance:
  from.poi = from.nb = F
  if (!inherits(poi.glmm, 'try-error')) {
    if (poi.glmm$converged) from.poi = T
  }
  if (!inherits(nb.glmm, 'try-error')) {
    if (nb.glmm$converged) from.nb = T
  }
  if (from.nb) {
    beta.init = nb.glmm$coefficients
    S.init = as.numeric(nb.glmm$D)
  }
  if (!from.nb) {
    if (from.poi) {
      beta.init = poi.glmm$coefficients
      S.init = as.numeric(poi.glmm$D)
    }
    else if (!from.poi){ # try lme4 Poisson GLMM or Poisson GLM
      requireNamespace('lme4')
      full.formula = update.formula(fixef.formula, '~ . + (1|id)')
      poi.lme4 = try(lme4::glmer(full.formula, data = data, family = 'poisson',
                       nAGQ = 21), silent = T)
      from.lme4 = F
      if (!inherits(poi.lme4, 'try-error')) {
        beta.init = lme4::fixef(poi.lme4)
        S.init = as.data.frame(lme4::VarCorr(poi.lme4))[1,'sdcor'] ^2
      }
      if (inherits(poi.lme4, 'try-error')) {
        poi.glm = glm(fixef.formula, family = poisson, data = data)
        beta.init = coef(poi.glm)
        S.init = 0.2
      }
    }
  }
  # deviance and power:
  D.init = harm.cond.dev(y, id, time)
  if (D.init <= 1) D.init = geom.cond.dev(y, id, time)
  if (D.init <= 1) D.init = conditional.deviance(y, id, time)
  if (D.init <= 1) D.init = 2
  a.init = a.moment.estimator(y)
  if (a.init >= 1) a.init = 0.5
  D.trasf.init = log(D.init-1)
  a.trasf.init = log(1-a.init)
  S.trasf.init = log(S.init)
  theta.init = c(beta.init, D.trasf.init, a.trasf.init, S.trasf.init)
  return(theta.init)
}

negbin.family <- function () {
  # taken from vignettes of GLMMadaptive package
  stats <- make.link(link = "log")
  log_dens <- function (y, eta, mu_fun, phis, eta_zi) {
    # the log density function
    phis <- exp(phis)
    mu <- mu_fun(eta)
    log_mu_phis <- log(mu + phis)
    comp1 <- lgamma(y + phis) - lgamma(phis) - lgamma(y + 1)
    comp2 <- phis * log(phis) - phis * log_mu_phis
    comp3 <- y * log(mu) - y * log_mu_phis
    out <- comp1 + comp2 + comp3
    attr(out, "mu_y") <- mu
    out
  }
  score_eta_fun <- function (y, mu, phis, eta_zi) {
    # the derivative of the log density w.r.t. mu
    phis <- exp(phis)
    mu_phis <- mu + phis
    comp2 <- - phis / mu_phis
    comp3 <- y / mu - y / mu_phis
    # the derivative of mu w.r.t. eta (this depends on the chosen link function)
    mu.eta <- mu
    (comp2 + comp3) * mu.eta
  }
  score_phis_fun <- function (y, mu, phis, eta_zi) {
    # the derivative of the log density w.r.t. phis
    phis <- exp(phis)
    mu_phis <- mu + phis
    comp1 <- digamma(y + phis) - digamma(phis)
    comp2 <- log(phis) + 1 - log(mu_phis) - phis / mu_phis
    comp3 <- - y / mu_phis
    (comp1 + comp2 + comp3) * phis
  }
  structure(list(family = "user Neg Binom", link = stats$name, linkfun = stats$linkfun,
                 linkinv = stats$linkinv, log_dens = log_dens,
                 score_eta_fun = score_eta_fun, score_phis_fun = score_phis_fun),
            class = "family")
}

conditional.deviance = function(x, id, time) {
  df.long = data.frame(x, id, 'time' = time+1)
  df.wide = reshape(df.long, idvar = 'id', timevar = 'time', direction = 'wide')[,-1]
  deviance.within = apply(df.wide, 1, function(x) var(x)/mean(x))
  return(mean(deviance.within, na.rm = T))
}

geom.cond.dev = function(x, id, time) {
  df.long = data.frame(x, id, time = time + 1)
  df.wide = reshape(df.long, idvar = "id", timevar = "time", 
                    direction = "wide")[, -1]
  deviance.within = apply(df.wide, 1, function(x) var(x)/mean(x))
  #return(mean(deviance.within, na.rm = T))
  return(exp(mean(log(deviance.within), na.rm = T)))
}

harm.cond.dev = function(x, id, time) {
  df.long = data.frame(x, id, time = time + 1)
  df.wide = reshape(df.long, idvar = "id", timevar = "time", 
                    direction = "wide")[, -1]
  deviance.within = apply(df.wide, 1, function(x) var(x)/mean(x))
  #return(mean(deviance.within, na.rm = T))
  return(1/(mean(1/deviance.within, na.rm = T)))
}

a.moment.estimator = function(x) {
  requireNamespace('moments')
  mean = mean(x)
  dev = var(x)/mean(x)
  skew = moments::skewness(x)
  r = skew * sqrt(mean * dev^3)
  s = (dev-1)^2
  t = 3*dev - 2
  a.est = 1 + s / (s+t-r)
  return(a.est)
}
