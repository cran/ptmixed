get.initial.theta = function(fixef.formula, data, y, id) {
  # NB: the offset should be added as offset(log.offset) in fixef.formula!
  # warning messages
  warn1 = 'A preliminary fit of a negative binomial mixed model indicates that a Poisson mixed model would fit the data better than the negative binomial one'
  warn2 = 'The response variable may be underdispersed. You may want to consider using either a Poisson mixed model or a model suitable for an underdispersed dependent response'
  warning.list = list()
  requireNamespace('GLMMadaptive')
  a.init = a.moment.estimator(y)
  if (a.init >= 0.8) {
    a.init = 0.5
  }
  # start with Poisson glmm:
  requireNamespace('lme4')
  full.formula = update.formula(fixef.formula, '~ . + (1|id)')
  poi.lme4 = try(lme4::glmer(full.formula, data = data, family = 'poisson',
                             nAGQ = 21), silent = T)
  if (!inherits(poi.lme4, 'try-error')) {
    nb.glmm = try(
      GLMMadaptive::mixed_model(fixed = fixef.formula, random = ~ 1 | id, 
                                data = data, family = GLMMadaptive::negative.binomial(), 
                                n_phis = 1, nAGQ=21,
                                initial_values = list("betas" = lme4::fixef(poi.lme4)),
                                control = list(iter_EM = 100, iter_qN = 100)),
      silent = T)
  }
  if (inherits(poi.lme4, 'try-error')) {
    poi.glmm = try(
      GLMMadaptive::mixed_model(fixed = fixef.formula, random = ~ 1 | id, 
                                data = data, family = poisson(), nAGQ=21,
                                initial_values = list(betas = poisson()),
                                control = list(iter_EM = 100, iter_qN = 100)),
      silent = T)
    nb.glmm = try(
      GLMMadaptive::mixed_model(fixed = fixef.formula, random = ~ 1 | id, 
                                data = data, family = GLMMadaptive::negative.binomial(), 
                                n_phis = 1, nAGQ=21,
                                initial_values = list("betas" = poi.glmm$coefficients),
                                control = list(iter_EM = 100, iter_qN = 100)),
      silent = T)
  }

  if (inherits(nb.glmm, 'try-error')) {
    if (nb.glmm == "Error in mixed_fit(y, X, Z, X_zi, Z_zi, id, offset, offset_zi, family,  : \n  A value greater than 22000 has been detected for the shape/size\n parameter of the negative binomial distribution. This typically\n indicates that the Poisson model would be better. Otherwise,\n adjust the 'max_phis_value' control argument.\n") {
      warning(warn1)
      warning.list = c(warning.list, 'A Poisson GLMM could be enough')
    }
  }
  
  # initial beta and variance:
  from.poi = from.nb = F
  if (!inherits(poi.lme4, 'try-error')) {
    if (poi.lme4@optinfo$conv$opt == 0) from.poi = T
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
      if (!inherits(poi.lme4, 'try-error')) { #lme4
        beta.init = lme4::fixef(poi.lme4)
        S.init = as.data.frame(lme4::VarCorr(poi.lme4))[1,'sdcor'] ^2
      }
      else if (!inherits(poi.glmm, 'try-error')) { #GLMMadaptive
        beta.init = poi.glmm$coefficients
        S.init = as.numeric(poi.glmm$D)
        }
      else { # try Poisson GLM
        poi.glm = glm(fixef.formula, family = poisson, data = data)
        beta.init = coef(poi.glm)
        S.init = 0.2
      }  
    }
  }
  # dispersion:
  temp = conditional.dispersion(x = y, id)
  D.init = temp$harmonic
  if (D.init <= 1.1) D.init = temp$geometric
  if (D.init <= 1.1) D.init = temp$arithmetic
  if (D.init <= 1.1) {
    warning(warn2)
    warning.list = c(warning.list, 'Data may be underdispersed')
    D.init = 1.5
  }
  # check on extreme values of D, s2 (from v 0.3.2)
  if (D.init > 20) D.init = 20
  if (S.init > 9) S.init = 5
  # transformation of D, a, S:
  D.trasf.init = log(D.init-1)
  a.trasf.init = log(1-a.init)
  S.trasf.init = log(S.init)
  theta.init = c(beta.init, D.trasf.init, a.trasf.init, S.trasf.init)
  out = list('theta.init' = theta.init, 'warnings' = warning.list)
  return(out)
}

a.moment.estimator = function(x) {
  requireNamespace('moments')
  mean = mean(x)
  disp = var(x)/mean(x)
  skew = moments::skewness(x)
  r = skew * sqrt(mean * disp^3)
  s = (disp-1)^2
  t = 3*disp - 2
  a.est = 1 + s / (s+t-r)
  return(a.est)
}

conditional.dispersion = function(x, id) {
  df.long = data.frame(x, id)
  df.long$seq = NA
  ids = unique(id)
  n = length(ids)
  for (i in 1:n) {
    rows = which(df.long$id == ids[i])
    mi = length(rows)
    df.long$seq[rows] = seq(1:mi)
  }
  df.wide = reshape(df.long, idvar = 'id', timevar = 'seq', direction = 'wide')[,-1]
  dispersion.within = apply(df.wide, 1, function(x) var(x)/mean(x))
  arit.mean = mean(dispersion.within, na.rm = T)
  geom.mean = exp(mean(log(dispersion.within), na.rm = T))
  harmon.mean = 1/(mean(1/dispersion.within, na.rm = T))
  out = list('arithmetic' = arit.mean, 'geometric' = geom.mean,
             'harmonic' = harmon.mean)
  return(out)
}
