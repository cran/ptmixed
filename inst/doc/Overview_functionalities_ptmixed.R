## ---- eval=FALSE, echo=TRUE, results='asis'-----------------------------------
#  install.packages('ptmixed')

## ---- eval=TRUE, echo=TRUE, results='asis'------------------------------------
library('ptmixed')

## ---- eval=T, echo=T, results='markup'----------------------------------------
example.df = simulate_ptglmm(n = 14, t = 4, seed = 1234,
                             beta = c(2.3, -0.9, -0.2, 0.5),
                             D = 1.5, a = -1,
                             sigma2 = 0.7)
data.long = example.df$data
head(data.long)

## ---- eval=T, echo=T, results='markup', fig.height=3.5, fig.width=5-----------
pmf(data.long$y, xlab = 'y', title = 'Distribution of y')
make.spaghetti(x = time, y = y, id = id,
  group = group, data = data.long,
  title = 'Trajectory ("spaghetti") plot',
  legend.title = 'GROUP')

## ---- eval=T, echo=T, results='markup'----------------------------------------
pt_glmm = ptmixed(fixef.formula = y ~ group*time, id = id,
                     data = data.long, npoints = 3, 
                     hessian = T, trace = F)

## ---- eval=T, echo=T, results='markup'----------------------------------------
summary(pt_glmm)

## ---- eval=T, echo=T, results='markup'----------------------------------------
L.group = matrix(0, nrow = 2, ncol = 4)
L.group[1, 2] = L.group[2, 4] = 1
L.group

## ---- eval=T, echo=T, results='markup'----------------------------------------
wald.test(pt_glmm, L = L.group, k = c(0, 0))

## ---- eval=T, echo=T, results='markup'----------------------------------------
null_model = ptmixed(fixef.formula = y ~ time, id = id,
                               data = data.long, npoints = 3, 
                               hessian = F, trace = F)

## ---- eval=T, echo=T, results='markup'----------------------------------------
lrt.stat = 2*(pt_glmm$logl - null_model$logl)
lrt.stat
p.lrt = pchisq(lrt.stat, df = 2, lower.tail = F)
p.lrt

## ---- eval=T, echo=T, results='markup'----------------------------------------
ranef(pt_glmm)

## ---- eval=T, echo=T, results='markup'----------------------------------------
nb_glmm = nbmixed(fixef.formula = y ~ group*time, id = id,
                     data = data.long, npoints = 3, 
                     hessian = T, trace = F)

## ---- eval=T, echo=T, results='markup'----------------------------------------
summary(nb_glmm)
ranef(nb_glmm)

## ---- eval=T, echo=T, results='markup'----------------------------------------
pt_glm = ptglm(formula = y ~ group*time, data = data.long, trace = F)
summary(pt_glm)

## ---- eval=T, echo=T, results='markup'----------------------------------------
nb_glm = nbglm(formula = y ~ group*time, data = data.long, trace = F)
summary(nb_glm)

