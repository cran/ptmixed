# Introduction

The aim of this document is to keep track of the changes made to the
different versions of the `R` package `ptmixed`.

The numbering of package versions follows the convention a.b.c, where a
and b are non-negative integers and c is a positive integer. When minor
changes are made to the package, a and b are kept fixed and only c is
increased. Major changes to the package, instead, are made apparent by
changing a or b.

Each section of this document corresponds to a major change in the
package - in other words, within a section you will find all those
package versions a.b.x where a and b are fixed whereas x = 1, 2, 3, …
Each subsection corresponds to a specific package version.

# 1.1.x

## ptmixed 1.1.2

-   Released: April 2022
-   Removed dependency on `aod` package (which is scheduled to be
    archived by `CRAN`)

## ptmixed 1.1.1

-   Released: December 2021
-   Fixed problem with `NEWS` file, which was not visible on `CRAN` any
    more

## ptmixed 1.1.0

-   Released: November 2021
-   Updated citation information with details of final (issued) version
    of the article published in *Statistical Modelling* in: package
    description, citation file, help pages and vignettes
-   Fixed bug in `make.spaghetti()` function (rows with `NA`s on either
    `x` or `y` do not cause problems any more)

# 1.0.x

## ptmixed 1.0.3

-   Released: February 2021
-   Added `.checkmle()` step in `ptmixed()` to flag as not converged
    problematic cases on the boundary of the parameter space

## ptmixed 1.0.2

-   Released: January 2021
-   Updated vignette
-   Changed displaying style for function arguments in the documentation
-   Corrected BibTeX formatting of author names in CITATION file
-   Added to `make.spaghetti()` code to restore `bty`, `mar` and `xpd`
    values as they were before the function call

## ptmixed 1.0.1

-   Released: October 2020
-   Updated URL in DESCRIPTION file
-   Added `na.rm = T` in computation of `ylim` within `make.spaghetti()`
-   Added warnings in `simulate_ptglmm`

## ptmixed 1.0.0

-   Released: August 2020
-   The article describing the Poisson-Tweedie GLMM is now published in
    *Statistical Modelling*! The article is published with Open Access,
    so anybody can freely download it [from the website of Statistical
    Modelling](https://doi.org/10.1177/1471082X20936017)
-   Updated package description, citation info, vignette and help pages
    to include references to the article
-   Added example on how to compute the likelihood ratio test in the
    vignette
-   Added `margins` and `legend.space` arguments to `make.spaghetti()`.
    Added automatic sorting of provided dataframe ( = no need to
    pre-sort it any more!)

# 0.5.x

## ptmixed 0.5.4

-   Released: June 2020
-   Added link to arXiv preprint in package description and vignette
-   Added file with citation info
-   Added further arguments to `make.spaghetti()`; `cex.lab` argument
    fixed

## ptmixed 0.5.2

-   Released: April 2020
-   Adapted code so that it runs both for balanced and unbalanced
    datasets (previously balanced design was assumed)
-   Fixed problem with visualization of vignette on CRAN page

## ptmixed 0.5.1

-   Released: April 2020
-   Added vignette to illustrate the functionalities of the package
-   Simplified syntax of `ptmixed()`, `ptglm()`, `nbmixed()` and
    `nbglm()` (wrt the `id` and `offset` arguments). `ranef()` function
    updated accordingly
-   Added example dataset `df1`, used in the `ptmixed()` and `nbmixed()`
    help pages. Examples in help pages revised
-   Added `simulate_ptglmm()` function, to be used for illustration
    purposes (in the vignettes)
-   Added `pmf()` function to visualize the pmf of a discrete variable
-   `make.spaghetti()`: fixed minor bug in that arose when the `col`
    argument was specified + added `legend.inset` argument

# 0.4.x

## ptmixed 0.4.3

-   Released: February 2020
-   Added the possibility to use Laplace approximation, which is a
    special case of the adaptive Gauss Hermite quadrature method where
    just 1 quadrature point is used (simply set `npoints = 1` in
    `ptmixed()` or `nbmixed()`). Note: use of the Laplace is not
    recommended, because it is less accurate than the adaptive GH,
    results in lower convergence rates and can yield biased parameter
    estimates! We recommend using a sufficient number of quadrature
    points (5 typically produces a good likelihood approximation)
-   Added `make.spaghetti()` function to create a spaghetti plot /
    trajectory plot to visualize longitudinal data
-   Added example dataset `df1`
-   Added `silent` argument to `summary.ptglmm()`. Furthermore, printed
    output table with parameter estimates and Wald test is now presented
    with at most 4 decimals
-   Fixed bug that caused `ptglm()` and `nbglm()`to print detailed
    optimization info also when `trace = T`

## ptmixed 0.4.2

-   Released: January 2020
-   Changed class check within `wald.test()` to prevent problems with
    future `R` release (4.0.0)
-   Fixed bug that occurred when `freq.updates = 1` was set in
    `ptmixed()` and `nbmixed()`

## ptmixed 0.4.1

-   Released: October 2019
-   Computation of starting values for `ptmixed()` and `nbmixed()`
    improved
-   Added `wald.test()` function for computation of the multivariate
    Wald test
-   Added checks on `maxit[1] == 0` within `ptglm()` and `nbglm()` so as
    to make it possible to skip BFGS optimization and go straight to
    Nelder-Mead
-   Help files revised and improved
-   Added extra checks in `summary.ptglmm()` and `summary.ptglm()` (to
    verify that the smallest eigenvalue is not too small)

# 0.3.x

## ptmixed 0.3.1

-   Released: September 2019
-   Added `ptglm()` function for the estimation of a Poisson-Tweedie GLM
-   Added `nbmixed()` and `nbglm()` functions for the estimation of
    negative binomial GLMM and GLM using the Poisson-Tweedie
    parametrization (negative binomial: a = 0)
-   The package now comprises two classes: `ptglmm` for objects obtained
    from `ptmixed()` and `nbmixed()`, and `ptglm` for objects obtained
    from `ptglm()` and `nbglm()`. Summary functions for objects of both
    classes have been implemented
-   `min.var.init` argument added to `ptmixed()`

# 0.2.x

## ptmixed 0.2.1

-   Released: July 2019
-   Added function to compute the empirical Bayes estimates of the
    random intercept
-   Class name of `ptmixed()` output changed from `ptmm` to `ptglmm`
-   Corrected typo in `summary.ptglmm()` function (the MLE of the
    dispersion parameter was wrongly called “deviance” instead of
    dispersion in the previous versions)
-   Added NEWS file

# 0.1.x

## ptmixed 0.1.2

-   Released: June 2019
-   This is a major update aimed at speeding up the maximization of the
    loglikelihood. When `ptmixed()` is called, it first attempts to
    maximize the loglikelihood with the Nelder-Mead algorithm and then,
    if this fails, with the BFGS algorithm. Until version 0.0.4 the
    quadrature points were updated at every iteration for both
    Nelder-Mead and BFGS. Starting from this version, when Nelder-Mead
    is called it is possible to update the positioning of the quadrature
    points every *n* iterations by setting the `freq.updates` argument
    equal to *n*. Default is set to `freq.updates = 200` (this typically
    makes the optimization about 10 times faster than when
    `freq.updates = 1`)

# 0.0.x

## ptmixed 0.0.4

-   Released: May 2019
-   `ptmixed()` now outputs extra information (number of quadrature
    points used, initial values, warnings)
-   A mistake in the computation of the GH quadrature points was
    introduced from version 0.0.2. This has been fixed in this version

## ptmixed 0.0.3

-   Released: May 2019
-   Fixed a typo in message on initial loglikelihood value that is
    displayed when `trace = T` in `ptmixed()` function
-   Added exceptions and warnings for the case that `maxit[1]` and/or
    `maxit[2]` are set = 0

## ptmixed 0.0.2

-   Released: Apr. 2019
-   Function `ptmixed()` does not require the specification of a `time`
    argument any more
-   `maxit` argument default value in function `ptmixed()` increased to
    c(1e4, 100)
-   Internal function that computes starting values improved
-   Added warning with indication that a simpler Poisson mixed model may
    fit the data sufficiently well
-   Added warning when initial estimate of the variance parameter is
    &lt; 0.001

## ptmixed 0.0.1

-   Released: Feb. 2019
-   First version of the package
