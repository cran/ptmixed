ptmixed 0.2.1
=============

-   Released: July 2019
-   Added function to compute the empirical Bayes estimates of the
    random intercept
-   Class name of `ptmixed` output changed from ptmm to ptglmm
-   Corrected typo in `summary` function (the MLE of the dispersion
    parameter was wrongly called "deviance" instead of dispersion in the
    previous versions)
-   Added NEWS file

ptmixed 0.1.2
=============

-   Released: June 2019
-   This is a major update aimed at speeding up the maximization of the
    loglikelihood. When `ptmixed` is called, it first attempts to
    maximize the loglikelihood with the Nelder-Mead algorithm and then,
    if this fails, with the BFGS algorithm. Until version 0.0.4 the
    quadrature points were updated at every iteration for both
    Nelder-Mead and BFGS. Starting from this version, when Nelder-Mead
    is called it is possible to update the positioning of the quadrature
    points every *n* iterations by setting the `freq.updates` argument
    equal to *n*. Default is set to `freq.updates = 200` (on average,
    this makes the optimization 2.5 times faster than when
    `freq.updates = 1`)

ptmixed 0.0.4
=============

-   Released: May 2019
-   `ptmixed` now outputs extra information (number of quadrature points
    used, initial values, warnings)
-   A mistake in the computation of the GH quadrature points was
    introduced from version 0.0.2. This has been fixed in this version

ptmixed 0.0.3
=============

-   Released: May 2019
-   Fixed a typo in message on initial loglikelihood value that is
    displayed when `trace = T` in `ptmixed` function
-   Added exceptions and warnings for the case that `maxit[1]` and/or
    `maxit[2]` are set = 0

ptmixed 0.0.2
=============

-   Released: Apr. 2019
-   Function `ptmixed` does not require the specification of a `time`
    argument any more
-   `maxit` argument default value in function `ptmixed` increased to
    c(1e4, 100)
-   Internal function that computes starting values improved
-   Added warning with indication that a simpler Poisson mixed model may
    fit the data sufficiently well
-   Added warning when initial estimate of the variance parameter is
    &lt; 0.001

ptmixed 0.0.1
=============

-   Released: Feb. 2019
-   First version of the package