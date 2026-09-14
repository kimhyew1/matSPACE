
<!-- README.md is generated from README.Rmd. Please edit that file -->

# matSPACE

<!-- badges: start -->

[![R-CMD-check](https://github.com/kimhyew1/matSPACE/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/kimhyew1/matSPACE/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

This is a README file of the R package *matSPACE*. Our package fits
sparse partial correlation networks for matrix-variate data,
i.e. observations that are themselves p x q matrices rather than a
single vector of q measurements. We extend the SPACE (Sparse PArtial
Correlation Estimation) joint estimation framework to a
Kronecker-product precision structure, `Omega = V %x% U`, where `V` (q x
q) is the column-wise precision matrix and `U` (p x p) is the row-wise
precision matrix. All partial correlations are estimated jointly via an
L1-penalized (lasso) shooting algorithm, which preserves symmetry of the
estimated network and avoids the tuning-parameter selection difficulties
of separate node-wise regressions.

## Installation of the package

To install our package, you may simply execute the following codes:

``` r
install.packages("matSPACE")

# development version
# install.packages("pak")
pak::pak("kimhyew1/matSPACE")
```

## A basic example of using the package

We give a toy example to apply the main function `matSPACE`, which
estimates the Kronecker-structured row/column precision matrices `U` and
`V`.

### Generate data

We first generate simulated matrix-variate data.

``` r
set.seed(1)
p = 5
q = 4
n = 20

library(matSPACE)
data = replicate(n, matrix(rnorm(p * q), p, q), simplify = FALSE)
```

### Model fitting

Applying `matSPACE` to the simulated data, a lasso penalty path is
searched for `V` and `U` separately, and the BIC-minimizing fit is
returned for each.

``` r
fit = matSPACE(data, K = 10)
fit$V   # q x q column-wise precision matrix
#>      [,1]     [,2]     [,3]     [,4]
#> [1,]    1 0.000000 0.000000 0.000000
#> [2,]    0 1.109352 0.000000 0.000000
#> [3,]    0 0.000000 1.194146 0.000000
#> [4,]    0 0.000000 0.000000 1.294862
fit$U   # p x p row-wise precision matrix
#>            [,1]       [,2]     [,3]     [,4]      [,5]
#> [1,] 0.80938713 0.00979279 0.000000 0.000000 0.0000000
#> [2,] 0.00979279 1.00514343 0.000000 0.000000 0.0000000
#> [3,] 0.00000000 0.00000000 1.333242 0.000000 0.0000000
#> [4,] 0.00000000 0.00000000 0.000000 1.051813 0.0000000
#> [5,] 0.00000000 0.00000000 0.000000 0.000000 0.9207915
```

### Visualization of the lambda path

The full lambda path behind this selection is attached as a `"path"`
attribute, which lets us plot BIC against lambda without refitting.

``` r
path = attr(fit, "path")
plot(log(path$V$lambda), path$V$bic, type = "b",
     xlab = "lambda", ylab = "BIC", main = "Column (V) lambda path")
```

![](man/figures/README-unnamed-chunk-5-1.png)<!-- -->

### Fitting with a single lasso penalty

If you already know the penalty you want, `space` fits a single lasso
penalty directly, without searching a path.

``` r
fit_single = space(data, lam = 0.5)
fit_single$ParCor
#>             [,1]        [,2]       [,3]        [,4]
#> [1,]  1.00000000 -0.02296079 -0.0122418  0.17788420
#> [2,] -0.02296079  1.00000000  0.1118811 -0.06916583
#> [3,] -0.01224180  0.11188114  1.0000000 -0.16751241
#> [4,]  0.17788420 -0.06916583 -0.1675124  1.00000000
```

## Issues

We are happy to troubleshoot any issue with the package;

- please contact to the maintainer by <kimhw4126@gmail.com>, or
- please open an issue in the github repository.
