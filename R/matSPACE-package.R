#' matSPACE: Sparse Partial Correlation Estimation for Matrix-Variate Data
#'
#' Fits a sparse network of partial correlations among the columns of
#' matrix-variate data using an L1-penalized (lasso) SPACE-style shooting
#' algorithm, with optional per-column reweighting and residual variance
#' re-estimation across outer iterations. See [space()] for the main
#' model-fitting function.
#'
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom Rcpp sourceCpp
#' @useDynLib matSPACE, .registration = TRUE
## usethis namespace: end
NULL
