#' Compute a BIC-type criterion for a fitted matSPACE model
#'
#' Computes a BIC-type model-selection score from the per-column residual
#' sums of squares of a fitted `beta`/`sig` pair, following the extended
#' BIC used to choose the lasso penalty `lam` in [space()].
#'
#' @param data list of n matrices, each p x q; the same replicated
#'   matrix-variate observations used to obtain `beta` and `sig`.
#' @param beta q x q matrix of regression coefficients (`beta[i, j]` is
#'   the effect of column `i` on column `j`), as returned in the `beta`
#'   component of [space()].
#' @param sig length-q vector of column residual precisions (1/variance),
#'   as returned in the `sig` component of [space()].
#' @param tol numeric tolerance below which a `beta` entry is treated as
#'   zero when counting degrees of freedom.
#' @param sf scaling factor applied to the `log(rss)` term of the BIC
#'   formula.
#' @return A single numeric BIC value; smaller indicates a better-fitting,
#'   more parsimonious model.
#' @noRd
compute_BIC = function(data, beta, sig, tol = 1e-6, sf = 2) {
  n = length(data)
  p = nrow(data[[1]])
  q = ncol(data[[1]])

  B = outer(sig, sig, "/")
  rss = numeric(q)
  for (j in seq_len(q)) {
    for (h in seq_len(n)) {
      Y = data[[h]]  # p x q
      yhat = Y %*% (beta[, j] * B[, j])  # p-vector
      yhat = yhat - Y[, j] * beta[j, j] * B[j, j]
      resid = Y[, j] - yhat
      rss[j] = rss[j] + sum(resid^2)
    }
  }

  beta_off = beta
  diag(beta_off) = 0
  deg = colSums(abs(beta_off) > tol)

  N = n*p
  BIC_j = N * sf * log(rss) + log(N) * deg / 2
  BIC = sum(BIC_j)
}

#' Flatten a list of matrices into a row-major numeric vector
#'
#' Internal helper that concatenates the n matrices in `data` (each
#' converted to row-major / C order) into a single flat numeric vector,
#' for passing to the compiled shooting function `space_shooting()`.
#'
#' @param data list of n matrices, each p x q.
#' @param n number of matrices in `data`. Currently unused; kept for
#'   interface symmetry since the dimensions are inferred from `data`
#'   itself.
#' @param p number of rows in each matrix. Currently unused; see `n`.
#' @param q number of columns in each matrix. Currently unused; see `n`.
#' @return A numeric vector of length `n * p * q`.
#' @noRd
tensor_to_vec = function(data, n, p, q) {
  unlist(lapply(data, function(mat) as.vector(t(mat))), use.names = FALSE)
}

#' Compute per-column weights for the SPACE weighted lasso
#'
#' @param f_type weighting scheme: `"equal"` for uniform weights,
#'   `"variance"` to weight columns by their estimated residual variance
#'   (`SIG`), or `"degree"` to weight columns by their number of nonzero
#'   partial correlations.
#' @param SIG length-q vector of current column residual precisions
#'   (1/variance); used when `f_type = "variance"`.
#' @param ParCor.fit q x q matrix of current partial correlation
#'   coefficients; used when `f_type = "degree"`.
#' @param q number of response columns.
#' @return A length-q numeric vector of weights, scaled to sum to `q`.
#' @noRd
compute_weight = function(f_type, SIG, ParCor.fit, q) {
  if (f_type == "equal")    return(rep(1, q))
  if (f_type == "variance") { w = SIG; return(w / sum(w) * q) }
  if (f_type == "degree") {
    deg <- apply(abs(ParCor.fit) > 1e-6, 1, sum)
    deg <- deg + max(deg)
    return(deg / sum(deg) * q)
  }
  stop("f_type must be 'equal', 'variance', or 'degree'")  # 예외 처리
}

#' Estimate per-column residual precision from a fitted beta matrix
#'
#' @param data list of n matrices, each p x q.
#' @param Beta q x q coefficient matrix; its diagonal is ignored (treated
#'   as zero) before computing residuals.
#' @param n number of matrices in `data`.
#' @param p number of rows in each matrix.
#' @param q number of columns in each matrix.
#' @return A length-q numeric vector of residual precisions
#'   (`1 / residual variance`) per column.
#' @noRd
estimate_sigma = function(data, Beta, n, p, q) {
  Beta.tmp = Beta
  diag(Beta.tmp) = 0
  rss = numeric(q)
  for (h in seq_len(n)) {
    Y   = data[[h]]
    res = Y - Y %*% Beta.tmp
    rss = rss + colSums(res^2)
  }
  1 / (rss / (n * p))
}

#' Convert partial correlation coefficients to a beta matrix
#'
#' @param coef numeric vector of upper-triangular partial correlation
#'   coefficients, of length `q * (q - 1) / 2`, in the order used by
#'   `matrix(0, q, q)[upper.tri(matrix(0, q, q))] <- coef`.
#' @param sig.fit length-q vector of column residual precisions
#'   (1/variance) used to rescale partial correlations into regression
#'   coefficients.
#' @return A q x q numeric matrix `beta` where
#'   `beta[i, j] = rho[i, j] * sqrt(sig.fit[i] / sig.fit[j])`; note this
#'   is generally not symmetric even though the underlying partial
#'   correlations `rho` are.
#' @noRd
rho_to_beta = function(coef, sig.fit) {
  q      = length(sig.fit)
  result = matrix(0, q, q)
  result[upper.tri(result)] = coef
  result = result + t(result)
  result * outer(sqrt(sig.fit), 1 / sqrt(sig.fit))
}

#' Generate a log-spaced lasso penalty sequence for [space()]
#'
#' Computes a decreasing, log-spaced sequence of `K` candidate lasso
#' penalties (`lam` values) for tuning [space()], analogous to the
#' `lambda_max`-based grids used in lasso path algorithms. The largest
#' value, `lambda_max`, is the largest off-diagonal entry of a
#' variance-rescaled Gram matrix — the smallest penalty above which every
#' off-diagonal partial correlation coefficient is driven to zero; the
#' sequence descends geometrically to `lambda_max * eps`.
#'
#' @param dt list of n matrices, each p x q, in the same format expected
#'   by the `data` argument of [space()].
#' @param eps ratio of the smallest to the largest penalty in the
#'   returned sequence, i.e. `lambda_min = lambda_max * eps`.
#' @param K number of penalty values to generate.
#' @return A numeric vector of length `K`, decreasing geometrically from
#'   `lambda_max` to `lambda_max * eps`, suitable to pass one at a time
#'   as the `lam` argument of [space()] (e.g. selecting among the fits
#'   with a BIC-type criterion).
#' @importFrom stats var
#' @export
#' @examples
#' set.seed(1)
#' p <- 5; q <- 4; n <- 3
#' data <- replicate(n, matrix(rnorm(p * q), p, q), simplify = FALSE)
#' lambda.bound(data, K = 10)
lambda.bound = function(dt, eps = 1e-6, K = 30) {

  dt_rbind = do.call(rbind, dt)

  sig = apply(dt_rbind, 2, var)
  W   = outer(sig, sig, "/")

  C   = crossprod(dt_rbind)
  YTX = t(W) * C

  G   = abs(YTX + t(YTX))

  lambda_max = max(G[upper.tri(G)])
  lambda_min = lambda_max * eps

  exp(seq(log(lambda_max), log(lambda_min), length.out = K))
}

# ── space ─────────────────────────────────────────
#' Fit a matrix-variate SPACE sparse partial correlation model
#'
#' Estimates a sparse network of partial correlations among the q columns
#' of matrix-variate data using an L1-penalized (lasso) SPACE-style
#' shooting algorithm, with optional per-column reweighting and residual
#' variance (`sig`) re-estimation across outer iterations.
#'
#' @param data      list of n matrices, each p × q. All matrices share
#'                  the same p × q shape; n is the number of independent
#'                  replicates and the q columns are the variables whose
#'                  pairwise partial correlations are estimated.
#' @param lam       lasso penalty applied to the off-diagonal partial
#'                  correlation coefficients.
#' @param sig       optional length-q sigma vector; NULL = iterative estimate
#' @param f_type    "equal" | "variance" | "degree"
#' @param iter      number of outer iterations
#' @param beta_init optional warm-start for beta: a q x q matrix/flat vector
#'                  (row-major, i*q+j layout) forwarded to the first internal
#'                  `space_shooting()` call instead of the default univariate
#'                  soft-threshold initialization. Only meaningful for the
#'                  first outer iteration, since that is the only one whose
#'                  sigma_sr (sqrt(SIG/WEIGHT)) does not depend on lam — the
#'                  first pass always starts from SIG = rep(1, q),
#'                  WEIGHT = rep(1, q) regardless of lam, so a beta fit at a
#'                  neighboring lambda is a valid warm start there. Requires
#'                  the compiled `space_shooting()` routine in
#'                  `src/matSPACE.cpp`, which accepts a `beta_init` argument.
#' @param sig_init  optional length-q warm-start for the STARTING value of
#'                  SIG when SIG.update = TRUE (i.e. sig = NULL). Unlike
#'                  `sig`, this does not pin SIG constant — it is still
#'                  re-estimated every outer iteration via estimate_sigma();
#'                  it only replaces the rep(1, q) cold start at i = 1 with a
#'                  neighboring lambda's converged SIG, so f_type =
#'                  "variance"/"degree" don't force a uniform-weight first
#'                  pass before ever seeing a realistic weighting.
#' @return A list with components:
#'   \item{ParCor}{q x q matrix of estimated partial correlations
#'     (diagonal = 1).}
#'   \item{beta}{q x q matrix of estimated regression coefficients.}
#'   \item{sig}{length-q vector of estimated column residual precisions
#'     (1/variance).}
#'   \item{f}{length-q vector of final per-column weights.}
#'   \item{total_iter}{total number of shooting iterations across all
#'     outer iterations.}
#'   \item{E}{list of n residual matrices (p x q each).}
#' @export
#' @examples
#' set.seed(1)
#' p <- 5; q <- 4; n <- 3
#' data <- replicate(n, matrix(rnorm(p * q), p, q), simplify = FALSE)
#' fit <- space(data, lam = 0.1, iter = 1)
#' fit$ParCor
space = function(data, lam, sig = NULL, f_type = "equal",
                  iter = 2, beta_init = NULL,
                  sig_init = NULL) {

  f_type = match.arg(f_type, c("equal", "variance", "degree"))

  n = length(data)
  p = nrow(data[[1]])
  q = ncol(data[[1]])

  SIG.update    = is.null(sig)
  SIG           = if (SIG.update) (if (is.null(sig_init)) rep(1, q) else sig_init) else sig

  WEIGHT.update = f_type %in% c("variance", "degree")
  ITER          = if (WEIGHT.update) max(2L, iter) else iter

  WEIGHT     = compute_weight(f_type, SIG, diag(q), q)
  ParCor.fit = NULL
  fit        = NULL

  total_iter_all = 0L

  beta_init_vec  = if (is.null(beta_init)) NULL else as.numeric(t(beta_init))
  prev_iter_beta = NULL

  for (i in seq_len(ITER)) {
    message("iter = ", i)

    w_mat  = matrix(sqrt(WEIGHT), nrow = p, ncol = q, byrow = TRUE)
    data.u = lapply(data, function(mat) mat * w_mat)
    sig.u  = SIG / WEIGHT

    Y_vec  = tensor_to_vec(data.u, n, p, q)

    fit = space_shooting(
      Y_data    = Y_vec,
      sigma_sr  = sqrt(sig.u),
      NN        = n,
      PP        = p,
      QQ        = q,
      L1        = lam,
      n_iter    = 500L,
      beta_init = if (i == 1L) beta_init_vec else prev_iter_beta
    )

    prev_iter_beta = fit$beta
    total_iter_all = total_iter_all + fit$iter_count

    ParCor.fit       = matrix(fit$beta, nrow = q, ncol = q, byrow = TRUE)
    diag(ParCor.fit) = 1

    coef     = ParCor.fit[upper.tri(ParCor.fit)]
    beta.cur = rho_to_beta(coef, SIG)

    if (!WEIGHT.update && !SIG.update) break

    if (SIG.update)    SIG    = estimate_sigma(data, beta.cur, n, p, q)
    if (WEIGHT.update) WEIGHT = compute_weight(f_type, SIG, ParCor.fit, q)
  }

  E_m_vec = fit$E_m
  E_list  = vector("list", n)
  for (h in seq_len(n)) {
    base        = (h - 1L) * p * q
    E_list[[h]] = matrix(
      E_m_vec[base + seq_len(p * q)], nrow = p, ncol = q, byrow = TRUE
    )
  }

  list(
    ParCor     = ParCor.fit,
    beta       = matrix(fit$beta, nrow = q, ncol = q, byrow = TRUE),
    sig        = SIG,
    f          = WEIGHT,
    total_iter = total_iter_all,
    E          = E_list
  )
}
