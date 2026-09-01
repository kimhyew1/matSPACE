#include <Rcpp.h>
using namespace Rcpp;

static inline double soft_thresh(double x, double lam) {
  double s = (x >= 0) ? 1.0 : -1.0;
  double a = (x >= 0) ? x : -x;
  double t = a - lam;
  return (t <= 0.0) ? 0.0 : s * t;
}

//' Active-set coordinate-descent shooting algorithm for matSPACE
//'
//' Compiled engine behind [space()]: fits the symmetric off-diagonal
//' partial-correlation coefficients via lasso-penalized coordinate
//' descent, restricting most sweeps to the current active set (nonzero
//' coefficients) for speed and falling back to a full global sweep to
//' detect newly active pairs. Called internally by [space()]; end users
//' should call [space()] instead.
//'
//' @param Y_data numeric vector of length `NN * PP * QQ`: the `NN`
//'   replicate `PP x QQ` matrices (weighted data), each flattened in
//'   row-major order and concatenated, as produced by R's
//'   `tensor_to_vec()`.
//' @param sigma_sr numeric vector of length `QQ`, the square root of
//'   `SIG / WEIGHT` per column; used to build the `B[i, j] = sigma_sr[i]
//'   / sigma_sr[j]` rescaling applied to each coefficient.
//' @param NN number of replicate matrices (`n`).
//' @param PP number of rows per matrix (`p`).
//' @param QQ number of columns/variables per matrix (`q`).
//' @param L1 lasso penalty (`lam`) applied to the off-diagonal
//'   coefficients via soft-thresholding.
//' @param n_iter maximum number of outer sweep iterations before giving
//'   up on convergence (each iteration is an active-set pass followed,
//'   when converged or empty, by a global sweep).
//' @param beta_init optional warm start: a length `QQ * QQ` numeric
//'   vector (row-major, `i * QQ + j` layout) used as the initial
//'   coefficient matrix instead of the default univariate
//'   soft-threshold initialization. Pass `R_NilValue` (the default) for
//'   a cold start.
//' @return A list with components:
//'   \item{beta}{length `QQ * QQ` numeric vector (row-major), the fitted
//'     symmetric coefficient matrix with a zero diagonal; in R this is
//'     reshaped and its diagonal set to 1 to form `ParCor`.}
//'   \item{Y_m}{length `NN * PP * QQ` numeric vector (row-major), `Y_data`
//'     after column-wise mean centering.}
//'   \item{E_m}{length `NN * PP * QQ` numeric vector (row-major), the
//'     residuals of `Y_m` after the fitted coefficients.}
//'   \item{iter_count}{total number of coordinate-descent coefficient
//'     updates performed across all active-set and global sweeps.}
//' @noRd
// [[Rcpp::export]]
List space_shooting(NumericVector Y_data,
                    NumericVector sigma_sr,
                    int NN, int PP, int QQ,
                    double L1,
                    int n_iter = 500,
                    Rcpp::Nullable<Rcpp::NumericVector> beta_init = R_NilValue)
{
  int n=NN, p=PP, q=QQ;
  double eps1 = 1e-6;

  bool warm = beta_init.isNotNull();
  NumericVector beta_init_vec;
  if (warm) {
    beta_init_vec = beta_init.get();
    if (beta_init_vec.size() != q*q)
      stop("beta_init length must equal QQ*QQ");
  }

  std::vector<double> Y_m(Y_data.begin(), Y_data.end());
  std::vector<double> beta_new(q*q, 0.0);
  std::vector<double> beta_old(q*q, 0.0);
  std::vector<double> beta_last(q*q, 0.0);
  std::vector<double> B(q*q), B_s(q*q);
  std::vector<double> normx(q, 0.0);

  double beta_change = 0.0;
  int change_i = 0, change_j = 1;
  int iter_count = 0;

  for (int j=0; j<q; j++) {
    double mn = 0.0;
    for (int h=0; h<n; h++)
      for (int i=0; i<p; i++)
        mn += Y_m[h*p*q + i*q + j];
    mn /= (double)(n*p);
    for (int h=0; h<n; h++)
      for (int i=0; i<p; i++) {
        Y_m[h*p*q + i*q + j] -= mn;
        normx[j] += Y_m[h*p*q + i*q + j] * Y_m[h*p*q + i*q + j];
      }
  }

  for (int i=0; i<q; i++)
    for (int j=0; j<q; j++)
      B[i*q+j] = sigma_sr[i] / sigma_sr[j];
  for (int i=0; i<q; i++)
    for (int j=0; j<=i; j++) {
      double v = B[i*q+j]*B[i*q+j]*normx[i] + B[j*q+i]*B[j*q+i]*normx[j];
      B_s[i*q+j] = B_s[j*q+i] = v;
    }

    std::vector<double> C(q*q, 0.0);
  for (int i=0; i<q; i++)
    for (int j=0; j<=i; j++) {
      double tmp = 0.0;
      for (int h=0; h<n; h++)
        for (int r=0; r<p; r++)
          tmp += Y_m[h*p*q + r*q + i] * Y_m[h*p*q + r*q + j];
      C[i*q+j] = C[j*q+i] = tmp;
    }

    // ── Step 0: 초기 beta (warm-start 없을 때) ───────
    if (!warm) {
      for (int i=0; i<q; i++)
        for (int j=0; j<=i; j++) {
          double t1 = C[i*q+j] * (B[j*q+i] + B[i*q+j]);
          double v  = soft_thresh(t1, L1) / B_s[i*q+j];
          beta_new[i*q+j] = beta_new[j*q+i] = v;
        }
    } else {
      for (int i=0; i<q*q; i++) beta_new[i] = beta_init_vec[i];
    }
    for (int i=0; i<q; i++) beta_new[i*q+i] = 0.0;

    std::vector<double> W(q*q, 0.0);
    for (int k=0; k<q; k++)
      for (int j=0; j<q; j++)
        W[k*q+j] = beta_new[k*q+j] * B[k*q+j];

    std::vector<double> S(q*q, 0.0);
    for (int i=0; i<q; i++)
      for (int j=0; j<q; j++) {
        double acc = 0.0;
        for (int k=0; k<q; k++) acc += C[i*q+k] * W[k*q+j];
        S[i*q+j] = C[i*q+j] - acc;
      }

      bool has_active = false;
    for (int jj=q-1; jj>=1 && !has_active; jj--)
      for (int ii=jj-1; ii>=0 && !has_active; ii--)
        if (beta_new[ii*q+jj] > eps1 || beta_new[ii*q+jj] < -eps1) {
          change_i = ii; change_j = jj;
          has_active = true;
        }

        if (has_active) {
          for (int i=0; i<q*q; i++) beta_old[i] = beta_new[i];
          int ci=change_i, cj=change_j;
          double Aij = S[ci*q+cj] * B[ci*q+cj];
          double Aji = S[cj*q+ci] * B[cj*q+ci];
          double beta_next = (Aij+Aji)/B_s[ci*q+cj] + beta_old[ci*q+cj];
          double new_val   = soft_thresh(beta_next, L1/B_s[ci*q+cj]);
          beta_change      = beta_old[ci*q+cj] - new_val;
          beta_new[ci*q+cj] = beta_new[cj*q+ci] = new_val;
        }

        std::vector<bool> in_active(q*q, false);
        std::vector<int>  active_pos(q*q, -1);
        std::vector<std::pair<int,int>> active_list;
        active_list.reserve(q*(q-1)/2);

        for (int jj=q-1; jj>=1; jj--)
          for (int ii=jj-1; ii>=0; ii--)
            if (beta_new[ii*q+jj] > eps1 || beta_new[ii*q+jj] < -eps1) {
              in_active[ii*q+jj] = true;
              active_pos[ii*q+jj] = (int)active_list.size();
              active_list.push_back({ii, jj});
            }

            std::vector<bool> touched(q, true);

            for (int iter=0; iter<n_iter; iter++) {

              for (int i=0; i<q*q; i++) beta_last[i] = beta_new[i];

              int nrow_pick = (int)active_list.size();
              double maxdif = -100.0;

              if (nrow_pick > 0) {
                for (int jrow=0; jrow<nrow_pick; jrow++) {
                  int ci = active_list[jrow].first;
                  int cj = active_list[jrow].second;

                  beta_old[change_i*q+change_j] = beta_new[change_i*q+change_j];
                  beta_old[change_j*q+change_i] = beta_new[change_j*q+change_i];

                                    double bc  = beta_change;
                  double Bji = B[change_j*q+change_i];
                  double Bij = B[change_i*q+change_j];
                  for (int i=0; i<q; i++) {
                    S[i*q+change_i] += bc * Bji * C[i*q+change_j];
                    S[i*q+change_j] += bc * Bij * C[i*q+change_i];
                  }

                  double Aij = S[ci*q+cj] * B[ci*q+cj];
                  double Aji = S[cj*q+ci] * B[cj*q+ci];

                  double beta_next = (Aij+Aji)/B_s[ci*q+cj] + beta_old[ci*q+cj];
                  double new_val   = soft_thresh(beta_next, L1/B_s[ci*q+cj]);
                  beta_change      = beta_old[ci*q+cj] - new_val;
                  beta_new[ci*q+cj] = beta_new[cj*q+ci] = new_val;
                  change_i=ci; change_j=cj;
                  iter_count++;
                  touched[ci] = true;
                  touched[cj] = true;
                }

                for (int i=0; i<q*q; i++) {
                  double d = beta_last[i]-beta_new[i]; if(d<0) d=-d;
                  if(d>maxdif) maxdif=d;
                }
              }

              if (maxdif < 1e-6 || nrow_pick < 1) {
                for (int i=0; i<q*q; i++) beta_last[i] = beta_new[i];

                for (int ci=0; ci<q-1; ci++)
                  for (int cj=ci+1; cj<q; cj++) {
                    if (!touched[ci] && !touched[cj]) continue;

                    beta_old[change_i*q+change_j] = beta_new[change_i*q+change_j];
                    beta_old[change_j*q+change_i] = beta_new[change_j*q+change_i];

                    double bc  = beta_change;
                    double Bji = B[change_j*q+change_i];
                    double Bij = B[change_i*q+change_j];
                    bool do_update = (bc > eps1 || bc < -eps1);
                    if (do_update) {
                      for (int i=0; i<q; i++) {
                        S[i*q+change_i] += bc * Bji * C[i*q+change_j];
                        S[i*q+change_j] += bc * Bij * C[i*q+change_i];
                      }
                    }

                    double Aij = S[ci*q+cj] * B[ci*q+cj];
                    double Aji = S[cj*q+ci] * B[cj*q+ci];

                    double beta_next = (Aij+Aji)/B_s[ci*q+cj] + beta_old[ci*q+cj];
                    double new_val   = soft_thresh(beta_next, L1/B_s[ci*q+cj]);
                    beta_change      = beta_old[ci*q+cj] - new_val;
                    beta_new[ci*q+cj] = beta_new[cj*q+ci] = new_val;
                    change_i=ci; change_j=cj;
                    iter_count++;

                    bool was_active = in_active[ci*q+cj];
                    bool now_active = (new_val > eps1 || new_val < -eps1);
                    if (!was_active && now_active) {
                      in_active[ci*q+cj] = true;
                      active_pos[ci*q+cj] = (int)active_list.size();
                      active_list.push_back({ci, cj});
                    } else if (was_active && !now_active) {
                      in_active[ci*q+cj] = false;
                      int pos  = active_pos[ci*q+cj];
                      int last = (int)active_list.size() - 1;
                      if (pos != last) {
                        active_list[pos] = active_list[last];
                        int mi = active_list[pos].first, mj = active_list[pos].second;
                        active_pos[mi*q+mj] = pos;
                      }
                      active_list.pop_back();
                      active_pos[ci*q+cj] = -1;
                    }
                  }

                  for (int t=0; t<q; t++) touched[t] = false;

                maxdif = -100.0;
                for (int i=0; i<q*q; i++) {
                  double d = beta_last[i]-beta_new[i]; if(d<0) d=-d;
                  if(d>maxdif) maxdif=d;
                }
                if (maxdif < 1e-6) break;
              }
            }

            std::vector<double> Em(n*p*q, 0.0);
            for (int h=0; h<n; h++)
              for (int r=0; r<p; r++)
                for (int j=0; j<q; j++) {
                  double yhat = 0.0;
                  for (int i=0; i<q; i++)
                    yhat += Y_m[h*p*q+r*q+i] * beta_new[i*q+j] * B[i*q+j];
                  Em[h*p*q+r*q+j] = Y_m[h*p*q+r*q+j] - yhat;
                }

                return List::create(
                  _["beta"]       = wrap(beta_new),
                  _["Y_m"]        = wrap(Y_m),
                  _["E_m"]        = wrap(Em),
                  _["iter_count"] = iter_count
                );
}
