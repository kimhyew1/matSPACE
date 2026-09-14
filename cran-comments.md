## R CMD check results

0 errors | 0 warnings | 1 note

* The "Days since last update" note is expected, as this is a follow-up
  release shortly after 0.1.0.

## Summary of changes since 0.1.0

* `matSPACE()`'s public API is simplified: the `sf_vec` argument (a
  vector of BIC scaling factors) is removed, and the return value is
  flattened from a list keyed by scaling factor (e.g. `fit$sf_1$V`) to a
  single result list (`fit$V`). This is a breaking change, hence the
  0.2.0 version bump.

* Internal cleanup: an experimental active-set coordinate-descent
  variant of the shooting algorithm was removed in favor of the
  original full-sweep implementation, with no user-facing effect.
