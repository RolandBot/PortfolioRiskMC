
simulationKernel <- function(pf, Z, r,
                             J,
                             K, mk = seq_len(K * nrow(Z)),
                             agg = factor(rep("PF", nrow(pf))),
                             seed,
                             timer = function(expr, ...) expr) {

  # check arguments
  stopifnot(all(c("j", "V0", "R", "PD") %in% colnames(pf)))
  stopifnot(all.equal(dim(Z), dim(r)))
  stopifnot(all(pf$j > 0 & pf$j <= ncol(Z)))
  stopifnot(all(mk > 0 & mk <= nrow(Z) * K))
  stopifnot(!is.unsorted(mk))
  stopifnot(is.function(timer))

  # aggregation levels
  agg <- as.factor(agg)
  agg_idx <- as.numeric(agg)
  agg_set <- as.character(levels(agg))

  # pre-allocate the aggregated output
  L_agg <- matrix(0.0, length(mk), max(agg_idx),
                  dimnames = list(mk = as.character(mk),
                                  agg = agg_set))

  # timer wrapping the core simulation only
  t <- timer(
    simulationKernel_C(pf = pf, Z = Z, r = r,
                       J = J,
                       K = K, mk = mk,
                       agg = agg_idx,
                       seed = seed,
                       L_agg = L_agg)
  )
  if (!is.null(t)) {
    print(t)
  }

  L_agg

}
