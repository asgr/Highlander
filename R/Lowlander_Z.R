Lowlander_kept = function(Lowlander = NULL, liketype = NULL, parm.names = NULL) {
  # The kept (plausible) parameter samples from a Lowlander() return, ranked
  # best-first. Returns NULL with a warning when the object is unusable.
  #
  # Lowlander restricts `keep` to finite-likelihood samples past its quantile
  # cut, so these rows are in-bounds by construction.

  if (is.null(Lowlander)) {
    return(NULL)
  }

  if (!is.list(Lowlander) || is.null(Lowlander$latin_mod) ||
      is.null(Lowlander$keep) || is.null(Lowlander$output)) {
    warning('Lowlander object is missing latin_mod/keep/output.')
    return(NULL)
  }

  latin_mod = as.matrix(Lowlander$latin_mod)
  keep = as.integer(Lowlander$keep)
  output = as.numeric(Lowlander$output)

  if (length(keep) == 0 || any(keep < 1L | keep > nrow(latin_mod))) {
    warning('Lowlander object has no valid kept samples.')
    return(NULL)
  }

  # Defensive: a hand-edited or older object should not be able to smuggle a
  # non-finite row through.
  outk = output[keep]
  ok = is.finite(outk)
  if (!any(ok)) {
    warning('All kept Lowlander samples have non-finite likelihood.')
    return(NULL)
  }
  kept = latin_mod[keep[ok], , drop = FALSE]
  outk = outk[ok]

  # Prefer the direction Lowlander actually used over anything passed in here.
  lt = Lowlander$liketype
  if (is.null(lt)) lt = liketype
  if (is.null(lt)) lt = 'min'

  ord = if (identical(lt, 'max')) order(outk, decreasing = TRUE) else order(outk)
  kept = kept[ord, , drop = FALSE]

  if (!is.null(parm.names) && length(parm.names) == ncol(kept)) {
    colnames(kept) = parm.names
  }

  return(kept)
}

Lowlander_Z = function(Lowlander = NULL, N = 1, liketype = NULL, parm.names = NULL) {
  # A design matrix of N plausible parameter samples drawn from a Lowlander()
  # return, to seed ensemble samplers (AIES Z) or anything else that needs
  # several starting points inside the sensible region. Returns NULL when
  # unavailable so callers can fall back to their own defaults.
  #
  # Rows sit at evenly-spaced quantile ranks of the kept samples, ranked
  # best-first. This guarantees distinct walkers and covers the plausible
  # region rather than piling up on the mode; AIES proposes
  # theta' = theta_s + z (theta_i - theta_s), so walkers that all sit on top of
  # each other can only ever propose where they already are. Choosing a
  # different rank spacing (e.g. denser near the best point) performed
  # similarly in testing, so nothing here should be read as an optimised
  # ranking - the benefit comes from using Lowlander's kept samples at all
  # rather than leaving Z = NULL.

  kept = Lowlander_kept(Lowlander, liketype = liketype, parm.names = parm.names)
  if (is.null(kept)) {
    return(NULL)
  }

  # One point cannot support an ensemble: the stretch move divides by the
  # difference between two distinct walkers.
  nk = nrow(kept)
  if (nk < 2L) {
    warning('Fewer than 2 usable Lowlander samples; cannot build a spread Z.')
    return(NULL)
  }

  N = max(1L, as.integer(N))

  if (N <= nk) {
    idx = floor(approx(seq(0, 1, length.out = nk), seq_len(nk),
                       xout = (seq_len(N) - 0.5) / N)$y)
    idx = pmin(pmax(idx, 1L), nk)
    Z = kept[idx, , drop = FALSE]
  } else {
    # More walkers than samples: recycle the kept set, then jitter the repeats
    # (exactly-duplicated rows degenerate the stretch move for any pair sharing
    # a position). Jitter is clamped back inside the plausible box.
    reps = rep(seq_len(nk), length.out = N)
    Z = kept[reps, , drop = FALSE]
    dup = which(duplicated(reps))
    if (length(dup) > 0) {
      span = apply(kept, 2, max) - apply(kept, 2, min)
      flat = !is.finite(span) | span <= 0
      span[flat] = pmax(abs(kept[1, ][flat]), 1)
      jit = matrix(rnorm(length(dup) * ncol(Z)), nrow = length(dup)) *
             matrix(rep(span * 1e-3, each = length(dup)), nrow = length(dup))
      Z[dup, ] = Z[dup, , drop = FALSE] + jit
      lo = apply(kept, 2, min)
      hi = apply(kept, 2, max)
      Z = matrix(pmin(pmax(as.vector(Z), rep(lo, each = nrow(Z))),
                      rep(hi, each = nrow(Z))), nrow = nrow(Z))
    }
  }

  Z = matrix(as.numeric(Z), nrow = N)
  if (!is.null(colnames(kept))) colnames(Z) = colnames(kept)
  return(Z)
}

Lowlander_best = function(Lowlander = NULL, parm.names = NULL) {
  # The single most plausible point from a Lowlander() return as a plain
  # numeric vector, for seeding Initial.Values / twalk SIV / IM mu.
  if (is.null(Lowlander)) {
    return(NULL)
  }
  best = Lowlander$best
  if (is.null(best)) {
    warning('Lowlander object has no $best element.')
    return(NULL)
  }
  best = as.numeric(unname(best))
  if (!is.null(parm.names) && length(parm.names) == length(best)) {
    names(best) = parm.names
  }
  return(best)
}
