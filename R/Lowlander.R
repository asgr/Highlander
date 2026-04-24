Lowlander = function(lower, upper, Nsamp = 100, pcut = 0.1,
                     likefunc, Data, liketype = 'min',
                     seed = NULL, latin = NULL, ncores = 1) {

  # Input validation
  lower = as.numeric(lower)
  upper = as.numeric(upper)

  if (length(lower) == 0 || length(lower) != length(upper) ||
      anyNA(lower) || anyNA(upper)) {
    stop("'lower' and 'upper' must be numeric vectors of the same length (> 0) with no NA values.")
  }

  Nsamp = as.integer(Nsamp)
  if (length(Nsamp) != 1L || is.na(Nsamp) || Nsamp < 1L) {
    stop("'Nsamp' must be a positive integer (>= 1).")
  }

  if (length(pcut) != 1L || is.na(as.numeric(pcut)) || as.numeric(pcut) <= 0 || as.numeric(pcut) >= 1) {
    stop("'pcut' must be a number strictly between 0 and 1.")
  }
  pcut = as.numeric(pcut)

  ncores = as.integer(ncores)
  if (length(ncores) != 1L || is.na(ncores) || ncores < 1L) {
    stop("'ncores' must be a positive integer.")
  }

  # Check lhs package is available
  if (!requireNamespace("lhs", quietly = TRUE)) {
    stop("The 'lhs' package is required but not installed. Install it with: install.packages('lhs')")
  }

  Npar = length(lower)

  # Set seed for reproducibility before generating LHS
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Generate LHS if not supplied
  if (is.null(latin)) {
    latin = lhs::randomLHS(Nsamp, Npar)
  }

  # Map unit-cube LHS samples to physical parameter bounds
  latin_mod = matrix(NA_real_, nrow = Nsamp, ncol = Npar)
  for (j in seq_len(Npar)) {
    latin_mod[, j] = lower[j] + latin[, j] * (upper[j] - lower[j])
  }

  # Evaluate likefunc in parallel (Windows-safe)
  if (.Platform$OS.type == "windows" && ncores > 1L) {
    # PSOCK cluster for Windows multi-core
    cl = parallel::makeCluster(ncores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterExport(cl, varlist = c("latin_mod", "Data", "likefunc"),
                            envir = environment())
    output = unlist(parallel::parLapply(cl, seq_len(Nsamp),
                                        function(i) likefunc(latin_mod[i, ], Data)))
  } else {
    output = unlist(parallel::mclapply(seq_len(Nsamp),
                                       function(i) likefunc(latin_mod[i, ], Data),
                                       mc.cores = ncores))
  }

  # Determine which samples to keep based on quantile cut
  if (liketype == 'min') {
    cutval = as.numeric(quantile(output, pcut, na.rm = TRUE))
    keep = which(output <= cutval)
  } else if (liketype == 'max') {
    cutval = as.numeric(quantile(output, 1 - pcut, na.rm = TRUE))
    keep = which(output >= cutval)
  } else {
    stop("'liketype' must be 'min' or 'max'.")
  }

  # Ensure enough kept points to shrink the parameter space safely.
  # Need at least Npar + 2 points so that min/max bounds differ from each other
  # across all Npar dimensions and there is room to shrink.
  if (length(keep) < Npar + 2L) {
    stop(paste0("Too few plausible points kept (", length(keep), ") to shrink safely; ",
                "need at least Npar + 2 = ", Npar + 2L, ". ",
                "Try increasing Nsamp or pcut."))
  }

  # Compute new parameter bounds from kept samples
  kept_pars = latin_mod[keep, , drop = FALSE]
  new_lower = apply(kept_pars, 2, min)
  new_upper = apply(kept_pars, 2, max)

  frac_keep = as.numeric(length(keep) / Nsamp)

  return(invisible(list(
    lower     = new_lower,
    upper     = new_upper,
    output    = output,
    keep      = keep,
    frac_keep = frac_keep,
    cutval    = cutval,
    latin     = latin,
    latin_mod = latin_mod,
    kept_pars = kept_pars
  )))
}
