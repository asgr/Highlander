Lowlander = function(lower, upper, Nsamp = 'auto', pcut = 0.1, sigcut = NULL, corcut = 0,
                     likefunc, Data, liketype = 'min', parm.names = NULL,
                     seed = 666, latin = NULL, ncores = 1L) {

  # Input validation
  lower = as.numeric(lower)
  upper = as.numeric(upper)

  if (length(lower) == 0 || length(lower) != length(upper) ||
      anyNA(lower) || anyNA(upper)) {
    stop("'lower' and 'upper' must be numeric vectors of the same length (> 0) with no NA values.")
  }

  if(is.character(Nsamp)){
    if(Nsamp[1] == 'auto'){
      Nsamp = min(50L * length(lower), 2000L)
    }else{
      stop("Nsamp must be a positive integer (>= 1) or 'auto'.")
    }
  }

  if (length(Nsamp) != 1L || is.na(Nsamp) || !is.finite(Nsamp)) {
    stop("Nsamp must be a positive integer (>= 1) or 'auto'.")
  }
  Nsamp = as.integer(Nsamp)

  if (length(pcut) != 1L || is.na(as.numeric(pcut)) || as.numeric(pcut) <= 0 || as.numeric(pcut) >= 1) {
    stop("'pcut' must be a number strictly between 0 and 1.")
  }
  pcut = as.numeric(pcut)

  ncores = as.integer(ncores)
  if (ncores < 1 ) {
    stop("'ncores' must be integer >= 1.")
  }

  Npar = length(lower)

  # Set seed for reproducibility before generating LHS
  if (!is.null(seed)) {
    set.seed(seed)
  }

  if (!is.null(parm.names) && length(parm.names) != Npar){
    stop("parm.names must be NULL or have length equal to number of parameters.")
  }


  # Generate LHS if not supplied; only require lhs package when needed
  if (is.null(latin)) {
    if (!requireNamespace("lhs", quietly = TRUE)) {
      stop("The 'lhs' package is required to generate LHS designs but is not installed. ",
           "Install it with: install.packages('lhs')")
    }
    latin = lhs::randomLHS(Nsamp, Npar)
  } else {
    # Validate user-supplied latin: must be a numeric matrix of the right shape with values in [0, 1]
    if (!is.matrix(latin) || !is.numeric(latin) || nrow(latin) < Nsamp || ncol(latin) < Npar) {
      stop(paste0("'latin' must be a numeric matrix with at least", Nsamp, " rows (Nsamp) ",
                  "andat least", Npar, " columns (length(lower))."))
    }
    if (anyNA(latin) || any(latin < 0 | latin > 1)) {
      stop("All values in 'latin' must be finite numbers in [0, 1] with no NAs.")
    }
    if(dim(latin)[2] > Npar){
      latin = latin[,1:Npar]
    }
  }

  colnames(latin) = parm.names

  # Map unit-cube LHS samples to physical parameter bounds
  latin_mod = matrix(NA_real_, nrow = Nsamp, ncol = Npar)

  colnames(latin_mod) = parm.names

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
  good = is.finite(output)
  output_good = output[good]
  latin_mod_good = latin_mod[good,]

  cor_data = suppressWarnings(
    cor(cbind(latin_mod_good, output_good), use = "pairwise.complete.obs")
  )


  if(!is.null(parm.names)){
    colnames(cor_data) = c(parm.names, 'Out')
    rownames(cor_data) = c(parm.names, 'Out')
  }else{
    colnames(cor_data) = NULL
    rownames(cor_data) = NULL
  }

  if (liketype == 'min') {
    cutval = as.numeric(quantile(output_good, pcut, na.rm = TRUE))
    keep = which(output <= cutval & good) #this is correct since we want IDs rel to original output rowID
    best = latin_mod_good[which.min(output_good),]
  } else if (liketype == 'max') {
    cutval = as.numeric(quantile(output_good, 1 - pcut, na.rm = TRUE))
    keep = which(output >= cutval & good) #this is correct since we want IDs rel to original output rowID
    best = latin_mod_good[which.max(output_good),]
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

  kept_pars = latin_mod[keep, , drop = FALSE]

  if(is.null(sigcut)){
    # Compute new parameter bounds from kept samples
    new_lower = apply(kept_pars, 2, min)
    new_upper = apply(kept_pars, 2, max)
  }else{
    sigcut_pval = pnorm(-sigcut)
    new_lower = apply(kept_pars, 2, quantile, sigcut_pval)
    new_upper = apply(kept_pars, 2, quantile, 1 - sigcut_pval)
  }

  if(corcut > 0){
    cor_sel = is.finite(cor_data[1:Npar, Npar + 1L]) & (abs(cor_data[1:Npar, (Npar + 1L)]) < corcut)
    new_lower[cor_sel] = lower[cor_sel]
    new_upper[cor_sel] = upper[cor_sel]
  }

  vol_red = prod(new_upper - new_lower) / prod(upper - lower)

  frac_keep = as.numeric(length(keep) / Nsamp)

  return(invisible(list(
    lower     = new_lower,
    upper     = new_upper,
    best      = best,
    output    = output,
    keep      = keep,
    bad       = which(!good),
    vol_red   = vol_red,
    frac_keep = frac_keep,
    cutval    = cutval,
    cor_data  = cor_data,
    latin     = latin,
    latin_mod = latin_mod
  )))
}
