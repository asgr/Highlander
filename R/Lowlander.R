Lowlander = function(lower, upper, Data, likefunc, Nsamp = 10000, pcut = 0.1, ncores = 4,
                     seed=666, liketype = 'min', latin = NULL) {

  Npar = length(lower)

  # 1. LHS in unit cube
  if(is.null(latin)){
    latin = randomLHS(Nsamp, Npar)
  }

  # 2. Map to physical parameters
  latin_mod <- sapply(seq_len(d), function(i)
    lower[i] + latin[, i] * (upper[i] - lower[i])
  )

  # 3. Evaluate likelihoods in parallel
  output = unlist(mclapply(1:Nsamp, FUN = function(i){likefunc(latin_mod[i, ], Data)}, mc.cores = ncores))

  # 4. Plausibility cut
  if(liketype == 'min'){
    keep = which(output <= quantile(output, pcut, na.rm=TRUE))
  }else if(liketype == 'max'){
    keep = which(output >= quantile(output, 1 - pcut, na.rm=TRUE))
  }else{
    stop('liketype must be min or max')
  }

  if (sum(keep) < d + 2)
    stop("Too few plausible points to shrink safely")

  # 5. New bounds
  new_lower = apply(latin_mod[keep, , drop = FALSE], 2, min, na.rm=TRUE)
  new_upper = apply(latin_mod[keep, , drop = FALSE], 2, max, na.rm=TRUE)

  list(
    lower = new_lower,
    upper = new_upper,
    output = output,
    keep = keep,
    frac_keep = length(keep) / Nsamp,
    latin = latin,
    latin_mod = latin_mod
  )
}
