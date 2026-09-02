Specs_help = function(Algorithm = 'CHARM', Data = NULL, Lowlander = NULL,
                      Iterations = NULL, Thinning = 1, ...) {
  # Returns the default Specs list for the given LaplacesDemon Algorithm name.
  # For algorithms that require user-supplied values (no built-in defaults),
  # a template list is returned with placeholder values and a message is
  # printed to guide the user.
  #
  # If a Lowlander() return object is supplied it is used to fill in the
  # parts that LaplacesDemon would otherwise guess blind (via GIV), which is
  # the usual cause of out-of-bounds failures on restricted likelihoods.

  dots = list(...)

  # Number of parameters, from whichever source is available.
  Npar = NULL
  if (!is.null(Data$parm.names)) {
    Npar = length(Data$parm.names)
  } else if (!is.null(Lowlander$latin_mod)) {
    Npar = ncol(as.matrix(Lowlander$latin_mod))
  }
  if (!is.null(Npar) && Npar < 1L) {
    stop('Cannot determine the number of parameters from Data or Lowlander.')
  }

  .even_Nc = function(Nc, CPUs = 1) {
    Nc = max(3L, as.integer(Nc))
    if (CPUs > 1 && Nc %% 2 != 0) {
      Nc = Nc + 1L
      message(Algorithm, ': Nc rounded up to ', Nc, ' (must be even when CPUs > 1).')
    }
    Nc
  }
  .Nc_ens = function(CPUs = 1) {
    .even_Nc(if (is.null(Npar)) 3L else max(4L, 2L * Npar), CPUs)
  }

  # LaplacesDemon rewrites Iterations/Thinning before it ever looks at Specs:
  # Iterations becomes round(abs(Iterations)) clamped up to a minimum of 11, and
  # Thinning becomes round(abs(Thinning)) but falls back to 1 unless it sits in
  # [1, Iterations]. DEMC then sizes Z as floor(Iterations/Thinning)+1 rows, so
  # an array built from the *raw* values can mismatch and be hard-stopped with
  # "The first dimension of Z is incorrect". Mirror the same rules here.
  .ld_iters = function(Iters) {
    Iters = round(abs(Iters))
    max(11L, as.integer(Iters))
  }
  .ld_thin = function(Thin, Iters_adj) {
    Thin = round(abs(Thin))
    if (is.na(Thin) || Thin < 1L || Thin > Iters_adj) 1L else as.integer(Thin)
  }
  .z_rows = function(Iters, Thin) {
    Iters = .ld_iters(Iters)
    floor(Iters / .ld_thin(Thin, Iters)) + 1L
  }
  # Non-finite Iterations cannot size anything (as.integer(Inf) is NA), so the
  # callers below must fall back to Z = NULL rather than build a broken array.
  .sizeable = !is.null(Iterations) && is.finite(Iterations) && !is.null(Npar)

  # Plausible starting points mined from the Lowlander run, if we have one.
  # With no Lowlander object these stay NULL and LaplacesDemon falls back to
  # GIV(), its own blind guesser; that is the pre-existing behaviour and is not
  # worth warning about.
  .Zmatrix = function(Nc) {
    if (is.null(Lowlander)) return(NULL)
    Z = Lowlander_Z(Lowlander, N = Nc, parm.names = Data$parm.names)
    if (is.null(Z)) {
      message(paste0(Algorithm, ': no usable Z could be built from the Lowlander ',
                     'object; LaplacesDemon will fall back to GIV(), which often ',
                     'fails on out-of-bounds proposals.'))
    }
    Z
  }
  .seed = function(what) {
    if (is.null(Lowlander)) return(NULL)
    b = Lowlander_best(Lowlander, parm.names = Data$parm.names)
    if (is.null(b)) {
      message(paste0(Algorithm, ': ', what, ' not auto-filled from the Lowlander object.'))
    }
    b
  }
  # Per-parameter evaluation grids for the Griddy-Gibbs family. Grid holds
  # OFFSETS from the current value (LaplacesDemon evaluates at
  # prop[j] + Grid[[j]]), so the built-in default of seq(-0.1, 0.1) is a fixed
  # +-0.1 probe that is far too coarse for parameters spanning several decades
  # and too fine for tightly-bounded ones. When a Lowlander object gives us the
  # plausible width of each parameter, scale each grid to that width.
  .grid = function(n = 7, frac = 0.5) {
    if (is.null(Npar)) return(NULL)
    lo = Lowlander$lower; hi = Lowlander$upper
    if (is.null(lo) || is.null(hi) || length(lo) != Npar || length(hi) != Npar ||
        any(!is.finite(lo)) || any(!is.finite(hi)) || any(lo >= hi)) {
      return(NULL)
    }
    n = max(3L, as.integer(n))
    half = frac * (hi - lo) / 2
    lapply(seq_len(Npar), function(j) seq(-half[j], half[j], length.out = n))
  }

  if (Algorithm == 'ADMG') {
    # Adaptive Directional Metropolis-within-Gibbs
    Specs_list = list(n = 0, Periodicity = 1)

  } else if (Algorithm == 'AFSS') {
    # Automated Factor Slice Sampler
    Specs_list = list(A = Inf, B = NULL, m = Inf, n = 0, w = 1)

  } else if (Algorithm == 'AGG') {
    # Adaptive Griddy-Gibbs (Specs required – no default)
    Grid = .grid()
    if (is.null(Grid)) {
      message('AGG requires a per-parameter Grid; pass a Lowlander object to auto-build one. Returning template.')
      Grid = NULL
    }
    # dparm must stay NULL: it names the *discrete* parameters, and passing 0
    # makes LaplacesDemon clamp it to a real index, silently treating one
    # continuous parameter as discrete. NULL is the "none" sentinel.
    Specs_list = list(Grid = Grid, dparm = NULL, smax = 0.1, CPUs = 1,
                Packages = NULL, Dyn.libs = NULL)

  } else if (Algorithm == 'AHMC') {
    # Adaptive Hamiltonian Monte Carlo
    # epsilon must be length(Initial.Values); LaplacesDemon recycles a scalar,
    # so we build the per-parameter vector it is asking for.
    npar = if (is.null(Npar)) 1L else Npar
    message('AHMC: epsilon length must equal length(Initial.Values).')
    # m = NULL makes LaplacesDemon use an identity mass matrix; a plain vector
    # is rejected with "x must be a square matrix".
    Specs_list = list(epsilon = rep(1 / npar, npar), L = 2, m = NULL,
                Periodicity = 1)

  } else if (Algorithm == 'AIES') {
    # Affine-Invariant Ensemble Sampler (Specs required – no default)
    Nc = if (is.null(dots$Nc)) .Nc_ens(if (is.null(dots$CPUs)) 1 else dots$CPUs)
         else .even_Nc(dots$Nc, if (is.null(dots$CPUs)) 1 else dots$CPUs)
    if (is.null(dots$Nc) && is.null(Npar)) {
      message('AIES requires user-supplied Nc or Data. Returning template default Nc = 3.')
    } else if (is.null(dots$Nc)) {
      message('AIES: using Nc = 2 * length(parm.names) = ', Nc, '.')
    }
    Specs_list = list(Nc = Nc, Z = .Zmatrix(Nc), beta = 2, CPUs = 1,
                Packages = NULL, Dyn.libs = NULL)

  } else if (Algorithm == 'AM') {
    # Adaptive Metropolis
    Specs_list = list(Adaptive = 1000, Periodicity = 1)

  } else if (Algorithm == 'AMM') {
    # Adaptive-Mixture Metropolis
    Specs_list = list(Adaptive = 1000, B = NULL, n = 0, Periodicity = 1, w = 0.05)

  } else if (Algorithm == 'AMWG') {
    # Adaptive Metropolis-within-Gibbs
    Specs_list = list(B = NULL, n = 0, Periodicity = 50)

  } else if (Algorithm == 'CHARM') {
    # Componentwise Hit-And-Run Metropolis
    if('alpha.star' %in% names(dots)){
      Specs_list = list(alpha.star = 0.234)
    }else{
      Specs_list = NULL
    }
  } else if (Algorithm == 'DEMC') {
    # Differential Evolution Markov Chain (Specs required – no default)
    # LaplacesDemon wants Z as an array: floor(Iterations/Thinning)+1 x
    # length(parm) x Nc. Without Iterations we cannot size it, so leave Z
    # NULL and let LaplacesDemon generate it (via GIV).
    Nc = .Nc_ens(1)
    if (!is.null(dots$Nc)) Nc = max(3L, as.integer(dots$Nc))
    Z = NULL
    if (is.null(Iterations) || !is.finite(Iterations)) {
      message("DEMC: pass a finite Iterations (and Thinning) to auto-build Z; ",
              "leaving Z = NULL so LaplacesDemon falls back to GIV().")
    } else if (is.null(Npar)) {
      message('DEMC: cannot determine the parameter count; leaving Z = NULL.')
    } else {
      pool = Lowlander_kept(Lowlander, parm.names = Data$parm.names)
      if (is.null(Lowlander)) {
        message('DEMC: pass a Lowlander object to seed Z with in-bounds samples; ',
                'leaving Z = NULL so LaplacesDemon falls back to GIV().')
      } else if (is.null(pool)) {
        message('DEMC: no usable Lowlander samples; Z = NULL, LaplacesDemon ',
                'will fall back to GIV().')
      } else {
        # Each chain gets its own history of plausible samples. LaplacesDemon
        # proposes gamma * (Z[r1,,s1] - Z[r2,,s2]) from these rows, so they all
        # need to be usable. Note it first overwrites Z[1,,1] with the value
        # returned at Initial.Values (so that one slot is ours to fill but not
        # to keep, exactly as AIES ignores Z[1,]), then evaluates Model() at
        # Z[1,,i] for i = 2..Nc with no try() wrapper -- an error there aborts
        # the run outright, so those rows must be in-bounds.
        nrow_z = .z_rows(Iterations, Thinning)
        Z = array(0, dim = c(nrow_z, Npar, Nc))
        for (i in seq_len(dim(Z)[3])) {
          Z[, , i] = pool[sample.int(nrow(pool), nrow_z, replace = TRUE), , drop = FALSE]
        }
      }
    }
    Specs_list = list(Nc = Nc, Z = Z, gamma = NULL, w = 0.1)

  } else if (Algorithm == 'DRAM') {
    # Delayed Rejection Adaptive Metropolis
    Specs_list = list(Adaptive = 1000, Periodicity = 1)

  } else if (Algorithm == 'DRM') {
    # Delayed Rejection Metropolis (no Specs needed)
    Specs_list = NULL

  } else if (Algorithm == 'ESS') {
    # Elliptical Slice Sampler
    Specs_list = list(B = NULL)

  } else if (Algorithm == 'GG') {
    # Griddy-Gibbs (Specs required – no default)
    Grid = .grid()
    if (is.null(Grid)) {
      message('GG requires a per-parameter Grid; pass a Lowlander object to auto-build one. Returning template.')
      Grid = NULL
    }
    Specs_list = list(Grid = Grid, dparm = NULL, CPUs = 1,
                Packages = NULL, Dyn.libs = NULL)

  } else if (Algorithm == 'Gibbs') {
    # Gibbs Sampler (FC must be a function; MWG is optional)
    message('Gibbs requires FC to be a full-conditional function. Returning template.')
    Specs_list = list(FC = NULL, MWG = NULL)

  } else if (Algorithm == 'HARM') {
    # Hit-And-Run Metropolis
    if(any(c('alpha.star', 'B') %in% names(dots))){
      Specs_list = list(alpha.star = 0.234, B = NULL)
    }else{
      Specs_list = NULL
    }
  } else if (Algorithm == 'HMC') {
    # Hamiltonian Monte Carlo
    npar = if (is.null(Npar)) 1L else Npar
    message('HMC: epsilon length must equal length(Initial.Values).')
    Specs_list = list(epsilon = rep(1 / npar, npar), L = 2, m = NULL)

  } else if (Algorithm == 'HMCDA') {
    # Hamiltonian Monte Carlo with Dual-Averaging (Specs required – no default)
    # A NULL epsilon lets dual-averaging tune the step size adaptively.
    message('HMCDA: A must be well below Iterations; epsilon = NULL lets dual-averaging tune it.')
    Specs_list = list(A = max(1L, if (is.null(Iterations)) floor(1000) else floor(Iterations / 2)),
                delta = 0.65, epsilon = NULL, Lmax = 1000, lambda = 0.1)

  } else if (Algorithm == 'IM') {
    # Independence Metropolis (Specs required – no default)
    # mu is the mean of the proposal distribution and must be exactly
    # length(Initial.Values) or LaplacesDemon stops. Left NULL rather than
    # guessed when there is no Lowlander object: an arbitrary mu (e.g. all
    # zeros) still "runs" but silently biases where proposals land, which is
    # worse than LaplacesDemon's own clear length error.
    mu = .seed('mu')
    if (is.null(mu)) {
      message('IM requires user-supplied mu of length equal to Initial.Values (or pass a Lowlander object).')
    }
    Specs_list = list(mu = mu)

  } else if (Algorithm == 'INCA') {
    # Interchain Adaptation
    Specs_list = list(Adaptive = 1000, Periodicity = 1)

  } else if (Algorithm == 'MALA') {
    # Metropolis-Adjusted Langevin Algorithm
    # Note: gamma is required by LaplacesDemon but absent from its built-in
    # default; a typical value of 0.6 is used here.
    Specs_list = list(A = 1e7, alpha.star = 0.574, delta = 1, gamma = 0.6,
                epsilon = c(1e-6, 1e-7))

  } else if (Algorithm == 'MCMCMC') {
    # Metropolis-Coupled Markov Chain Monte Carlo
    Specs_list = list(lambda = 1, CPUs = 1, Packages = NULL, Dyn.libs = NULL)

  } else if (Algorithm == 'MTM') {
    # Multiple-Try Metropolis
    Specs_list = list(K = 4, CPUs = 1, Packages = NULL, Dyn.libs = NULL)

  } else if (Algorithm == 'MWG') {
    # Metropolis-within-Gibbs (LaplacesDemon default algorithm)
    Specs_list = list(B = NULL)

  } else if (Algorithm == 'NUTS') {
    # No-U-Turn Sampler (Specs required – no default)
    message('NUTS requires user-supplied Specs. Returning template.')
    Specs_list = list(A = 1000, delta = 0.6, epsilon = NULL, Lmax = 1000)

  } else if (Algorithm == 'OHSS') {
    # Oblique Hyperrectangle Slice Sampler
    Specs_list = list(A = Inf, n = 0)

  } else if (Algorithm == 'pCN') {
    # Preconditioned Crank-Nicolson
    Specs_list = list(beta = 0.01)

  } else if (Algorithm == 'RAM') {
    # Robust Adaptive Metropolis
    Specs_list = list(alpha.star = 0.234, B = NULL, Dist = 'N', gamma = 0.66,
                n = 0)

  } else if (Algorithm == 'RDMH') {
    # Random Dive Metropolis-Hastings (no Specs needed)
    Specs_list = NULL

  } else if (Algorithm == 'Refractive') {
    # Refractive Sampler
    Specs_list = list(Adaptive = 1, m = 2, w = 0.1, r = 1.3)

  } else if (Algorithm == 'RJ') {
    # Reversible-Jump (Specs required – no default)
    message('RJ requires user-supplied Specs. Returning template.')
    Specs_list = list(bin.n = 1, bin.p = 0.5, parm.p = 0.5,
                selectable = NULL, selected = NULL)

  } else if (Algorithm == 'RSS') {
    # Reflective Slice Sampler (Specs required – no default)
    message('RSS requires user-supplied Specs. Returning template.')
    Specs_list = list(m = 10, w = 1)

  } else if (Algorithm == 'RWM') {
    # Random-Walk Metropolis
    Specs_list = list(B = list())

  } else if (Algorithm == 'SAMWG') {
    # Sequential Adaptive Metropolis-within-Gibbs (Specs required – no default)
    message('SAMWG requires user-supplied Specs. Dyn must be a matrix.')
    Specs_list = list(Dyn = NULL, Periodicity = 50)

  } else if (Algorithm == 'SGLD') {
    # Stochastic Gradient Langevin Dynamics (Specs required – no default)
    message('SGLD requires user-supplied Specs. Returning template.')
    Specs_list = list(epsilon = NULL, file = NULL, Nr = NULL, Nc = NULL,
                size = NULL)

  } else if (Algorithm == 'Slice') {
    # Slice Sampler
    Specs_list = list(B = NULL, Bounds = c(-Inf, Inf), m = Inf,
                Type = 'Continuous', w = 1)

  } else if (Algorithm == 'SMWG') {
    # Sequential Metropolis-within-Gibbs (Specs required – no default)
    message('SMWG requires user-supplied Specs. Dyn must be a matrix.')
    Specs_list = list(Dyn = NULL)

  } else if (Algorithm == 'THMC') {
    # Tempered Hamiltonian Monte Carlo (Specs required – no default)
    message('THMC requires user-supplied Specs. epsilon and m length must equal length(Initial.Values).')
    Specs_list = list(epsilon = 0.1, L = 2, m = NULL, Temperature = 1)

  } else if (Algorithm == 'twalk') {
    # t-walk. LaplacesDemon stops if SIV gives a non-finite posterior, or if
    # SIV equals Initial.Values after the model update.
    SIV = NULL
    if (!is.null(Lowlander)) {
      # Take a plausible point that is *not* the best one, to reduce the chance
      # of colliding with an Initial.Values that is already near the mode.
      Zt = Lowlander_Z(Lowlander, N = 2, parm.names = Data$parm.names)
      if (!is.null(Zt)) SIV = as.numeric(Zt[nrow(Zt), ])
    }
    if (is.null(SIV)) {
      message("twalk: pass a Lowlander object to auto-fill SIV; it must give a finite posterior and must differ from Initial.Values.")
    }
    Specs_list = list(SIV = SIV, n1 = 4, at = 6, aw = 1.5)

  } else if (Algorithm == 'UESS') {
    # Univariate Eigenvector Slice Sampler
    Specs_list = list(A = Inf, B = NULL, m = 100, n = 0)

  } else if (Algorithm == 'USAMWG') {
    # Updating Sequential Adaptive Metropolis-within-Gibbs (Specs required)
    message('USAMWG requires user-supplied Specs. Dyn must be a matrix.')
    Specs_list = list(Dyn = NULL, Periodicity = 1, Fit = NULL, Begin = NULL)

  } else if (Algorithm == 'USMWG') {
    # Updating Sequential Metropolis-within-Gibbs (Specs required)
    message('USMWG requires user-supplied Specs. Dyn must be a matrix.')
    Specs_list = list(Dyn = NULL, Fit = NULL, Begin = NULL)

  } else {
    stop(paste0('Unknown Algorithm: "', Algorithm, '". See LaplacesDemon documentation for supported algorithm names.'))
  }

  if (!is.null(Specs_list) && length(dots) > 0) {
    hit = names(dots) %in% names(Specs_list)
    Specs_list[names(dots)[hit]] = dots[hit]
    # Names that match no known element of this algorithm's Specs are silently
    # dropped by LaplacesDemon, so surface them here rather than let a typo
    # (e.g. NC = 8) pass unnoticed.
    if (any(!hit)) {
      warning(paste0("Specs_help: ignoring unknown argument(s) for ", Algorithm, ": ",
                     paste(names(dots)[!hit], collapse = ", "),
                     ". Valid names are: ", paste(names(Specs_list), collapse = ", "), "."))
    }
  } else if (is.null(Specs_list) && length(dots) > 0) {
    # DRM/RDMH take no Specs; a bare NULL keeps LaplacesDemon's own default path.
    warning(paste0('Specs_help: ', Algorithm, ' takes no Specs; the following ',
                   'argument(s) are ignored: ', paste(names(dots), collapse = ', '), '.'))
  }

  # User overrides above can desynchronise Z from Nc (LaplacesDemon stops on a
  # mismatch), so repair the pair here.
  if (Algorithm %in% c('AIES', 'DEMC') && !is.null(Specs_list$Z)) {
    Nc = max(3L, abs(round(Specs_list$Nc)))
    CPUs = if (is.null(Specs_list$CPUs)) 1 else Specs_list$CPUs
    if (Algorithm == 'AIES' && CPUs > 1 && Nc %% 2 != 0) {
      Nc = Nc + 1L
      message('AIES: Nc rounded up to ', Nc, ' (must be even when CPUs > 1).')
    }
    Specs_list$Nc = Nc
    if (Algorithm == 'AIES') {
      Z = Specs_list$Z
      if (is.matrix(Z) && nrow(Z) != Nc) {
        # Only rebuild from Lowlander when Z is ours; a user-supplied Z is
        # tiled/trimmed to match Nc rather than discarded.
        if ('Z' %in% names(dots)) {
          Specs_list$Z = matrix(as.numeric(Z[rep(seq_len(nrow(Z)), length.out = Nc), , drop = FALSE]),
                               nrow = Nc, ncol = ncol(Z))
        } else {
          Zn = Lowlander_Z(Lowlander, N = Nc, parm.names = Data$parm.names)
          if (!is.null(Zn)) Specs_list$Z = Zn
        }
      }
    } else if (!is.array(Specs_list$Z) || length(dim(Specs_list$Z)) != 3 ||
               dim(Specs_list$Z)[3] != Nc ||
               # identical(), not !=: Npar can be NULL here (no Data, no
               # Lowlander), and `dim(Z)[2] != NULL` is logical(0), which turns
               # the whole condition into NA and errors inside if().
               !identical(dim(Specs_list$Z)[2], Npar) ||
               (.sizeable && dim(Specs_list$Z)[1] != .z_rows(Iterations, Thinning))) {
      # DEMC wants a 3-D array: floor(Iterations/Thinning)+1 x Npar x Nc, with
      # Iterations/Thinning as LaplacesDemon *adjusts* them. A bare matrix
      # silently hits LaplacesDemon's rbind-recycle and dies with "subscript out
      # of bounds", so rebuild from whatever we have.
      Z = Specs_list$Z
      Zn = NULL
      if (.sizeable) {
        nrow_z = .z_rows(Iterations, Thinning)
        src = NULL
        if ('Z' %in% names(dots)) {
          if (is.matrix(Z) && ncol(Z) == Npar) src = Z
          if (is.null(src) && is.array(Z) && length(dim(Z)) == 3 &&
              dim(Z)[2] == Npar) src = Z[, , 1]
        } else {
          src = Lowlander_kept(Lowlander, parm.names = Data$parm.names)
        }
        if (!is.null(src)) {
          Zn = array(0, dim = c(nrow_z, Npar, Nc))
          for (i in seq_len(Nc)) {
            Zn[, , i] = if ('Z' %in% names(dots)) {
              src[rep(seq_len(nrow(src)), length.out = nrow_z), , drop = FALSE]
            } else {
              # Draw independently per chain: DEMC proposes gamma *
              # (Z[r1,,s1] - Z[r2,,s2]), so chains sharing an identical history
              # produce degenerate (zero) difference vectors.
              src[sample.int(nrow(src), nrow_z, replace = TRUE), , drop = FALSE]
            }
          }
        }
      }
      if (is.null(Zn)) {
        warning(paste0('DEMC: Z must be a 3-D array of dim ',
                       if (.sizeable) .z_rows(Iterations, Thinning) else 'floor(Iterations/Thinning)+1',
                       ' x ', Npar, ' x ', Nc, '. Set to NULL so LaplacesDemon ',
                       'regenerates it. Pass Iterations and Thinning to Specs_help ',
                       'to have the array built for you.'))
        Specs_list$Z = NULL
      } else {
        Specs_list$Z = Zn
      }
    }
  }

  return(Specs_list)
}
