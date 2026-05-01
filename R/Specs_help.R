#' Return the default Specs list for a LaplacesDemon algorithm
#'
#' @description
#' Given a LaplacesDemon algorithm name, \code{Specs_help} returns the
#' corresponding default \code{Specs} list that can be passed directly to
#' \code{\link[LaplacesDemon]{LaplacesDemon}}.  For algorithms that have no
#' built-in defaults (i.e. the user must supply all values), a template list
#' is returned with placeholder values and an informative \code{message} is
#' printed to guide the user.
#'
#' @param Algorithm Character scalar. The short algorithm name recognised by
#'   \code{LaplacesDemon}, e.g. \code{'CHARM'}, \code{'AIES'}, \code{'MWG'}.
#'   Defaults to \code{'CHARM'}.
#'
#' @return A named list suitable for passing as the \code{Specs} argument of
#'   \code{LaplacesDemon}, or \code{NULL} for algorithms that do not use
#'   \code{Specs} (e.g. \code{'DRM'} and \code{'RDMH'}).
#'
#' @examples
#' Specs_help('CHARM')
#' Specs_help('AIES')
#' Specs_help('MWG')
#'
#' @export
Specs_help = function(Algorithm = 'CHARM') {
  # Returns the default Specs list for the given LaplacesDemon Algorithm name.
  # For algorithms that require user-supplied values (no built-in defaults),
  # a template list is returned with placeholder values and a message is
  # printed to guide the user.

  if (Algorithm == 'ADMG') {
    # Adaptive Directional Metropolis-within-Gibbs
    return(list(n = 0, Periodicity = 1))

  } else if (Algorithm == 'AFSS') {
    # Automated Factor Slice Sampler
    return(list(A = Inf, B = NULL, m = Inf, n = 0, w = 1))

  } else if (Algorithm == 'AGG') {
    # Adaptive Griddy-Gibbs (Specs required – no default)
    message('AGG requires user-supplied Specs. Returning template.')
    return(list(Grid = NULL, dparm = NULL, smax = 0.1, CPUs = 1,
                Packages = NULL, Dyn.libs = NULL))

  } else if (Algorithm == 'AHMC') {
    # Adaptive Hamiltonian Monte Carlo
    # epsilon and m depend on the number of parameters; placeholders shown here
    message('AHMC: epsilon and m length must equal length(Initial.Values).')
    return(list(epsilon = 0.1, L = 2, m = NULL, Periodicity = 100))

  } else if (Algorithm == 'AIES') {
    # Affine-Invariant Ensemble Sampler (Specs required – no default)
    message('AIES requires user-supplied Specs. Returning template.')
    return(list(Nc = 3, Z = NULL, beta = 2, CPUs = 1,
                Packages = NULL, Dyn.libs = NULL))

  } else if (Algorithm == 'AM') {
    # Adaptive Metropolis
    return(list(Adaptive = 1000, Periodicity = 1))

  } else if (Algorithm == 'AMM') {
    # Adaptive-Mixture Metropolis
    return(list(Adaptive = 1000, B = NULL, n = 0, Periodicity = 1, w = 0.05))

  } else if (Algorithm == 'AMWG') {
    # Adaptive Metropolis-within-Gibbs
    return(list(B = NULL, n = 0, Periodicity = 50))

  } else if (Algorithm == 'CHARM') {
    # Componentwise Hit-And-Run Metropolis
    return(list(alpha.star = NA))

  } else if (Algorithm == 'DEMC') {
    # Differential Evolution Markov Chain (Specs required – no default)
    message('DEMC requires user-supplied Specs. Returning template.')
    return(list(Nc = 3, Z = NULL, gamma = NULL, w = 0.1))

  } else if (Algorithm == 'DRAM') {
    # Delayed Rejection Adaptive Metropolis
    return(list(Adaptive = 1000, Periodicity = 1))

  } else if (Algorithm == 'DRM') {
    # Delayed Rejection Metropolis (no Specs needed)
    return(NULL)

  } else if (Algorithm == 'ESS') {
    # Elliptical Slice Sampler
    return(list(B = NULL))

  } else if (Algorithm == 'GG') {
    # Griddy-Gibbs (Specs required – no default)
    message('GG requires user-supplied Specs. Returning template.')
    return(list(Grid = NULL, dparm = NULL, CPUs = 1,
                Packages = NULL, Dyn.libs = NULL))

  } else if (Algorithm == 'Gibbs') {
    # Gibbs Sampler (FC must be a function; MWG is optional)
    message('Gibbs requires FC to be a full-conditional function. Returning template.')
    return(list(FC = NULL, MWG = NULL))

  } else if (Algorithm == 'HARM') {
    # Hit-And-Run Metropolis
    return(list(alpha.star = NA, B = NULL))

  } else if (Algorithm == 'HMC') {
    # Hamiltonian Monte Carlo
    # epsilon and m depend on the number of parameters; placeholders shown here
    message('HMC: epsilon and m length must equal length(Initial.Values).')
    return(list(epsilon = 0.1, L = 2, m = NULL))

  } else if (Algorithm == 'HMCDA') {
    # Hamiltonian Monte Carlo with Dual-Averaging (Specs required – no default)
    message('HMCDA requires user-supplied Specs. Returning template.')
    return(list(A = 1000, delta = 0.65, epsilon = NULL, Lmax = 1000,
                lambda = 0.1))

  } else if (Algorithm == 'IM') {
    # Independence Metropolis (Specs required – no default)
    message('IM requires user-supplied Specs. mu must equal Initial.Values in length.')
    return(list(mu = NULL))

  } else if (Algorithm == 'INCA') {
    # Interchain Adaptation
    return(list(Adaptive = 1000, Periodicity = 1))

  } else if (Algorithm == 'MALA') {
    # Metropolis-Adjusted Langevin Algorithm
    # Note: gamma is required by LaplacesDemon but absent from its built-in
    # default; a typical value of 0.6 is used here.
    return(list(A = 1e7, alpha.star = 0.574, delta = 1, gamma = 0.6,
                epsilon = c(1e-6, 1e-7)))

  } else if (Algorithm == 'MCMCMC') {
    # Metropolis-Coupled Markov Chain Monte Carlo
    return(list(lambda = 1, CPUs = 1, Packages = NULL, Dyn.libs = NULL))

  } else if (Algorithm == 'MTM') {
    # Multiple-Try Metropolis
    return(list(K = 4, CPUs = 1, Packages = NULL, Dyn.libs = NULL))

  } else if (Algorithm == 'MWG') {
    # Metropolis-within-Gibbs (LaplacesDemon default algorithm)
    return(list(B = NULL))

  } else if (Algorithm == 'NUTS') {
    # No-U-Turn Sampler (Specs required – no default)
    message('NUTS requires user-supplied Specs. Returning template.')
    return(list(A = 1000, delta = 0.6, epsilon = NULL, Lmax = 1000))

  } else if (Algorithm == 'OHSS') {
    # Oblique Hyperrectangle Slice Sampler
    return(list(A = Inf, n = 0))

  } else if (Algorithm == 'pCN') {
    # Preconditioned Crank-Nicolson
    return(list(beta = 0.01))

  } else if (Algorithm == 'RAM') {
    # Robust Adaptive Metropolis
    return(list(alpha.star = 0.234, B = NULL, Dist = 'N', gamma = 0.66,
                n = 0))

  } else if (Algorithm == 'RDMH') {
    # Random Dive Metropolis-Hastings (no Specs needed)
    return(NULL)

  } else if (Algorithm == 'Refractive') {
    # Refractive Sampler
    return(list(Adaptive = 1, m = 2, w = 0.1, r = 1.3))

  } else if (Algorithm == 'RJ') {
    # Reversible-Jump (Specs required – no default)
    message('RJ requires user-supplied Specs. Returning template.')
    return(list(bin.n = 1, bin.p = 0.5, parm.p = 0.5,
                selectable = NULL, selected = NULL))

  } else if (Algorithm == 'RSS') {
    # Reflective Slice Sampler (Specs required – no default)
    message('RSS requires user-supplied Specs. Returning template.')
    return(list(m = 10, w = 1))

  } else if (Algorithm == 'RWM') {
    # Random-Walk Metropolis
    return(list(B = list()))

  } else if (Algorithm == 'SAMWG') {
    # Sequential Adaptive Metropolis-within-Gibbs (Specs required – no default)
    message('SAMWG requires user-supplied Specs. Dyn must be a matrix.')
    return(list(Dyn = NULL, Periodicity = 1))

  } else if (Algorithm == 'SGLD') {
    # Stochastic Gradient Langevin Dynamics (Specs required – no default)
    message('SGLD requires user-supplied Specs. Returning template.')
    return(list(epsilon = NULL, file = NULL, Nr = NULL, Nc = NULL,
                size = NULL))

  } else if (Algorithm == 'Slice') {
    # Slice Sampler
    return(list(B = NULL, Bounds = c(-Inf, Inf), m = Inf,
                Type = 'Continuous', w = 1))

  } else if (Algorithm == 'SMWG') {
    # Sequential Metropolis-within-Gibbs (Specs required – no default)
    message('SMWG requires user-supplied Specs. Dyn must be a matrix.')
    return(list(Dyn = NULL))

  } else if (Algorithm == 'THMC') {
    # Tempered Hamiltonian Monte Carlo (Specs required – no default)
    message('THMC requires user-supplied Specs. epsilon and m length must equal length(Initial.Values).')
    return(list(epsilon = 0.1, L = 2, m = NULL, Temperature = 1))

  } else if (Algorithm == 'twalk') {
    # t-walk
    return(list(SIV = NULL, n1 = 4, at = 6, aw = 1.5))

  } else if (Algorithm == 'UESS') {
    # Univariate Eigenvector Slice Sampler
    return(list(A = Inf, B = NULL, m = 100, n = 0))

  } else if (Algorithm == 'USAMWG') {
    # Updating Sequential Adaptive Metropolis-within-Gibbs (Specs required)
    message('USAMWG requires user-supplied Specs. Dyn must be a matrix.')
    return(list(Dyn = NULL, Periodicity = 1, Fit = NULL, Begin = NULL))

  } else if (Algorithm == 'USMWG') {
    # Updating Sequential Metropolis-within-Gibbs (Specs required)
    message('USMWG requires user-supplied Specs. Dyn must be a matrix.')
    return(list(Dyn = NULL, Fit = NULL, Begin = NULL))

  } else {
    stop(paste0('Unknown Algorithm: "', Algorithm, '". See LaplacesDemon documentation for supported algorithm names.'))
  }
}
