\name{Specs_help}
\alias{Specs_help}
\title{
Return the default Specs list for a LaplacesDemon algorithm
}
\description{
Given a \code{LaplacesDemon} algorithm name, \code{Specs_help} returns the
corresponding default \code{Specs} list that can be passed directly to
\code{\link[LaplacesDemon]{LaplacesDemon}}. For algorithms that have no
built-in defaults (i.e. the user must supply all values), a template list
is returned with placeholder values and an informative message is printed
to guide the user.
}
\usage{
Specs_help(Algorithm = 'CHARM', Data = NULL, Lowlander = NULL,
  Iterations = NULL, Thinning = 1, ...)
}
\arguments{
  \item{Algorithm}{
Character scalar; the short algorithm name recognised by
\code{\link[LaplacesDemon]{LaplacesDemon}}, e.g. \option{CHARM},
\option{AIES}, \option{MWG}. Defaults to \option{CHARM}. Supported values
are: ADMG, AFSS, AGG, AHMC, AIES, AM, AMM, AMWG, CHARM, DEMC, DRAM, DRM,
ESS, GG, Gibbs, HARM, HMC, HMCDA, IM, INCA, MALA, MCMCMC, MTM, MWG, NUTS,
OHSS, pCN, RAM, RDMH, Refractive, RJ, RSS, RWM, SAMWG, SGLD, Slice, SMWG,
THMC, twalk, UESS, USAMWG, USMWG.
}
  \item{Data}{
List; used to determine the number of parameters from \code{Data$parm.names},
which drives \option{Nc} for the ensemble samplers (\option{AIES}, \option{DEMC})
and the length of \option{epsilon}/\option{m} for the Hamiltonian family. If
\code{Data$parm.names} is absent, the parameter count is taken from the
\option{Lowlander} object instead.
}
  \item{Lowlander}{
List; the return value of a \code{\link{Lowlander}} call. When supplied, the
parts of \code{Specs} that \code{LaplacesDemon} would otherwise guess blind
(via its \code{GIV} initial-value generator) are filled from Lowlander's kept,
in-bounds samples. This is the recommended way to run \option{AIES},
\option{DEMC}, \option{twalk}, \option{IM}, \option{GG} and \option{AGG} on
likelihoods that are only well-behaved inside a restricted region, since
\code{GIV} routinely proposes wildly out-of-bounds points and the samplers then
either abort or stall. See Details.
}
  \item{Iterations}{
Integer scalar; the number of MCMC iterations that will be passed to
\code{\link[LaplacesDemon]{LaplacesDemon}}. Required to size the \option{Z}
array for \option{DEMC}, which must have
\code{floor(Iterations/Thinning) + 1} rows. If \code{NULL}, DEMC's \option{Z} is
left as \code{NULL} (and \code{LaplacesDemon} generates it) with a message.
}
  \item{Thinning}{
Integer scalar; the thinning interval corresponding to \option{Iterations}, used
only for sizing DEMC's \option{Z} array. Defaults to 1.
}
  \item{\dots}{
Named values used to override the defaults for a particular algorithm. See
Examples. Names that do not correspond to an element of that algorithm's
\code{Specs} produce a warning rather than being silently ignored.
}
}
\details{
For algorithms where \code{LaplacesDemon} has a built-in default, the
returned list mirrors that default exactly. For algorithms that require
the user to supply all \code{Specs} values (e.g. \option{AIES},
\option{NUTS}, \option{SGLD}), a template list with sensible placeholder
values is returned and a message is printed advising the user of what must
be filled in. The algorithms \option{DRM} and \option{RDMH} do not use a
\code{Specs} argument, so \code{NULL} is returned.

\subsection{Seeding from Lowlander}{
When \code{Specs$Z} (or \code{SIV}, \code{mu}, \code{Grid}) is left \code{NULL},
\code{LaplacesDemon} falls back to \code{GIV()}, which searches for usable
initial values by evaluating the model at random Normal draws with means in
\code{[-100, 100]} and standard deviations in \code{[0.1, 1000]}. On
astrophysical likelihoods that are undefined (or return non-finite values)
outside a much narrower region, \code{GIV} frequently gives up and returns
\code{NA}, which leaves ensemble walkers stacked at \code{NA} and the chain
either aborts or never approaches the mode.

Passing the output of a \code{\link{Lowlander}} call avoids this. Lowlander has
already identified a set of in-bounds, plausible parameter samples, so:
\itemize{
\item \option{AIES} \code{Z} becomes an \code{Nc} x \code{Npar} matrix of kept
samples taken at evenly-spaced quantile ranks (\code{\link{Lowlander_Z}}).
\item \option{DEMC} \code{Z} becomes the required 3-D array, with each chain
filled from an independent draw of the kept samples.
\item \option{twalk} \code{SIV} becomes a plausible point that is deliberately
not the best one, since \code{LaplacesDemon} stops if \code{SIV} coincides with
\code{Initial.Values}.
\item \option{IM} \code{mu} becomes the best kept point.
\item \option{GG}/\option{AGG} \code{Grid} becomes a per-parameter list of
\emph{offset} grids scaled to each parameter's own plausible width, rather than
\code{LaplacesDemon}'s fixed \code{seq(-0.1, 0.1)}.
}
}

\subsection{A good \code{Z} cannot rescue a bad \code{Initial.Values}}{
\code{LaplacesDemon} evaluates the model at \code{Initial.Values} before it ever
looks at \code{Z}, and stops outright with \code{"The posterior is infinite!"} if
that value is not finite. For \option{AIES} an explicitly supplied \code{Z} is
taken \emph{whole} -- rows 2 to \code{Nc} become walkers 2 to \code{Nc}, but
walker 1 starts from \code{Initial.Values} and \code{Z[1,]} is ignored. The
reported thinned chain is walker 1's trajectory, so it is \code{Initial.Values},
not \code{Z[1,]}, that anchors what you get back. (Walker 1 is not frozen at
\code{Initial.Values}: it proposes and moves like any other walker. Only the
\emph{starting} value is taken from there.)

\code{Highlander} passes a \code{Lowlander}-seeded \code{Z} through as-is, so if
\code{parm} is outside the region where the likelihood is finite the run can
abort regardless of how good \code{Z} is. When using Lowlander to bound the
problem, take \code{Initial.Values} from \code{Lowlander(...)$best} as well --
which is what \code{Highlander} does automatically if \code{parm} is left
\code{NULL}.
}
}
\value{
A named list suitable for passing as the \code{Specs} argument of
\code{\link[LaplacesDemon]{LaplacesDemon}}, or \code{NULL} for algorithms
that do not use \code{Specs} (i.e. \option{DRM} and \option{RDMH}).
}
\author{
Aaron Robotham
}
\seealso{
\code{\link[LaplacesDemon]{LaplacesDemon}}, \code{\link{Highlander}}
}
\examples{
# Default algorithm (CHARM):
Specs_help('CHARM')
# Use a non-default NA alpha.star (0.234 is suggested for multi-variate sampling):
Specs_help('CHARM', alpha.star=0.234)

# Affine-Invariant Ensemble Sampler:
Specs_help('AIES')

# Adaptive Metropolis:
Specs_help('AM')

# Metropolis-within-Gibbs:
Specs_help('MWG')

# Algorithms with no Specs return NULL:
Specs_help('DRM')

\dontrun{
# The interesting case: seed the ensemble samplers from a Lowlander run so that
# LaplacesDemon never has to guess initial values blind.
lower = c(-4, -2)
upper = c(3, 1)
Data = list(parm.names = c('mu', 'sigma'), mon.names = 'LP', N = 100)

# Note the two likefunc contracts: Highlander wants the LaplacesDemon-style
# list (LP, Dev, Monitor), while Lowlander collapses its likefunc to one
# numeric scalar per sample. Wrap it and pull out the LP, as here (assuming
# mylikefunc is the LD-style list-returning function; 'max' matches how
# Highlander's own liketype defaults to 'max' for LD-style outputs:
intervals = Lowlander(lower, upper, Data = Data, liketype = 'max',
                      likefunc = function(parm, Data) mylikefunc(parm, Data)$LP)

# Z now holds Nc in-bounds walker positions instead of NULL:
Specs_help('AIES', Data, Lowlander = intervals)

# DEMC needs Iterations/Thinning to size its 3-D Z array:
Specs_help('DEMC', Data, Lowlander = intervals, Iterations = 1000, Thinning = 1)

# And take Initial.Values from Lowlander's best point, because LaplacesDemon
# evaluates the model there before it ever looks at Z:
Highlander(intervals$best, Data, mylikefunc, lower = intervals$lower,
           upper = intervals$upper, Algorithm = 'AIES',
           Specs = Specs_help('AIES', Data, Lowlander = intervals))

# Simpler: Highlander does all of the above when given the Lowlander result,
# seeding Specs and filling parm/lower/upper wherever they are left NULL.
Highlander(Data = Data, likefunc = mylikefunc, Algorithm = 'AIES',
           Lowlander = intervals)
}
}
\concept{ MCMC }
