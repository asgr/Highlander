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
Specs_help(Algorithm = 'CHARM', Data = NULL, ...)
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
List; AIES can generate sensible default \option{Nc} if provided with
LD style \option{Data}. This defaults to 2 x length(\option{parm.names}).
}
  \item{\dots}{
List; Named list where we wish to over-ride the defaults for a particular algorithm. See Examples.
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

# Affine-Invariant Ensemble Sampler (requires user values; prints a message):
Specs_help('AIES')

# Algorithms with no Specs return NULL:
Specs_help('DRM')
}
\concept{ MCMC }
