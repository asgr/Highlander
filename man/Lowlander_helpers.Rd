\name{Lowlander_helpers}
\alias{Lowlander_helpers}
\alias{Lowlander_kept}
\alias{Lowlander_Z}
\alias{Lowlander_best}
\title{
Extract usable starting values from a Lowlander result
}
\description{
Helpers that pull in-bounds, plausible parameter values out of the return value
of \code{\link{Lowlander}}, for seeding algorithms that need them. These are
used by \code{\link{Specs_help}} when a \option{Lowlander} object is supplied,
and are exported so they can be used directly.

\code{Lowlander_kept} returns the kept samples ranked best-first, \code{Lowlander_Z}
returns a design matrix of \option{N} such samples spread across the plausible
region, and \code{Lowlander_best} returns the single most plausible point.
}
\usage{
Lowlander_kept(Lowlander = NULL, liketype = NULL, parm.names = NULL)

Lowlander_Z(Lowlander = NULL, N = 1, liketype = NULL, parm.names = NULL)

Lowlander_best(Lowlander = NULL, parm.names = NULL)
}
\arguments{
  \item{Lowlander}{
List; the (invisible) return value of a \code{\link{Lowlander}} call.
\code{Lowlander_kept} and \code{Lowlander_Z} require the \code{latin_mod},
\code{keep} and \code{output} elements; \code{Lowlander_best} requires
\code{best}. If \code{NULL}, \code{NULL} is returned so callers can fall back to
their own defaults.
}
  \item{N}{
Integer scalar; number of rows (walker positions) requested from
\code{Lowlander_Z}. If \option{N} exceeds the number of kept samples the samples
are recycled and the repeats jittered, since exactly-duplicated walker positions
degenerate the affine-invariant stretch move.
}
  \item{liketype}{
Character scalar; \option{'min'} or \option{'max'}, the direction in which
\code{output} should be read. Only used as a fallback: when the
\option{Lowlander} object records the \code{liketype} it was run with (the case
for \code{Lowlander} versions that return it), that value takes precedence so the
ranking cannot disagree with how the samples were selected.
}
  \item{parm.names}{
Character vector; optional names for the columns of the result (and for the
returned vector from \code{Lowlander_best}). Only applied when its length matches
the number of parameters.
}
}
\details{
Lowlander restricts \code{keep} to samples with a finite likelihood that also
beat its quantile cut, so every row returned here lies inside the region
Lowlander judged plausible. This is what makes them safe as sampler seeds, in
contrast to \code{LaplacesDemon}'s own \code{GIV()} initial-value generator,
which proposes from broad Normal distributions and frequently returns
\code{NA} for likelihoods that are undefined outside a narrow region.

The rows of \code{Lowlander_Z} are taken at evenly-spaced quantile ranks of the
kept samples rather than as the best \option{N} points, so the ensemble covers
the width of the plausible region instead of collapsing onto the mode. In
testing the choice of rank spacing made little difference; what matters is
supplying in-bounds seeds at all rather than leaving \code{Z = NULL}.

\code{Lowlander_Z} requires at least two usable kept samples and warns (returning
\code{NULL}) when fewer are available.
}
\value{
\code{Lowlander_kept} and \code{Lowlander_Z} return a numeric matrix with one
column per parameter, or \code{NULL} if nothing usable is available.
\code{Lowlander_best} returns a numeric vector, or \code{NULL}.
}
\author{
Aaron Robotham
}
\seealso{
\code{\link{Lowlander}}, \code{\link{Specs_help}}, \code{\link{Highlander}}
}
\examples{
\dontrun{
likefunc = function(parm, Data) {
  sum((parm - c(2, 5))^2)
}
Data = list(parm.names = c('a', 'b'), N = 2)

LL = Lowlander(c(0, 0), c(10, 10), Data = Data, likefunc = likefunc,
               liketype = 'min', Nsamp = 500, pcut = 0.1)

# The kept samples, ranked best-first:
Lowlander_kept(LL, parm.names = Data$parm.names)

# 8 spread walker positions, suitable as Specs$Z for AIES:
Lowlander_Z(LL, N = 8, parm.names = Data$parm.names)

# The single best point, suitable as Initial.Values:
Lowlander_best(LL, parm.names = Data$parm.names)
}
}
\concept{ optimisation }
\concept{ MCMC }
