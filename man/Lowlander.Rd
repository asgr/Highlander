\name{Lowlander}
\alias{Lowlander}
\title{
Lowlander: Latin-Hypercube Parameter-Space Shrink
}
\description{
Performs one iteration of a Latin-hypercube-based parameter-space shrink step.
Samples \option{Nsamp} parameter vectors uniformly across the box
[\option{lower}, \option{upper}] using a Latin Hypercube design, evaluates
\option{likefunc} at each sample (in parallel where requested), and returns
tightened bounds equal to the column-wise min/max of the most plausible
\option{pcut} fraction of samples.
}
\usage{
Lowlander(lower, upper, Nsamp = 100, pcut = 0.1,
          likefunc, Data, liketype = 'min',
          seed = NULL, latin = NULL, ncores = 1)
}
\arguments{
  \item{lower}{
Numeric vector; lower limits of the parameter box. Must be the same length as
\option{upper} and contain no NA values.
}
  \item{upper}{
Numeric vector; upper limits of the parameter box. Must be the same length as
\option{lower} and contain no NA values.
}
  \item{Nsamp}{
Positive integer scalar; number of Latin-hypercube samples to draw (default 100).
}
  \item{pcut}{
Numeric scalar strictly in (0, 1); fraction of samples to \emph{keep}. The
most plausible \option{pcut} fraction is retained (lowest values for
\code{liketype='min'}, highest for \code{liketype='max'}).
}
  \item{likefunc}{
Function; takes arguments \option{parm} (numeric vector) and \option{Data}
(list) and returns a single numeric scalar.
}
  \item{Data}{
List; passed as the second argument to \option{likefunc} at every evaluation.
}
  \item{liketype}{
Character scalar; \code{'min'} (default) keeps the \option{pcut} fraction of
samples with the \emph{smallest} \option{likefunc} values (Chi^2 / cost
convention). \code{'max'} keeps the fraction with the \emph{largest} values
(log-likelihood convention).
}
  \item{seed}{
Integer scalar or NULL (default); if not NULL, \code{set.seed(seed)} is called
immediately before generating the LHS design for reproducibility.
}
  \item{latin}{
Numeric matrix or NULL (default); a pre-computed \option{Nsamp} x \option{Npar}
Latin-hypercube design in the unit cube [0, 1]. If NULL a new design is
generated using \code{lhs::randomLHS} (requires the \pkg{lhs} package). When
supplied, the matrix is validated for correct dimensions and value range.
}
  \item{ncores}{
Positive integer scalar; number of CPU cores to use for parallel evaluation of
\option{likefunc} (default 1). On Windows a PSOCK cluster is used; on
Unix-like systems \code{parallel::mclapply} is used.
}
}
\details{
The workflow is:
\enumerate{
  \item Optionally set the RNG seed.
  \item Generate (or validate) an \option{Nsamp} x \option{Npar} Latin
        hypercube design in [0, 1].
  \item Map each column \emph{j} to physical coordinates:
        \code{lower[j] + latin[,j] * (upper[j] - lower[j])}.
  \item Evaluate \option{likefunc} at each row in parallel.
  \item Compute the quantile threshold at \option{pcut} (or \code{1-pcut}
        for \code{liketype='max'}) with \code{na.rm = TRUE}.
  \item Keep indices whose \option{likefunc} value is at or below (min) /
        above (max) the threshold.
  \item Require at least \code{Npar + 2} kept points so that meaningful
        column-wise min/max bounds can be computed.
  \item Return new bounds and bookkeeping information.
}
}
\value{
An invisible list containing:

\item{lower}{Numeric vector; new lower bounds (column-wise min of kept samples).}
\item{upper}{Numeric vector; new upper bounds (column-wise max of kept samples).}
\item{output}{Numeric vector; \option{likefunc} values for all \option{Nsamp} samples.}
\item{keep}{Integer vector; indices of the kept samples.}
\item{frac_keep}{Numeric scalar; \code{length(keep) / Nsamp}.}
\item{cutval}{Numeric scalar; quantile threshold used to determine \code{keep}.}
\item{latin}{Numeric matrix; the \option{Nsamp} x \option{Npar} unit-cube LHS design used.}
\item{latin_mod}{Numeric matrix; the \option{Nsamp} x \option{Npar} parameter matrix mapped to physical bounds.}
\item{kept_pars}{Numeric matrix; \code{latin_mod[keep, , drop=FALSE]} — the kept parameter vectors.}
}
\author{
Aaron Robotham
}
\seealso{
\code{\link{Highlander}}, \code{\link[lhs]{randomLHS}}
}
\examples{
# Simple 2-parameter example: find bounds around the minimum of a quadratic
likefunc = function(parm, Data) {
  sum((parm - Data$target)^2)
}

Data = list(target = c(2, 5))

result = Lowlander(
  lower = c(0, 0), upper = c(10, 10),
  Nsamp = 200, pcut = 0.1,
  likefunc = likefunc, Data = Data,
  liketype = 'min', seed = 42
)

# New bounds should be tighter and centred near c(2, 5):
result$lower
result$upper
result$frac_keep
}
\concept{ optimisation }
