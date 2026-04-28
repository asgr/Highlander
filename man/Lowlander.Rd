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
Lowlander(lower, upper, Nsamp = 'auto', pcut = 0.1, corcut = 0, likefunc, Data,
  liketype = 'min', seed = 666, latin = NULL, ncores = 1)
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
Positive integer scalar; number of Latin-hypercube samples to draw. Default 'auto' will
scale the number of paramters (i.e. length of \option{lower}) by x50 or use 2000
(whichever is smaller). You probably do not want to go a lot lower than this threshold practically.
}
  \item{pcut}{
Numeric scalar strictly in [0-1]; fraction of samples to \emph{keep}. The
most plausible \option{pcut} fraction is retained (lowest values for
\code{liketype='min'}, highest for \code{liketype='max'}).
}
  \item{corcut}{
Numeric scalar; the threshold of absolute correlation required between parameters and the
\option{output} [0-1]. Larger values means more correlation is necessary before any interval
shrinkage is allowed. If used, sensible values are in the range 0.1-0.3 (much larger than this
and parameters will often not be allowed to shrink that could reasonably be). Default (0)
means this test is not required, but the output \option{cor_data} still contains the relevant
correlation matrix that allows users to inspect the covariance post running \code{Lowlander}.
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
Latin-hypercube design in the unit cube [0-1]. If NULL a new design is
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
        hypercube design in [0-1].
  \item Map each column \emph{j} to physical coordinates:
        \code{lower[j] + latin[,j] * (upper[j] - lower[j])}.
  \item Evaluate \option{likefunc} at each row in parallel.
  \item Compute the quantile threshold at \option{pcut} (or \code{1-pcut}
        for \code{liketype='max'}) with \code{na.rm = TRUE}.
  \item Keep indices whose \option{likefunc} value is at or below (min) /
        above (max) the threshold.
  \item Require at least \code{Npar + 2} kept points so that meaningful
        column-wise min/max bounds can be computed.
  \item Optional require the correlation meets the threshold required by \option{corcut}.
  \item Return new \option{lower} and \option{upper} bounds and bookkeeping information.
}

Conceptually this function is used to create reasonable bounds for \code{\link{Highlander}} if the user really has no sense of plausible parameter regimes.
}
\value{
An invisible list containing:

\item{lower}{Numeric vector; new lower bounds (column-wise min of kept samples).}
\item{upper}{Numeric vector; new upper bounds (column-wise max of kept samples).}
\item{best}{Numeric vector; current 'best' parameters in \option{latin_mod}.}
\item{output}{Numeric vector; \option{likefunc} values for all \option{Nsamp} samples.}
\item{keep}{Integer vector; indices of the kept samples.}
\item{bad}{Integer vector; indices of the non-finite \option{output}.}
\item{vol_red}{Numeric scalar; fraction of original lower-upper bounds new volume covers [0-1]. Lower values indicate a greater reduction factor.}
\item{frac_keep}{Numeric scalar; \code{length(keep) / Nsamp}.}
\item{cutval}{Numeric scalar; quantile threshold used to determine \code{keep}.}
\item{cor_data}{Numeric matric; correlation matrix between all parameters and \option{output} (latter is the last row and column).}
\item{latin}{Numeric matrix; the \option{Nsamp} x \option{Npar} unit-cube LHS design used.}
\item{latin_mod}{Numeric matrix; the \option{Nsamp} x \option{Npar} parameter matrix mapped to physical bounds.}
}
\author{
Aaron Robotham
}
\seealso{
\code{\link{Highlander}}, \code{\link[lhs]{randomLHS}}
}
\examples{
\dontrun{
library(lhs)

# Simple 2-parameter example: find bounds around the minimum of a quadratic
likefunc = function(parm, Data) {
  sum((parm - Data)^2)
}

Data = c(2, 5) #Try with more numbers and see how rapidly this becomes 'hard'
Ndata = length(Data)

result = Lowlander(
  lower = rep(0, Ndata), upper = rep(10, Ndata), Nsamp = 1000, pcut = 0.05,
  likefunc = likefunc, Data = Data, liketype = 'min'
)

# New bounds should be tighter and centred near c(2, 5):
result$lower
result$upper
result$best

# It is interesting to see how much convergence it is possible to get for more parameters:

Data = 1:9
Ndata = length(Data)

lower = rep(0, Ndata)
upper = rep(10, Ndata)

for(i in 1:8){
    message('Iteration ', i)
    result = Lowlander(
        lower = lower, upper = upper, Nsamp = 1000, pcut = 0.05, corcut = 0.3,
        likefunc = likefunc, Data = Data, liketype = 'min'
    )
    lower = result$lower
    upper = result$upper
    message("Volume reduction: ", result$vol_red, "\nNew Lower: ", paste(round(lower,2), collapse=' '), '\nNew Upper: ', paste(round(upper,2), collapse=' '))
}

result$lower
result$upper
result$best

# So after only 4 iterations we cannot reliably shrink the limits safely!
# Notice how even after one iteration we have reduced our search volume by 75\%

# Total volume reduction

prod((upper - lower)) / 10^Ndata
}
}
\concept{ optimisation }
