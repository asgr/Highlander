\name{Highlander}
\alias{Highlander}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Highlander
}
\description{
There Can be Only One! Finds the highest points in log-likelihood space by switching between CMA and MCMC phases. Defaults should work pretty well.
}
\usage{
Highlander(parm = NULL, Data, likefunc, likefunctype = NULL, liketype = NULL,
  seed = 666, lower = NULL, upper = NULL, dynlim = 2, ablim = 0,
  optim_iters = 2, Niters = c(100,100), NfinalMCMC = Niters[2], walltime = Inf,
  CMAargs = list(control=list(maxit=Niters[1])),
  LDargs = list(control=list(abstol=0.1), Iterations = Niters[2],
  Algorithm = 'CHARM', Thinning = 1)
  )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{parm}{
Numeric vector; usual parameter vector- input to \option{likefunc}
}
  \item{Data}{
List; usual data input to \option{likefunc} This is where you should put all additional arguments needs since \code{\link{LaplacesDemon}} cannot be passed additional likelihood function arguments.
}
  \item{likefunc}{
Function; takes arguments \option{parm} and \option{Data}; can output CMA scalar or LD list classes of outputs (this needs to be communicated by setting \option{likefunctype} as appropriate.
}
  \item{likefunctype}{
Character scalar; if \option{likefunc} outputs just the \code{abs(LP)} then 'CMA', if the list for LD then 'LD' (see \code{\link{LaplacesDemon}}). If left as NULL then it will attempt to guess the format based on the length of the \option{likefunc} output: if it is length 1 then it assumes \option{likefunctype} = 'CMA', otherwise it will assume 'LD'. This is usually correct.
}
  \item{liketype}{
Character scalar; is the likelihood output by \option{likefunc} defined in a way that means we want to minimise ('min'; Chi^2 type) or maximise ('max'; log-likelihood type) the function? If left as NULL the default assumes the standard minimisation setup of CMA (so 'min'), and maximisation setup of LD (so 'max').
}
  \item{seed}{
Integer scalar; random seed to start the \code{Highlander} function with.
}
  \item{lower}{
Numeric vector; lower limit allowed for \option{parm}; if NULL it becomes the minimum of \option{parm}*(1/\option{dynlim})  - \option{ablim} or \option{parm}*\option{dynlim} - \option{ablim}.
}
  \item{upper}{
Numeric vector; upper limit allowed for \option{parm}; if NULL it becomes the maximum of \option{parm}*(1/\option{dynlim})  + \option{ablim} or \option{parm}*\option{dynlim} + \option{ablim}.
}
  \item{dynlim}{
Numeric vector; dynamic limit range for automatically defined limits (effectively ignored if \option{dynlim} = 1). See \option{lower} and \option{upper} to see examples of how \option{dynlim} is used in practice.
}
  \item{ablim}{
Numeric vector; absolute limit range for automatically defined limits (effectively ignored if \option{ablim} = 0). See \option{lower} and \option{upper} to see examples of how \option{dynlim} is used in practice.
}
  \item{optim_iters}{
Integer scalar; number of CMA / LD loops.
}
  \item{Niters}{
Integer vector; number of iterations per CMA and LD (so can be vector length 2, otherwise value is repeated).
}
  \item{NfinalMCMC}{
Integer scalar; number of iterations to run for the final MCMC run. You might want to make this large if the idea is you are using the previous optimisation steps to find a good starting point, but you want to the final MCMC run using \code{\link{LaplacesDemon}} to provide detailed posterior inference.
}
  \item{walltime}{
Numeric scalar; the maximum allowed CPU time for \code{Highlander}. The function will stop with the best solution to date once the walltime has been exceeded (although it does not stop mid CMA or MCMC, so the \option{walltime} will usually be exceeded a bit).
}
  \item{CMAargs}{
List; arguments to be passed into \code{\link{cmaeshpc}}.
}
  \item{LDargs}{
List; arguments to be passed into \code{\link{LaplacesDemon}}.
}
}
\details{
\code{Highlander} is designed to find good global solutions to highly multi-modal likelihood problems. It does this by oscillating between genetic algorithm stages, and classical MCMC stages. The latter produces
}
\value{
A list output containing:

\item{parm }{Numeric vector; best \option{parm} of all iterations.}
\item{LP }{Numeric scalar; best log-posterior of all iterations.}
\item{diff }{Numeric scalar; LP difference between current best LP and last best LP (if large, might need more \option{optim_iters} and/or \option{Niters}).}
\item{best}{Character scalar; optimisation type of the best solution found, one of 'CMA' / 'LD_Median' / 'LD_Mean'.}
\item{iteration}{Integer scalar; iteration number of best solution.}
\item{CMA_last}{List; output of last CMA optimisation (see \code{\link{cmaeshpc}}).}
\item{LD_last}{List; output of last MCMC optimisation (see \code{\link{LaplacesDemon}}).}
\item{call}{List; the full call made to \code{Highlander}.}
\item{date}{Date item; the date and time of the call.}
\item{time}{Numeric scalar; run time in minutes.}
}
\author{
Aaron Robotham
}
\seealso{
\code{\link{cmaeshpc}}, \code{\link{LaplacesDemon}}
}
\examples{
# Nothing to see here yet!
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\concept{ optimisation }% use one of  RShowDoc("KEYWORDS")

