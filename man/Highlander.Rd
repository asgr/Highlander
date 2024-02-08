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
  seed = 666, lower = NULL, upper = NULL, applyintervals = TRUE, applyconstraints = TRUE,
  dynlim = 2, ablim = 0, optim_iters = 2, Niters = c(100,100), NfinalMCMC = Niters[2],
  walltime = Inf, CMAargs = list(control=list(maxit=Niters[1])),
  LDargs = list(control=list(abstol=0.1), Iterations = Niters[2],
  Algorithm = 'CHARM', Thinning = 1), parm.names = NULL, keepall = FALSE
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
Numeric vector; lower limit allowed for \option{parm}; if \option{lower} = NULL and Data$intervals$lo = NULL it becomes the minimum of \option{parm}*(1/\option{dynlim})  - \option{ablim} or \option{parm}*\option{dynlim} - \option{ablim}; else if Data$intervals$lo is not NULL then this will be used instead.
}
  \item{upper}{
Numeric vector; upper limit allowed for \option{parm}; if \option{upper} = NULL and Data$intervals$hi = NULL it becomes the maximum of \option{parm}*(1/\option{dynlim})  + \option{ablim} or \option{parm}*\option{dynlim} + \option{ablim}; else if Data$intervals$hi is not NULL then this will be used instead.
}
  \item{applyintervals}{
Logical; should the lower and upper intervals be applied by \code{Highlander} directly? This means it will create the list items Data$intervals$lo = lower and Data$intervals$hi = upper. You might want to set this to FALSE when your likelihood function deals with intervals itself using a different structure to the Data$intervals and creating such structures will break functionality inside \option{likefunc} (e.g. this is the case when running with \code{ProFit}).
}
  \item{applyconstraints}{
Logical; should constraints be applied by \code{Highlander} directly? This means it will pass the current \option{parm} through Data$contraints. You might want to set this to FALSE when your likelihood function deals with constraints itself and applying such structures will break functionality inside \option{likefunc} (e.g. this is the case when running with \code{ProFit}).
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
Integer vector; number of iterations per CMA and LD respectively (must be length 2).
}
  \item{NfinalMCMC}{
Integer scalar; number of iterations to run for the final MCMC run. You might want to make this large if the idea is you are using the previous optimisation steps to find a good starting point, but you want to the final MCMC run using \code{\link{LaplacesDemon}} to provide detailed posterior inference.
}
  \item{walltime}{
Numeric scalar; the maximum allowed CPU time for \code{Highlander} in minutes. The function will stop with the best solution to date once the walltime has been exceeded (although it does not stop mid CMA or MCMC, so the \option{walltime} will usually be exceeded a bit).
}
  \item{CMAargs}{
List; arguments to be passed into \code{\link{cmaeshpc}}.
}
  \item{LDargs}{
List; arguments to be passed into \code{\link{LaplacesDemon}}.
}
  \item{parm.names}{
Character vector; optional vector of parameter names. This will over-write the current names in \option{parm} if these also exist. If provided it ensures the \option{parm} and \option{CMA_last$par} outputs are named.
}
  \item{keepall}{
Logical; if FALSE (default) then does nothing. If TRUE then output objects \option{CMAall} and \option{LDall} will be concatenated lists of all CMA and LD iterations respectively.
}
}
\details{
\code{Highlander} is designed to find good global solutions to highly multi-modal likelihood problems. It does this by oscillating between genetic algorithm stages, and classical MCMC stages. The latter produces
}
\value{
A list output containing:

\item{parm}{Numeric vector; best \option{parm} of all iterations.}
\item{LP}{Numeric scalar; best log-posterior of all iterations.}
\item{diff}{Numeric scalar; LP difference between current best LP and last best LP (if large, might need more \option{optim_iters} and/or \option{Niters}).}
\item{best}{Character scalar; optimisation type of the best solution found, one of 'CMA' / 'LD_Median' / 'LD_Mean'.}
\item{iteration}{Integer scalar; iteration number of best solution.}
\item{CMA_last}{List; output of last CMA optimisation (see \code{\link{cmaeshpc}}).}
\item{LD_last}{List; output of last MCMC optimisation (see \code{\link{LaplacesDemon}}).}
\item{N}{Integer scalar; number of observations (approximately the degrees of freedom).}
\item{RedChi2}{Numeric scalar; the approximate reduced Chi^2 of the fit (more for reference than use).}
\item{call}{List; the full call made to \code{Highlander}.}
\item{date}{Date item; the date and time of the call.}
\item{time}{Numeric scalar; run time in minutes.}
\item{CMA_all}{If \option{keepall} is TRUE then a list of all the CMA runs.}
\item{LD_all}{If \option{keepall} is TRUE then a list of all the LD runs.}
}
\author{
Aaron Robotham
}
\seealso{
\code{\link{cmaeshpc}}, \code{\link{LaplacesDemon}}
}
\examples{
# Create some fake data (note it needs to be a list, but data need not be called 'x'):

mock = list(x=rnorm(1e3, mean=20, sd=100))

# Make a likelihood evaluator (this is the CMA style simple likelihood):

likefunc = function(parm, Data){
 sum(dnorm(Data$x, parm[1], parm[2], log=TRUE), na.rm=TRUE)
}

# Fit it!

highfit = Highlander(parm=c(-1,2), Data=mock, likefunc=likefunc, liketype='max',
  lower=c(-1e3,1), upper=c(1e3,1e3), parm.names=c('Mean', 'SD'))

# Should be near to 20 (mean) 100 (SD):

highfit$parm

# And we can see the impact of tighter intervals:

highfit = Highlander(parm=c(-1,2), Data=mock, likefunc=likefunc, liketype='max',
  lower=c(19.9,99), upper=c(20.1,101), parm.names=c('Mean', 'SD'))

highfit$CMA_last$par

highfit$LD_last$Summary1
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\concept{ optimisation }% use one of  RShowDoc("KEYWORDS")

