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
  Algorithm = "CHARM", seed = 666, lower = NULL, upper = NULL,
  applyintervals = TRUE, updateintervals = FALSE, applyconstraints = TRUE,
  dynlim = 2, ablim = 0, optim_iters = 2, Niters = c(100,100), NfinalMCMC = Niters[2],
  walltime = Inf, Specs = Specs_help(Algorithm, Data),
  CMAargs = list(control=list(maxit=Niters[1])),
  LDargs = list(control = list(abstol = 0.1), Iterations = Niters[2],
    Algorithm = Algorithm, Thinning = 1, Specs = Specs),
  parm.names = NULL, keepall = FALSE, cores = 1L,
  jitter = NULL, jitter_init = NULL, jitter_lower = -30, jitter_upper = 5
  )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{parm}{
Numeric vector; usual parameter vector- input to \option{likefunc}
}
  \item{Data}{
List; usual data input to \option{likefunc} This is where you should put all additional arguments needs since \code{\link[LaplacesDemon]{LaplacesDemon}} cannot be passed additional likelihood function arguments.
}
  \item{likefunc}{
Function; takes arguments \option{parm} and \option{Data}; can output CMA scalar or LD list classes of outputs (this needs to be communicated by setting \option{likefunctype} as appropriate.
}
  \item{likefunctype}{
Character scalar; if \option{likefunc} outputs just the \code{abs(LP)} then 'CMA', if the list for LD then 'LD' (see \code{\link[LaplacesDemon]{LaplacesDemon}}). If left as NULL then it will attempt to guess the format based on the length of the \option{likefunc} output: if it is length 1 then it assumes \option{likefunctype} = 'CMA', otherwise it will assume 'LD'. This is usually correct.
}
  \item{liketype}{
Character scalar; is the likelihood output by \option{likefunc} defined in a way that means we want to minimise ('min'; Chi^2 type) or maximise ('max'; log-likelihood type) the function? If left as NULL the default assumes the standard minimisation setup of CMA (so 'min'), and maximisation setup of LD (so 'max').
}
  \item{Algorithm}{
Character scalar; the short algorithm name recognised by
\code{\link[LaplacesDemon]{LaplacesDemon}}, e.g. \option{CHARM},
\option{AIES}, \option{MWG}. Defaults to \option{CHARM}. Supported values
are: ADMG, AFSS, AGG, AHMC, AIES, AM, AMM, AMWG, CHARM, DEMC, DRAM, DRM,
ESS, GG, Gibbs, HARM, HMC, HMCDA, IM, INCA, MALA, MCMCMC, MTM, MWG, NUTS,
OHSS, pCN, RAM, RDMH, Refractive, RJ, RSS, RWM, SAMWG, SGLD, Slice, SMWG,
THMC, twalk, UESS, USAMWG, USMWG.
}
  \item{seed}{
Integer scalar/vector; random seed to start the \code{Highlander} function with. If using multiple cores and the length of \option{seed} equals the number of cores, then those values will be used, else the seeds will become \option{seed} + 0:(\option{cores} - 1), i.e. the first job will use the \option{seed} provided, and the second job becomes \option{seed} + 1 etc.
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
  \item{updateintervals}{
Logical; can intervals be updated based on results of computing the Hessian after each CMA fit? You might never want the \option{lower} and \option{upper} intervals (and the related intervals structure inside \option{data}) to be modified, in which case leave as the default of FALSE. If you do want to update then the marginalised standard deviation is computed from the Hessian matrix, and x3 each value become the new limit if it sits inside the old limit (otherwise the old limit is used instead). When this is TRUE, a matrix of relevant update metrics is printed to screen so users can check it is behaving reasonably.
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
Integer scalar; number of iterations to run for the final MCMC run. You might want to make this large if the idea is you are using the previous optimisation steps to find a good starting point, but you want to the final MCMC run using \code{\link[LaplacesDemon]{LaplacesDemon}} to provide detailed posterior inference.
}
  \item{walltime}{
Numeric scalar; the maximum allowed CPU time for \code{Highlander} in minutes. The function will stop with the best solution to date once the walltime has been exceeded (although it does not stop mid CMA or MCMC, so the \option{walltime} will usually be exceeded a bit).
}
  \item{Specs}{
List; \option{Specs} to be passed into \option{LDargs} as a convenience. \option{Specs} is by default populated with the settings obtained from \code{\link{Specs_help}}, if non-defaults are desired then extra arguments can be passed into \code{\link{Specs_help}} via \dots or a fully user-defined list can be constructed obeying the requirements of \code{\link[LaplacesDemon]{LaplacesDemon}} for the \option{Algorithm} being used.
}
  \item{CMAargs}{
List; arguments to be passed into \code{\link[cmaeshpc]{cmaeshpc}}.
}
  \item{LDargs}{
List; arguments to be passed into \code{\link[LaplacesDemon]{LaplacesDemon}}. Note that \option{Algorithm} and \option{Specs} are passed into the \option{LDargs} by default.
}
  \item{parm.names}{
Character vector; optional vector of parameter names. This will over-write the current names in \option{parm} if these also exist. If provided it ensures the \option{parm} and \option{CMA_last$par} outputs are named.
}
  \item{keepall}{
Logical; if FALSE (default) then does nothing extra. If TRUE then: in serial mode (\option{cores} = 1) output objects \option{CMA_all} and \option{LD_all} will be concatenated lists of all CMA and LD iterations respectively; in parallel mode (\option{cores} > 1) the full list of all chain outputs is also returned as \option{High_jobs}, making it easy to inspect or compare individual chains.
}
  \item{cores}{
Positive integer scalar; number of CPU cores to use (default 1). When \option{cores} > 1, \code{Highlander} runs \option{cores} independent optimisation chains in parallel, each starting from a different random seed (\option{seed}, \option{seed}+1, ..., \option{seed}+\option{cores}-1), and returns the chain that achieved the highest LP. On Windows a PSOCK cluster is used (requires the \pkg{Highlander} package to be installed); on Unix-like systems (Linux, macOS) \code{parallel::mclapply} is used. Setting \option{cores} = 1 (the default) gives identical serial behaviour to previous versions and requires no special setup.
}
  \item{jitter}{
Character scalar/vector or NULL (default); when not NULL, names the field of \option{Data} that contains the per-observation sigma (error) values, e.g. \option{"sigma"} for \code{Data$sigma}. If a vector then this is the list path to the field. When set, \code{Highlander} appends an extra parameter \code{log_jitter} to \option{parm} and, before every likelihood evaluation, inflates the named sigma field by \code{sqrt(sigma^2 + exp(log_jitter)^2)}. This introduces the correct additive-variance jitter term to the log-posterior without modifying \option{likefunc}. Both the CMA and MCMC phases jointly optimise \code{log_jitter} together with the other parameters; the best-fit jitter in natural units is \code{exp(result$parm["log_jitter"])}.
}
  \item{jitter_init}{
Numeric scalar or NULL (default); initial value for the \code{log_jitter} parameter. When NULL the value is set to \code{log(median(abs(Data[[jitter]])))} so that the jitter starts at the median sigma scale.
}
  \item{jitter_lower}{
Numeric scalar; lower bound for \code{log_jitter} (default -30, i.e. a jitter near zero).
}
  \item{jitter_upper}{
Numeric scalar; upper bound for \code{log_jitter} (default 5, i.e. a jitter of about 150 in natural units).
}
}
\details{
\code{Highlander} is designed to find good global solutions to highly multi-modal likelihood problems. It does this by oscillating between genetic algorithm stages, and classical MCMC stages.
}
\value{
A list output containing:

\item{parm}{Numeric vector; best \option{parm} of all iterations.}
\item{LP}{Numeric scalar; best log-posterior of all iterations.}
\item{diff}{Numeric scalar; LP difference between current best LP and last best LP (if large, might need more \option{optim_iters} and/or \option{Niters}).}
\item{best}{Character scalar; optimisation type of the best solution found, one of 'CMA' / 'LD_Median' / 'LD_Mean'.}
\item{iteration}{Integer scalar; iteration number of best solution.}
\item{CMA_last}{List; output of last CMA optimisation (see \code{\link[cmaeshpc]{cmaeshpc}}).}
\item{LD_last}{List; output of last MCMC optimisation (see \code{\link[LaplacesDemon]{LaplacesDemon}}).}
\item{N}{Integer scalar; number of observations (approximately the degrees of freedom).}
\item{RedChi2}{Numeric scalar; the approximate reduced Chi^2 of the fit (more for reference than use).}
\item{call}{List; the full call made to \code{Highlander}.}
\item{date}{Date item; the date and time of the call.}
\item{time}{Numeric scalar; run time in minutes.}
\item{CMA_all}{If \option{keepall} is TRUE then a list of all the CMA runs.}
\item{LD_all}{If \option{keepall} is TRUE then a list of all the LD runs.}
\item{LD_last_comb}{List; present when \option{cores} > 1 and \option{keepall} = TRUE. The combined outputs of all \option{LD_last} outputs from the various Highlander jobs run. This is probably what you want to use for exploring combined posteriors etc.}
\item{best_job}{Integer scalar; present when \option{cores} > 1 and \option{keepall} = TRUE. Index (1-based) of the parallel chain that achieved the highest LP.}
\item{High_jobs}{List; present when \option{cores} > 1 and \option{keepall} = TRUE. Contains the full \code{Highlander} output list for every parallel chain, in chain order. Use \option{best_job} to select the winning chain: \code{High_out$High_jobs[[High_out$best_job]]}.}

When run in multicore mode \option{parm} through to \option{LD_all} inclusive are the results of the 'best' run. \option{LD_last_comb}, \option{best_job and \option{High_jobs} are only present when \option{cores} > 1 and \option{keepall} = TRUE.}
}
\author{
Aaron Robotham
}
\seealso{
\code{\link[cmaeshpc]{cmaeshpc}}, \code{\link[LaplacesDemon]{LaplacesDemon}}
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

