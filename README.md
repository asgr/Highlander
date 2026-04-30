# Highlander (R package)

<!-- badges: start -->
![R-CMD-check](https://github.com/asgr/Highlander/workflows/R-CMD-check/badge.svg)
<!-- badges: end -->

## Synopsis

There Can be Only One! Finds the highest points in log-likelihood space by switching between CMA (Covariance Matrix Adaptation) and MCMC (Markov Chain Monte Carlo) phases. Defaults should work pretty well.

In recent work fitting spectral energy distributions with **ProSpect** we nearly always use **Highlander** for rapid and well explored fitting.

## MCMC Algorithm Selection

**Highlander** delegates MCMC sampling to the **LaplacesDemon** package, which provides dozens of algorithms (see `?LaplacesDemon::LaplacesDemon`). The default is CHARM (Componentwise Hit-And-Run Metropolis), which is a reliable general-purpose choice. Drawing on MCMC literature and the LaplacesDemon tutorial, the following guidelines can help when selecting an alternative via the `Algorithm` field of the `LDargs` argument:

**4–8 parameters:** CHARM (default) performs well at low dimensions. NUTS (No-U-Turn Sampler) and HMCDA are very efficient when the log-posterior is differentiable and smooth. RAM (Robust Adaptive Metropolis) and Slice sampling are also strong choices, and the t-walk (twalk) requires no gradient information and is easy to apply. AIES also works well in this domain, and is popular in astronomy.

**8–20 parameters:** Gradient-based samplers become increasingly valuable as dimension grows. NUTS and HMCDA are recommended when a differentiable log-posterior is available. DEMC (Differential Evolution Markov Chain) handles correlated parameters well in this range. AFSS (Automated Factor Slice Sampler) is the LaplacesDemon-recommended general-purpose algorithm and adapts automatically. CHARM remains viable but may converge more slowly than gradient-based alternatives. AIES also works well in this domain, and is popular in astronomy.

**20–40 parameters:** NUTS continues to perform well when gradients are accessible. DEMC and AFSS (used with parameter blocks) are strong gradient-free options. AMWG (Adaptive Metropolis-within-Gibbs) with carefully selected blocks scales efficiently at this dimensionality. CHARM can still be used but is likely to mix slowly; blockwise versions such as AMM or AMWG are preferable.

**40-100 parameters:** fitting will be getting very hard, especially if likelihoods are not auto-differentialable (if they are, then probably you should be using **Stan** or similar at this point). The best approch will take some user experimentation, but CHARM is often not a bad starting point.

To switch algorithm, pass `Algorithm` via `LDargs`:

```R
Highlander(..., LDargs=list(Algorithm='**SEE BELOW**', Specs = list(**SEE BELOW**), Iterations=1000, Thinning=1))
```

The `LaplacesDemon::Juxtapose` function can compare algorithm efficiencies on a simplified version of your problem before committing to a full run.

## Installation

### Getting Highlander

Source installation from GitHub should be easy:

```R
install.packages('remotes')
remotes::install_github("asgr/Highlander")
library(Highlander)
```

A few Mac people seem to have issues with the above due to the backend used to download files. A work around seems to be to either use devtools (which I do not use as the default since it has a low more dependencies, and is tricky to install on HPCs):

```R
install.packages('devtools')
devtools::install_github("asgr/Highlander")
library(Highlander)
```

Or try the following:

```R
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS="true")
remotes::install_github("asgr/Highlander")
```

I also have these options set by default in my .Rprofile, which seems to help with some of the remote install issues some people face:

```R
options(download.file.method = "libcurl")
options(repos="http://cran.rstudio.com/")
options(rpubs.upload.method = "internal")
```

If all of these do not work than the nuclear option is to download (or clone) the GitHub repo, cd to where the tar.gz file is and run in the **console** (or **Terminal** on Mac):

```console
R CMD install Highlander_X.Y.Z.tar.gz
```

where X, Y and Z should be set as appropriate for the version downloaded (check the name of the file basically).

If none of the above works then you should consider burning your computer in sacrifice to the IO Gods. Then buy a newer *better* computer, and try all the above steps again.

Failing all of the above, please email me for help (or perhaps raise an Issue here, if it really does not seem like a local issue).

#### Package Dependencies

The above should also install the required packages. If you have trouble with this you can try installing the required packages manually first and then retry the installation for **Highlander**:

```R
install.packages(c('knitr', 'magicaxis', 'ProSpect', 'celestial', 'foreach', 'rmarkdown')) # CRAN packages
install.packages('remotes')
remotes::install_github("asgr/LaplacesDemon")
remotes::install_github("asgr/cmaeshpc")
remotes::install_github("asgr/Highlander")
```

Assuming you have installed all of the packages that you need/want, you should now be able to load **Highlander** within **R** with the usual:

```R
library(Highlander)
```

## Code Example

```R
library(Highlander)

# Simple example: fit a Gaussian
Data = list(
  x = seq(-3, 3, len=100),
  y = dnorm(seq(-3, 3, len=100), mean=0.5, sd=1.2) + rnorm(100, sd=0.01),
  intervals = list(lo=c(-5, 0.01), hi=c(5, 5))
)

likefunc = function(parm, Data){
  model = dnorm(Data$x, mean=parm[1], sd=parm[2])
  LP = sum(dnorm(Data$y, mean=model, sd=0.01, log=TRUE))
  return(list(LP=LP, Dev=-2*LP, Monitor=parm, yhat=model))
}

result = Highlander(parm=c(0, 1), Data=Data, likefunc=likefunc,
                    lower=c(-5, 0.01), upper=c(5, 5))
```

To find more long-form examples, including complicated fitting use-cases, please check the vignettes provided. You can browse these with:

```R
browseVignettes('Highlander')
```

Or if that does not work they are all hosted externally at <http://rpubs.com/asgr/>

## Motivation

**Highlander** was developed to provide a robust and efficient optimisation strategy for fitting complex models, in particular spectral energy distributions with **ProSpect**. It combines the global search power of CMA-ES with the statistical sampling of MCMC, switching between phases to find the highest posterior probability solution quickly and reliably.

## Contributors

Aaron Robotham

## License

LGPL-3+
