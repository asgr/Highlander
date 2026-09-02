Highlander=function(parm=NULL, Data, likefunc, likefunctype=NULL, liketype=NULL,
                    prior=NULL, Algorithm='CHARM', seed=666, lower=NULL, upper=NULL,
                    Lowlander=NULL,
                    applyintervals=TRUE, updateintervals=FALSE,
                    applyconstraints=TRUE, dynlim=2, ablim=0, optim_iters=2,
                    Niters=c(100,100), NfinalMCMC=Niters[2], walltime = Inf,
                    Specs=Specs_help(Algorithm, Data),
                    CMAargs=list(control=list(maxit=Niters[1])),
                    LDargs=list(control=list(abstol=0.1), Iterations=Niters[2], Algorithm=Algorithm,
                    Thinning=1, Specs=Specs), parm.names=NULL, keepall=FALSE, cores=1L
                    ){

  timestart = proc.time()[3] # start timer
  date = date()
  call = match.call(expand.dots=TRUE)

  # Capture what the user actually supplied *before* anything is forced. The
  # Specs and LDargs defaults reference each other (LDargs embeds Specs), so
  # evaluating one lazily evaluates the other and destroys the ability to tell
  # "absent" from "supplied as-is".
  parm_given = !missing(parm)
  specs_given = !missing(Specs)
  ldargs_given = !missing(LDargs)

  # Validate the Lowlander object up front, before the parallel dispatch below.
  # Each worker wraps its own call in try(), so a malformed object would
  # otherwise surface only as "every chain failed" rather than one clear error.
  if(!is.null(Lowlander)){
    if(!is.list(Lowlander) || is.null(Lowlander$lower) || is.null(Lowlander$upper) ||
       is.null(Lowlander$latin_mod) || is.null(Lowlander$keep) || is.null(Lowlander$output)){
      stop("'Lowlander' must be the return value of a Lowlander() call (needs $lower, $upper, $latin_mod, $keep, $output).")
    }
    Npar_ll = length(Lowlander$lower)
    arglens = list(parm=parm, lower=lower, upper=upper)
    for (nm in names(arglens)) {
      v = arglens[[nm]]
      if(!is.null(v) && length(v) != Npar_ll){
        stop(paste0("'", nm, "' has length ", length(v), " but the Lowlander object has ",
                    Npar_ll, " parameters."))
      }
    }
  }

  # Validate cores
  cores = as.integer(cores)
  if(is.na(cores) || cores < 1L){
    stop("'cores' must be an integer >= 1.")
  }

  # Parallel dispatch: run `cores` independent Highlander chains with staggered seeds
  # and return the best result (highest LP). Uses mclapply on Unix/macOS and a PSOCK
  # cluster on Windows so it is safe across all platforms.
  if(cores > 1L){
    if(length(seed) == cores){
      seeds = seed
    }else{
      # Take just the first to be safe.
      seeds = seed[1] + seq_len(cores) - 1L
    }

    run_args = list(
      parm=parm, Data=Data, likefunc=likefunc, likefunctype=likefunctype,
      liketype=liketype, prior=prior, lower=lower, upper=upper, Lowlander=Lowlander,
      applyintervals=applyintervals,
      updateintervals=updateintervals, applyconstraints=applyconstraints, dynlim=dynlim,
      ablim=ablim, optim_iters=optim_iters, Niters=Niters, NfinalMCMC=NfinalMCMC,
      walltime=walltime, CMAargs=CMAargs, parm.names=parm.names,
      keepall=keepall, cores=1L
    )
    # Algorithm must be forwarded explicitly now. Previously the whole default
    # LDargs was passed down and carried Algorithm (and Specs) inside it, which
    # is how the worker ended up running the requested algorithm; but forwarding
    # an already-evaluated LDargs would mark Specs as user-supplied inside the
    # worker and disable the Lowlander seeding, including the per-call DEMC
    # resize. So forward the algorithm itself, and only those of Specs/LDargs
    # the user actually wrote out.
    run_args$Algorithm = Algorithm
    if(specs_given){run_args$Specs = Specs}
    if(ldargs_given){run_args$LDargs = LDargs}

    if(.Platform$OS.type == "windows"){
      cl = parallel::makeCluster(cores)
      on.exit(parallel::stopCluster(cl), add=TRUE)
      parallel::clusterExport(cl, varlist="run_args", envir=environment())
      parallel::clusterEvalQ(cl, library(Highlander))
      results = parallel::parLapply(cl, seeds, function(seed_in){
        run_args$seed = seed_in
        try(do.call(Highlander, run_args), silent=TRUE)
      })
    } else {
      results = parallel::mclapply(seeds, function(seed_in){
        run_args$seed = seed_in
        try(do.call(Highlander, run_args), silent=TRUE)
      }, mc.cores=cores)
    }

    LP_vals = sapply(results, function(r){
      if(is.null(r) || inherits(r, "try-error") || is.null(r$LP) || is.na(r$LP) || !is.finite(r$LP)) {
        -Inf
      } else {
        r$LP
      }
    })
    best_idx = which.max(LP_vals)
    best_result = results[[best_idx]]

    # Patch the call, date and elapsed time from this top-level invocation
    best_result$call = call
    best_result$date = date
    best_result$time = (proc.time()[3] - timestart) / 60

    if(keepall){
      #Make sure the mon.names will match up with what we do internally
      Data[['mon.names']] = c("LP", Data[['mon.names']][! Data[['mon.names']] == 'LP'])
      best_result$LD_last_comb = try(
        LaplacesDemon::Combine(lapply(results, function(x) x$LD_last), Data=Data)
      )
      best_result$best_job = best_idx
      best_result$High_jobs = results
    }

    return(invisible(best_result))
  }

  # Inputs:

  # parm: usual parameter vector; input to likefunc
  # Data: usual data; input to likefunc
  # likefunc: likelihood function that takes in parm and Data as arguments; can output CMA scalar of LD list outputs
  # likefunctype: if likefunc outputs just the abs(LP) then 'CMA', if the list for LD then 'LD'
  # prior: optional function that takes parm (and Data) as arguments and returns the scalar
  #        log-likelihood of the prior, which is added to the log-posterior (LP) computed
  #        from likefunc, correctly propagated through both the CMA and LD stages
  # seed: random seed to start with
  # lower: lower limit vector
  # upper: upper limit vector
  # dynlim: dynamic range for auto limits (ignored if 1)
  # ablim: additional absolute range for auto limits (ignored if 0)
  # optim_iters: number of CMA / LD loops
  # Niters: iters per CMA and LD (so can be vector length 2, otherwise value is repeated)

  # Lowlander: optional result of a Lowlander() run, used to seed the parts of
  # Specs that LaplacesDemon would otherwise guess blind (its GIV() generator
  # probes from very wide Normals and returns NA on restricted likelihoods).
  # It also supplies parm/lower/upper, but *only* where the user left them NULL.
  if(!is.null(Lowlander)){
    if(!parm_given || is.null(parm)){
      # Prefer explicit parm.names, then the names Lowlander itself recorded,
      # then Data$parm.names (Lowlander ignores the latter, so an object built
      # without parm.names would otherwise hand back an unnamed start).
      ll_names = if(!is.null(parm.names)) parm.names
                 else if(!is.null(colnames(Lowlander$latin_mod))) colnames(Lowlander$latin_mod)
                 else Data$parm.names
      parm = Lowlander_best(Lowlander, parm.names = ll_names)
      if(is.null(parm)){
        stop("The Lowlander object has no $best to start from.")
      }
      parm_given = TRUE
    }else{
      # LaplacesDemon evaluates the model at Initial.Values before it considers
      # Z at all, and stops outright with "The posterior is infinite!" if that
      # is not finite. For AIES the reported thinned chain is walker 1's
      # trajectory and Z[1,] is ignored, so an out-of-region parm can abort the
      # run no matter how good the seeded Z is.
      if(any(parm < Lowlander$lower | parm > Lowlander$upper, na.rm = TRUE)){
        message('Highlander: parm sits outside the Lowlander plausible region; LaplacesDemon evaluates the model there first and stops if the posterior is not finite, so this can abort the run regardless of the seeded Specs. Consider parm = Lowlander(...)$best.')
      }
    }
    if(is.null(lower)){lower = Lowlander$lower}
    if(is.null(upper)){upper = Lowlander$upper}
  }

  if(is.null(parm) & !is.null(lower) & !is.null(upper)){
    parm = (lower + upper)/2
  }

  if(is.null(parm)){
    stop('parm is NULL!')
  }

  if(!is.null(parm.names)){
    names(parm) = parm.names
  }else if(!is.null(names(parm))){
    parm.names = names(parm)
  }

  Data[['applyintervals']] = applyintervals
  Data[['applyconstraints']] = applyconstraints

  if(is.null(lower)){
    if(!is.null(Data[['intervals']]$lo)){
      lower = Data[['intervals']]$lo
    }else{
      lower = parm*(1/dynlim)
      lower[(parm*dynlim)<lower] = (parm*dynlim)[(parm*dynlim)<lower]
      lower = lower - abs(ablim)
      if(applyintervals){
        Data[['intervals']]$lo = lower
      }else{
        lower[lower == 0 & ablim==0] = -Inf
      }
    }
  }else{
    if(applyintervals){Data[['intervals']]$lo = lower}
  }
  if(is.null(upper)){
    if(!is.null(Data[['intervals']]$hi)){
      upper = Data[['intervals']]$hi
    }else{
      upper = parm*dynlim
      upper[(parm*(1/dynlim))>upper] = (upper*(1/dynlim))[(upper*(1/dynlim))>upper]
      upper = upper + abs(ablim)
      if(applyintervals){
        Data[['intervals']]$hi = upper
      }else{
        upper[upper == 0 & ablim==0] = Inf
      }
    }
  }else{
    if(applyintervals){Data[['intervals']]$hi = upper}
  }

  if(any(lower == upper)){
    stop('lower and upper cannot have the same values!')
  }

  # likefunctype is detected by calling likefunc once at parm and looking at the
  # length of what comes back. Two problems with doing that raw: this is the
  # first evaluation in the whole run, so any error surfaces with no hint that
  # type detection is what tripped (and in multicore mode not one but `cores`
  # identical failures appear); and Data has not had mon.names / parm.names / N
  # filled in yet, which is exactly what an LD-style likefunc will be handed
  # later, so a function can fail here and then work fine when fitted. Retry
  # once against the augmented Data before giving up.
  if(is.null(likefunctype)){
    .probe = function(D) try(likefunc(parm, D), silent=TRUE)

    probe_out = .probe(Data)
    probe_dat = 'Data exactly as supplied'

    if(inherits(probe_out, 'try-error')){
      DataProbe = Data
      if(is.null(DataProbe[['mon.names']])){
        DataProbe[['mon.names']] = 'LP'
      }else{
        DataProbe[['mon.names']] = c('LP', DataProbe[['mon.names']][!DataProbe[['mon.names']] == 'LP'])
      }
      if(is.null(DataProbe[['parm.names']])){
        if(!is.null(parm.names) && length(parm.names) == length(parm)){
          DataProbe[['parm.names']] = parm.names
        }else{
          DataProbe[['parm.names']] = letters[seq_along(parm)]
        }
      }
      if(is.null(DataProbe[['N']])){DataProbe[['N']] = 1}

      probe_out2 = .probe(DataProbe)
      if(inherits(probe_out2, 'try-error')){
        stop(paste0(
          "likefunc() failed when Highlander called it once at the start parm to work out whether it ",
          "returns a CMA scalar or an LD list, so the fit has not started. Both the Data as supplied and ",
          "Data with mon.names/parm.names/N filled in were tried. likefunc must be able to run at parm = ",
          paste(if(is.null(names(parm))) sprintf('%g', parm)
                else sprintf('%s=%g', names(parm), parm), collapse = ', '),
          " (an out-of-region start is the usual cause: LaplacesDemon likewise evaluates the model at ",
          "Initial.Values before it does anything else, so it would stop there too). If likefunc legitimately ",
          "cannot be evaluated at this parm, set likefunctype explicitly ('CMA' for a single numeric return, ",
          "'LD' for the list with LP/Dev/Monitor) to skip the probe. Underlying error: ",
          conditionMessage(attr(probe_out2, 'condition'))))
      }
      probe_out = probe_out2
      probe_dat = 'Data with mon.names/parm.names/N filled in'
    }

    probe_len = length(probe_out)
    if(probe_len == 0){
      stop(paste0("likefunc(", probe_dat, ") returned nothing at the start parm, so Highlander cannot ",
                  "tell whether it is a CMA scalar or an LD list. It must return either a single numeric ",
                  "value or the LaplacesDemon list (LP, Dev, Monitor, ...); check that it returns a value ",
                  "on every path, rather than falling off the end or returning NULL."))
    }
    if(probe_dat != 'Data exactly as supplied'){
      message("Highlander: likefunc failed with the Data as supplied but worked once mon.names/parm.names/N ",
              "were present, so likefunctype = '", if(probe_len == 1) 'CMA' else 'LD', "' was detected from ",
              "the second call. Pass likefunctype explicitly to avoid the extra probe.")
    }

    if(probe_len == 1){
      likefunctype = 'CMA'
    }else{
      likefunctype = 'LD'
    }
  }

  if(is.null(liketype)){
    if(likefunctype == 'CMA'){liketype = 'min'}
    if(likefunctype == 'LD'){liketype = 'max'}
  }

  if(!is.null(Data$prior) && !is.null(prior)){
    stop('prior is provided in input Data and as an argument to prior! Resolve conflict before running Highlander.')
  }

  if(is.null(prior) && !is.null(Data$prior)){
    prior = Data$prior
    #Just to be safe, always return 0 if user likelihood expects to process prior function
    Data$prior = function(...) 0
  }

  DataCMA = Data

  if(likefunctype == 'CMA'){
    CMAfunc = function(parm, Data, inlikefunc=likefunc, inliketype=liketype, inprior=prior){
      .convert_CMA2CMA(parm=parm, Data=Data, likefunc=inlikefunc, liketype=inliketype, prior=inprior)
    }
  }else{
    DataCMA[['mon.names']] = ''
    CMAfunc = function(parm, Data, inlikefunc=likefunc, inliketype=liketype, inprior=prior){
      .convert_LD2CMA(parm=parm, Data=Data, likefunc=inlikefunc, liketype=inliketype, prior=inprior)
    }
  }

  DataLD = Data

  if(likefunctype == 'LD'){
    if(is.null(DataLD[['mon.names']])){
      DataLD[['mon.names']] = "LP"
    }else{
      DataLD[['mon.names']] = c("LP", DataLD[['mon.names']][! DataLD[['mon.names']] == 'LP'])
    }

    if(is.null(DataLD[['parm.names']])){
      DataLD[['parm.names']] = letters[1:length(parm)]
    }
    if(is.null(DataLD[['N']])){
      DataLD[['N']] = 1
    }
    LDfunc = function(parm, Data, inlikefunc=likefunc, inliketype=liketype, inprior=prior){
      .convert_LD2LD(parm=parm, Data=Data, likefunc=inlikefunc, liketype=inliketype, prior=inprior)
    }
  }else{
    DataLD[['mon.names']] = "LP"
    if(is.null(DataLD[['parm.names']])){
      if(is.null(parm.names)){
        DataLD[['parm.names']] = letters[1:length(parm)]
      }else if(length(parm) == length(parm.names)){
        DataLD[['parm.names']] = parm.names
      }else{
        message('parm.names does not match length of parm!')
        DataLD[['parm.names']] = letters[1:length(parm)]
      }
    }
    if(is.null(DataLD[['N']])){
      DataLD[['N']] = 1
    }
    LDfunc = function(parm, Data, inlikefunc=likefunc, inliketype=liketype, inprior=prior){
      .convert_CMA2LD(parm=parm, Data=Data, likefunc=inlikefunc, liketype=inliketype, prior=inprior)
    }
  }

  # Rebuild the default LDargs by hand rather than forcing the promise: the
  # signature's default embeds Specs, and evaluating it would print a second
  # round of Specs_help guidance when we re-seed below. Keep the raw Data for
  # that default -- DataLD has synthesized parm.names, which would change the
  # auto-derived Nc (and epsilon length) versus the Specs promise default.
  # When the user supplied LDargs we deliberately leave Specs absent, as before,
  # so LaplacesDemon still falls back to its own per-algorithm defaults.
  if(!ldargs_given){
    LDargs = list(control=list(abstol=0.1), Iterations=Niters[2],
                  Algorithm=Algorithm, Thinning=1)
    LDargs[['Specs']] = if(specs_given) Specs else Specs_help(Algorithm, Data)
  }

  # Normalize the pieces of LDargs that Highlander itself needs to know about.
  # A user-supplied LDargs may omit any of these.
  if(is.null(LDargs[['Iterations']])){LDargs[['Iterations']] = Niters[2]}
  if(is.null(LDargs[['Thinning']])){LDargs[['Thinning']] = 1}
  if(is.null(LDargs[['Algorithm']])){LDargs[['Algorithm']] = Algorithm}

  # Specs Highlander generated itself get seeded from the Lowlander samples.
  # Anything the user passed by hand is left exactly as-is, including a Specs
  # element embedded in their own LDargs.
  user_spec = specs_given || (ldargs_given && !is.null(LDargs[['Specs']]))
  seed_specs = (!user_spec && !is.null(Lowlander))
  preloop_iters = NULL
  if(seed_specs){
    preloop_iters = LDargs[['Iterations']]
    LDargs[['Specs']] = Specs_help(Algorithm, DataLD, Lowlander = Lowlander,
                                   Iterations = LDargs[['Iterations']],
                                   Thinning = LDargs[['Thinning']])
  }

  parm_out = parm

  CMA_out = NULL
  LD_out = list()

  CMA_all = list()
  LD_all = list()

  LP_out = -Inf
  diff = NA
  best = NA
  iteration = NA

  set.seed(seed)

  for(i in 1:ceiling(optim_iters)){
    message('Iteration ',i)

    time=(proc.time()[3] - timestart)/60
    if(time > walltime){break}

    if(Niters[1] > 0){
      tempsafe = try(
        do.call('cmaeshpc', c(list(par=parm_out, fn=CMAfunc, Data=quote(DataCMA), lower=lower,
             upper=upper), CMAargs))
      )

      if(inherits(tempsafe, "try-error")){
        message('CMA failed!')
        CMA_out = list(
          value = Inf,
          par = parm_out
        )
      }else{
        CMA_out = tempsafe
        if(is.null(CMA_out[['par']])){ #Catch bad initial starting positions and jitter
          tempsafe = try(
            do.call('cmaeshpc', c(list(par=jitter(parm_out), fn=CMAfunc, Data=quote(DataCMA), lower=lower,
                                       upper=upper), CMAargs))
          )
          if(inherits(tempsafe, "try-error")){
            message('CMA failed!')
            CMA_out = list(
              value = Inf,
              par = parm_out
            )
          }else{
            CMA_out = tempsafe
            if(is.null(CMA_out[['par']])){
              message('CMA is failing- something must be badly wrong!')
              return(NULL)
            }
          }
        }
      }

      if(is.finite(CMA_out[['value']])){
        if(updateintervals){
          # create new limits based on CMA
          hess = numDeriv::hessian(CMAfunc, x=CMA_out[['par']], method.args=list(eps=(upper-lower)/100, d=1), Data=Data)
          errors = sqrt(abs(diag(solve(hess))))

          CMA_out$hess = hess
          CMA_out$errors = errors

          lower_old = lower
          upper_old = upper

          lower = pmax(lower, CMA_out[['par']] - 5*errors)
          upper = pmin(upper, CMA_out[['par']] + 5*errors)

          if(applyintervals){
            DataLD[['intervals']]$lo = lower
            DataLD[['intervals']]$hi = upper
          }

          out_print = rbind(round(CMA_out[['par']],2), round(errors,2), round(lower_old,2), round(lower,2), round(upper_old,2), round(upper,2))
          colnames(out_print) = Data$parm.names
          rownames(out_print) = c('Best', 'Error', 'Low_old', 'Low_new', 'High_old', 'High_new')
          print(out_print)
        }
      }

      if(keepall){
        CMA_all = c(CMA_all, list(CMA_out))
      }

      if(i==1){
        diff = NA
        LP_out = -CMA_out[['value']]
        parm_out = CMA_out[['par']]
        best = 'CMA'
        iteration = i
        message('CMA ',i,': ',round(LP_out,3), ' ', paste(round(parm_out,3),collapse = ' '))
      }else{
        if(LP_out < -CMA_out[['value']]){ #this means new CMA is larger LP and preferred
          diff = abs(LP_out - -CMA_out[['value']])
          LP_out = -CMA_out[['value']]
          parm_out = CMA_out[['par']]
          best = 'CMA'
          iteration = i
          message('CMA ',i,': ',round(LP_out,3), ' ', paste(round(parm_out,3),collapse = ' '))
        }
      }
    }

    time = (proc.time()[3] - timestart)/60
    if(time > walltime){break}
    if(i > optim_iters){break}

    if(i == optim_iters){LDargs[['Iterations']] = NfinalMCMC}

    # The run length of this particular LaplacesDemon call is only known here,
    # and DEMC sizes Z as floor(Iterations/Thinning)+1, so Specs seeded above
    # (which used the pre-loop Iterations) are rebuilt at the true length.
    # Iterations/Thinning do not otherwise feed into Specs, so this is a no-op
    # in shape for every other algorithm.
    if(seed_specs && !identical(LDargs[['Iterations']], preloop_iters)){
      LDargs[['Specs']] = suppressMessages(Specs_help(Algorithm, DataLD, Lowlander = Lowlander,
                                     Iterations = LDargs[['Iterations']],
                                     Thinning = LDargs[['Thinning']]))
    }

    if(LDargs[['Iterations']] > 0){
      LD_out = do.call('LaplacesDemon', c(list(Model=LDfunc, Data=quote(DataLD),  Initial.Values=parm_out),
                            LDargs))

      LD_out$Model = NULL #don't want this in case it is big!
      LD_out$Call = NULL #don't want this in case it is big!

      if(keepall){
        LD_all = c(LD_all, list(LD_out))
      }

      if(LP_out < max(LD_out[['Monitor']][,'LP'])){ #this means new LD_Monitor is larger LP and preferred
        diff = abs(LP_out - max(LD_out[['Monitor']][,'LP']))
        LP_out = max(LD_out[['Monitor']][,'LP'])
        parm_out = LD_out$Posterior1[which.max(LD_out[['Monitor']][,'LP']),]
        best = 'LD_Mode'
        iteration = i
        message('LD Mode ',i,': ',round(LP_out,3), ' ', paste(round(parm_out,3),collapse = ' '))
      }

      if(LP_out < LD_out$Summary1['LP','Median']){ #this means new LD_Median is larger LP and preferred
        diff = abs(LP_out - LD_out$Summary1['LP','Median'])
        LP_out = LD_out$Summary1['LP','Median']
        parm_out = LD_out$Summary1[1:length(parm_out),'Median']
        best = 'LD_Median'
        iteration = i
        message('LD Median ',i,': ',round(LP_out,3), ' ', paste(round(parm_out,3),collapse = ' '))
      }

      if(LP_out < LD_out$Summary1['LP','Mean']){ #this means new LD_Mean is larger LP and preferred
        diff = abs(LP_out - LD_out$Summary1['LP','Mean'])
        LP_out = LD_out$Summary1['LP','Mean']
        parm_out = LD_out$Summary1[1:length(parm_out),'Mean']
        best = 'LD_Mean'
        iteration = i
        message('LD Mean ',i,': ',round(LP_out,3), ' ', paste(round(parm_out,3),collapse = ' '))
      }
    }
  }

  # Outputs:

  # parm: best parm of all iters
  # LP: best LP of all iters
  # diff: LP difference between current best LP and last best LP (if large, might need more optim_iters and/or Niters)
  # best: optim type of best solution, one of CMA / LD_Median / LD_Mean
  # iteration: iteration number of best solution (if last, might need more optim_iters and/or Niters)
  # CMA_last: last CMA output
  # LD_last: Last LD output

  if(applyconstraints & !is.null(Data[['constraints']])){
    parm_out = Data[['constraints']](parm_out)
  }

  if(applyintervals & !is.null(Data[['intervals']])){
    parm_out[parm_out < Data[['intervals']]$lo] = Data[['intervals']]$lo[parm_out < Data[['intervals']]$lo]
    parm_out[parm_out > Data[['intervals']]$hi] = Data[['intervals']]$hi[parm_out > Data[['intervals']]$hi]
  }

  RedChi2 = LP_out/(-1.418939 * DataLD[['N']])

  if(!is.null(parm.names)){
    names(parm_out) = parm.names
    if(!is.null(CMA_out)){
      names(CMA_out$par) = parm.names
    }
  }

  time = (proc.time()[3]-timestart)/60

  return(invisible(list(parm=parm_out, LP=LP_out, diff=diff, best=best, iteration=iteration,
                        CMA_last=CMA_out, LD_last=LD_out, N = DataLD[['N']], RedChi2 = RedChi2, call=call, date=date,
                        time=time, CMA_all=CMA_all, LD_all=LD_all)))
}

.convert_CMA2CMA=function(parm, Data, likefunc, liketype='min', prior=NULL){
  # Convert CMA type output to LD
  output = likefunc(parm, Data)
  if(liketype=='min'){
    fnscale = 1
  }else if(liketype=='max'){
    fnscale = -1
  }
  V = fnscale*output # value to be minimised (higher LP => lower V)
  if(!is.null(prior)){
    V = V - prior(parm, Data)
  }
  return(V)
}

.convert_CMA2LD=function(parm, Data, likefunc, liketype='min', prior=NULL){
  # Convert CMA type output to LD
  if(Data[['applyconstraints']] & !is.null(Data[['constraints']])){
    parm = Data[['constraints']](parm)
  }

  if(Data[['applyintervals']] & !is.null(Data[['intervals']]$lo) & !is.null(Data[['intervals']]$hi)){
    parm[parm<Data[['intervals']]$lo] = Data[['intervals']]$lo[parm<Data[['intervals']]$lo]
    parm[parm>Data[['intervals']]$hi] = Data[['intervals']]$hi[parm>Data[['intervals']]$hi]
  }
  output = likefunc(parm, Data)
  if(liketype=='min'){
    fnscale = -1
  }else if(liketype=='max'){
    fnscale = 1
  }
  LL = fnscale * output

  LP = LL
  if(!is.null(prior)){
    LP = LP + prior(parm, Data)
  }

  return(list(LP = LP, Dev = -2 * LL, Monitor = LP, yhat = 1,parm = parm))

}

.convert_LD2CMA=function(parm, Data, likefunc, liketype='min', prior=NULL){
  # Convert LD type output to CMA
  output = likefunc(parm, Data)
  if(liketype=='min'){
    fnscale = 1
  }else if(liketype=='max'){
    fnscale = -1
  }
  V = fnscale*output$LP # value to be minimised (higher LP => lower V)
  if(!is.null(prior)){
    V = V - prior(parm, Data)
  }
  return(V)
}

.convert_LD2LD=function(parm, Data, likefunc, liketype='min', prior=NULL){
  # Convert LD type output to LD
  if(Data[['applyconstraints']] & !is.null(Data[['constraints']])){
    parm = Data[['constraints']](parm)
  }

  if(Data[['applyintervals']] & !is.null(Data[['intervals']]$lo) & !is.null(Data[['intervals']]$hi)){
    parm[parm<Data[['intervals']]$lo] = Data[['intervals']]$lo[parm<Data[['intervals']]$lo]
    parm[parm>Data[['intervals']]$hi] = Data[['intervals']]$hi[parm>Data[['intervals']]$hi]
  }

  if(length(Data[['mon.names']]) > 1){
    # This is so we remove the leading LP internally since some
    # likefunc actually use the contents of mon.names to determine outputs
    Data[['mon.names']] = Data[['mon.names']][2:length(Data[['mon.names']])]
    useful_mon = TRUE
  }else{
    Data[['mon.names']] == ""
    useful_mon = FALSE
  }

  output = likefunc(parm, Data)

  if(liketype=='min'){
    fnscale = -1
  }else if(liketype=='max'){
    fnscale = 1
  }

  LP = fnscale*output$LP
  Dev = output$Dev #should not need scaling
  if(!is.null(prior)){
    LP = LP + prior(parm, Data)
  }

  if(useful_mon){
    # Need to check we don't also return LP elsewhere
    Monitor = c(LP, output$Monitor[!names(output$Monitor) == 'LP'])
  }else{
    Monitor = output$Monitor
  }

  if(!is.null(output[['parm']])){
    parm = output[['parm']]
  }

  # We add the expected LP back to the front of Monitor output
  return(list(LP=LP, Dev=Dev, Monitor=Monitor, yhat=output$yhat, parm=parm))
}
