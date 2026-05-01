Highlander=function(parm=NULL, Data, likefunc, likefunctype=NULL, liketype=NULL,
                    seed=666, lower=NULL, upper=NULL, applyintervals=TRUE, applyconstraints=TRUE, dynlim=2, ablim=0, optim_iters=2,
                    Niters=c(100,100), NfinalMCMC=Niters[2], walltime = Inf,
                    CMAargs=list(control=list(maxit=Niters[1])),
                    LDargs=list(control=list(abstol=0.1), Iterations=Niters[2], Algorithm='CHARM',
                    Thinning=1), parm.names=NULL, keepall=FALSE, cores=1L
                    ){

  timestart = proc.time()[3] # start timer
  date = date()
  call = match.call(expand.dots=TRUE)

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
      liketype=liketype, lower=lower, upper=upper, applyintervals=applyintervals,
      applyconstraints=applyconstraints, dynlim=dynlim, ablim=ablim,
      optim_iters=optim_iters, Niters=Niters, NfinalMCMC=NfinalMCMC,
      walltime=walltime, CMAargs=CMAargs, LDargs=LDargs, parm.names=parm.names,
      keepall=keepall, cores=1L
    )

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
      best_result$LD_last_comb = LaplacesDemon::Combine(lapply(results,\(x) x$LD_last), Data=Data)
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
  # seed: random seed to start with
  # lower: lower limit vector
  # upper: upper limit vector
  # dynlim: dynamic range for auto limits (ignored if 1)
  # ablim: additional absolute range for auto limits (ignored if 0)
  # optim_iters: number of CMA / LD loops
  # Niters: iters per CMA and LD (so can be vector length 2, otherwise value is repeated)

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

  if(is.null(likefunctype)){
    if(length(likefunc(parm, Data)) == 1){
      likefunctype = 'CMA'
    }else{
      likefunctype = 'LD'
    }
  }

  if(is.null(liketype)){
    if(likefunctype == 'CMA'){liketype = 'min'}
    if(likefunctype == 'LD'){liketype = 'max'}
  }

  DataCMA = Data

  if(likefunctype == 'CMA'){
    CMAfunc = function(parm, Data, inlikefunc=likefunc, inliketype=liketype){
      .convert_CMA2CMA(parm=parm, Data=Data, likefunc=inlikefunc, liketype=inliketype)
    }
  }else{
    DataCMA[['mon.names']] = ''
    CMAfunc = function(parm, Data, inlikefunc=likefunc, inliketype=liketype){
      .convert_LD2CMA(parm=parm, Data=Data, likefunc=inlikefunc, liketype=inliketype)
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
    LDfunc = function(parm, Data, inlikefunc=likefunc, inliketype=liketype){
      .convert_LD2LD(parm=parm, Data=Data, likefunc=inlikefunc, liketype=inliketype)
    }
  }else{
    DataLD[['mon.names']] = "LP"
    if(is.null(DataLD[['parm.names']])){
      DataLD[['parm.names']] = letters[1:length(parm)]
    }
    if(is.null(DataLD[['N']])){
      DataLD[['N']] = 1
    }
    LDfunc = function(parm, Data, inlikefunc=likefunc, inliketype=liketype){
      .convert_CMA2LD(parm=parm, Data=Data, likefunc=inlikefunc, liketype=inliketype)
    }
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
        if(LP_out < -CMA_out[['value']]){
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

    if(LDargs[['Iterations']] > 0){
      LD_out = do.call('LaplacesDemon', c(list(Model=LDfunc, Data=quote(DataLD),  Initial.Values=parm_out),
                            LDargs))

      LD_out$Model = NULL #don't want this in case it is big!
      LD_out$Call = NULL #don't want this in case it is big!

      if(keepall){
        LD_all = c(LD_all, list(LD_out))
      }

      if(LP_out < max(LD_out[['Monitor']][,1])){
        diff = abs(LP_out - max(LD_out[['Monitor']][,1]))
        LP_out = max(LD_out[['Monitor']][,1])
        parm_out = LD_out$Posterior1[which.max(LD_out[['Monitor']][,1]),]
        best = 'LD_Mode'
        iteration = i
        message('LD Mode ',i,': ',round(LP_out,3), ' ', paste(round(parm_out,3),collapse = ' '))
      }

      if(LP_out < LD_out$Summary1['LP','Median']){
        diff = abs(LP_out - LD_out$Summary1['LP','Median'])
        LP_out = LD_out$Summary1['LP','Median']
        parm_out = LD_out$Summary1[1:length(parm_out),'Median']
        best = 'LD_Median'
        iteration = i
        message('LD Median ',i,': ',round(LP_out,3), ' ', paste(round(parm_out,3),collapse = ' '))
      }

      if(LP_out < LD_out$Summary1['LP','Mean']){
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

.convert_CMA2CMA=function(parm, Data, likefunc, liketype='min'){
  # Convert CMA type output to LD
  output = likefunc(parm, Data)
  if(liketype=='min'){
    fnscale = 1
  }else if(liketype=='max'){
    fnscale = -1
  }
  return(fnscale*output)
}

.convert_CMA2LD=function(parm, Data, likefunc, liketype='min'){
  # Convert CMA type output to LD
  if(Data[['applyconstraints']] & !is.null(Data[['constraints']])){
    parm = Data[['constraints']](parm)
  }

  if(Data[['applyconstraints']] & !is.null(Data[['intervals']]$lo) & !is.null(Data[['intervals']]$hi)){
    parm[parm<Data[['intervals']]$lo] = Data[['intervals']]$lo[parm<Data[['intervals']]$lo]
    parm[parm>Data[['intervals']]$hi] = Data[['intervals']]$hi[parm>Data[['intervals']]$hi]
  }
  output = likefunc(parm, Data)
  if(liketype=='min'){
    fnscale = -1
  }else if(liketype=='max'){
    fnscale = 1
  }
  return(list(LP=fnscale*output, Dev=-fnscale*2*output, Monitor=fnscale*output, yhat=1, parm=parm))
}

.convert_LD2CMA=function(parm, Data, likefunc, liketype='min'){
  # Convert LD type output to CMA
  output = likefunc(parm, Data)
  if(liketype=='min'){
    fnscale = 1
  }else if(liketype=='max'){
    fnscale = -1
  }
  return(fnscale*output$LP)
}

.convert_LD2LD=function(parm, Data, likefunc, liketype='min'){
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

  if(useful_mon){
    Monitor = c(fnscale*output$LP, output$Monitor)
  }else{
    Monitor = output$Monitor
  }

  if(!is.null(output[['parm']])){
    parm = output[['parm']]
  }

  # We add the expected LP back to the front of Monitor output
  return(list(LP=fnscale*output$LP, Dev=fnscale*output$Dev, Monitor=Monitor, yhat=output$yhat, parm=parm))
}
