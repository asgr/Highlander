Highlander=function(parm=NULL, Data, likefunc, likefunctype=NULL, liketype=NULL,
                    Algorithm='CHARM', seed=666, lower=NULL, upper=NULL,
                    applyintervals=TRUE, updateintervals=FALSE,
                    applyconstraints=TRUE, dynlim=2, ablim=0, optim_iters=2,
                    Niters=c(100,100), NfinalMCMC=Niters[2], walltime = Inf,
                    Specs=Specs_help(Algorithm, Data),
                    CMAargs=list(control=list(maxit=Niters[1])),
                    LDargs=list(control=list(abstol=0.1), Iterations=Niters[2], Algorithm=Algorithm,
                    Thinning=1, Specs=Specs), parm.names=NULL, keepall=FALSE, cores=1L,
                    jitter=NULL, jitter_init=NULL, jitter_lower=-30, jitter_upper=5
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
      updateintervals=updateintervals, applyconstraints=applyconstraints, dynlim=dynlim,
      ablim=ablim, optim_iters=optim_iters, Niters=Niters, NfinalMCMC=NfinalMCMC,
      walltime=walltime, CMAargs=CMAargs, LDargs=LDargs, parm.names=parm.names,
      keepall=keepall, cores=1L,
      jitter=jitter, jitter_init=jitter_init, jitter_lower=jitter_lower, jitter_upper=jitter_upper
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
      best_result$LD_last_comb = try(
        LaplacesDemon::Combine(lapply(results,\(x) x$LD_last), Data=Data)
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
    LDfunc = function(parm, Data, inlikefunc=likefunc, inliketype=liketype){
      .convert_CMA2LD(parm=parm, Data=Data, likefunc=inlikefunc, liketype=inliketype)
    }
  }

  # Set up jitter parameter if requested.
  # When jitter is not NULL it names the field in Data containing per-observation
  # sigma values. An extra parameter log_jitter is appended to parm so that both
  # CMA and LD jointly optimise it. Before each likefunc call the wrapper inflates
  # Data[[jitter]] by sqrt(sigma^2 + exp(log_jitter)^2), producing the correct
  # additive-variance jitter term without modifying likefunc itself.
  if (!is.null(jitter)) {
    if (!is.character(jitter) || length(jitter) != 1L) {
      stop("'jitter' must be NULL or a single character string naming the sigma field in Data.")
    }
    if (is.null(Data[[jitter]])) {
      stop(paste0("'Data$", jitter, "' not found. 'jitter' must name a field in Data ",
                  "containing per-observation sigma (error) values."))
    }
    if (is.null(jitter_init)) {
      jitter_init = log(median(abs(Data[[jitter]]), na.rm = TRUE))
    }

    # Append log_jitter to parm and update bounds (CMA uses explicit lower/upper).
    # DataLD$intervals is intentionally left at N_orig length so the inner wrapper
    # functions clip only the original parameters.
    parm = c(parm, jitter_init)
    names(parm)[length(parm)] = 'log_jitter'
    if (!is.null(parm.names)) {
      parm.names = c(parm.names, 'log_jitter')
    }
    lower = c(lower, jitter_lower)
    upper = c(upper, jitter_upper)

    # Update DataLD parm.names to include log_jitter
    DataLD[['parm.names']] = c(DataLD[['parm.names']], 'log_jitter')

    # Capture sigma field name in closures
    .jitter_field = jitter

    # Wrap CMAfunc: strip log_jitter, inflate sigma, call inner, return scalar
    .inner_CMAfunc = CMAfunc
    CMAfunc = function(parm, Data, ...) {
      log_jitter_val = parm[length(parm)]
      parm = parm[-length(parm)]
      s = exp(log_jitter_val)
      Data[[.jitter_field]] = sqrt(Data[[.jitter_field]]^2 + s^2)
      .inner_CMAfunc(parm, Data, ...)
    }

    # Wrap LDfunc: strip log_jitter, inflate sigma, call inner, restore log_jitter to parm
    .inner_LDfunc = LDfunc
    LDfunc = function(parm, Data, ...) {
      log_jitter_val = parm[length(parm)]
      parm = parm[-length(parm)]
      s = exp(log_jitter_val)
      Data[[.jitter_field]] = sqrt(Data[[.jitter_field]]^2 + s^2)
      result = .inner_LDfunc(parm, Data, ...)
      result$parm = c(result$parm, log_jitter_val)
      result
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
            # When jitter is active, DataLD$intervals must stay at N_orig length so
            # the inner wrapper functions can clip original parameters without a
            # dimension mismatch.  Jitter bounds are kept only in lower/upper (CMA).
            if (!is.null(jitter)) {
              n_orig = length(lower) - 1L
              DataLD[['intervals']]$lo = lower[seq_len(n_orig)]
              DataLD[['intervals']]$hi = upper[seq_len(n_orig)]
            } else {
              DataLD[['intervals']]$lo = lower
              DataLD[['intervals']]$hi = upper
            }
          }

          out_print = rbind(round(CMA_out[['par']],2), round(errors,2), round(lower_old,2), round(lower,2), round(upper_old,2), round(upper,2))
          colnames(out_print) = parm.names
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
    if (!is.null(jitter)) {
      # User constraints do not know about log_jitter; apply only to original params
      log_jitter_end = parm_out[length(parm_out)]
      parm_out = c(Data[['constraints']](parm_out[-length(parm_out)]), log_jitter_end)
    } else {
      parm_out = Data[['constraints']](parm_out)
    }
  }

  if(applyintervals & !is.null(Data[['intervals']])){
    # Data$intervals has N_orig elements; clip only those to avoid dimension mismatch
    # when jitter has been appended to parm_out.
    n_clip = length(Data[['intervals']]$lo)
    parm_clip = parm_out[seq_len(n_clip)]
    parm_clip[parm_clip < Data[['intervals']]$lo] = Data[['intervals']]$lo[parm_clip < Data[['intervals']]$lo]
    parm_clip[parm_clip > Data[['intervals']]$hi] = Data[['intervals']]$hi[parm_clip > Data[['intervals']]$hi]
    parm_out[seq_len(n_clip)] = parm_clip
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
    # Need to check we don't also return LP elsewhere
    Monitor = c(fnscale*output$LP, output$Monitor[!names(output$Monitor) == 'LP'])
  }else{
    Monitor = output$Monitor
  }

  if(!is.null(output[['parm']])){
    parm = output[['parm']]
  }

  # We add the expected LP back to the front of Monitor output
  return(list(LP=fnscale*output$LP, Dev=fnscale*output$Dev, Monitor=Monitor, yhat=output$yhat, parm=parm))
}
