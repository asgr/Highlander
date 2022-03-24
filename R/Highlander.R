Highlander=function(parm=NULL, Data, likefunc, likefunctype=NULL, liketype=NULL,
                    seed=666, lower=NULL, upper=NULL, applyintervals=TRUE, applyconstraints=TRUE, dynlim=2, ablim=0, optim_iters=2,
                    Niters=c(100,100), NfinalMCMC=Niters[2], walltime = Inf,
                    CMAargs=list(control=list(maxit=Niters[1])),
                    LDargs=list(control=list(abstol=0.1), Iterations=Niters[2], Algorithm='CHARM',
                    Thinning=1), parm.names=NULL, keepall=FALSE
                    ){

  timestart = proc.time()[3] # start timer
  date = date()
  call = match.call(expand.dots=TRUE)

  # Inputs:

  # parm: usual parameter vector; input to likefunc
  # Data: usual data; input to likefunc
  # likefunc: likelihood function that takes in parm and Data as arguments; can output CMA scalar of LD list outputs
  # likefunctype: if likefunc outputs just the abs(LP) then 'CMA', if the list for LD then 'LD'
  # seed: random seed to start with
  # lower: lower limit vector; if NULL it becomes the minimum of parm*(1/dynlim) | parm*dynlim - ablim
  # upper: upper limit vector; if NULL it becomes the maximum of parm*(1/dynlim) | parm*dynlim + ablim
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
    DataLD[['mon.names']] = c('LP_mon', DataLD[['mon.names']])
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
    DataLD[['mon.names']] = 'LP'
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

  CMAout = list()
  LDout = list()

  CMA_all = list()
  LD_all = list()

  set.seed(seed)

  for(i in 1:ceiling(optim_iters)){
    message('Iteration ',i)

    time=(proc.time()[3]-timestart)/60
    if(time > walltime){break}

    tempsafe = try(
      do.call('cmaeshpc', c(list(par=parm_out, fn=CMAfunc, Data=quote(DataCMA), lower=lower,
           upper=upper), CMAargs))
    )

    if(class(tempsafe)=="try-error"){
      message('CMA failed!')
      CMAout = list(
        value = -Inf,
        par = parm_out
      )
    }else{
      CMAout = tempsafe
      if(is.null(CMAout[['par']])){ #Catch bad initial starting positions and jitter
        tempsafe = try(
          do.call('cmaeshpc', c(list(par=jitter(parm_out), fn=CMAfunc, Data=quote(DataCMA), lower=lower,
                                     upper=upper), CMAargs))
        )
        if(class(tempsafe)=="try-error"){
          message('CMA failed!')
          CMAout = list(
            value = -Inf,
            par = parm_out
          )
        }else{
          CMAout = tempsafe
          if(is.null(CMAout[['par']])){
            message('CMA is failing- something must be badly wrong!')
            return(NULL)
          }
        }
      }
    }

    if(keepall){
      CMA_all = c(CMA_all, list(CMAout))
    }

    if(i==1){
      diff=NA
      LP_out = -CMAout[['value']]
      parm_out = CMAout[['par']]
      best = 'CMA'
      iteration = i
      message('CMA ',i,': ',round(LP_out,3), ' ', paste(round(parm_out,3),collapse = ' '))
    }else{
      if(LP_out < -CMAout[['value']]){
        diff = abs(LP_out - -CMAout[['value']])
        LP_out = -CMAout[['value']]
        parm_out = CMAout[['par']]
        best = 'CMA'
        iteration = i
        message('CMA ',i,': ',round(LP_out,3), ' ', paste(round(parm_out,3),collapse = ' '))
      }
    }

    time=(proc.time()[3]-timestart)/60
    if(time > walltime){break}
    if(i > optim_iters){break}

    if(i == optim_iters){LDargs[['Iterations']] = NfinalMCMC}

    if(LDargs[['Iterations']] > 0){
      LDout = do.call('LaplacesDemon', c(list(Model=LDfunc, Data=quote(DataLD),  Initial.Values=parm_out),
                            LDargs))

      if(keepall){
        LD_all = c(LD_all, list(LDout))
      }

      if(LP_out < max(LDout[['Monitor']][,1])){
        diff = abs(LP_out - max(LDout[['Monitor']][,1]))
        LP_out = max(LDout[['Monitor']][,1])
        parm_out = LDout$Posterior1[which.max(LDout[['Monitor']][,1]),]
        best = 'LD_Mode'
        iteration = i
        message('LD Mode ',i,': ',round(LP_out,3), ' ', paste(round(parm_out,3),collapse = ' '))
      }

      if(LP_out < LDout$Summary1['LP','Median']){
        diff = abs(LP_out - LDout$Summary1['LP','Median'])
        LP_out = LDout$Summary1['LP','Median']
        parm_out = LDout$Summary1[1:length(parm_out),'Median']
        best = 'LD_Median'
        iteration = i
        message('LD Median ',i,': ',round(LP_out,3), ' ', paste(round(parm_out,3),collapse = ' '))
      }

      if(LP_out < LDout$Summary1['LP','Mean']){
        diff = abs(LP_out - LDout$Summary1['LP','Mean'])
        LP_out = LDout$Summary1['LP','Mean']
        parm_out = LDout$Summary1[1:length(parm_out),'Mean']
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
    names(CMAout$par) = parm.names
  }

  time = (proc.time()[3]-timestart)/60

  LDout$Model = NULL #don't want this in case it is big!
  LDout$Call = NULL #don't want this in case it is big!

  return(invisible(list(parm=parm_out, LP=LP_out, diff=diff, best=best, iteration=iteration,
                        CMA_last=CMAout, LD_last=LDout, N = DataLD[['N']], RedChi2 = RedChi2, call=call, date=date,
                        time=time, CMA_all=CMA_all, LD_all=LD_all)))
}

.convert_CMA2CMA=function(parm, Data, likefunc, liketype='min'){
  #Convert CMA type output to LD
  output = likefunc(parm, Data)
  if(liketype=='min'){
    fnscale = 1
  }else if(liketype=='max'){
    fnscale = -1
  }
  return(fnscale*output)
}

.convert_CMA2LD=function(parm, Data, likefunc, liketype='min'){
  #Convert CMA type output to LD
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
  #Convert LD type output to CMA
  output = likefunc(parm, Data)
  if(liketype=='min'){
    fnscale = 1
  }else if(liketype=='max'){
    fnscale = -1
  }
  return(fnscale*output$LP)
}

.convert_LD2LD=function(parm, Data, likefunc, liketype='min'){
  #Convert LD type output to LD
  if(Data[['applyconstraints']] & !is.null(Data[['constraints']])){
    parm = Data[['constraints']](parm)
  }

  if(Data[['applyintervals']] & !is.null(Data[['intervals']]$lo) & !is.null(Data[['intervals']]$hi)){
    parm[parm<Data[['intervals']]$lo] = Data[['intervals']]$lo[parm<Data[['intervals']]$lo]
    parm[parm>Data[['intervals']]$hi] = Data[['intervals']]$hi[parm>Data[['intervals']]$hi]
  }
  if(length(Data[['mon.names']]) > 1){
    Data[['mon.names']] = Data[['mon.names']][2:length(Data[['mon.names']])]
  }else{
    Data[['mon.names']] = ''
  }
  output = likefunc(parm, Data)
  if(liketype=='min'){
    fnscale = -1
  }else if(liketype=='max'){
    fnscale = 1
  }
  if(!is.null(output[['parm']])){
    parm = output[['parm']]
  }
  return(list(LP=fnscale*output$LP, Dev=fnscale*output$Dev, Monitor=c(fnscale*output$LP,output$Monitor), yhat=output$yhat, parm=parm))
}
