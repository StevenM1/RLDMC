library(snowfall)

getDataPpBPIC <- function(model, samples, dataName, do.plot=FALSE, BPIConly=FALSE) {
  dat <- loadData(dataName, removeBlock = NULL)[['dat']]
  # Load, generate posterior preds -------------------------------------------
  data <- lapply(samples, function(x) x$data)
  if(do.plot) plot.dmc(samples, hyper=TRUE, density=TRUE, layout=c(4,4))
  if(!BPIConly) {
    pp = h.post.predict.dmc(samples = samples, adapt=TRUE,save.simulation = TRUE, cores=30)
    ppNoSim <- h.pp.summary(pp, samples=samples)
    
    #### Append stimulus set info to data & model --------
    nBins <- 10
    pp2 <- lapply(1:length(pp), addStimSetInfo, input=pp, orig_dat=dat)
    data2 <- lapply(1:length(data), addStimSetInfo, input=data, orig_dat=dat)
    if(!sfIsRunning()) sfInit(parallel=TRUE, cpus =30); sfLibrary(moments)
    pp3 <- sfLapply(pp2, calculateByBin)
    data3 <- lapply(data2, calculateByBin)
    
    bpics <- h.IC.dmc(samples)
    return(list('pp3'=pp3, 'data3'=data3, 'BPIC'=bpics))
  } else {
    return(list('BPIC'=h.IC.dmc(samples)))
  }
}

getqRTs <- function(data3, pp3) {
  q10RTsOverTime <- list(getDescriptives(data3, dep.var='RT.10.', attr.name='qRTsCorrect', id.var1='~bin*s', id.var2="~bin"),
                         getDescriptives(pp3, dep.var='RT.10.', attr.name='qRTsCorrect'))
  q50RTsOverTime <- list(getDescriptives(data3, dep.var='RT.50.', attr.name='qRTsCorrect', id.var1='~bin*s', id.var2="~bin"),
                         getDescriptives(pp3, dep.var='RT.50.', attr.name='qRTsCorrect'))
  q90RTsOverTime <- list(getDescriptives(data3, dep.var='RT.90.', attr.name='qRTsCorrect', id.var1='~bin*s', id.var2="~bin"),
                         getDescriptives(pp3, dep.var='RT.90.', attr.name='qRTsCorrect'))
  q10RTsOverTimeE <- list(getDescriptives(data3, dep.var='RT.10.', attr.name='qRTsError', id.var1='~bin*s', id.var2="~bin"),
                          getDescriptives(pp3, dep.var='RT.10.', attr.name='qRTsError'))
  q50RTsOverTimeE <- list(getDescriptives(data3, dep.var='RT.50.', attr.name='qRTsError', id.var1='~bin*s', id.var2="~bin"),
                          getDescriptives(pp3, dep.var='RT.50.', attr.name='qRTsError'))
  q90RTsOverTimeE <- list(getDescriptives(data3, dep.var='RT.90.', attr.name='qRTsError', id.var1='~bin*s', id.var2="~bin"),
                          getDescriptives(pp3, dep.var='RT.90.', attr.name='qRTsError'))
  
  meanAccOverTime <- list(getDescriptives(data3, dep.var='acc', attr.name='AccOverBins', id.var1='~bin*s', id.var2="~bin"),
                          getDescriptives(pp3, dep.var='acc', attr.name='AccOverBins'))
  return(list('q10RTsOverTime'=q10RTsOverTime,
              'q50RTsOverTime'=q50RTsOverTime,
              'q90RTsOverTime'=q90RTsOverTime,
              'q10RTsOverTimeE'=q10RTsOverTimeE,
              'q50RTsOverTimeE'=q50RTsOverTimeE,
              'q90RTsOverTimeE'=q90RTsOverTimeE,
              'meanAccOverTime'=meanAccOverTime))
}