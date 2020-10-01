## load data
rm(list=ls())
library(snowfall)
source ("dmc/dmc.R")
source('utils.R')
source('models.R')
samplesDir <- 'samples'
savePlot <- FALSE

calculateByBin <- function(df) {
  df$acc <- as.integer(df$R)==2
  
  attr(df, 'qRTs') <- do.call(data.frame, aggregate(RT~reps*bin, df, quantile, probs=seq(.1, .9, .4)))
  attr(df, 'qRTsCorrect') <- do.call(data.frame, aggregate(RT~reps*bin, df[df$acc==1,], quantile, probs=seq(.1, .9, .4)))
  attr(df, 'qRTsError') <- do.call(data.frame, aggregate(RT~reps*bin, df[df$acc==0,], quantile, probs=seq(.1, .9, .4)))
  
  attr(df, 'RTsOverBins') <- aggregate(RT~reps*bin, df, mean)
  attr(df, 'AccOverBins') <- aggregate(acc~reps*bin, df, mean)
  
  if(!is.null(attr(df, 'adapt'))) {
    adapted <- attr(df, 'adapt')
    adapted$ease <- df$ease
    adapted$bin <- df$bin
    for(colName in colnames(adapted)) {
      form <- as.formula(paste0(colName, '~reps*ease*bin'))
      attr(df, paste0(colName, 'OverBins')) <- aggregate(form, adapted, mean)
    }
  }
  df
}

getDataPpBPIC <- function(modelName, dataName, do.plot=FALSE, BPIConly=FALSE) {
  modelName <- 'arw-RL-risk-mag-Niek'
  dataName <- 'exp2'
  model <- setupModel(modelName)  # calls load_model(), which loads transform.dmc() and transform2.dmc()
  dat <- loadData(dataName, removeBlock = NULL)[['dat']]
  dat$bin <- dat$trialNreversal   ## important for the reversals!
  #fn <- paste0('model-', modelName, '_data-', dataName)
  #fn <- paste0('model-', modelName, ".R", '_data-', dataName)
  # Load, generate posterior preds -------------------------------------------
  samples <- loadSamples("parameterRecoveries/samples/model-arw-RL-risk-mag-Niek_data-parameterRecovery-exp2-Risk")
  data <- lapply(samples, function(x) x$data)
  if(do.plot) plot.dmc(samples, hyper=TRUE, density=TRUE, layout=c(4,4))
  if(!BPIConly) {
    pp = h.post.predict.dmc(samples = samples, adapt=TRUE,save.simulation = TRUE, cores=30)
    ppNoSim <- h.pp.summary(pp, samples=samples)
    
    #### Append stimulus set info to data & model --------
    nBins <- 10
    pp2 <- lapply(1:length(pp), addStimSetInfo, input=pp, orig_dat=dat, addColumns='bin')
    data2 <- lapply(1:length(data), addStimSetInfo, input=data, orig_dat=dat, addColumns='bin')
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
  q10RTsOverTime <- list(getDescriptives(data3, dep.var='RT.10.', attr.name='qRTs', id.var1='~bin*s', id.var2="~bin"),
                         getDescriptives(pp3, dep.var='RT.10.', attr.name='qRTs'))
  q50RTsOverTime <- list(getDescriptives(data3, dep.var='RT.50.', attr.name='qRTs', id.var1='~bin*s', id.var2="~bin"),
                         getDescriptives(pp3, dep.var='RT.50.', attr.name='qRTs'))
  q90RTsOverTime <- list(getDescriptives(data3, dep.var='RT.90.', attr.name='qRTs', id.var1='~bin*s', id.var2="~bin"),
                         getDescriptives(pp3, dep.var='RT.90.', attr.name='qRTs'))

  meanAccOverTime <- list(getDescriptives(data3, dep.var='acc', attr.name='AccOverBins', id.var1='~bin*s', id.var2="~bin"),
                          getDescriptives(pp3, dep.var='acc', attr.name='AccOverBins'))
  return(list('q10RTsOverTime'=q10RTsOverTime,
              'q50RTsOverTime'=q50RTsOverTime,
              'q90RTsOverTime'=q90RTsOverTime,
              'meanAccOverTime'=meanAccOverTime))
}


# Load --------------------------------------------------------------------
# DDM
tmp <- getDataPpBPIC('arw-RL-risk-mag-Niek', 'exp2')
BPICDDM <- tmp$BPIC
qRTsDDM <- getqRTs(tmp[['data3']], tmp[['pp3']])

allqRTs <- list(qRTsDDM)



# Plot --------------------------------------------------------------------
layoutM <- matrix(1:15, nrow=3, byrow=TRUE)
layoutM[c(1, 3),1:2] <- 1
layoutM[c(1, 3),4:5] <- 5
layoutM[2,1:2] <- 2:3 #matrix(c(2:13), nrow=3, byrow=TRUE)
layoutM[2,4:5] <- 6:7 #matrix(c(2:13), nrow=3, byrow=TRUE)
layoutM[,3] <- 4
layoutM

if(savePlot) pdf('./figures/exp2-reversals-ddm-st0.pdf', width=4, height=2.5)
layout(layoutM, heights = c(0.01, 1, 0.01), widths=c(1,1,0,0,0)) #.1,1,1))
par(oma=c(3,2,2,0), mar=c(0, 2, 1, 0.5) + 0.1, #mfcol=c(3,4), 
    mgp=c(2.75,.75,0), las=1, bty='l')
i <- 0
data.cex=1.5
for(qRTs in allqRTs[1]) {
  par(mar=c(0,0,0,0))
  plot.new()
  if(i == 0) mtext('RL-DDM A3', side=3, cex=.66*1.2, font=2, line=1)
  if(i == 1) {plot.new(); mtext('RL-fARD', side=3, cex=.66*1.2, font=2, line=1)}
  i <- i+1
  par(mar=c(0,2,1,.5)+.1)
  plotDataPPBins(data=qRTs$meanAccOverTime[[1]], pp=qRTs$meanAccOverTime[[2]],
                 xaxt='n', xlim=c(-60, 40), plot.model.points=FALSE,
                 dep.var='acc', ylab=expression('Proportion choice A'),
                 xlab = '', data.lwd=1.5, data.cex=data.cex,
                 draw.legend=i==1,
                 legend.pos='bottomleft', ylim=c(0.25, 0.85), hline.by=0.05, axvlines=seq(-50, 50, 10))
  
  axis(1, at=seq(-50, 50, 10), lwd=2)
  abline(v=0, lty=2, lwd=2)
  if(i == 1) {
    mtext('Proportion choice A', side=2, cex=.66, line=3, las=0, font=1)
    axis(2, at=seq(.1, .9, .1), lwd=1.5)
  } else {
    axis(2, at=seq(.1, .9, .1), labels=rep(NA, 5), lwd=1.5)
  }
  mtext('Trial (relative to reversal)', side=1, cex=.66, line=2)
  title('Choices')
  
  ##
  par(mar=c(0,3,1,.5)+.1)
  plotDataPPBins(data=qRTs$q10RTsOverTime[[1]], pp=qRTs$q10RTsOverTime[[2]],
                 dep.var='RT.10.', ylim=c(.4, 1.1), ylab='',
                 xaxt='n', xlim=c(-60, 40), plot.model.points = FALSE, data.lwd=1.5, data.cex=data.cex,
                 xlab='', hline.by=0.05, axvlines=seq(-50, 50, 10), legend.pos = FALSE)
  plotDataPPBins(data=qRTs$q50RTsOverTime[[1]], pp=qRTs$q50RTsOverTime[[2]], dep.var='RT.50.', plot.new = FALSE, draw.legend=FALSE,
                 plot.model.points=FALSE,data.cex = data.cex)
  plotDataPPBins(data=qRTs$q90RTsOverTime[[1]], pp=qRTs$q90RTsOverTime[[2]], dep.var='RT.90.', plot.new = FALSE, draw.legend=FALSE, 
                 plot.model.points=FALSE, data.cex = data.cex)
  axis(1, at=seq(-50, 50, 10), lwd=2)
  abline(v=0, lty=2, lwd=2)
  if(i == 1) {
    mtext('RTs (s)', side=2, cex=.66, line=2, las=0, font=1)
    axis(2, seq(.4, 1.2, .2), lwd=1.5)
  } else {
    axis(2, seq(.4, 1.2, .2), lwd=1.5)
  }
  mtext('Trial (relative to reversal)', side=1, cex=.66, line=2)
  title('RTs')
}
if(savePlot) dev.off()

