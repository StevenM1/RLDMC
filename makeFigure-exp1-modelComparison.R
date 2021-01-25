rm(list=ls())
library(snowfall)
source ("dmc/dmc.R")
source('utils.R')
source('models.R')
samplesDir <- 'samples'
savePlot <- TRUE

getDataPpBPIC <- function(modelName, dataName, do.plot=FALSE, BPIConly=FALSE) {
  model <- setupModel(modelName)  # calls load_model(), which loads transform.dmc() and transform2.dmc()
  dat <- loadData(dataName, removeBlock = NULL)[['dat']]
  fn <- paste0('model-', modelName, '_data-', dataName)
  
  # Load, generate posterior preds -------------------------------------------
  samples <- loadSamples(fn, samplesDir)
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


# Load BPICs & quantiles per bin ---------------------------------------------------------------
# DDM
tmp <- getDataPpBPIC('ddm-RL', 'exp1')
BPICDDM <- tmp$BPIC
qRTsDDM <- getqRTs(tmp[['data3']], tmp[['pp3']])

# RD
tmp <- getDataPpBPIC('rw-RL', 'exp1')
BPICRD <- tmp$BPIC
qRTsRD <- getqRTs(tmp[['data3']], tmp[['pp3']])

# ARDl
tmp <- getDataPpBPIC('arw-RL', 'exp1')
BPICARDl <- tmp$BPIC
qRTsARDl <- getqRTs(tmp[['data3']], tmp[['pp3']])

# ARDf
tmp <- getDataPpBPIC('arw-RL-mag', 'exp1')
BPICARD <- tmp$BPIC
qRTsARD <- getqRTs(tmp[['data3']], tmp[['pp3']])

# Combine accuracy per bin, quantile per bin
allqRTs <- list(qRTsDDM, qRTsRD, qRTsARDl, qRTsARD)

# Main text model comparison for experiment 1
allBPICs <- cbind(BPICDDM[,2], BPICRD[,2], BPICARDl[,2], BPICARD[,2])
apply(allBPICs, 2, sum)

# Plot posterior predictives
if(savePlot) pdf('./figures/modelcomparison-exp1.pdf', width=7, height=7/4*3)
par(oma=c(3,4,2,0), mar=c(0, 0, 1, 0.5) + 0.1, mfcol=c(3,4), mgp=c(2.75,.75,0), las=1, bty='l')
i <- 0
data.cex=1.5
corrRTylim <- errRTylim <- c(.45, 1.1)
for(qRTs in allqRTs) {
  i <- i+1
  
  plotDataPPBins(data=qRTs$meanAccOverTime[[1]], pp=qRTs$meanAccOverTime[[2]],
                 xaxt='n', draw.legend = i==1, data.cex = data.cex,
                 dep.var='acc', ylab='', xlab = '', yaxt='n',
                 legend.pos='topleft', ylim=c(0.5, 0.85), hline.by=0.1)
  axis(1, at=seq(2, 10, 2), labels=rep(NA, 5), lwd=2)
  if(i == 1) {
    mtext('Accuracy', side=2, cex=.66, line=3, las=0, font=1)
    axis(2, at=seq(.5, .9, .1), lwd=1.5)
  } else {
    axis(2, at=seq(.5, .9, .1), labels=rep(NA, 5), lwd=1.5)
  }
  par(xpd=NA)
  if(i == 1) title('RL-DDM', line=1.2)
  if(i == 2) title('RL-RD', line=1.2)
  if(i == 3) title('RL-lARD', line=1.2)
  if(i == 4) title('RL-ARD', line=1.2)
  mtext(paste0('BPIC = ', round(apply(allBPICs, 2, sum)[i])), cex=.66) 
  par(xpd=FALSE)
  
  ##
  plotDataPPBins(data=qRTs$q10RTsOverTime[[1]], pp=qRTs$q10RTsOverTime[[2]], dep.var='RT.10.', 
                 ylim=corrRTylim, xaxt='n', ylab='', yaxt='n', draw.legend = FALSE, data.cex = data.cex)
  plotDataPPBins(data=qRTs$q50RTsOverTime[[1]], pp=qRTs$q50RTsOverTime[[2]], dep.var='RT.50.', plot.new = FALSE, draw.legend=FALSE, data.cex = data.cex)
  plotDataPPBins(data=qRTs$q90RTsOverTime[[1]], pp=qRTs$q90RTsOverTime[[2]], dep.var='RT.90.', plot.new = FALSE, draw.legend=FALSE, data.cex = data.cex)
  axis(1, at=seq(2, 10, 2), labels=NA, lwd=1.5)
  if(i == 1) {
    mtext('Correct RTs (s)', side=2, cex=.66, line=3, las=0, font=1)
    axis(2, seq(.4, 1.2, .2), lwd=1.5)
  } else {
    axis(2, seq(.4, 1.2, .2), labels=NA, lwd=1.5)
  }
  
  ##
  plotDataPPBins(data=qRTs$q10RTsOverTimeE[[1]], pp=qRTs$q10RTsOverTimeE[[2]], dep.var='RT.10.', ylim=errRTylim, xaxt='n', ylab='', yaxt='n', draw.legend = FALSE, data.cex = data.cex)
  plotDataPPBins(data=qRTs$q50RTsOverTimeE[[1]], pp=qRTs$q50RTsOverTimeE[[2]], dep.var='RT.50.', plot.new = FALSE, draw.legend=FALSE, data.cex = data.cex)
  plotDataPPBins(data=qRTs$q90RTsOverTimeE[[1]], pp=qRTs$q90RTsOverTimeE[[2]], dep.var='RT.90.', plot.new = FALSE, draw.legend=FALSE, data.cex = data.cex)
  if(i == 1) {
    mtext('Error RTs (s)', side=2, cex=.66, line=3, las=0, font=1)
    axis(2, seq(.4, 1.2, .2), lwd=1.5)
  } else {
    axis(2, seq(.4, 1.2, .2), labels=NA, lwd=1.5)
  }
  axis(1, at=seq(2, 10, 2), lwd=1.5)
  mtext('Trial bin', side=1, cex=.66, line=2)
}
if(savePlot) dev.off()
