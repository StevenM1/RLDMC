rm(list=ls())
library(snowfall)
source ("dmc/dmc.R")
source('utils.R')
source('models.R')
samplesDir <- 'samples'
savePlot <- FALSE

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

# Load BPICs & quantiles per bin ---------------------------------------------------------------
modelName <- 'arw-RL-mag'  #'ddm-RL-st0'
tmp <- getDataPpBPIC(modelName, 'exp1')
BPIC <- tmp$BPIC
qRTs <- getqRTs(tmp[['data3']], tmp[['pp3']])
data3 <- tmp[['data3']]
pp3 <- tmp[['pp3']]

q10RTsByEase <- list(getDescriptives(data3, dep.var='RT.10.', attr.name='qRTsCorrectByEase', id.var1='~bin*ease', id.var2=NULL),
                     getDescriptives(pp3, dep.var='RT.10.', attr.name='qRTsCorrectByEase', id.var1='~reps*bin*ease', id.var2=NULL))
q50RTsByEase <- list(getDescriptives(data3, dep.var='RT.50.', attr.name='qRTsCorrectByEase', id.var1='~bin*ease', id.var2=NULL),
                     getDescriptives(pp3, dep.var='RT.50.', attr.name='qRTsCorrectByEase', id.var1='~reps*bin*ease', id.var2=NULL))
q90RTsByEase <- list(getDescriptives(data3, dep.var='RT.90.', attr.name='qRTsCorrectByEase', id.var1='~bin*ease', id.var2=NULL),
                     getDescriptives(pp3, dep.var='RT.90.', attr.name='qRTsCorrectByEase', id.var1='~reps*bin*ease', id.var2=NULL))
q10RTsByEaseE <- list(getDescriptives(data3, dep.var='RT.10.', attr.name='qRTsErrorByEase', id.var1='~bin*ease', id.var2=NULL),
                      getDescriptives(pp3, dep.var='RT.10.', attr.name='qRTsErrorByEase', id.var1='~reps*bin*ease', id.var2=NULL))
q50RTsByEaseE <- list(getDescriptives(data3, dep.var='RT.50.', attr.name='qRTsErrorByEase', id.var1='~bin*ease', id.var2=NULL),
                      getDescriptives(pp3, dep.var='RT.50.', attr.name='qRTsErrorByEase', id.var1='~reps*bin*ease', id.var2=NULL))
q90RTsByEaseE <- list(getDescriptives(data3, dep.var='RT.90.', attr.name='qRTsErrorByEase', id.var1='~bin*ease', id.var2=NULL),
                      getDescriptives(pp3, dep.var='RT.90.', attr.name='qRTsErrorByEase', id.var1='~reps*bin*ease', id.var2=NULL))
meanAccByEase <- list(getDescriptives(data3, dep.var='acc', attr.name='AccByEase', id.var1='~bin*ease', id.var2=NULL),
                      getDescriptives(pp3, dep.var='acc', attr.name='AccByEase', id.var1='~reps*bin*ease', id.var2=NULL))


# Plot posterior predictives
if(savePlot) pdf(file=paste0('./figures/exp1_difficulty_', modelName, '-QQ-horizontal.pdf'), width=7, height=7/4*3)
par(oma=c(3,3,1,0), mar=c(0, 1, 1, 0) + 0.1, mfcol=c(3,4), mgp=c(2.75,.75,0), las=1, bty='l')
i <- 0
corrRTylim <- c(0.45, 1.1)
errRTylim <- c(0.45, 1.1)
data.cex=1.5
for(ease in unique(meanAccByEase[[1]]$ease)) {
  i <- i+1
  idxD = meanAccByEase[[1]]$ease == ease
  idxM = meanAccByEase[[2]]$ease == ease
  
  plotDataPPBins(data=meanAccByEase[[1]][idxD,], pp=meanAccByEase[[2]][idxM,],
                 xaxt='n', draw.legend = i==1, data.cex = data.cex,
                 dep.var='acc', ylab='', xlab = '', yaxt='n',
                 legend.pos='topleft', ylim=c(0.5, 0.95), hline.by=0.1)
  axis(1, at=seq(2, 10, 2), labels=rep(NA, 5), lwd=2)
  if(i == 1) {
    mtext('Accuracy', side=2, cex=.66, line=3, las=0, font=1)
    axis(2, at=seq(.5, .9, .1), lwd=1.5)
  } else {
    axis(2, at=seq(.5, .9, .1), labels=rep(NA, 5), lwd=1.5)
  }
  if(i == 1) title('0.6/0.4 (Hardest)')
  if(i == 2) title('0.65/0.35')
  if(i == 3) title('0.7/0.3')
  if(i == 4) title('0.8/0.2 (Easiest)')
  
  ##
  plotDataPPBins(data=q10RTsByEase[[1]][idxD,], pp=q10RTsByEase[[2]][idxM,], dep.var='RT.10.', 
                 ylim=corrRTylim, xaxt='n', ylab='', yaxt='n', draw.legend = FALSE, data.cex = data.cex)
  plotDataPPBins(data=q50RTsByEase[[1]][idxD,], pp=q50RTsByEase[[2]][idxM,], dep.var='RT.50.', plot.new = FALSE, draw.legend=FALSE, data.cex = data.cex)
  plotDataPPBins(data=q90RTsByEase[[1]][idxD,], pp=q90RTsByEase[[2]][idxM,], dep.var='RT.90.', plot.new = FALSE, draw.legend=FALSE, data.cex = data.cex)
  axis(1, at=seq(2, 10, 2), labels=NA, lwd=1.5)
  if(i == 1) {
    mtext('Correct RTs (s)', side=2, cex=.66, line=3, las=0, font=1)
    axis(2, seq(.4, 1.2, .2), lwd=1.5)
  } else {
    axis(2, seq(.4, 1.2, .2), labels=NA, lwd=1.5)
  }
  
  ##
  plotDataPPBins(data=q10RTsByEaseE[[1]][idxD,], pp=q10RTsByEaseE[[2]][idxM,], dep.var='RT.10.', ylim=errRTylim, xaxt='n', ylab='', yaxt='n', draw.legend = FALSE, data.cex = data.cex)
  plotDataPPBins(data=q50RTsByEaseE[[1]][idxD,], pp=q50RTsByEaseE[[2]][idxM,], dep.var='RT.50.', plot.new = FALSE, draw.legend=FALSE, data.cex = data.cex)
  plotDataPPBins(data=q90RTsByEaseE[[1]][idxD,], pp=q90RTsByEaseE[[2]][idxM,], dep.var='RT.90.', plot.new = FALSE, draw.legend=FALSE, data.cex = data.cex)
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



# Q-values, drift rates ---------------------------------------------------
get.color <- function(ease) {
  if(ease == "0.6") return(1)
  if(ease == "0.4") return(2)
  if(ease == "0.3") return(3)
  if(ease == "0.2") return(4)
}


draw.polygon <- function(pp, dep.var, xaxis='bin', colorM='blue', plot.model.points=TRUE) {
  lowerQ <- aggregate(as.formula(paste0(dep.var, '~bin')), pp, quantile, .025)
  upperQ <- aggregate(as.formula(paste0(dep.var, '~bin')), pp, quantile, .975)
  xs <- c(lowerQ[,xaxis], rev(lowerQ[,xaxis]))
  ys <- c(lowerQ[,dep.var], rev(upperQ[,dep.var]))
  polygon(xs, ys, col=rgb(col2rgb(colorM)[1]/255, col2rgb(colorM)[2]/255, col2rgb(colorM)[3]/255, alpha=.3), lty = NULL, border=NA)
  if(plot.model.points) points(pp[,xaxis], pp[,dep.var], pch=20, col=colorM, cex=.01)
}

allQ1OverTimeM <- do.call(rbind, (lapply(1:length(pp3), function(x) {tmp <- attr(pp3[[x]], 'SR.r1OverBins'); tmp$s <- x; tmp})))
allQ2OverTimeM <- do.call(rbind, (lapply(1:length(pp3), function(x) {tmp <- attr(pp3[[x]], 'SR.r2OverBins'); tmp$s <- x; tmp})))
meanQ1OverTimeM <- aggregate(SR.r1~reps*bin*ease, allQ1OverTimeM, mean)
meanQ2OverTimeM <- aggregate(SR.r2~reps*bin*ease, allQ2OverTimeM, mean)

# differences
deltaQ <- allQ1OverTimeM
deltaQ$SR.r2 <- allQ2OverTimeM$SR.r2
deltaQ$deltaQ <- deltaQ$SR.r1 - deltaQ$SR.r2

# sums
deltaQ$sumQ <- deltaQ$SR.r1 + deltaQ$SR.r2

meanDeltaQOverTime <- aggregate(deltaQ~reps*bin*ease, deltaQ, mean)
meanSumQOverTime <- aggregate(sumQ~reps*bin*ease, deltaQ, mean)

# drift rates
allV1OverTimeM <- do.call(rbind, (lapply(1:length(pp3), function(x) {tmp <- attr(pp3[[x]], 'mean_v.r1OverBins'); tmp$s <- x; tmp})))
allV2OverTimeM <- do.call(rbind, (lapply(1:length(pp3), function(x) {tmp <- attr(pp3[[x]], 'mean_v.r2OverBins'); tmp$s <- x; tmp})))
meanV1OverTimeM <- aggregate(mean_v.r1~reps*bin*ease, allV1OverTimeM, mean)
meanV2OverTimeM <- aggregate(mean_v.r2~reps*bin*ease, allV2OverTimeM, mean)


# Plot --------------------------------------------------------------------
pdf(file='./figures/q-values.pdf', width=7, height=2.5)
par(mfrow=c(1,4), las=1, bty='l', oma=c(0,1,1,0), mar=c(4, 3, 2, 0.5) + 0.1, mgp=c(2.25,.75,0))
# Q-values
plot(0,0, type='n', xlim=range(meanQ1OverTimeM$bin)+c(-.5, .5), ylim=c(0, .85), xlab='Trial bin', ylab='Q-values', main='Q-values')
abline(h=seq(0, 1, .1), col='grey')
abline(v=seq(0, 10, 2), col='grey')
for(ease in unique(meanQ1OverTimeM$ease)) {
  draw.polygon(meanQ1OverTimeM[meanQ1OverTimeM$ease==ease,], dep.var='SR.r1', colorM=get.color(ease)) #ifelse(ease=="0.6", 1, ifelse(ease=="0.4", 2, ifelse(ease=="0.3", 3, 4))))
  draw.polygon(meanQ2OverTimeM[meanQ2OverTimeM$ease==ease,], dep.var='SR.r2', colorM=get.color(ease))
}

# delta here
plot(0,0, type='n', xlim=range(meanDeltaQOverTime$bin)+c(-.5, .5), ylim=c(0, .85), xlab='Trial bin', ylab=expression(paste(Delta, 'Q-values')), main=expression(bold(paste(Delta, 'Q-values'))))
abline(h=seq(0, 1, .1), col='grey')
abline(v=seq(0, 10, 2), col='grey')
for(ease in unique(meanDeltaQOverTime$ease)) {
  draw.polygon(meanDeltaQOverTime[meanDeltaQOverTime$ease==ease,], dep.var='deltaQ', colorM=get.color(ease))
}

# sum here
plot(0,0, type='n', xlim=range(meanDeltaQOverTime$bin)+c(-.5, .5), ylim=c(0, .85), xlab='Trial bin', ylab=expression(paste(Sigma, 'Q-values')), main=expression(bold(paste(Sigma, 'Q-values'))))
abline(h=seq(0, 1, .1), col='grey')
abline(v=seq(0, 10, 2), col='grey')
for(ease in unique(meanDeltaQOverTime$ease)) {
  draw.polygon(meanSumQOverTime[meanSumQOverTime$ease==ease,], dep.var='sumQ', colorM=get.color(ease))
}
legend('bottomright', legend=c('0.8/0.2', '0.7/0.3', '0.65/0.35', '0.6/0.4'), col=1:4, pch=15, bty='n', title='Difficulty')

# drift rates
plot(0,0, type='n', xlim=range(meanV1OverTimeM$bin)+c(-.5, .5), ylim=c(1.0, 4.5), xlab='Trial bin', ylab='Drift rates', main='Drift rates')
abline(h=seq(0, 5, .25), col='grey')
abline(v=seq(0, 10, 2), col='grey')
for(ease in unique(allV1OverTimeM$ease)) {
  draw.polygon(meanV1OverTimeM[meanV1OverTimeM$ease==ease,], dep.var='mean_v.r1', colorM=get.color(ease))
  draw.polygon(meanV2OverTimeM[meanV2OverTimeM$ease==ease,], dep.var='mean_v.r2', colorM=get.color(ease))
}

dev.off()

