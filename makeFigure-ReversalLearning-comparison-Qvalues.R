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
  model <- setupModel(modelName)  # calls load_model(), which loads transform.dmc() and transform2.dmc()
  dat <- loadData(dataName, removeBlock = NULL)[['dat']]
  dat$bin <- dat$trialNreversal   ## important for the reversals!
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
tmp <- getDataPpBPIC('ddm-RL', 'exp2')
BPICDDM <- tmp$BPIC
qRTsDDM <- getqRTs(tmp[['data3']], tmp[['pp3']])

# ARD
tmp <- getDataPpBPIC('arw-RL-mag', 'exp2')
BPICARD <- tmp$BPIC
qRTsARD <- getqRTs(tmp[['data3']], tmp[['pp3']])

# Combine
allqRTs <- list(qRTsDDM, qRTsARD)
allBPICs <- cbind(BPICDDM[,2],  BPICARD[,2])
apply(allBPICs, 2, sum)


# # Plot --------------------------------------------------------------------
pdf('./figures/exp2-reversals-2x2.pdf', width=5, height=4)
par(oma=c(3,4,2,0), mar=c(1, 2, 1, 0.5) + 0.1, mfcol=c(2,2),
    mgp=c(2.75,.75,0), las=1, bty='l', cex=.66)
i <- 0
data.cex=1.
lwd.axis=1.5
for(qRTs in allqRTs) {
  i <- i+1
  plotDataPPBins(data=qRTs$meanAccOverTime[[1]], pp=qRTs$meanAccOverTime[[2]],
                 xaxt='n', yaxt='n', xlim=c(-60, 40), plot.model.points=FALSE,
                 dep.var='acc', ylab=expression('Proportion choice A'),
                 xlab = '', data.lwd=1.5, data.cex=data.cex,
                 draw.legend=i==1,
                 legend.pos='bottomleft', ylim=c(0.25, 0.85), hline.by=0.05, axvlines=seq(-50, 50, 10))
  abline(v=0, lty=2, lwd=2)
  abline(h=0.5, lty=2, lwd=2)
  axis(1, at=seq(-60, 50, 10), labels=NA, lwd=lwd.axis)
  if(i == 1) {
    mtext('Proportion choice A', side=2, cex=.66, line=3, las=0, font=1)
    mtext('RL-DDM', side=3, cex=.66*1.2, font=2, line=1)
    axis(2, at=seq(.3, .9, .1), lwd=lwd.axis)
  } else {
    mtext('RL-ARD', side=3, cex=.66*1.2, font=2, line=1)
    axis(2, at=seq(.3, .9, .1), labels=rep(NA, 5), lwd=lwd.axis)
  }

  ##
  plotDataPPBins(data=qRTs$q10RTsOverTime[[1]], pp=qRTs$q10RTsOverTime[[2]],
                 dep.var='RT.10.', ylim=c(.4, 1.1), ylab='RT (s)',
                 xaxt='n',  yaxt='n', xlim=c(-60, 40), plot.model.points = FALSE, data.lwd=1.5, data.cex=data.cex,
                 xlab='', hline.by=0.05, axvlines=seq(-50, 50, 10), legend.pos = FALSE)
  plotDataPPBins(data=qRTs$q50RTsOverTime[[1]], pp=qRTs$q50RTsOverTime[[2]], dep.var='RT.50.', plot.new = FALSE, draw.legend=FALSE,
                 plot.model.points=FALSE,data.cex = data.cex)
  plotDataPPBins(data=qRTs$q90RTsOverTime[[1]], pp=qRTs$q90RTsOverTime[[2]], dep.var='RT.90.', plot.new = FALSE, draw.legend=FALSE,
                 plot.model.points=FALSE, data.cex = data.cex)
  axis(1, at=seq(-60, 50, 10), lwd=lwd.axis)
  abline(v=0, lty=2, lwd=2)
  if(i == 1) {
    mtext('RTs (s)', side=2, cex=.66, line=3, las=0, font=1)
    axis(2, seq(.4, 1.2, .2), lwd=lwd.axis)
  } else {
    axis(2, seq(.4, 1.2, .2), labels=NA, lwd=lwd.axis)
  }
  mtext('Trial (relative to reversal)', side=1, cex=.66, line=2)
}
dev.off()


# Q-values over time ------------------------------------------------------
pp3 <- tmp[['pp3']]

# Q-values, drift rates ---------------------------------------------------
allQ1OverTimeM <- do.call(rbind, (lapply(1:length(pp3), function(x) {tmp <- attr(pp3[[x]], 'SR.r1OverBins'); tmp$s <- x; tmp})))
allQ2OverTimeM <- do.call(rbind, (lapply(1:length(pp3), function(x) {tmp <- attr(pp3[[x]], 'SR.r2OverBins'); tmp$s <- x; tmp})))
meanQ1OverTimeM <- aggregate(SR.r1~bin*ease, allQ1OverTimeM, mean)
meanQ2OverTimeM <- aggregate(SR.r2~bin*ease, allQ2OverTimeM, mean)

# differences
deltaQ <- allQ1OverTimeM
deltaQ$SR.r2 <- allQ2OverTimeM$SR.r2
deltaQ$deltaQ <- deltaQ$SR.r1 - deltaQ$SR.r2

# sums
deltaQ$sumQ <- deltaQ$SR.r1 + deltaQ$SR.r2

meanDeltaQOverTime <- aggregate(deltaQ~bin*ease, deltaQ, mean)
meanSumQOverTime <- aggregate(sumQ~bin*ease, deltaQ, mean)

# drift rates
allV1OverTimeM <- do.call(rbind, (lapply(1:length(pp3), function(x) {tmp <- attr(pp3[[x]], 'mean_v.r1OverBins'); tmp$s <- x; tmp})))
allV2OverTimeM <- do.call(rbind, (lapply(1:length(pp3), function(x) {tmp <- attr(pp3[[x]], 'mean_v.r2OverBins'); tmp$s <- x; tmp})))
meanV1OverTimeM <- aggregate(mean_v.r1~bin*ease, allV1OverTimeM, mean)
meanV2OverTimeM <- aggregate(mean_v.r2~bin*ease, allV2OverTimeM, mean)

library(RColorBrewer)
dark2 <- brewer.pal(n=8, name='Dark2')
set2 <- brewer.pal(n=8, name='Set2')
palette(c(set2[1], dark2[1], set2[2], dark2[2]))

get.color <- function(ease, choice='A') {
  if(ease == "0.6") i <- 1 #return(1)
  if(ease == "0.4") i <- 3 #return(2)
  return(i+(choice=='B'))
}


# plot
pdf(file='./figures/exp3-q-values.pdf', width=7, height=2.5)
par(mfrow=c(1,4), las=1, bty='l', oma=c(0,1,1,0), mar=c(4, 3, 4, 0.5) + 0.1, mgp=c(2.25,.75,0))
# Q-values
plot(0,0, type='n', ylim=c(0, .85), xlim=c(-60, 40), xlab='', xaxt='n', ylab='Q-values', main='')
title(expression(bold('A. Q-values')), line=2.5)
mtext('Trial (relative to reversal)', side=1, cex=.66, line=2)
abline(h=seq(0, 1, .1), col='grey')
abline(v=seq(-60, 50, 10), col='grey')
abline(v=0, lty=2, col='black', lwd=2)
axis(1, at=seq(-60, 50, 10), lwd=lwd.axis)

for(ease in unique(meanQ1OverTimeM$ease)) {
  lines(meanQ1OverTimeM[meanQ1OverTimeM$ease==ease,'bin'], meanQ1OverTimeM[meanQ1OverTimeM$ease==ease,'SR.r1'], col=get.color(ease, 'B'), lwd=1.5, lty=1)
  lines(meanQ2OverTimeM[meanQ2OverTimeM$ease==ease,'bin'], meanQ2OverTimeM[meanQ2OverTimeM$ease==ease,'SR.r2'], col=get.color(ease, 'B'), lwd=1.5, lty=2)
}

par(xpd=TRUE)
legend('top', legend=c('0.8/0.2 A',
                       '0.8/0.2 B',
                       '0.7/0.3 A',
                       '0.7/0.3 B'), col=c(2,2,4,4), lty=c(1,2), lwd=1.5, ncol = 2, bg='white', inset=c(0, -.25), bty='n', cex=.8)
par(xpd=FALSE)

# delta here
plot(0,0, type='n', xlim=c(-60, 40), ylim=c(-.6, .7), xlab='',  xaxt='n',
     ylab=expression(paste(Delta, 'Q-values')), main='')# expression(bold(paste(Delta, 'Q-values'))))
title(expression(bold(paste('B. ', Delta, 'Q-values'))), line=2.5)
mtext('Trial (relative to reversal)', side=1, cex=.66, line=2)
abline(h=seq(-1, 1, .1), col='grey')
abline(v=seq(-60, 50, 10), col='grey')
abline(v=0, lty=2, col='black', lwd=2)
abline(h=0, lty=2, col='black', lwd=2)
axis(1, at=seq(-60, 50, 10), lwd=lwd.axis)
for(ease in unique(meanDeltaQOverTime$ease)) {
  lines(meanDeltaQOverTime[meanDeltaQOverTime$ease==ease,'bin'], meanDeltaQOverTime[meanDeltaQOverTime$ease==ease,'deltaQ'], col=get.color(ease, 'B'), lwd=1.5)
}
par(xpd=TRUE)
legend('top', legend=c('0.8/0.2', '0.7/0.3'), col=c(2,4), lty=1, lwd=1.5, ncol = 2, bg='white', inset=c(0, -.175), bty='n', cex=.8)
par(xpd=FALSE)

# sum here
plot(0,0, type='n', xlim=c(-60, 40), ylim=c(0.1, 1), xlab='',  xaxt='n',ylab=expression(paste(Sigma, 'Q-values')), main='')
title(expression(bold(paste('C. ', Sigma, 'Q-values'))), line=2.5)
mtext('Trial (relative to reversal)', side=1, cex=.66, line=2)
abline(h=seq(0, 1, .1), col='grey')
abline(v=seq(-60, 50, 10), col='grey')
abline(v=0, lty=2, col='black', lwd=2)
axis(1, at=seq(-60, 50, 10), lwd=lwd.axis)
for(ease in unique(meanSumQOverTime$ease)) {
  lines(meanSumQOverTime[meanSumQOverTime$ease==ease,'bin'], meanSumQOverTime[meanSumQOverTime$ease==ease,'sumQ'], col=get.color(ease, 'B'), lwd=1.5)
}
par(xpd=TRUE)
legend('top', legend=c('0.8/0.2', '0.7/0.3'), col=c(2,4), lty=1, lwd=1.5, ncol = 2, bg='white', inset=c(0, -.175), bty='n', cex=.8)
par(xpd=FALSE)

# drift rates
plot(0,0, type='n', xlim=c(-60, 40), ylim=c(0.8, 3.5), xlab='', xaxt='n', ylab='Drift rates', main='')
mtext('Trial (relative to reversal)', side=1, cex=.66, line=2)
title(expression(bold('D. Drift rates')), line=2.5)
abline(h=seq(0, 10, .5), col='grey')
abline(v=seq(-60, 50, 10), col='grey')
abline(v=0, lty=2, col='black', lwd=2)
axis(1, at=seq(-60, 50, 10), lwd=lwd.axis)
for(ease in unique(allV1OverTimeM$ease)) {
  lines(meanV1OverTimeM[meanV1OverTimeM$ease==ease,'bin'], meanV1OverTimeM[meanV1OverTimeM$ease==ease,'mean_v.r1'], col=get.color(ease, 'B'), lwd=1.5, lty=2)
  lines(meanV2OverTimeM[meanV2OverTimeM$ease==ease,'bin'], meanV2OverTimeM[meanV2OverTimeM$ease==ease,'mean_v.r2'], col=get.color(ease, 'B'), lwd=1.5, lty=1)
}
par(xpd=TRUE)
legend('top', legend=c('0.8/0.2 A',
                       '0.8/0.2 B',
                       '0.7/0.3 A',
                       '0.7/0.3 B'), col=c(2,2,4,4), lty=c(1,2,1,2),lwd=1.5, ncol = 2, bg='white', inset=c(0, -.25), bty='n', cex=.8)
par(xpd=FALSE)
dev.off()




# Old stuff ---------------------------------------------------------------
# layoutM <- matrix(1:15, nrow=3, byrow=TRUE)
# layoutM[c(1, 3),1:2] <- 1
# layoutM[c(1, 3),4:5] <- 5
# layoutM[2,1:2] <- 2:3 #matrix(c(2:13), nrow=3, byrow=TRUE)
# layoutM[2,4:5] <- 6:7 #matrix(c(2:13), nrow=3, byrow=TRUE)
# layoutM[,3] <- 4
# layoutM
# 
# if(savePlot) pdf('./figures/exp2-reversals.pdf', width=7, height=2.25)
# layout(layoutM, heights = c(0.01, 1, 0.01), widths=c(1,1,.1,1,1))
# par(oma=c(3,4,2,0), mar=c(0, 2, 1, 0.5) + 0.1, #mfcol=c(3,4), 
#     mgp=c(2.75,.75,0), las=1, bty='l')
# i <- 0
# data.cex=1.5
# for(qRTs in allqRTs) {
#   par(mar=c(0,0,0,0))
#   plot.new()
#   if(i == 0) mtext('RL-DDM', side=3, cex=.66*1.2, font=2, line=1)
#   if(i == 1) {plot.new(); mtext('RL-ARD', side=3, cex=.66*1.2, font=2, line=1)}
#   i <- i+1
#   par(mar=c(0,2,1,.5)+.1)
#   plotDataPPBins(data=qRTs$meanAccOverTime[[1]], pp=qRTs$meanAccOverTime[[2]],
#                  xaxt='s', xlim=c(-60, 40), plot.model.points=FALSE,
#                  dep.var='acc', ylab=expression('Proportion choice A'),
#                  xlab = '', data.lwd=1.5, data.cex=data.cex,
#                  draw.legend=i==1,
#                  legend.pos='bottomleft', ylim=c(0.25, 0.85), hline.by=0.05, axvlines=seq(-50, 50, 10))
#   
#   axis(1, at=seq(-50, 50, 10), labels=rep(NA, 5), lwd=2)
#   abline(v=0, lty=2, lwd=2)
#   if(i == 1) {
#     mtext('Proportion choice A', side=2, cex=.66, line=3, las=0, font=1)
#     axis(2, at=seq(.5, .9, .1), lwd=1.5)
#   } else {
#     axis(2, at=seq(.5, .9, .1), labels=rep(NA, 5), lwd=1.5)
#   }
#   mtext('Trial (relative to reversal)', side=1, cex=.66, line=2)
#   title('Choices')
#   
#   ##
#   plotDataPPBins(data=qRTs$q10RTsOverTime[[1]], pp=qRTs$q10RTsOverTime[[2]],
#                  dep.var='RT.10.', ylim=c(.4, 1.1), ylab='RT (s)',
#                  xaxt='s', xlim=c(-60, 40), plot.model.points = FALSE, data.lwd=1.5, data.cex=data.cex,
#                  xlab='', hline.by=0.05, axvlines=seq(-50, 50, 10), legend.pos = FALSE)
#   plotDataPPBins(data=qRTs$q50RTsOverTime[[1]], pp=qRTs$q50RTsOverTime[[2]], dep.var='RT.50.', plot.new = FALSE, draw.legend=FALSE,
#                  plot.model.points=FALSE,data.cex = data.cex)
#   plotDataPPBins(data=qRTs$q90RTsOverTime[[1]], pp=qRTs$q90RTsOverTime[[2]], dep.var='RT.90.', plot.new = FALSE, draw.legend=FALSE, 
#                  plot.model.points=FALSE, data.cex = data.cex)
#   axis(1, at=seq(-50, 50, 10), labels=rep(NA, 5), lwd=2)
#   abline(v=0, lty=2, lwd=2)
#   if(i == 1) {
#     mtext('RTs (s)', side=2, cex=.66, line=3, las=0, font=1)
#     axis(2, seq(.4, 1.2, .2), lwd=1.5)
#   } else {
#     axis(2, seq(.4, 1.2, .2), labels=NA, lwd=1.5)
#   }
#   mtext('Trial (relative to reversal)', side=1, cex=.66, line=2)
#   title('RTs')
# }
# if(savePlot) dev.off()

# 
# # Plot --------------------------------------------------------------------
# if(savePlot) pdf(file='./figures/exp3-q-values2.pdf', width=7, height=2.5)
# lwd.axis=1.5
# par(mfrow=c(1,4), las=1, bty='l', oma=c(0,1,1,0), mar=c(4, 3, 4, 0.5) + 0.1, mgp=c(2.25,.75,0))
# # Q-values
# plot(0,0, type='n', ylim=c(0, .85), xlim=c(-60, 40), xlab='', xaxt='n', ylab='Q-values', main='')
# title(expression(bold('Q-values')), line=2.5)
# mtext('Trial (relative to reversal)', side=1, cex=.66, line=2)
# abline(h=seq(0, 1, .1), col='grey')
# abline(v=seq(-60, 50, 10), col='grey')
# abline(v=0, lty=2, col='black', lwd=2)
# axis(1, at=seq(-60, 50, 10), lwd=lwd.axis)
# 
# #if(byEase) {}
# # for(ease in unique(meanQ1OverTimeM$ease)) {
# # draw.polygon(meanQ1OverTimeM[meanQ1OverTimeM$ease==ease,], dep.var='SR.r1', colorM=get.color(ease, 'B'), plot.model.points = FALSE) #ifelse(ease=="0.6", 1, ifelse(ease=="0.4", 2, ifelse(ease=="0.3", 3, 4))))
# # draw.polygon(meanQ2OverTimeM[meanQ2OverTimeM$ease==ease,], dep.var='SR.r2', colorM=get.color(ease, 'A'), plot.model.points = FALSE)
# # }
# draw.polygon(meanQ1OverTimeM, dep.var='SR.r1', colorM=get.color(0.4, 'B'), plot.model.points = FALSE) #ifelse(ease=="0.6", 1, ifelse(ease=="0.4", 2, ifelse(ease=="0.3", 3, 4))))
# draw.polygon(meanQ2OverTimeM, dep.var='SR.r2', colorM=get.color(0.4, 'A'), plot.model.points = FALSE)
# # par(xpd=TRUE)
# # legend('top', legend=c('0.8/0.2 A',
# #                        '0.8/0.2 B',
# #                        '0.7/0.3 A',
# #                        '0.7/0.3 B'), col=c(2,1,4,3), pch=15, ncol = 2, bg='white', inset=c(0, -.25), bty='n')
# # par(xpd=FALSE)
# par(xpd=TRUE)
# legend('top', legend=c('A', 'B'), col=c(2,4), pch=15, ncol = 2, bg='white', inset=c(0, -.175), bty='n')
# par(xpd=FALSE)
# 
# #legend('topleft', legend=c('0.8/0.2', '0.7/0.3'), col=1:2, pch=15, bty='n', title='Difficulty')
# 
# # delta here
# plot(0,0, type='n', xlim=c(-60, 40), ylim=c(-.6, .7), xlab='',  xaxt='n',
#      ylab=expression(paste(Delta, 'Q-values')), main='')# expression(bold(paste(Delta, 'Q-values'))))
# title(expression(bold(paste(Delta, 'Q-values'))), line=2.5)
# mtext('Trial (relative to reversal)', side=1, cex=.66, line=2)
# abline(h=seq(-1, 1, .1), col='grey')
# abline(v=seq(-60, 50, 10), col='grey')
# abline(v=0, lty=2, col='black', lwd=2)
# abline(h=0, lty=2, col='black', lwd=2)
# axis(1, at=seq(-60, 50, 10), lwd=lwd.axis)
# draw.polygon(meanDeltaQOverTime, dep.var='deltaQ', colorM='blue', alpha=0.8, plot.model.points = FALSE)
# # for(ease in unique(meanDeltaQOverTime$ease)) {
# #   draw.polygon(meanDeltaQOverTime[meanDeltaQOverTime$ease==ease,], dep.var='deltaQ', colorM=get.color(ease, 'B'), plot.model.points = FALSE)
# # }
# par(xpd=TRUE)
# # legend('top', legend=c('0.8/0.2', '0.7/0.3'), col=c(2,4), pch=15, ncol = 2, bg='white', inset=c(0, -.175), bty='n')
# par(xpd=FALSE)
# 
# # sum here
# plot(0,0, type='n', xlim=c(-60, 40), ylim=c(0.1, 1), xlab='',  xaxt='n',ylab=expression(paste(Sigma, 'Q-values')), main='')
# title(expression(bold(paste(Sigma, 'Q-values'))), line=2.5)
# mtext('Trial (relative to reversal)', side=1, cex=.66, line=2)
# abline(h=seq(0, 1, .1), col='grey')
# abline(v=seq(-60, 50, 10), col='grey')
# abline(v=0, lty=2, col='black', lwd=2)
# axis(1, at=seq(-60, 50, 10), lwd=lwd.axis)
# draw.polygon(meanSumQOverTime, dep.var='sumQ', colorM='blue', alpha=.8, plot.model.points = FALSE)
# # for(ease in unique(meanDeltaQOverTime$ease)) {
# #   draw.polygon(meanSumQOverTime[meanSumQOverTime$ease==ease,], dep.var='sumQ', colorM=get.color(ease, 'B'), plot.model.points = FALSE)
# # }
# # par(xpd=TRUE)
# # legend('top', legend=c('0.8/0.2', '0.7/0.3'), col=c(2,4), pch=15, ncol = 2, bg='white', inset=c(0, -.175), bty='n')
# # par(xpd=FALSE)
# 
# # drift rates
# plot(0,0, type='n', xlim=c(-60, 40), ylim=c(0.8, 3.5), xlab='', xaxt='n', ylab='Drift rates', main='')
# mtext('Trial (relative to reversal)', side=1, cex=.66, line=2)
# title(expression(bold('Drift rates')), line=2.5)
# abline(h=seq(0, 10, .5), col='grey')
# abline(v=seq(-60, 50, 10), col='grey')
# abline(v=0, lty=2, col='black', lwd=2)
# axis(1, at=seq(-60, 50, 10), lwd=lwd.axis)
# draw.polygon(meanV1OverTimeM, dep.var='mean_v.r1', colorM=get.color(0.4, 'A'), plot.model.points = FALSE)
# draw.polygon(meanV2OverTimeM, dep.var='mean_v.r2', colorM=get.color(0.4, 'B'), plot.model.points = FALSE)
# # for(ease in unique(allV1OverTimeM$ease)) {
# #   draw.polygon(meanV1OverTimeM[meanV1OverTimeM$ease==ease,], dep.var='mean_v.r1', colorM=get.color(ease, 'A'), plot.model.points = FALSE)
# #   draw.polygon(meanV2OverTimeM[meanV2OverTimeM$ease==ease,], dep.var='mean_v.r2', colorM=get.color(ease, 'B'), plot.model.points = FALSE)
# # }
# par(xpd=TRUE)
# # legend('top', legend=c('0.8/0.2 A',
# #                        '0.8/0.2 B',
# #                        '0.7/0.3 A',
# #                        '0.7/0.3 B'), col=c(2,1,4,3), pch=15, ncol = 2, bg='white', inset=c(0, -.25), bty='n')
# legend('top', legend=c('A', 'B'), col=c(2,4), pch=15, ncol = 2, bg='white', inset=c(0, -.175), bty='n')
# par(xpd=FALSE)
# 
# if(savePlot) dev.off()
# 


# get.color <- function(ease, choice='A') {
# if(ease == "0.6") i <- 1 #return(1)
# if(ease == "0.4") i <- 3 #return(2)
# return(i+(choice=='B'))
# }
# 
# draw.polygon <- function(pp, dep.var, xaxis='bin', colorM='blue', plot.model.points=TRUE, alpha=1) {
#   lowerQ <- aggregate(as.formula(paste0(dep.var, '~bin')), pp, quantile, .025)
#   upperQ <- aggregate(as.formula(paste0(dep.var, '~bin')), pp, quantile, .975)
#   xs <- c(lowerQ[,xaxis], rev(lowerQ[,xaxis]))
#   ys <- c(lowerQ[,dep.var], rev(upperQ[,dep.var]))
#   polygon(xs, ys, col=rgb(col2rgb(colorM)[1]/255, col2rgb(colorM)[2]/255, col2rgb(colorM)[3]/255, alpha=alpha), lty = NULL, border=NA)
#   if(plot.model.points) points(pp[,xaxis], pp[,dep.var], pch=20, col=colorM, cex=.01)
# }
# 
# draw.polygon <- function(pp, dep.var, xaxis='bin', colorM='blue', plot.model.points=TRUE, alpha=1) {
#   lowerQ <- aggregate(as.formula(paste0(dep.var, '~bin')), pp, quantile, .025)
#   upperQ <- aggregate(as.formula(paste0(dep.var, '~bin')), pp, quantile, .975)
#   xs <- c(lowerQ[,xaxis], rev(lowerQ[,xaxis]))
#   ys <- c(lowerQ[,dep.var], rev(upperQ[,dep.var]))
#   polygon(xs, ys, col=rgb(col2rgb(colorM)[1]/255, col2rgb(colorM)[2]/255, col2rgb(colorM)[3]/255, alpha=alpha), lty = NULL, border=NA)
#   if(plot.model.points) points(pp[,xaxis], pp[,dep.var], pch=20, col=colorM, cex=.01)
# }
# 

# 
# # Q-values, drift rates ---------------------------------------------------
# pp3 <- tmp[['pp3']]   # make sure this is based on the RL-ARD
# allQ1OverTimeM <- do.call(rbind, (lapply(1:length(pp3), function(x) {tmp <- attr(pp3[[x]], 'SR.r1OverBins'); tmp$s <- x; tmp})))
# allQ2OverTimeM <- do.call(rbind, (lapply(1:length(pp3), function(x) {tmp <- attr(pp3[[x]], 'SR.r2OverBins'); tmp$s <- x; tmp})))
# meanQ1OverTimeM <- aggregate(SR.r1~reps*bin, allQ1OverTimeM, mean)
# meanQ2OverTimeM <- aggregate(SR.r2~reps*bin, allQ2OverTimeM, mean)
# 
# # differences
# deltaQ <- allQ1OverTimeM
# deltaQ$SR.r2 <- allQ2OverTimeM$SR.r2
# deltaQ$deltaQ <- deltaQ$SR.r1 - deltaQ$SR.r2
# 
# # sums
# deltaQ$sumQ <- deltaQ$SR.r1 + deltaQ$SR.r2
# 
# meanDeltaQOverTime <- aggregate(deltaQ~reps*bin, deltaQ, mean)
# meanSumQOverTime <- aggregate(sumQ~reps*bin, deltaQ, mean)
# 
# # drift rates
# allV1OverTimeM <- do.call(rbind, (lapply(1:length(pp3), function(x) {tmp <- attr(pp3[[x]], 'mean_v.r1OverBins'); tmp$s <- x; tmp})))
# allV2OverTimeM <- do.call(rbind, (lapply(1:length(pp3), function(x) {tmp <- attr(pp3[[x]], 'mean_v.r2OverBins'); tmp$s <- x; tmp})))
# meanV1OverTimeM <- aggregate(mean_v.r1~reps*bin, allV1OverTimeM, mean)
# meanV2OverTimeM <- aggregate(mean_v.r2~reps*bin, allV2OverTimeM, mean)
# 
# library(RColorBrewer)
# dark2 <- brewer.pal(n=8, name='Dark2')
# set2 <- brewer.pal(n=8, name='Set2')
# palette(c(set2[1], dark2[1], set2[2], dark2[2]))
# # Plot --------------------------------------------------------------------
# pdf(file='./figures/exp3-q-values2.pdf', width=7, height=2.5)
# par(mfrow=c(1,4), las=1, bty='l', oma=c(0,1,1,0), mar=c(4, 3, 4, 0.5) + 0.1, mgp=c(2.25,.75,0))
# # Q-values
# plot(0,0, type='n', ylim=c(0, .85), xlim=c(-60, 40), xlab='', xaxt='n', ylab='Q-values', main='')
# title(expression(bold('Q-values')), line=2.5)
# mtext('Trial (relative to reversal)', side=1, cex=.66, line=2)
# abline(h=seq(0, 1, .1), col='grey')
# abline(v=seq(-60, 50, 10), col='grey')
# abline(v=0, lty=2, col='black', lwd=2)
# axis(1, at=seq(-60, 50, 10), lwd=lwd.axis)
# 
# # for(ease in unique(meanQ1OverTimeM$ease)) {
#   # draw.polygon(meanQ1OverTimeM[meanQ1OverTimeM$ease==ease,], dep.var='SR.r1', colorM=get.color(ease, 'B'), plot.model.points = FALSE) #ifelse(ease=="0.6", 1, ifelse(ease=="0.4", 2, ifelse(ease=="0.3", 3, 4))))
#   # draw.polygon(meanQ2OverTimeM[meanQ2OverTimeM$ease==ease,], dep.var='SR.r2', colorM=get.color(ease, 'A'), plot.model.points = FALSE)
# # }
# draw.polygon(meanQ1OverTimeM, dep.var='SR.r1', colorM=get.color(0.4, 'B'), plot.model.points = FALSE) #ifelse(ease=="0.6", 1, ifelse(ease=="0.4", 2, ifelse(ease=="0.3", 3, 4))))
# draw.polygon(meanQ2OverTimeM, dep.var='SR.r2', colorM=get.color(0.4, 'A'), plot.model.points = FALSE)
# # par(xpd=TRUE)
# # legend('top', legend=c('0.8/0.2 A',
# #                        '0.8/0.2 B',
# #                        '0.7/0.3 A',
# #                        '0.7/0.3 B'), col=c(2,1,4,3), pch=15, ncol = 2, bg='white', inset=c(0, -.25), bty='n')
# # par(xpd=FALSE)
# par(xpd=TRUE)
# legend('top', legend=c('A', 'B'), col=c(2,4), pch=15, ncol = 2, bg='white', inset=c(0, -.175), bty='n')
# par(xpd=FALSE)
# 
# #legend('topleft', legend=c('0.8/0.2', '0.7/0.3'), col=1:2, pch=15, bty='n', title='Difficulty')
# 
# # delta here
# plot(0,0, type='n', xlim=c(-60, 40), ylim=c(-.6, .7), xlab='',  xaxt='n',
#      ylab=expression(paste(Delta, 'Q-values')), main='')# expression(bold(paste(Delta, 'Q-values'))))
# title(expression(bold(paste(Delta, 'Q-values'))), line=2.5)
# mtext('Trial (relative to reversal)', side=1, cex=.66, line=2)
# abline(h=seq(-1, 1, .1), col='grey')
# abline(v=seq(-60, 50, 10), col='grey')
# abline(v=0, lty=2, col='black', lwd=2)
# abline(h=0, lty=2, col='black', lwd=2)
# axis(1, at=seq(-60, 50, 10), lwd=lwd.axis)
# draw.polygon(meanDeltaQOverTime, dep.var='deltaQ', colorM='blue', alpha=0.8, plot.model.points = FALSE)
# # for(ease in unique(meanDeltaQOverTime$ease)) {
# #   draw.polygon(meanDeltaQOverTime[meanDeltaQOverTime$ease==ease,], dep.var='deltaQ', colorM=get.color(ease, 'B'), plot.model.points = FALSE)
# # }
# par(xpd=TRUE)
# # legend('top', legend=c('0.8/0.2', '0.7/0.3'), col=c(2,4), pch=15, ncol = 2, bg='white', inset=c(0, -.175), bty='n')
# par(xpd=FALSE)
# 
# # sum here
# plot(0,0, type='n', xlim=c(-60, 40), ylim=c(0.1, 1), xlab='',  xaxt='n',ylab=expression(paste(Sigma, 'Q-values')), main='')
# title(expression(bold(paste(Sigma, 'Q-values'))), line=2.5)
# mtext('Trial (relative to reversal)', side=1, cex=.66, line=2)
# abline(h=seq(0, 1, .1), col='grey')
# abline(v=seq(-60, 50, 10), col='grey')
# abline(v=0, lty=2, col='black', lwd=2)
# axis(1, at=seq(-60, 50, 10), lwd=lwd.axis)
# draw.polygon(meanSumQOverTime, dep.var='sumQ', colorM='blue', alpha=.8, plot.model.points = FALSE)
# # for(ease in unique(meanDeltaQOverTime$ease)) {
# #   draw.polygon(meanSumQOverTime[meanSumQOverTime$ease==ease,], dep.var='sumQ', colorM=get.color(ease, 'B'), plot.model.points = FALSE)
# # }
# # par(xpd=TRUE)
# # legend('top', legend=c('0.8/0.2', '0.7/0.3'), col=c(2,4), pch=15, ncol = 2, bg='white', inset=c(0, -.175), bty='n')
# # par(xpd=FALSE)
# 
# # drift rates
# plot(0,0, type='n', xlim=c(-60, 40), ylim=c(0.8, 3.5), xlab='', xaxt='n', ylab='Drift rates', main='')
# mtext('Trial (relative to reversal)', side=1, cex=.66, line=2)
# title(expression(bold('Drift rates')), line=2.5)
# abline(h=seq(0, 10, .5), col='grey')
# abline(v=seq(-60, 50, 10), col='grey')
# abline(v=0, lty=2, col='black', lwd=2)
# axis(1, at=seq(-60, 50, 10), lwd=lwd.axis)
# draw.polygon(meanV1OverTimeM, dep.var='mean_v.r1', colorM=get.color(0.4, 'A'), plot.model.points = FALSE)
# draw.polygon(meanV2OverTimeM, dep.var='mean_v.r2', colorM=get.color(0.4, 'B'), plot.model.points = FALSE)
# # for(ease in unique(allV1OverTimeM$ease)) {
# #   draw.polygon(meanV1OverTimeM[meanV1OverTimeM$ease==ease,], dep.var='mean_v.r1', colorM=get.color(ease, 'A'), plot.model.points = FALSE)
# #   draw.polygon(meanV2OverTimeM[meanV2OverTimeM$ease==ease,], dep.var='mean_v.r2', colorM=get.color(ease, 'B'), plot.model.points = FALSE)
# # }
# par(xpd=TRUE)
# # legend('top', legend=c('0.8/0.2 A',
# #                        '0.8/0.2 B',
# #                        '0.7/0.3 A',
# #                        '0.7/0.3 B'), col=c(2,1,4,3), pch=15, ncol = 2, bg='white', inset=c(0, -.25), bty='n')
# legend('top', legend=c('A', 'B'), col=c(2,4), pch=15, ncol = 2, bg='white', inset=c(0, -.175), bty='n')
# par(xpd=FALSE)
# 
# dev.off()


#### small
# pdf(file=paste0('./figures/exp2_reversal_', modelName, '-QQsmall.pdf'), width=6, height=2.5)
# par(oma=c(0,1,0,0), mar=c(3, 4, 2, 0.5) + 0.1, mfrow=c(1,2), mgp=c(2.75,.75,0), las=1, cex=.66, bty='l')
# plotDataPPBins(data=q10RTsOverTrials[[1]], pp=q10RTsOverTrials[[2]], 
#                dep.var='RT.10.', ylim=c(.4, 1.1), ylab='RT (s)',
#                xaxt='s', xlim=c(-60, 40), plot.model.points = FALSE, data.lwd=1.5, data.cex=1.5,
#                xlab='', hline.by=0.05, axvlines=seq(-50, 50, 10))
# 
# plotDataPPBins(data=q50RTsOverTrials[[1]], pp=q50RTsOverTrials[[2]], dep.var='RT.50.', plot.new = FALSE, draw.legend=FALSE, plot.model.points = FALSE, data.lwd=1.5, data.cex=1.5)
# plotDataPPBins(data=q90RTsOverTrials[[1]], pp=q90RTsOverTrials[[2]], dep.var='RT.90.', plot.new = FALSE, draw.legend=FALSE, plot.model.points = FALSE, data.lwd=1.5, data.cex=1.5)
# abline(v=0, lty=2)
# title('RTs')
# mtext('Trial (relative to reversal point)', side=1, cex=.66, line=2)
# 
# plotDataPPBins(data=meanAccOverTrials[[1]], pp=meanAccOverTrials[[2]],
               # xaxt='s', xlim=c(-60, 40), plot.model.points=FALSE,
               # dep.var='acc', ylab=expression('Proportion choice A'),
               # xlab = '', data.lwd=1.5, data.cex=1.5,
               # legend.pos='topright', ylim=c(0.25, 0.85), hline.by=0.05, axvlines=seq(-50, 50, 10))
# title('Choices')
# axis(1, at=seq(-50, 50, 10), labels=rep(NA, 5))
# abline(v=0, lty=2)
# mtext('Trial (relative to reversal point)', side=1, cex=.66, line=2)
# dev.off()

