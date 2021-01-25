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
modelName <- 'arw-RL-WA' #'arw-RL-mag-SAT-BV02'  #'ddm-RL-st0'
dataName <- 'expMA'
tmp <- getDataPpBPIC(modelName, dataName)
BPIC <- tmp$BPIC
excludePerfectAcc <- FALSE
if(excludePerfectAcc) {
  tmp2 <- tmp
  tmp2[['data3']] <- tmp[['data3']][sapply(tmp[['data3']], function(x) mean(x[x$bin==1&x$ease=='0.6','acc'])<1)]
  tmp2[['pp3']] <- tmp[['pp3']][sapply(tmp[['data3']], function(x) mean(x[x$bin==1&x$ease=='0.6','acc'])<1)]
  data3 <- tmp2[['data3']]
  pp3 <- tmp2[['pp3']]
  modelName <- paste0(modelName, '-exclperfacc')
} else {
  data3 <- tmp[['data3']]
  pp3 <- tmp[['pp3']]
}

data3 <- do.call(rbind, tmp[['data3']])
pp3 <- do.call(rbind, tmp[['pp3']])
getDensities <- function(dat) {
  d1 <- density(dat$RT[dat$acc==1])
  d2 <- density(dat$RT[dat$acc==0])
  acc <- mean(dat$acc)
  d1$y <- d1$y*acc
  d2$y <- d2$y*(1-acc)
  return(list(d1=d1, d2=d2))
}

plotdPDFs <- function(data, pp, legend.pos='topleft', colorM='cornflowerblue', colorD='black', alpha=0.3, 
                      yaxis.mult=1.3, ylab='', xlab='') {
  densities <- getDensities(data)
#  print(max(c(densities[[1]]$y, densities[[2]]$y)))
  plot(0,0, xlab=xlab, ylab=ylab, type='n',
       xlim=c(-max(data$RT), max(data$RT)), ylim=c(0, max(c(densities[[1]]$y, densities[[2]]$y)))*c(0,yaxis.mult))
  abline(h=seq(0, 10, 0.25), col='grey')
  abline(v=seq(-3, 3, 0.5), col='grey')
  
  for(rep in 1:max(pp$reps)) {
    densities <- getDensities(pp[pp$reps==rep,])
    lines(densities[[1]]$x, densities[[1]]$y, type='l', lwd=1, col=rgb(col2rgb(colorM)[1]/255, col2rgb(colorM)[2]/255, col2rgb(colorM)[3]/255, alpha=alpha))
    lines(-densities[[2]]$x, densities[[2]]$y, type='l', lwd=1, col=rgb(col2rgb(colorM)[1]/255, col2rgb(colorM)[2]/255, col2rgb(colorM)[3]/255, alpha=alpha))
  }
  densities <- getDensities(data)
  lines(densities[[1]]$x, densities[[1]]$y, type='l', lwd=2, col=colorD)
  lines(-densities[[2]]$x, densities[[2]]$y, type='l', lwd=2, col=colorD)
  legend(legend.pos, c('Data', 'Model'), lty=c(1, NA), 
         fill=c(NA, rgb(col2rgb(colorM)[1]/255, col2rgb(colorM)[2]/255, col2rgb(colorM)[3]/255, alpha=1)),
         border=c(NA, 'black'),
         lwd=c(2, NA), col=c(colorD,colorM), pch=c(NA, NA), bty='n')
}

## all subjects
if(savePlot) pdf(file=paste0('./figures/', dataName, '-defpdfs-allsubs.pdf'), width=7, height=7)
par(oma=c(3,3,0,0), mar=c(0, 2, 3, 0) + 0.1, mfrow=c(3,4), mgp=c(2,.75,0), las=1, bty='l')
layout(matrix(c(1,1:11), nrow=3, byrow=TRUE))
plotdPDFs(data=data3, pp=pp3, alpha=1, yaxis.mult=1)
mtext('Grand average', side=3, line=0, font=2, cex=.66*1.2)
for(subN in 1:length(tmp[['data3']])) {
  plotdPDFs(data=tmp[['data3']][[subN]], pp=tmp[['pp3']][[subN]])
  mtext(paste0('Subject ', subN), side=3, line=0, font=2, cex=.66*1.2)
  if(subN == 10) {
    layout(matrix(c(1:12), nrow=3, byrow=TRUE))
  }
  if(subN %in% seq(10, 100, 12)) {
    mtext('Defective probability density', side=2, outer=TRUE, cex=.66*1.2, font=1, las=0, line=1)
    mtext('RT (s)', side=1, outer=TRUE, cex=.66*1.2, font=1, line=2)
  }
}
mtext('Defective probability density', side=2, outer=TRUE, cex=.66*1.2, font=1, las=0, line=1)
mtext('RT (s)', side=1, outer=TRUE, cex=.66*1.2, font=1, line=2)
if(savePlot) dev.off()

## Plot for paper?
if(savePlot) pdf(file=paste0('./figures/', dataName, '-defpdfs-paper.pdf'), width=7, height=7)
par(oma=c(3,3,0,0), mar=c(0, 2, 3, 0) + 0.1, mfrow=c(3,4), mgp=c(2,.75,0), las=1, bty='l')
layout(matrix(c(1,1:11), nrow=3, byrow=TRUE))
plotdPDFs(data=data3, pp=pp3, alpha=1, yaxis.mult=1)
mtext('Grand average', side=3, line=0, font=2, cex=.66*1.2)
#title('Grand average')
for(subN in 1:10) {
  plotdPDFs(data=tmp[['data3']][[subN]], pp=tmp[['pp3']][[subN]])
  mtext(paste0('Subject ', subN), side=3, line=0, font=2, cex=.66*1.2)
}
mtext('Defective probability density', side=2, outer=TRUE, cex=.66*1.2, font=1, las=0, line=1)
mtext('RT (s)', side=1, outer=TRUE, cex=.66*1.2, font=1, line=2)
if(savePlot) dev.off()




# Experiment 1: standard errors ------------------------------------------
# Okay, now get data + standard errors
getDescriptivesSE <- function(x, dep.var='RT', attr.name='RTsOverBins', id.var1='~reps*bin*s', id.var2='~reps*bin') {
  #  print('huh')
  allOverTime <- do.call(rbind, (lapply(1:length(x), function(y) {tmp <- attr(x[[y]], attr.name); tmp$s <- y; tmp})))
  #  print(head(allOverTime))
  
  form1 <- as.formula(paste0(dep.var, id.var1))
  res <- aggregate(form1, allOverTime, function(x) sd(x)/sqrt(length(x)))
  
  if(!is.null(id.var2)) {
    form2 <- as.formula(paste0(dep.var, id.var2))
    res <- aggregate(form2, aggregate(form1, res, mean), mean)
  }
  return(res)
}


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

SEq10RTsByEase <- list(getDescriptivesSE(data3, dep.var='RT.10.', attr.name='qRTsCorrectByEase', id.var1='~bin*ease', id.var2=NULL)) 
SEq50RTsByEase <- list(getDescriptivesSE(data3, dep.var='RT.50.', attr.name='qRTsCorrectByEase', id.var1='~bin*ease', id.var2=NULL))
SEq90RTsByEase <- list(getDescriptivesSE(data3, dep.var='RT.90.', attr.name='qRTsCorrectByEase', id.var1='~bin*ease', id.var2=NULL))
SEq10RTsByEaseE <- list(getDescriptivesSE(data3, dep.var='RT.10.', attr.name='qRTsErrorByEase', id.var1='~bin*ease', id.var2=NULL))
SEq50RTsByEaseE <- list(getDescriptivesSE(data3, dep.var='RT.50.', attr.name='qRTsErrorByEase', id.var1='~bin*ease', id.var2=NULL))
SEq90RTsByEaseE <- list(getDescriptivesSE(data3, dep.var='RT.90.', attr.name='qRTsErrorByEase', id.var1='~bin*ease', id.var2=NULL))
SEmeanAccByEase <- list(getDescriptivesSE(data3, dep.var='acc', attr.name='AccByEase', id.var1='~bin*ease', id.var2=NULL))


# Plot posterior predictives
if(savePlot) pdf(file=paste0('./figures/exp1_difficulty_', modelName, '-QQ-horizontal-SE.pdf'), width=7, height=7/4*3)
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
  arrows(x0=SEmeanAccByEase[[1]][idxD,'bin'], x1=SEmeanAccByEase[[1]][idxD,'bin'],
         y0=meanAccByEase[[1]][idxD,'acc']-SEmeanAccByEase[[1]][idxD,'acc'], 
         y1=meanAccByEase[[1]][idxD,'acc']+SEmeanAccByEase[[1]][idxD,'acc'], code=3, angle=90, lwd=1.5, length = 0.05)
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
  arrows(x0=q10RTsByEase[[1]][idxD,'bin'], x1=q10RTsByEase[[1]][idxD,'bin'],
         y0=q10RTsByEase[[1]][idxD,'RT.10.']-SEq10RTsByEase[[1]][idxD,'RT.10.'], 
         y1=q10RTsByEase[[1]][idxD,'RT.10.']+SEq10RTsByEase[[1]][idxD,'RT.10.'], code=3, angle=90, lwd=1.5, length = 0.05)
  arrows(x0=q50RTsByEase[[1]][idxD,'bin'], x1=q50RTsByEase[[1]][idxD,'bin'],
         y0=q50RTsByEase[[1]][idxD,'RT.50.']-SEq50RTsByEase[[1]][idxD,'RT.50.'], 
         y1=q50RTsByEase[[1]][idxD,'RT.50.']+SEq50RTsByEase[[1]][idxD,'RT.50.'], code=3, angle=90, lwd=1.5, length = 0.05)
  arrows(x0=q90RTsByEase[[1]][idxD,'bin'], x1=q90RTsByEase[[1]][idxD,'bin'],
         y0=q90RTsByEase[[1]][idxD,'RT.90.']-SEq90RTsByEase[[1]][idxD,'RT.90.'], 
         y1=q90RTsByEase[[1]][idxD,'RT.90.']+SEq90RTsByEase[[1]][idxD,'RT.90.'], code=3, angle=90, lwd=1.5, length = 0.05)
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
  arrows(x0=q10RTsByEaseE[[1]][idxD,'bin'], x1=q10RTsByEaseE[[1]][idxD,'bin'],
         y0=q10RTsByEaseE[[1]][idxD,'RT.10.']-SEq10RTsByEaseE[[1]][idxD,'RT.10.'], 
         y1=q10RTsByEaseE[[1]][idxD,'RT.10.']+SEq10RTsByEaseE[[1]][idxD,'RT.10.'], code=3, angle=90, lwd=1.5, length = 0.05)
  arrows(x0=q50RTsByEaseE[[1]][idxD,'bin'], x1=q50RTsByEaseE[[1]][idxD,'bin'],
         y0=q50RTsByEaseE[[1]][idxD,'RT.50.']-SEq50RTsByEaseE[[1]][idxD,'RT.50.'], 
         y1=q50RTsByEaseE[[1]][idxD,'RT.50.']+SEq50RTsByEaseE[[1]][idxD,'RT.50.'], code=3, angle=90, lwd=1.5, length = 0.05)
  arrows(x0=q90RTsByEaseE[[1]][idxD,'bin'], x1=q90RTsByEaseE[[1]][idxD,'bin'],
         y0=q90RTsByEaseE[[1]][idxD,'RT.90.']-SEq90RTsByEaseE[[1]][idxD,'RT.90.'], 
         y1=q90RTsByEaseE[[1]][idxD,'RT.90.']+SEq90RTsByEaseE[[1]][idxD,'RT.90.'], code=3, angle=90, lwd=1.5, length = 0.05)
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





# Experiment 2: Standard errors -------------------------------------------


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

data3 <- tmp[['data3']]
SEq10RTs <- list(getDescriptivesSE(data3, dep.var='RT.10.', attr.name='qRTs', id.var1='~bin', id.var2=NULL)) 
SEq50RTs <- list(getDescriptivesSE(data3, dep.var='RT.50.', attr.name='qRTs', id.var1='~bin', id.var2=NULL))
SEq90RTs <- list(getDescriptivesSE(data3, dep.var='RT.90.', attr.name='qRTs', id.var1='~bin', id.var2=NULL))
SEq10RTsE <- list(getDescriptivesSE(data3, dep.var='RT.10.', attr.name='qRTs', id.var1='~bin', id.var2=NULL))
SEq50RTsE <- list(getDescriptivesSE(data3, dep.var='RT.50.', attr.name='qRTs', id.var1='~bin', id.var2=NULL))
SEq90RTsE <- list(getDescriptivesSE(data3, dep.var='RT.90.', attr.name='qRTs', id.var1='~bin', id.var2=NULL))
SEmeanAcc <- list(getDescriptivesSE(data3, dep.var='acc', attr.name='AccOverBins', id.var1='~bin', id.var2=NULL))


# # Plot --------------------------------------------------------------------
if(savePlot) pdf('./figures/exp2-reversals-2x2-SE.pdf', width=5, height=4)
par(oma=c(3,4,2,0), mar=c(1, 2, 1, 0.5) + 0.1, mfcol=c(2,2),
    mgp=c(2.75,.75,0), las=1, bty='l', cex=.66)
i <- 0
data.cex=0.75
lwd.axis=1.5
for(qRTs in allqRTs) {
  i <- i+1
  plotDataPPBins(data=qRTs$meanAccOverTime[[1]], pp=qRTs$meanAccOverTime[[2]],
                 xaxt='n', yaxt='n', xlim=c(-60, 40), plot.model.points=FALSE,
                 dep.var='acc', ylab=expression('Proportion choice A'),
                 xlab = '', data.lwd=1.5, data.cex=data.cex,
                 draw.legend=i==1,
                 legend.pos='bottomleft', ylim=c(0.25, 0.85), hline.by=0.05, axvlines=seq(-50, 50, 10))
  arrows(x0=SEmeanAcc[[1]][,'bin'], x1=SEmeanAcc[[1]][,'bin'],
         y0=qRTs$meanAccOverTime[[1]][,'acc']-SEmeanAcc[[1]][,'acc'], 
         y1=qRTs$meanAccOverTime[[1]][,'acc']+SEmeanAcc[[1]][,'acc'], code=3, angle=90, lwd=1, length = 0.05, col=rgb(0,0,0,.5))
  abline(v=0, lty=2, lwd=2)
  abline(h=0.5, lty=2, lwd=2)
  axis(1, at=seq(-60, 50, 10), labels=NA, lwd=lwd.axis)
  if(i == 1) {
    mtext('Proportion choice A', side=2, cex=.66, line=3, las=0, font=1)
    mtext('RL-DDM', side=3, cex=.66*1.2, font=2, line=2)
    axis(2, at=seq(.3, .9, .1), lwd=lwd.axis)
  } else {
    mtext('RL-ARD', side=3, cex=.66*1.2, font=2, line=2)
    axis(2, at=seq(.3, .9, .1), labels=rep(NA, 5), lwd=lwd.axis)
  }
  mtext(paste0('BPIC = ', round(apply(allBPICs, 2, sum)[i])), line=1, cex=.66)
  
  ##
  plotDataPPBins(data=qRTs$q10RTsOverTime[[1]], pp=qRTs$q10RTsOverTime[[2]],
                 dep.var='RT.10.', ylim=c(.4, 1.1), ylab='RT (s)',
                 xaxt='n',  yaxt='n', xlim=c(-60, 40), plot.model.points = FALSE, data.lwd=1.5, data.cex=data.cex,
                 xlab='', hline.by=0.05, axvlines=seq(-50, 50, 10), legend.pos = FALSE)
  plotDataPPBins(data=qRTs$q50RTsOverTime[[1]], pp=qRTs$q50RTsOverTime[[2]], dep.var='RT.50.', plot.new = FALSE, draw.legend=FALSE,
                 plot.model.points=FALSE,data.cex = data.cex)
  plotDataPPBins(data=qRTs$q90RTsOverTime[[1]], pp=qRTs$q90RTsOverTime[[2]], dep.var='RT.90.', plot.new = FALSE, draw.legend=FALSE,
                 plot.model.points=FALSE, data.cex = data.cex)
  arrows(x0=qRTs$q10RTsOverTime[[1]][,'bin'], x1=qRTs$q10RTsOverTime[[1]][,'bin'],
         y0=qRTs$q10RTsOverTime[[1]][,'RT.10.']-SEq10RTs[[1]][,'RT.10.'], 
         y1=qRTs$q10RTsOverTime[[1]][,'RT.10.']+SEq10RTs[[1]][,'RT.10.'], code=3, angle=90, lwd=1, length = 0.05, col=rgb(0,0,0,.5))
  arrows(x0=qRTs$q50RTsOverTime[[1]][,'bin'], x1=qRTs$q50RTsOverTime[[1]][,'bin'],
         y0=qRTs$q50RTsOverTime[[1]][,'RT.50.']-SEq50RTs[[1]][,'RT.50.'], 
         y1=qRTs$q50RTsOverTime[[1]][,'RT.50.']+SEq50RTs[[1]][,'RT.50.'], code=3, angle=90, lwd=1, length = 0.05, col=rgb(0,0,0,.5))
  arrows(x0=qRTs$q90RTsOverTime[[1]][,'bin'], x1=qRTs$q90RTsOverTime[[1]][,'bin'],
         y0=qRTs$q90RTsOverTime[[1]][,'RT.90.']-SEq90RTs[[1]][,'RT.90.'], 
         y1=qRTs$q90RTsOverTime[[1]][,'RT.90.']+SEq90RTs[[1]][,'RT.90.'], code=3, angle=90, lwd=1, length = 0.05, col=rgb(0,0,0,.5))
  
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
if(savePlot) dev.off()


# Experiment 3 (SAT) ------------------------------------------------------
calculateByBin <- function(df) {
  df$acc <- as.integer(df$R)==2
  
  attr(df, 'qRTs') <- do.call(data.frame, aggregate(RT~reps*bin, df, quantile, probs=seq(.1, .9, .4))) #cbind(quants[quants$acc==1,c('reps', 'bin')], quants[quants$acc==1,'RT'][,1])
  attr(df, 'qRTsCorrect') <- do.call(data.frame, aggregate(RT~reps*bin, df[df$acc==1,], quantile, probs=seq(.1, .9, .4))) #cbind(quants[quants$acc==1,c('reps', 'bin')], quants[quants$acc==1,'RT'][,1])
  attr(df, 'qRTsError') <- do.call(data.frame, aggregate(RT~reps*bin, df[df$acc==0,], quantile, probs=seq(.1, .9, .4))) #cbind(quants[quants$acc==1,c('reps', 'bin')], quants[quants$acc==1,'RT'][,1])
  attr(df, 'qRTsCorrectByCue') <- do.call(data.frame, aggregate(RT~reps*bin*cue, df[df$acc==1,], quantile, probs=seq(.1, .9, .4))) #cbind(quants[quants$acc==1,c('reps', 'bin')], quants[quants$acc==1,'RT'][,1])
  attr(df, 'qRTsErrorByCue') <- do.call(data.frame, aggregate(RT~reps*bin*cue, df[df$acc==0,], quantile, probs=seq(.1, .9, .4))) #cbind(quants[quants$acc==1,c('reps', 'bin')], quants[quants$acc==1,'RT'][,1])
  
  attr(df, 'RTsOverBins') <- aggregate(RT~reps*bin, df, mean)
  attr(df, 'AccOverBins') <- aggregate(acc~reps*bin, df, mean)
  attr(df, 'SkewOverBins') <- aggregate(RT~reps*bin, df, skewness)
  attr(df, 'RTsByChoiceByEase') <- aggregate(RT~reps*bin*ease*R, df, mean)
  attr(df, 'RTsByEase') <- aggregate(RT~reps*ease*bin, df, mean)
  attr(df, 'AccByEase') <- aggregate(acc~reps*ease*bin, df, mean)
  attr(df, 'SkewByEase') <- aggregate(RT~reps*ease*bin, df, skewness)
  
  attr(df, 'RTsByCue') <- aggregate(RT~reps*cue*bin, df, mean)
  attr(df, 'AccByCue') <- aggregate(acc~reps*cue*bin, df, mean)
  attr(df, 'SkewByCue') <- aggregate(RT~reps*cue*bin, df, skewness)
  
  
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

getqRTsByCue <- function(data3, pp3) {
  q10RTsByCue <- list(getDescriptives(data3, dep.var='RT.10.', attr.name='qRTsCorrectByCue', id.var1='~bin*cue', id.var2=NULL),
                      getDescriptives(pp3, dep.var='RT.10.', attr.name='qRTsCorrectByCue', id.var1='~reps*bin*cue', id.var2=NULL))
  q50RTsByCue <- list(getDescriptives(data3, dep.var='RT.50.', attr.name='qRTsCorrectByCue', id.var1='~bin*cue', id.var2=NULL),
                      getDescriptives(pp3, dep.var='RT.50.', attr.name='qRTsCorrectByCue', id.var1='~reps*bin*cue', id.var2=NULL))
  q90RTsByCue <- list(getDescriptives(data3, dep.var='RT.90.', attr.name='qRTsCorrectByCue', id.var1='~bin*cue', id.var2=NULL),
                      getDescriptives(pp3, dep.var='RT.90.', attr.name='qRTsCorrectByCue', id.var1='~reps*bin*cue', id.var2=NULL))
  q10RTsByCueE <- list(getDescriptives(data3, dep.var='RT.10.', attr.name='qRTsErrorByCue', id.var1='~bin*cue', id.var2=NULL),
                       getDescriptives(pp3, dep.var='RT.10.', attr.name='qRTsErrorByCue', id.var1='~reps*bin*cue', id.var2=NULL))
  q50RTsByCueE <- list(getDescriptives(data3, dep.var='RT.50.', attr.name='qRTsErrorByCue', id.var1='~bin*cue', id.var2=NULL),
                       getDescriptives(pp3, dep.var='RT.50.', attr.name='qRTsErrorByCue', id.var1='~reps*bin*cue', id.var2=NULL))
  q90RTsByCueE <- list(getDescriptives(data3, dep.var='RT.90.', attr.name='qRTsErrorByCue', id.var1='~bin*cue', id.var2=NULL),
                       getDescriptives(pp3, dep.var='RT.90.', attr.name='qRTsErrorByCue', id.var1='~reps*bin*cue', id.var2=NULL))
  meanAccByCue <- list(getDescriptives(data3, dep.var='acc', attr.name='AccByCue', id.var1='~bin*cue', id.var2=NULL),
                       getDescriptives(pp3, dep.var='acc', attr.name='AccByCue', id.var1='~reps*bin*cue', id.var2=NULL))
  return(list('q10RTsByCue'=q10RTsByCue,
              'q50RTsByCue'=q50RTsByCue,
              'q90RTsByCue'=q90RTsByCue,
              'q10RTsByCueE'=q10RTsByCueE,
              'q50RTsByCueE'=q50RTsByCueE,
              'q90RTsByCueE'=q90RTsByCueE,
              'meanAccByCue'=meanAccByCue))
}

# Load quantiles of winning model & RL-fARD, BPICs ---------------------------------------------------
# DDM
tmp <- getDataPpBPIC('ddm-RL-SAT-a', 'exp3')
BPICDDM <- getDataPpBPIC('ddm-RL-SAT-a', 'exp3', BPIConly = TRUE)$BPIC
qRTsDDM <- getqRTsByCue(tmp[['data3']], tmp[['pp3']])

# ARDf
tmp <- getDataPpBPIC('arw-RL-mag-SAT-BV02', 'exp3')
BPICARD <- tmp$BPIC
qRTsARD <- getqRTsByCue(tmp[['data3']], tmp[['pp3']])

# Combine -----------------------------------------------------------------
allqRTs <- list(qRTsDDM, qRTsARD)
allBPICs <- cbind(BPICDDM[,2], BPICARD[,2])
sBPICs <- apply(allBPICs, 2, sum)
sBPICs - sBPICs[2]

data3 <- tmp[['data3']]
SEq10RTsByCue <- list(getDescriptivesSE(data3, dep.var='RT.10.', attr.name='qRTsCorrectByCue', id.var1='~bin*cue', id.var2=NULL))
SEq50RTsByCue <- list(getDescriptivesSE(data3, dep.var='RT.50.', attr.name='qRTsCorrectByCue', id.var1='~bin*cue', id.var2=NULL))
SEq90RTsByCue <- list(getDescriptivesSE(data3, dep.var='RT.90.', attr.name='qRTsCorrectByCue', id.var1='~bin*cue', id.var2=NULL))
SEq10RTsByCueE <- list(getDescriptivesSE(data3, dep.var='RT.10.', attr.name='qRTsErrorByCue', id.var1='~bin*cue', id.var2=NULL))
SEq50RTsByCueE <- list(getDescriptivesSE(data3, dep.var='RT.50.', attr.name='qRTsErrorByCue', id.var1='~bin*cue', id.var2=NULL))
SEq90RTsByCueE <- list(getDescriptivesSE(data3, dep.var='RT.90.', attr.name='qRTsErrorByCue', id.var1='~bin*cue', id.var2=NULL))
SEmeanAccByCue <- list(getDescriptivesSE(data3, dep.var='acc', attr.name='AccByCue', id.var1='~bin*cue', id.var2=NULL))



# Plot --------------------------------------------------------------------
layoutM <- matrix(1:25, nrow=5, byrow=TRUE)
layoutM[c(1, 5),1:2] <- 1
layoutM[c(1, 5),4:5] <- 9
layoutM[2:4,1:2] <- 2:7 #matrix(c(2:13), nrow=3, byrow=TRUE)
layoutM[2:4,4:5] <- 10:15 #matrix(c(2:13), nrow=3, byrow=TRUE)
layoutM[,3] <- 8
layoutM

if(savePlot) pdf('./figures/exp3-SAT-SE.pdf', width=7, height=7/4*3)
layout(layoutM, heights = c(0.01, .8, 1, 1, 0.01), widths=c(1,1,.1,1,1))
par(oma=c(3,4,2,0), mar=c(0, 0, 1, 0.5) + 0.1, #mfcol=c(3,4), 
    mgp=c(2.75,.75,0), las=1, bty='l')
i <- 0
data.cex=1.5
corrRTylim <- errRTylim <- c(.35,1.1)
for(qRTs in allqRTs) {
  plot.new()
  if(i == 0) { mtext('RL-DDM', side=3, cex=.66*1.2, font=2, line=2); mtext(paste0('BPIC = ', round(apply(allBPICs, 2, sum)[1])), cex=.66, line=1)}
  if(i == 2) {plot.new(); mtext('RL-ARD', side=3, cex=.66*1.2, font=2, line=2); mtext(paste0('BPIC = ', round(apply(allBPICs, 2, sum)[2])), cex=.66, line=1)}
  for(cue in c('SPD', 'ACC')) {
    i <- i+1
    idxD <- qRTs$meanAccByCue[[1]]$cue==cue
    idxM <- qRTs$meanAccByCue[[2]]$cue==cue
    
    plotDataPPBins(data=qRTs$meanAccByCue[[1]][idxD,], pp=qRTs$meanAccByCue[[2]][idxM,],
                   xaxt='n', draw.legend = i==1, data.cex = data.cex,
                   dep.var='acc', ylab='', xlab = '', yaxt='n',
                   legend.pos='bottomright', ylim=c(0.5, 0.9), hline.by=0.1)
    arrows(x0=qRTs$meanAccByCue[[1]][idxD,'bin'], x1=qRTs$meanAccByCue[[1]][idxD,'bin'],
           y0=qRTs$meanAccByCue[[1]][idxD,'acc']-SEmeanAccByCue[[1]][idxD,'acc'], 
           y1=qRTs$meanAccByCue[[1]][idxD,'acc']+SEmeanAccByCue[[1]][idxD,'acc'], code=3, angle=90, lwd=1.5, length = 0.05)
    
    axis(1, at=seq(2, 10, 2), labels=rep(NA, 5), lwd=2)
    if(i == 1) {
      mtext('Accuracy', side=2, cex=.66, line=3, las=0, font=1)
      axis(2, at=seq(.5, .9, .1), lwd=1.5)
    } else {
      axis(2, at=seq(.5, .9, .1), labels=rep(NA, 5), lwd=1.5)
    }
    if(i == 1) title('Speed')
    if(i == 2) title('Accuracy')
    if(i == 3) title('Speed')
    if(i == 4) title('Accuracy')
    
    ##
    plotDataPPBins(data=qRTs$q10RTsByCue[[1]][idxD,], pp=qRTs$q10RTsByCue[[2]][idxM,], dep.var='RT.10.', 
                   ylim=corrRTylim, xaxt='n', ylab='', yaxt='n', draw.legend = FALSE, data.cex = data.cex)
    plotDataPPBins(data=qRTs$q50RTsByCue[[1]][idxD,], pp=qRTs$q50RTsByCue[[2]][idxM,], dep.var='RT.50.', plot.new = FALSE, draw.legend=FALSE, data.cex = data.cex)
    plotDataPPBins(data=qRTs$q90RTsByCue[[1]][idxD,], pp=qRTs$q90RTsByCue[[2]][idxM,], dep.var='RT.90.', plot.new = FALSE, draw.legend=FALSE, data.cex = data.cex)
    arrows(x0=qRTs$q10RTsByCue[[1]][idxD,'bin'], x1=qRTs$q10RTsByCue[[1]][idxD,'bin'],
           y0=qRTs$q10RTsByCue[[1]][idxD,'RT.10.']-SEq10RTsByCue[[1]][idxD,'RT.10.'], 
           y1=qRTs$q10RTsByCue[[1]][idxD,'RT.10.']+SEq10RTsByCue[[1]][idxD,'RT.10.'], code=3, angle=90, lwd=1.5, length = 0.05)
    arrows(x0=qRTs$q10RTsByCue[[1]][idxD,'bin'], x1=qRTs$q10RTsByCue[[1]][idxD,'bin'],
           y0=qRTs$q50RTsByCue[[1]][idxD,'RT.50.']-SEq50RTsByCue[[1]][idxD,'RT.50.'], 
           y1=qRTs$q50RTsByCue[[1]][idxD,'RT.50.']+SEq50RTsByCue[[1]][idxD,'RT.50.'], code=3, angle=90, lwd=1.5, length = 0.05)
    arrows(x0=qRTs$q90RTsByCue[[1]][idxD,'bin'], x1=qRTs$q10RTsByCue[[1]][idxD,'bin'],
           y0=qRTs$q90RTsByCue[[1]][idxD,'RT.90.']-SEq90RTsByCue[[1]][idxD,'RT.90.'], 
           y1=qRTs$q90RTsByCue[[1]][idxD,'RT.90.']+SEq90RTsByCue[[1]][idxD,'RT.90.'], code=3, angle=90, lwd=1.5, length = 0.05)
    axis(1, at=seq(2, 10, 2), labels=NA, lwd=1.5)
    if(i == 1) {
      mtext('Correct RTs (s)', side=2, cex=.66, line=3, las=0, font=1)
      axis(2, seq(.4, 1.2, .2), lwd=1.5)
    } else {
      axis(2, seq(.4, 1.2, .2), labels=NA, lwd=1.5)
    }
    
    ##
    plotDataPPBins(data=qRTs$q10RTsByCueE[[1]][idxD,], pp=qRTs$q10RTsByCueE[[2]][idxM,], dep.var='RT.10.', ylim=errRTylim, xaxt='n', ylab='', yaxt='n', draw.legend = FALSE, data.cex = data.cex)
    plotDataPPBins(data=qRTs$q50RTsByCueE[[1]][idxD,], pp=qRTs$q50RTsByCueE[[2]][idxM,], dep.var='RT.50.', plot.new = FALSE, draw.legend=FALSE, data.cex = data.cex)
    plotDataPPBins(data=qRTs$q90RTsByCueE[[1]][idxD,], pp=qRTs$q90RTsByCueE[[2]][idxM,], dep.var='RT.90.', plot.new = FALSE, draw.legend=FALSE, data.cex = data.cex)
    arrows(x0=qRTs$q10RTsByCueE[[1]][idxD,'bin'], x1=qRTs$q10RTsByCueE[[1]][idxD,'bin'],
           y0=qRTs$q10RTsByCueE[[1]][idxD,'RT.10.']-SEq10RTsByCueE[[1]][idxD,'RT.10.'], 
           y1=qRTs$q10RTsByCueE[[1]][idxD,'RT.10.']+SEq10RTsByCueE[[1]][idxD,'RT.10.'], code=3, angle=90, lwd=1.5, length = 0.05)
    arrows(x0=qRTs$q10RTsByCueE[[1]][idxD,'bin'], x1=qRTs$q10RTsByCueE[[1]][idxD,'bin'],
           y0=qRTs$q50RTsByCueE[[1]][idxD,'RT.50.']-SEq50RTsByCueE[[1]][idxD,'RT.50.'], 
           y1=qRTs$q50RTsByCueE[[1]][idxD,'RT.50.']+SEq50RTsByCueE[[1]][idxD,'RT.50.'], code=3, angle=90, lwd=1.5, length = 0.05)
    arrows(x0=qRTs$q90RTsByCueE[[1]][idxD,'bin'], x1=qRTs$q10RTsByCueE[[1]][idxD,'bin'],
           y0=qRTs$q90RTsByCueE[[1]][idxD,'RT.90.']-SEq90RTsByCueE[[1]][idxD,'RT.90.'], 
           y1=qRTs$q90RTsByCueE[[1]][idxD,'RT.90.']+SEq90RTsByCueE[[1]][idxD,'RT.90.'], code=3, angle=90, lwd=1.5, length = 0.05)
    if(i == 1) {
      mtext('Error RTs (s)', side=2, cex=.66, line=3, las=0, font=1)
      axis(2, seq(.4, 1.2, .2), lwd=1.5)
    } else {
      axis(2, seq(.4, 1.2, .2), labels=NA, lwd=1.5)
    }
    axis(1, at=seq(2, 10, 2), lwd=1.5)
    mtext('Trial bin', side=1, cex=.66, line=2)
  }
}
if(savePlot) dev.off()
