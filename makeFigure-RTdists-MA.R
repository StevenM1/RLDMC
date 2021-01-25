rm(list=ls())
library(snowfall)
source ("dmc/dmc.R")
source('utils.R')
source('models.R')
samplesDir <- 'samples'
savePlot <- FALSE

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
calculateByBin <- function(df) {
  df$acc <- as.integer(df$R)==3
  
  attr(df, 'qRTs') <- do.call(data.frame, aggregate(RT~reps*bin, df, quantile, probs=seq(.1, .9, .4))) #cbind(quants[quants$acc==1,c('reps', 'bin')], quants[quants$acc==1,'RT'][,1])
  attr(df, 'qRTsTarget') <- do.call(data.frame, aggregate(RT~reps*bin, df[df$R=='r3',], quantile, probs=seq(.1, .9, .4))) #cbind(quants[quants$acc==1,c('reps', 'bin')], quants[quants$acc==1,'RT'][,1])
  attr(df, 'qRTsDistractors') <- do.call(data.frame, aggregate(RT~reps*bin, df[df$R!='r3',], quantile, probs=seq(.1, .9, .4))) #cbind(quants[quants$acc==1,c('reps', 'bin')], quants[quants$acc==1,'RT'][,1])
  attr(df, 'qRTsDistractor1') <- do.call(data.frame, aggregate(RT~reps*bin, df[df$R=='r1',], quantile, probs=seq(.1, .9, .4))) #cbind(quants[quants$acc==1,c('reps', 'bin')], quants[quants$acc==1,'RT'][,1])
  attr(df, 'qRTsDistractor2') <- do.call(data.frame, aggregate(RT~reps*bin, df[df$R=='r2',], quantile, probs=seq(.1, .9, .4))) #cbind(quants[quants$acc==1,c('reps', 'bin')], quants[quants$acc==1,'RT'][,1])
  
  attr(df, 'qRTsTargetByStimType') <- do.call(data.frame, aggregate(RT~reps*bin*stimulus_type, df[df$R=='r3',], quantile, probs=seq(.1, .9, .4))) #cbind(quants[quants$acc==1,c('reps', 'bin')], quants[quants$acc==1,'RT'][,1])
  attr(df, 'qRTsDistractorsByStimType') <- do.call(data.frame, aggregate(RT~reps*bin*stimulus_type, df[df$R!='r3',], quantile, probs=seq(.1, .9, .4))) #cbind(quants[quants$acc==1,c('reps', 'bin')], quants[quants$acc==1,'RT'][,1])
  attr(df, 'qRTsDistrator1ByStimType') <- do.call(data.frame, aggregate(RT~reps*bin*stimulus_type, df[df$R=='r1',], quantile, probs=seq(.1, .9, .4))) #cbind(quants[quants$acc==1,c('reps', 'bin')], quants[quants$acc==1,'RT'][,1])
  attr(df, 'qRTsDistrator2ByStimType') <- do.call(data.frame, aggregate(RT~reps*bin*stimulus_type, df[df$R=='r2',], quantile, probs=seq(.1, .9, .4))) #cbind(quants[quants$acc==1,c('reps', 'bin')], quants[quants$acc==1,'RT'][,1])
  
  attr(df, 'accByStimType') <- do.call(data.frame, aggregate(acc~reps*bin*stimulus_type, df, mean)) #cbind(quants[quants$acc==1,c('reps', 'bin')], quants[quants$acc==1,'RT'][,1])
  
  attr(df, 'RTsOverBins') <- aggregate(RT~reps*bin, df, mean)
  attr(df, 'AccOverBins') <- aggregate(acc~reps*bin, df, mean)
  attr(df, 'SkewOverBins') <- aggregate(RT~reps*bin, df, skewness)
  attr(df, 'RTsByChoiceByEase') <- aggregate(RT~reps*bin*ease*R, df, mean)
  attr(df, 'RTsByEase') <- aggregate(RT~reps*ease*bin, df, mean)
  attr(df, 'AccByEase') <- aggregate(acc~reps*ease*bin, df, mean)
  attr(df, 'SkewByEase') <- aggregate(RT~reps*ease*bin, df, skewness)
  
  if(!is.null(attr(df, 'adapt'))) {
    adapted <- attr(df, 'adapt')
    adapted$ease <- df$ease
    adapted$magnitude <- df$magnitude
    adapted$bin <- df$bin
    for(colName in colnames(adapted)) {
      if(colName %in% c('ease', 'reps', 'magnitude', 'bin')) next
      form <- as.formula(paste0('`', colName, '`~reps*ease*magnitude*bin'))
      attr(df, paste0(colName, 'OverBins')) <- aggregate(form, adapted, mean)
    }
  }
  df
}

getDataPpBPIC <- function(modelName, dataName, do.plot=FALSE, BPIConly=FALSE, n.post=100, nCores=30) {
  model <- setupModel(modelName)  # calls load_model(), which loads transform.dmc() and transform2.dmc()
  # debug(random.dmc)
  #  debug(WArnd)
  
  dat <- loadData('expMA')[['dat']]
  fn <- paste0('model-', modelName, '_data-', dataName)
  
  # Load, generate posterior preds -------------------------------------------
  samples <- loadSamples(fn, samplesDir)
  #  samples <- h.samples.dmc(samples=samples, add=TRUE, nmc=0, remove=1:10)
  
  # we need to add VVchoiceIdx to cvs to samples so we can simulate random.dmc...
  tmp <- loadData(dataName, removeBlock=NULL)
  data <- data.model.dmc(tmp[['data']], model)
  dat <- tmp[['dat']]
  # dat$ease <- round(dat$ease,3)
  # dat$magnitude <- round(dat$magnitude,3)
  dat$stimulus_type <- NA
  dat$stimulus_type[dat$ease==0.55 & dat$magnitude==1] <- 'easy_low'
  dat$stimulus_type[dat$ease==0.55 & dat$magnitude==1.3] <- 'easy_high'
  dat$stimulus_type[dat$ease==0.4 & dat$magnitude==1] <- 'hard_low'
  dat$stimulus_type[dat$ease==0.4 & dat$magnitude==1.3] <- 'hard_high'
  # 
  # # append cvs by sub as well
  # for(sub in unique(dat$sub)) {
  #   d <- prepareForFitting(dat[dat$sub==sub,])
  #   attr(data[[sub]], 'cvs') <- d$outcomes
  #   attr(attr(data[[sub]], 'cvs'), 'VVchoiceIdx') <- d$VVchoiceIdx    # Trick for random.dmc
  #   attr(data[[sub]], 'VVchoiceIdx') <- d$VVchoiceIdx
  #   attr(data[[sub]], 'startingValues') <- d$values
  #   attr(data[[sub]], 'trialsToIgnore') <- d$trialsToIgnore
  # }
  #   samples <- lapply(1:length(samples), function(x) {
  #     samples[[x]]$data <- data[[x]];
  #     samples[[x]]$theta[,'wS',] <- 0  # fix to 0, what does the result look like?
  #     return(samples[[x]])
  #   })
  data <- lapply(samples, function(x) x$data)
  if(do.plot) plot.dmc(samples, hyper=TRUE, density=TRUE, layout=c(4,4))
  if(!BPIConly) {
    pp = h.post.predict.dmc(samples = samples, adapt=TRUE,save.simulation = TRUE, n.post=n.post, cores=nCores)
    ppNoSim <- h.pp.summary(pp, samples=samples)
    
    #### Append stimulus set info to data & model --------
    nBins <- 10
    pp2 <- lapply(1:length(pp), addStimSetInfo, input=pp, orig_dat=dat, addColumns=c('stimulus_type', 'ease', 'magnitude'))
    data2 <- lapply(1:length(data), addStimSetInfo, input=data, orig_dat=dat, addColumns=c('stimulus_type', 'ease', 'magnitude'))
    if(!sfIsRunning()) sfInit(parallel=TRUE, cpus =30); sfLibrary(moments)
    pp3 <- sfLapply(pp2, calculateByBin)
    data3 <- lapply(data2, calculateByBin)
    
    bpics <- h.IC.dmc(samples)
    return(list('pp3'=pp3, 'data3'=data3, 'BPIC'=bpics))
  } else {
    return(list('BPIC'=h.IC.dmc(samples)))
  }
}


#  load -------------------------------------------------------------------
modelName = 'arw-RL-WA'
dataName <- 'expMA'
#debugonce(getDataPpBPIC)
tmp <- getDataPpBPIC(modelName, dataName, n.post=100)
BPIC <- tmp$BPIC

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






# With SEs ----------------------------------------------------------------

data3 <- tmp[['data3']]
pp3 <- tmp[['pp3']]

q10RTsTargetBySS <- list(getDescriptives(data3, dep.var='RT.10.', attr.name='qRTsTargetByStimType', id.var1='~bin*stimulus_type', id.var2=NULL),
                         getDescriptives(pp3, dep.var='RT.10.', attr.name='qRTsTargetByStimType', id.var1='~reps*bin*stimulus_type', id.var2=NULL))
q50RTsTargetBySS <- list(getDescriptives(data3, dep.var='RT.50.', attr.name='qRTsTargetByStimType', id.var1='~bin*stimulus_type', id.var2=NULL),
                         getDescriptives(pp3, dep.var='RT.50.', attr.name='qRTsTargetByStimType', id.var1='~reps*bin*stimulus_type', id.var2=NULL))
q90RTsTargetBySS <- list(getDescriptives(data3, dep.var='RT.90.', attr.name='qRTsTargetByStimType', id.var1='~bin*stimulus_type', id.var2=NULL),
                         getDescriptives(pp3, dep.var='RT.90.', attr.name='qRTsTargetByStimType', id.var1='~reps*bin*stimulus_type', id.var2=NULL))

q10RTsDBySS <- list(getDescriptives(data3, dep.var='RT.10.', attr.name='qRTsDistractorsByStimType', id.var1='~bin*stimulus_type', id.var2=NULL),
                    getDescriptives(pp3, dep.var='RT.10.', attr.name='qRTsDistractorsByStimType', id.var1='~reps*bin*stimulus_type', id.var2=NULL))
q50RTsDBySS <- list(getDescriptives(data3, dep.var='RT.50.', attr.name='qRTsDistractorsByStimType', id.var1='~bin*stimulus_type', id.var2=NULL),
                    getDescriptives(pp3, dep.var='RT.50.', attr.name='qRTsDistractorsByStimType', id.var1='~reps*bin*stimulus_type', id.var2=NULL))
q90RTsDBySS <- list(getDescriptives(data3, dep.var='RT.90.', attr.name='qRTsDistractorsByStimType', id.var1='~bin*stimulus_type', id.var2=NULL),
                    getDescriptives(pp3, dep.var='RT.90.', attr.name='qRTsDistractorsByStimType', id.var1='~reps*bin*stimulus_type', id.var2=NULL))

q10RTsD1BySS <- list(getDescriptives(data3, dep.var='RT.10.', attr.name='qRTsDistrator1ByStimType', id.var1='~bin*stimulus_type', id.var2=NULL),
                     getDescriptives(pp3, dep.var='RT.10.', attr.name='qRTsDistrator1ByStimType', id.var1='~reps*bin*stimulus_type', id.var2=NULL))
q50RTsD1BySS <- list(getDescriptives(data3, dep.var='RT.50.', attr.name='qRTsDistrator1ByStimType', id.var1='~bin*stimulus_type', id.var2=NULL),
                     getDescriptives(pp3, dep.var='RT.50.', attr.name='qRTsDistrator1ByStimType', id.var1='~reps*bin*stimulus_type', id.var2=NULL))
q90RTsD1BySS <- list(getDescriptives(data3, dep.var='RT.90.', attr.name='qRTsDistrator1ByStimType', id.var1='~bin*stimulus_type', id.var2=NULL),
                     getDescriptives(pp3, dep.var='RT.90.', attr.name='qRTsDistrator1ByStimType', id.var1='~reps*bin*stimulus_type', id.var2=NULL))

q10RTsD2BySS <- list(getDescriptives(data3, dep.var='RT.10.', attr.name='qRTsDistrator2ByStimType', id.var1='~bin*stimulus_type', id.var2=NULL),
                     getDescriptives(pp3, dep.var='RT.10.', attr.name='qRTsDistrator2ByStimType', id.var1='~reps*bin*stimulus_type', id.var2=NULL))
q50RTsD2BySS <- list(getDescriptives(data3, dep.var='RT.50.', attr.name='qRTsDistrator2ByStimType', id.var1='~bin*stimulus_type', id.var2=NULL),
                     getDescriptives(pp3, dep.var='RT.50.', attr.name='qRTsDistrator2ByStimType', id.var1='~reps*bin*stimulus_type', id.var2=NULL))
q90RTsD2BySS <- list(getDescriptives(data3, dep.var='RT.90.', attr.name='qRTsDistrator2ByStimType', id.var1='~bin*stimulus_type', id.var2=NULL),
                     getDescriptives(pp3, dep.var='RT.90.', attr.name='qRTsDistrator2ByStimType', id.var1='~reps*bin*stimulus_type', id.var2=NULL))

# q10RTsTargetBySSE <- list(getDescriptives(data3, dep.var='RT.10.', attr.name='qRTsErrorBySS', id.var1='~bin*stimulus_set_ABCD', id.var2=NULL),
#                     getDescriptives(pp3, dep.var='RT.10.', attr.name='qRTsErrorBySS', id.var1='~reps*bin*stimulus_set_ABCD', id.var2=NULL))
# q50RTsTargetBySSE <- list(getDescriptives(data3, dep.var='RT.50.', attr.name='qRTsErrorBySS', id.var1='~bin*stimulus_set_ABCD', id.var2=NULL),
#                     getDescriptives(pp3, dep.var='RT.50.', attr.name='qRTsErrorBySS', id.var1='~reps*bin*stimulus_set_ABCD', id.var2=NULL))
# q90RTsTargetBySSE <- list(getDescriptives(data3, dep.var='RT.90.', attr.name='qRTsErrorBySS', id.var1='~bin*stimulus_set_ABCD', id.var2=NULL),
#                     getDescriptives(pp3, dep.var='RT.90.', attr.name='qRTsErrorBySS', id.var1='~reps*bin*stimulus_set_ABCD', id.var2=NULL))
meanAccBySS <- list(getDescriptives(data3, dep.var='acc', attr.name='accByStimType', id.var1='~bin*stimulus_type', id.var2=NULL),
                    getDescriptives(pp3, dep.var='acc', attr.name='accByStimType', id.var1='~reps*bin*stimulus_type', id.var2=NULL))


SEq10RTsTargetBySS <- list(getDescriptivesSE(data3, dep.var='RT.10.', attr.name='qRTsTargetByStimType', id.var1='~bin*stimulus_type', id.var2=NULL))
SEq50RTsTargetBySS <- list(getDescriptivesSE(data3, dep.var='RT.50.', attr.name='qRTsTargetByStimType', id.var1='~bin*stimulus_type', id.var2=NULL))
SEq90RTsTargetBySS <- list(getDescriptivesSE(data3, dep.var='RT.90.', attr.name='qRTsTargetByStimType', id.var1='~bin*stimulus_type', id.var2=NULL))
SEq10RTsDBySS <- list(getDescriptivesSE(data3, dep.var='RT.10.', attr.name='qRTsDistractorsByStimType', id.var1='~bin*stimulus_type', id.var2=NULL))
SEq50RTsDBySS <- list(getDescriptivesSE(data3, dep.var='RT.50.', attr.name='qRTsDistractorsByStimType', id.var1='~bin*stimulus_type', id.var2=NULL))
SEq90RTsDBySS <- list(getDescriptivesSE(data3, dep.var='RT.90.', attr.name='qRTsDistractorsByStimType', id.var1='~bin*stimulus_type', id.var2=NULL))
SEmeanAccBySS <- list(getDescriptivesSE(data3, dep.var='acc', attr.name='accByStimType', id.var1='~bin*stimulus_type', id.var2=NULL))



# Plot posterior predictives
if(savePlot) pdf(file=paste0('./figures/', modelName, 'SE.pdf'), width=7, height=7/4*3)
par(oma=c(3,3,1,0), mar=c(0, 1, 1, 0) + 0.1, mfcol=c(3,4), mgp=c(2.75,.75,0), las=1, bty='l')
i <- 0
corrRTylim <- c(.6, 1.6)
errRTylim <- c(.6, 1.6)
data.cex=1.5
stimulus_type = unique(meanAccBySS[[1]]$stimulus_type)[1]
for(stimulus_type in unique(meanAccBySS[[1]]$stimulus_type)) {
  i <- i+1
  idxD = meanAccBySS[[1]]$stimulus_type == stimulus_type
  idxM = meanAccBySS[[2]]$stimulus_type == stimulus_type
  
  plotDataPPBins(data=meanAccBySS[[1]][idxD,], pp=meanAccBySS[[2]][idxM,],
                 xaxt='n', draw.legend = i==1, data.cex = data.cex,
                 dep.var='acc', ylab='', xlab = '', yaxt='n',
                 legend.pos='topleft', ylim=c(0.25, 0.8), hline.by=0.1)
  arrows(x0=meanAccBySS[[1]][idxD,'bin'], x1=meanAccBySS[[1]][idxD,'bin'],
         y0=meanAccBySS[[1]][idxD,'acc']-SEmeanAccBySS[[1]][idxD,'acc'],
         y1=meanAccBySS[[1]][idxD,'acc']+SEmeanAccBySS[[1]][idxD,'acc'], code=3, angle=90, lwd=1.5, length = 0.05)
  axis(1, at=seq(2, 10, 2), labels=rep(NA, 5), lwd=2)
  if(i == 1) {
    mtext('Accuracy', side=2, cex=.66, line=3, las=0, font=1)
    axis(2, at=seq(.3, .8, .1), lwd=1.5)
  } else {
    axis(2, at=seq(.3, .8, .1), labels=rep(NA, 5), lwd=1.5)
  }
  # title(ifelse(stimulus_type == 'easy_high', 'Easy, High',
  #              ifelse(stimulus_type)))
  if(i == 1) title('Easy, high')
  if(i == 2) title('Easy, low')
  if(i == 3) title('Hard, high')
  if(i == 4) title('Hard, low')
  # 
  
  # Target RTs
  plotDataPPBins(data=q10RTsTargetBySS[[1]][idxD,], pp=q10RTsTargetBySS[[2]][idxM,], dep.var='RT.10.',
                 ylim=corrRTylim, xaxt='n', ylab='', yaxt='n', draw.legend = FALSE, data.cex = data.cex)
  plotDataPPBins(data=q50RTsTargetBySS[[1]][idxD,], pp=q50RTsTargetBySS[[2]][idxM,], dep.var='RT.50.', plot.new = FALSE, draw.legend=FALSE, data.cex = data.cex)
  plotDataPPBins(data=q90RTsTargetBySS[[1]][idxD,], pp=q90RTsTargetBySS[[2]][idxM,], dep.var='RT.90.', plot.new = FALSE, draw.legend=FALSE, data.cex = data.cex)
  arrows(x0=q10RTsTargetBySS[[1]][idxD,'bin'], x1=q10RTsTargetBySS[[1]][idxD,'bin'],
         y0=q10RTsTargetBySS[[1]][idxD,'RT.10.']-SEq10RTsTargetBySS[[1]][idxD,'RT.10.'],
         y1=q10RTsTargetBySS[[1]][idxD,'RT.10.']+SEq10RTsTargetBySS[[1]][idxD,'RT.10.'], code=3, angle=90, lwd=1.5, length = 0.05)
  arrows(x0=q50RTsTargetBySS[[1]][idxD,'bin'], x1=q50RTsTargetBySS[[1]][idxD,'bin'],
         y0=q50RTsTargetBySS[[1]][idxD,'RT.50.']-SEq50RTsTargetBySS[[1]][idxD,'RT.50.'],
         y1=q50RTsTargetBySS[[1]][idxD,'RT.50.']+SEq50RTsTargetBySS[[1]][idxD,'RT.50.'], code=3, angle=90, lwd=1.5, length = 0.05)
  arrows(x0=q90RTsTargetBySS[[1]][idxD,'bin'], x1=q90RTsTargetBySS[[1]][idxD,'bin'],
         y0=q90RTsTargetBySS[[1]][idxD,'RT.90.']-SEq90RTsTargetBySS[[1]][idxD,'RT.90.'],
         y1=q90RTsTargetBySS[[1]][idxD,'RT.90.']+SEq90RTsTargetBySS[[1]][idxD,'RT.90.'], code=3, angle=90, lwd=1.5, length = 0.05)
  axis(1, at=seq(2, 10, 2), labels=NA, lwd=1.5)
  if(i == 1) {
    mtext('Correct RTs (s)', side=2, cex=.66, line=3, las=0, font=1)
    axis(2, seq(.4, 2.5, .2), lwd=1.5)
  } else {
    axis(2, seq(.4, 2.5, .2), labels=NA, lwd=1.5)
  }
  
  ##
  # plotDataPPBins(data=q10RTsD1BySS[[1]][idxD,], pp=q10RTsD1BySS[[2]][idxM,], dep.var='RT.10.', ylim=errRTylim, xaxt='n', ylab='', yaxt='n', draw.legend = FALSE, data.cex = data.cex)
  # plotDataPPBins(data=q50RTsD1BySS[[1]][idxD,], pp=q50RTsD1BySS[[2]][idxM,], dep.var='RT.50.', plot.new = FALSE, draw.legend=FALSE, data.cex = data.cex)
  # plotDataPPBins(data=q90RTsD1BySS[[1]][idxD,], pp=q90RTsD1BySS[[2]][idxM,], dep.var='RT.90.', plot.new = FALSE, draw.legend=FALSE, data.cex = data.cex)
  plotDataPPBins(data=q10RTsDBySS[[1]][idxD,], pp=q10RTsDBySS[[2]][idxM,], dep.var='RT.10.', ylim=errRTylim, xaxt='n', ylab='', yaxt='n', draw.legend = FALSE, data.cex = data.cex)
  plotDataPPBins(data=q50RTsDBySS[[1]][idxD,], pp=q50RTsDBySS[[2]][idxM,], dep.var='RT.50.', plot.new = FALSE, draw.legend=FALSE, data.cex = data.cex)
  plotDataPPBins(data=q90RTsDBySS[[1]][idxD,], pp=q90RTsDBySS[[2]][idxM,], dep.var='RT.90.', plot.new = FALSE, draw.legend=FALSE, data.cex = data.cex)
  arrows(x0=q10RTsDBySS[[1]][idxD,'bin'], x1=q10RTsDBySS[[1]][idxD,'bin'],
         y0=q10RTsDBySS[[1]][idxD,'RT.10.']-SEq10RTsDBySS[[1]][idxD,'RT.10.'],
         y1=q10RTsDBySS[[1]][idxD,'RT.10.']+SEq10RTsDBySS[[1]][idxD,'RT.10.'], code=3, angle=90, lwd=1.5, length = 0.05)
  arrows(x0=q50RTsDBySS[[1]][idxD,'bin'], x1=q50RTsDBySS[[1]][idxD,'bin'],
         y0=q50RTsDBySS[[1]][idxD,'RT.50.']-SEq50RTsDBySS[[1]][idxD,'RT.50.'],
         y1=q50RTsDBySS[[1]][idxD,'RT.50.']+SEq50RTsDBySS[[1]][idxD,'RT.50.'], code=3, angle=90, lwd=1.5, length = 0.05)
  arrows(x0=q90RTsDBySS[[1]][idxD,'bin'], x1=q90RTsDBySS[[1]][idxD,'bin'],
         y0=q90RTsDBySS[[1]][idxD,'RT.90.']-SEq90RTsDBySS[[1]][idxD,'RT.90.'],
         y1=q90RTsDBySS[[1]][idxD,'RT.90.']+SEq90RTsDBySS[[1]][idxD,'RT.90.'], code=3, angle=90, lwd=1.5, length = 0.05)
  if(i == 1) {
    mtext('Error RTs (s)', side=2, cex=.66, line=3, las=0, font=1)
    axis(2, seq(.4, 2.5, .2), lwd=1.5)
  } else {
    axis(2, seq(.4, 2.5, .2), labels=NA, lwd=1.5)
  }
  axis(1, at=seq(2, 10, 2), lwd=1.5)
  mtext('Trial bin', side=1, cex=.66, line=2)
}
if(savePlot) dev.off()
