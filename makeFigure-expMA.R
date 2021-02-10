rm(list=ls())
library(snowfall)
source ("dmc/dmc.R")
source('utils.R')
source('models.R')
samplesDir <- 'samples'
savePlot <- FALSE

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
#BPIC <- tmp$BPIC
#qRTs <- getqRTs(tmp[['data3']], tmp[['pp3']])
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

# q10RTsBySSE <- list(getDescriptives(data3, dep.var='RT.10.', attr.name='qRTsErrorBySS', id.var1='~bin*stimulus_set_ABCD', id.var2=NULL),
#                     getDescriptives(pp3, dep.var='RT.10.', attr.name='qRTsErrorBySS', id.var1='~reps*bin*stimulus_set_ABCD', id.var2=NULL))
# q50RTsBySSE <- list(getDescriptives(data3, dep.var='RT.50.', attr.name='qRTsErrorBySS', id.var1='~bin*stimulus_set_ABCD', id.var2=NULL),
#                     getDescriptives(pp3, dep.var='RT.50.', attr.name='qRTsErrorBySS', id.var1='~reps*bin*stimulus_set_ABCD', id.var2=NULL))
# q90RTsBySSE <- list(getDescriptives(data3, dep.var='RT.90.', attr.name='qRTsErrorBySS', id.var1='~bin*stimulus_set_ABCD', id.var2=NULL),
#                     getDescriptives(pp3, dep.var='RT.90.', attr.name='qRTsErrorBySS', id.var1='~reps*bin*stimulus_set_ABCD', id.var2=NULL))
meanAccBySS <- list(getDescriptives(data3, dep.var='acc', attr.name='accByStimType', id.var1='~bin*stimulus_type', id.var2=NULL),
                    getDescriptives(pp3, dep.var='acc', attr.name='accByStimType', id.var1='~reps*bin*stimulus_type', id.var2=NULL))

# Plot posterior predictives
if(savePlot) pdf(file=paste0('./figures/', modelName, '.pdf'), width=7, height=7/4*3)
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


# #
fn <- paste0('model-', modelName, '_data-', dataName)

# Load, generate posterior preds -------------------------------------------
samples <- loadSamples(fn, samplesDir)
plot.dmc(samples, hyper=TRUE, location=TRUE)

summ <- summary.dmc(samples)
summdf <- data.frame(do.call(rbind, lapply(summ, function(x) x$quantiles[,3])))
summdf$alpha <- pnorm(summdf[,4])
apply(summdf, 2, mean)
apply(summdf, 2, sd)


# Q-values, drift rates ---------------------------------------------------
get.color <- function(ease) {
  if(ease == "1") return (1)
  if(ease == "1.3") return(2)
  if(ease == "0.55") return(1)  # Easy
  if(ease == "0.4") return(2)   # Hard
}


draw.polygon <- function(pp, dep.var, xaxis='bin', colorM='blue', plot.model.points=FALSE) {
  lowerQ <- aggregate(as.formula(paste0(dep.var, '~bin')), pp, quantile, .025)
  upperQ <- aggregate(as.formula(paste0(dep.var, '~bin')), pp, quantile, .975)
  xs <- c(lowerQ[,xaxis], rev(lowerQ[,xaxis]))
  ys <- c(lowerQ[,dep.var], rev(upperQ[,dep.var]))
  polygon(xs, ys, col=rgb(col2rgb(colorM)[1]/255, col2rgb(colorM)[2]/255, col2rgb(colorM)[3]/255, alpha=.3), lty = NULL, border=NA)
  if(plot.model.points) points(pp[,xaxis], pp[,dep.var], pch=20, col=colorM, cex=.01)
}



# Low vs High magnitude rows ----------------------------------------------
savePlot <- TRUE
if(savePlot) pdf(file=paste0('./figures/MA-Q-values.pdf'), width=7, height=7/4*2)
par(mfrow=c(2,4), las=1, bty='l', oma=c(1,2,1,0), mar=c(2, 4, 2, 0.5) + 0.1, mgp=c(2.25,.75,0))

# By Difficulty -----------------------------------------------------------
allQ1OverTimeM <- do.call(rbind, (lapply(1:length(pp3), function(x) {tmp <- attr(pp3[[x]], 'SR.r1OverBins'); tmp$s <- x; tmp})))
allQ2OverTimeM <- do.call(rbind, (lapply(1:length(pp3), function(x) {tmp <- attr(pp3[[x]], 'SR.r2OverBins'); tmp$s <- x; tmp})))
allQ3OverTimeM <- do.call(rbind, (lapply(1:length(pp3), function(x) {tmp <- attr(pp3[[x]], 'SR.r3OverBins'); tmp$s <- x; tmp})))
# drifts
all12OverTimeM <- do.call(rbind, (lapply(1:length(pp3), function(x) {tmp <- attr(pp3[[x]], '1-2OverBins'); tmp$s <- x; tmp})))
all13OverTimeM <- do.call(rbind, (lapply(1:length(pp3), function(x) {tmp <- attr(pp3[[x]], '1-3OverBins'); tmp$s <- x; tmp})))
all21OverTimeM <- do.call(rbind, (lapply(1:length(pp3), function(x) {tmp <- attr(pp3[[x]], '2-1OverBins'); tmp$s <- x; tmp})))
all23OverTimeM <- do.call(rbind, (lapply(1:length(pp3), function(x) {tmp <- attr(pp3[[x]], '2-3OverBins'); tmp$s <- x; tmp})))
all31OverTimeM <- do.call(rbind, (lapply(1:length(pp3), function(x) {tmp <- attr(pp3[[x]], '3-1OverBins'); tmp$s <- x; tmp})))
all32OverTimeM <- do.call(rbind, (lapply(1:length(pp3), function(x) {tmp <- attr(pp3[[x]], '3-2OverBins'); tmp$s <- x; tmp})))
# change colnames because annoying
all12OverTimeM$Q12 <- all12OverTimeM$`1-2`
all13OverTimeM$Q13 <- all13OverTimeM$`1-3`
all21OverTimeM$Q21 <- all21OverTimeM$`2-1`
all23OverTimeM$Q23 <- all23OverTimeM$`2-3`
all31OverTimeM$Q31 <- all31OverTimeM$`3-1`
all32OverTimeM$Q32 <- all32OverTimeM$`3-2`

for(magnitude in c(1, 1.3)) {
  if(magnitude == 1.3) par(mar=c(3, 4, 0, 0.5) + 0.1) else par(mar=c(1, 4, 2, 0.5) + 0.1)
  meanQ1OverTimeM <- aggregate(SR.r1~reps*bin*ease, allQ1OverTimeM[allQ1OverTimeM$magnitude==magnitude,], mean)
  meanQ2OverTimeM <- aggregate(SR.r2~reps*bin*ease, allQ2OverTimeM[allQ2OverTimeM$magnitude==magnitude,], mean)
  meanQ3OverTimeM <- aggregate(SR.r3~reps*bin*ease, allQ3OverTimeM[allQ3OverTimeM$magnitude==magnitude,], mean)
  
  # differences
  deltaQ <- allQ1OverTimeM[allQ1OverTimeM$magnitude==magnitude,]
  deltaQ$SR.r2 <- allQ2OverTimeM[allQ1OverTimeM$magnitude==magnitude,]$SR.r2
  deltaQ$SR.r3 <- allQ3OverTimeM[allQ1OverTimeM$magnitude==magnitude,]$SR.r3
  deltaQ$delta1_2 <- deltaQ$SR.r1 - deltaQ$SR.r2
  deltaQ$delta1_3 <- deltaQ$SR.r1 - deltaQ$SR.r3
  deltaQ$delta2_1 <- deltaQ$SR.r2 - deltaQ$SR.r1
  deltaQ$delta2_3 <- deltaQ$SR.r2 - deltaQ$SR.r3
  deltaQ$delta3_1 <- deltaQ$SR.r3 - deltaQ$SR.r2
  deltaQ$delta3_2 <- deltaQ$SR.r3 - deltaQ$SR.r2
  
  # sums
  deltaQ$sumQ1_2 <- deltaQ$SR.r1 + deltaQ$SR.r2
  deltaQ$sumQ1_3 <- deltaQ$SR.r1 + deltaQ$SR.r3
  deltaQ$sumQ2_3 <- deltaQ$SR.r2 + deltaQ$SR.r3
  
  meanDeltaQ12OverTime <- aggregate(delta1_2~reps*bin*ease, deltaQ, mean)
  meanDeltaQ13OverTime <- aggregate(delta1_3~reps*bin*ease, deltaQ, mean)
  meanDeltaQ21OverTime <- aggregate(delta2_1~reps*bin*ease, deltaQ, mean)
  meanDeltaQ23OverTime <- aggregate(delta2_3~reps*bin*ease, deltaQ, mean)
  meanDeltaQ31OverTime <- aggregate(delta3_1~reps*bin*ease, deltaQ, mean)
  meanDeltaQ32OverTime <- aggregate(delta3_2~reps*bin*ease, deltaQ, mean)
  meanSumQ12OverTime <- aggregate(sumQ1_2~reps*bin*ease, deltaQ, mean)
  meanSumQ23OverTime <- aggregate(sumQ2_3~reps*bin*ease, deltaQ, mean)
  meanSumQ13OverTime <- aggregate(sumQ1_3~reps*bin*ease, deltaQ, mean)
  
  # drift rates
  meanV12OverTimeM <- aggregate(Q12~reps*bin*ease, all12OverTimeM[all12OverTimeM$magnitude==magnitude,], mean)
  meanV13OverTimeM <- aggregate(Q13~reps*bin*ease, all13OverTimeM[all13OverTimeM$magnitude==magnitude,], mean)
  meanV21OverTimeM <- aggregate(Q21~reps*bin*ease, all21OverTimeM[all21OverTimeM$magnitude==magnitude,], mean)
  meanV23OverTimeM <- aggregate(Q23~reps*bin*ease, all23OverTimeM[all23OverTimeM$magnitude==magnitude,], mean)
  meanV31OverTimeM <- aggregate(Q31~reps*bin*ease, all31OverTimeM[all31OverTimeM$magnitude==magnitude,], mean)
  meanV32OverTimeM <- aggregate(Q32~reps*bin*ease, all32OverTimeM[all32OverTimeM$magnitude==magnitude,], mean)
  
  # Plot --------------------------------------------------------------------
  #if(savePlot) pdf(file='./figures/q-values.pdf', width=7, height=2.5)
  # Q-values
  plot(0,0, type='n', xlim=range(meanQ1OverTimeM$bin)+c(-.5, .5), ylim=c(0, .7),
       xlab=ifelse(magnitude==1, '', 'Trial bin'), 
       ylab='Q-values', xaxt='n',
       main=ifelse(magnitude==1, 'A. Q-values', ''))
  if(magnitude == 1) {
    axis(at=seq(2, 10, 2), side=1, labels=rep(NA, 5))# ifelse(magnitude==1, rep(c(''), 5), seq(2, 10, 2)))
  } else {
    axis(at=seq(2, 10, 2), side=1)# ifelse(magnitude==1, rep(c(''), 5), seq(2, 10, 2)))
  }
#  axis(at=seq(2, 10, 2), side=1, labels=ifelse(magnitude==1, rep(c(''), 5), seq(2, 10, 2)))
  #mtext(paste0('Magnitude = ', magnitude), 2, las=0, line=4, cex=0.83, font=2)
  mtext(ifelse(magnitude==1, 'Low magnitude', 'High magnitude'), 2, las=0, line=4, cex=0.83, font=2)
  abline(h=seq(0, 1, .1), col='grey')
  abline(v=seq(0, 10, 2), col='grey')
  for(ease in unique(meanQ1OverTimeM$ease)) {
    draw.polygon(meanQ1OverTimeM[meanQ1OverTimeM$ease==ease,], dep.var='SR.r1', colorM=get.color(ease)) #ifelse(ease=="0.6", 1, ifelse(ease=="0.4", 2, ifelse(ease=="0.3", 3, 4))))
    draw.polygon(meanQ2OverTimeM[meanQ2OverTimeM$ease==ease,], dep.var='SR.r2', colorM=get.color(ease))
    draw.polygon(meanQ3OverTimeM[meanQ3OverTimeM$ease==ease,], dep.var='SR.r3', colorM=get.color(ease))
  }
  legend('topleft', legend=c('Easy', 'Hard'), bg='white', col=1:2, pch=15, title='Difficulty')
  text(8, ifelse(magnitude==1, 0.3, 0.4), labels = expression('Q'[T]))
  text(8, ifelse(magnitude==1, 0.12, 0.14), labels = expression('Q'[D]))
  
  # delta here
  plot(0,0, type='n', xlim=range(meanDeltaQ12OverTime$bin)+c(-.5, .5), ylim=c(0, 0.7), xaxt='n',
       xlab=ifelse(magnitude==1, '', 'Trial bin'),  
       ylab=expression(paste(Delta, 'Q-values')), 
       main=ifelse(magnitude==1, expression(bold(paste('B. ', Delta, 'Q-values'))), ''))
  if(magnitude == 1) {
    axis(at=seq(2, 10, 2), side=1, labels=rep(NA, 5))# ifelse(magnitude==1, rep(c(''), 5), seq(2, 10, 2)))
  } else {
    axis(at=seq(2, 10, 2), side=1)# ifelse(magnitude==1, rep(c(''), 5), seq(2, 10, 2)))
  }
  abline(h=seq(-1, 1, .1), col='grey')
  abline(v=seq(0, 10, 2), col='grey')
  for(ease in unique(meanDeltaQ12OverTime$ease)) {
    draw.polygon(meanDeltaQ12OverTime[meanDeltaQ12OverTime$ease==ease,], dep.var='delta1_2', colorM=get.color(ease))
    # draw.polygon(meanDeltaQ13OverTime[meanDeltaQ13OverTime$ease==ease,], dep.var='delta1_3', colorM=get.color(ease))
    draw.polygon(meanDeltaQ21OverTime[meanDeltaQ21OverTime$ease==ease,], dep.var='delta2_1', colorM=get.color(ease))
    # draw.polygon(meanDeltaQ23OverTime[meanDeltaQ23OverTime$ease==ease,], dep.var='delta2_3', colorM=get.color(ease))
    draw.polygon(meanDeltaQ31OverTime[meanDeltaQ31OverTime$ease==ease,], dep.var='delta3_1', colorM=get.color(ease))
    draw.polygon(meanDeltaQ32OverTime[meanDeltaQ32OverTime$ease==ease,], dep.var='delta3_2', colorM=get.color(ease))
  }
  text(8, ifelse(magnitude==1, 0.50, 0.54), labels = expression(paste(Delta, 'Q'[T-D])))
  text(8, 0.05, labels = expression(paste(Delta, 'Q'[D-D])))
  
  # sum here
  plot(0,0, type='n', xlim=range(meanDeltaQ12OverTime$bin)+c(-.5, .5), ylim=c(0, .7), 
       xlab=ifelse(magnitude==1, '', 'Trial bin'), 
       xaxt='n',
       ylab=expression(paste(Sigma, 'Q-values')), main=ifelse(magnitude==1, expression(bold(paste('C. ', Sigma, 'Q-values'))), ''))
  if(magnitude == 1) {
    axis(at=seq(2, 10, 2), side=1, labels=rep(NA, 5))# ifelse(magnitude==1, rep(c(''), 5), seq(2, 10, 2)))
  } else {
    axis(at=seq(2, 10, 2), side=1)# ifelse(magnitude==1, rep(c(''), 5), seq(2, 10, 2)))
  }
  abline(h=seq(0, 1, .1), col='grey')
  abline(v=seq(0, 10, 2), col='grey')
  for(ease in unique(meanDeltaQ12OverTime$ease)) {
    draw.polygon(meanSumQ12OverTime[meanSumQ12OverTime$ease==ease,], dep.var='sumQ1_2', colorM=get.color(ease))
    draw.polygon(meanSumQ13OverTime[meanSumQ13OverTime$ease==ease,], dep.var='sumQ1_3', colorM=get.color(ease))
    draw.polygon(meanSumQ23OverTime[meanSumQ23OverTime$ease==ease,], dep.var='sumQ2_3', colorM=get.color(ease))
  }
  text(9, ifelse(magnitude==1, 0.375, 0.475), labels = expression(paste(Sigma, 'Q'[T+D])))
  text(9, ifelse(magnitude==1, 0.2, 0.25), labels = expression(paste(Sigma, 'Q'[D+D])))
  
  plot(0,0, type='n', xlim=range(meanV12OverTimeM$bin)+c(-.5, .5), ylim=c(0.4, 2.1), 
       xlab=ifelse(magnitude==1, '', 'Trial bin'),xaxt='n',
       ylab='Drift rates', main=ifelse(magnitude==1, 'D. Drift rates', ''))
  if(magnitude == 1) {
    axis(at=seq(2, 10, 2), side=1, labels=rep(NA, 5))# ifelse(magnitude==1, rep(c(''), 5), seq(2, 10, 2)))
  } else {
    axis(at=seq(2, 10, 2), side=1)# ifelse(magnitude==1, rep(c(''), 5), seq(2, 10, 2)))
  }
  abline(h=seq(0, 5, .25), col='grey')
  abline(v=seq(0, 10, 2), col='grey')
  for(ease in unique(meanV12OverTimeM$ease)) {
    draw.polygon(meanV12OverTimeM[meanV12OverTimeM$ease==ease,], dep.var='Q12', colorM=get.color(ease))
    draw.polygon(meanV13OverTimeM[meanV13OverTimeM$ease==ease,], dep.var='Q13', colorM=get.color(ease))
    draw.polygon(meanV21OverTimeM[meanV21OverTimeM$ease==ease,], dep.var='Q21', colorM=get.color(ease))
    draw.polygon(meanV23OverTimeM[meanV23OverTimeM$ease==ease,], dep.var='Q23', colorM=get.color(ease))
    draw.polygon(meanV31OverTimeM[meanV31OverTimeM$ease==ease,], dep.var='Q31', colorM=get.color(ease))
    draw.polygon(meanV32OverTimeM[meanV32OverTimeM$ease==ease,], dep.var='Q32', colorM=get.color(ease))
  }
  text(6, ifelse(magnitude==1, 1.9, 2), labels = expression(paste('v'[T-D])))
  text(8, ifelse(magnitude==1, 1.3, 1.3), labels = expression(paste('v'[D-D])))
  text(6, ifelse(magnitude==1, 0.5, 0.5), labels = expression(paste('v'[D-T])))
}
if(savePlot) dev.off()
# drift rates
# plot(0,0, type='n', xlim=range(meanV1OverTimeM$bin)+c(-.5, .5), ylim=c(1.0, 4.5), xlab='Trial bin', ylab='Drift rates', main='D. Drift rates')
# abline(h=seq(0, 5, .25), col='grey')
# abline(v=seq(0, 10, 2), col='grey')
# for(ease in unique(allV1OverTimeM$ease)) {
#   draw.polygon(meanV1OverTimeM[meanV1OverTimeM$ease==ease,], dep.var='mean_v.r1', colorM=get.color(ease))
#   draw.polygon(meanV2OverTimeM[meanV2OverTimeM$ease==ease,], dep.var='mean_v.r2', colorM=get.color(ease))
# }
#if(savePlot) dev.off()



# By magnitude ------------------------------------------------------------
allQ1OverTimeM <- do.call(rbind, (lapply(1:length(pp3), function(x) {tmp <- attr(pp3[[x]], 'SR.r1OverBins'); tmp$s <- x; tmp})))
allQ2OverTimeM <- do.call(rbind, (lapply(1:length(pp3), function(x) {tmp <- attr(pp3[[x]], 'SR.r2OverBins'); tmp$s <- x; tmp})))
allQ3OverTimeM <- do.call(rbind, (lapply(1:length(pp3), function(x) {tmp <- attr(pp3[[x]], 'SR.r3OverBins'); tmp$s <- x; tmp})))
meanQ1OverTimeM <- aggregate(SR.r1~reps*bin*magnitude, allQ1OverTimeM, mean)
meanQ2OverTimeM <- aggregate(SR.r2~reps*bin*magnitude, allQ2OverTimeM, mean)
meanQ3OverTimeM <- aggregate(SR.r3~reps*bin*magnitude, allQ3OverTimeM, mean)

# differences
deltaQ <- allQ1OverTimeM
deltaQ$SR.r2 <- allQ2OverTimeM$SR.r2
deltaQ$SR.r3 <- allQ3OverTimeM$SR.r3
deltaQ$delta1_2 <- deltaQ$SR.r1 - deltaQ$SR.r2
deltaQ$delta1_3 <- deltaQ$SR.r1 - deltaQ$SR.r3
deltaQ$delta2_1 <- deltaQ$SR.r2 - deltaQ$SR.r1
deltaQ$delta2_3 <- deltaQ$SR.r2 - deltaQ$SR.r3
deltaQ$delta3_1 <- deltaQ$SR.r3 - deltaQ$SR.r2
deltaQ$delta3_2 <- deltaQ$SR.r3 - deltaQ$SR.r2

# sums
deltaQ$sumQ1_2 <- deltaQ$SR.r1 + deltaQ$SR.r2
deltaQ$sumQ1_3 <- deltaQ$SR.r1 + deltaQ$SR.r3
deltaQ$sumQ2_3 <- deltaQ$SR.r2 + deltaQ$SR.r3

meanDeltaQ12OverTime <- aggregate(delta1_2~reps*bin*magnitude, deltaQ, mean)
meanDeltaQ13OverTime <- aggregate(delta1_3~reps*bin*magnitude, deltaQ, mean)
meanDeltaQ21OverTime <- aggregate(delta2_1~reps*bin*magnitude, deltaQ, mean)
meanDeltaQ23OverTime <- aggregate(delta2_3~reps*bin*magnitude, deltaQ, mean)
meanDeltaQ31OverTime <- aggregate(delta3_1~reps*bin*magnitude, deltaQ, mean)
meanDeltaQ32OverTime <- aggregate(delta3_2~reps*bin*magnitude, deltaQ, mean)
meanSumQ12OverTime <- aggregate(sumQ1_2~reps*bin*magnitude, deltaQ, mean)
meanSumQ23OverTime <- aggregate(sumQ2_3~reps*bin*magnitude, deltaQ, mean)
meanSumQ13OverTime <- aggregate(sumQ1_3~reps*bin*magnitude, deltaQ, mean)

# drift rates
allV1OverTimeM <- do.call(rbind, (lapply(1:length(pp3), function(x) {tmp <- attr(pp3[[x]], 'mean_v.r1OverBins'); tmp$s <- x; tmp})))
allV2OverTimeM <- do.call(rbind, (lapply(1:length(pp3), function(x) {tmp <- attr(pp3[[x]], 'mean_v.r2OverBins'); tmp$s <- x; tmp})))
allV3OverTimeM <- do.call(rbind, (lapply(1:length(pp3), function(x) {tmp <- attr(pp3[[x]], 'mean_v.r3OverBins'); tmp$s <- x; tmp})))
meanV1OverTimeM <- aggregate(mean_v.r1~reps*bin*magnitude, allV1OverTimeM, mean)
meanV2OverTimeM <- aggregate(mean_v.r2~reps*bin*magnitude, allV2OverTimeM, mean)
meanV3OverTimeM <- aggregate(mean_v.r3~reps*bin*magnitude, allV3OverTimeM, mean)


# Plot --------------------------------------------------------------------
#if(savePlot) pdf(file='./figures/q-values.pdf', width=7, height=2.5)
#par(mfrow=c(1,4), las=1, bty='l', oma=c(0,1,1,0), mar=c(4, 3, 2, 0.5) + 0.1, mgp=c(2.25,.75,0))
# Q-values
plot(0,0, type='n', xlim=range(meanQ1OverTimeM$bin)+c(-.5, .5), ylim=c(0, .85), xlab='Trial bin', ylab='Q-values', main='A. Q-values')
abline(h=seq(0, 1, .1), col='grey')
abline(v=seq(0, 10, 2), col='grey')
for(magnitude in unique(meanQ1OverTimeM$magnitude)) {
  draw.polygon(meanQ1OverTimeM[meanQ1OverTimeM$magnitude==magnitude,], dep.var='SR.r1', colorM=get.color(magnitude)) #ifelse(ease=="0.6", 1, ifelse(ease=="0.4", 2, ifelse(ease=="0.3", 3, 4))))
  draw.polygon(meanQ2OverTimeM[meanQ2OverTimeM$magnitude==magnitude,], dep.var='SR.r2', colorM=get.color(magnitude))
  draw.polygon(meanQ3OverTimeM[meanQ3OverTimeM$magnitude==magnitude,], dep.var='SR.r3', colorM=get.color(magnitude))
}

# delta here
plot(0,0, type='n', xlim=range(meanDeltaQ12OverTime$bin)+c(-.5, .5), ylim=c(-0.6, 0.6), xlab='Trial bin', ylab=expression(paste(Delta, 'Q-values')), main=expression(bold(paste('B. ', Delta, 'Q-values'))))
abline(h=seq(-1, 1, .1), col='grey')
abline(v=seq(0, 10, 2), col='grey')
for(magnitude in unique(meanDeltaQ12OverTime$magnitude)) {
  draw.polygon(meanDeltaQ12OverTime[meanDeltaQ12OverTime$magnitude==magnitude,], dep.var='delta1_2', colorM=get.color(magnitude))
  draw.polygon(meanDeltaQ13OverTime[meanDeltaQ13OverTime$magnitude==magnitude,], dep.var='delta1_3', colorM=get.color(magnitude))
  draw.polygon(meanDeltaQ21OverTime[meanDeltaQ21OverTime$magnitude==magnitude,], dep.var='delta2_1', colorM=get.color(magnitude))
  draw.polygon(meanDeltaQ23OverTime[meanDeltaQ23OverTime$magnitude==magnitude,], dep.var='delta2_3', colorM=get.color(magnitude))
  draw.polygon(meanDeltaQ31OverTime[meanDeltaQ31OverTime$magnitude==magnitude,], dep.var='delta3_1', colorM=get.color(magnitude))
  draw.polygon(meanDeltaQ32OverTime[meanDeltaQ32OverTime$magnitude==magnitude,], dep.var='delta3_2', colorM=get.color(magnitude))
}

# sum here
plot(0,0, type='n', xlim=range(meanDeltaQ12OverTime$bin)+c(-.5, .5), ylim=c(0, .85), xlab='Trial bin', ylab=expression(paste(Sigma, 'Q-values')), main=expression(bold(paste('C. ', Sigma, 'Q-values'))))
abline(h=seq(0, 1, .1), col='grey')
abline(v=seq(0, 10, 2), col='grey')
for(magnitude in unique(meanDeltaQ12OverTime$magnitude)) {
  draw.polygon(meanSumQ12OverTime[meanSumQ12OverTime$magnitude==magnitude,], dep.var='sumQ1_2', colorM=get.color(magnitude))
  draw.polygon(meanSumQ13OverTime[meanSumQ13OverTime$magnitude==magnitude,], dep.var='sumQ1_3', colorM=get.color(magnitude))
  draw.polygon(meanSumQ23OverTime[meanSumQ23OverTime$magnitude==magnitude,], dep.var='sumQ2_3', colorM=get.color(magnitude))
}
legend('bottomright', legend=c('Low', 'High'), bg='white', col=1:2, pch=15, title='Magnitude')

#legend('bottomright', legend=c('0.8/0.2', '0.7/0.3', '0.65/0.35', '0.6/0.4'), bg='white', col=1:4, pch=15, title='Difficulty')

# drift rates
plot(0,0, type='n', xlim=range(meanV1OverTimeM$bin)+c(-.5, .5), ylim=c(1.0, 4.5), xlab='Trial bin', ylab='Drift rates', main='D. Drift rates')
abline(h=seq(0, 5, .25), col='grey')
abline(v=seq(0, 10, 2), col='grey')
for(ease in unique(allV1OverTimeM$ease)) {
  draw.polygon(meanV1OverTimeM[meanV1OverTimeM$ease==ease,], dep.var='mean_v.r1', colorM=get.color(ease))
  draw.polygon(meanV2OverTimeM[meanV2OverTimeM$ease==ease,], dep.var='mean_v.r2', colorM=get.color(ease))
}
#if(savePlot) dev.off()

