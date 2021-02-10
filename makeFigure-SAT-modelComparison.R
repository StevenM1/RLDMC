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



# Which is the winning model? -----------------------------------------------
BPICDDMa <- getDataPpBPIC('ddm-RL-SAT-a', 'exp3', BPIConly=TRUE)$BPIC
BPICDDMm <- getDataPpBPIC('ddm-RL-SAT-m', 'exp3', BPIConly=TRUE)$BPIC
BPICDDMam <- getDataPpBPIC('ddm-RL-SAT-am', 'exp3', BPIConly=TRUE)$BPIC
BPICARD1 <- getDataPpBPIC('arw-RL-mag-SAT-B2', 'exp3', BPIConly=TRUE)$BPIC
BPICARD2 <- getDataPpBPIC('arw-RL-mag-SAT-W2', 'exp3', BPIConly=TRUE)$BPIC
BPICARD3 <- getDataPpBPIC('arw-RL-mag-SAT-V02', 'exp3', BPIConly=TRUE)$BPIC
BPICARD4 <- getDataPpBPIC('arw-RL-mag-SAT-BW2', 'exp3', BPIConly=TRUE)$BPIC
BPICARD5 <- getDataPpBPIC('arw-RL-mag-SAT-BV02', 'exp3', BPIConly=TRUE)$BPIC
BPICARD6 <- getDataPpBPIC('arw-RL-mag-SAT-V0W2', 'exp3', BPIConly=TRUE)$BPIC
BPICARD7 <- getDataPpBPIC('arw-RL-mag-SAT-BV0W2', 'exp3', BPIConly=TRUE)$BPIC

allBPICs <- cbind('DDM-a'=BPICDDMa[,2],  
                  'DDM-m'=BPICDDMm[,2], 
                  'DDM-am'=BPICDDMam[,2], 
                  'ARD-B'=BPICARD1[,2], 
                  'ARD-W'=BPICARD2[,2], 
                  'ARD-V0'=BPICARD3[,2], 
                  'ARD-BW'=BPICARD4[,2], 
                  'ARD-BV0'=BPICARD5[,2], 
                  'ARD-V0W'=BPICARD6[,2], 
                  'ARD-BV0W'=BPICARD7[,2])
sBPICs <- apply(allBPICs, 2, sum)
sBPICs-min(sBPICs)

# Load quantiles of winning model & RL-fARD, BPICs ---------------------------------------------------
# DDM
tmp <- getDataPpBPIC('ddm-RL-SAT-a', 'exp3')
qRTsDDM <- getqRTsByCue(tmp[['data3']], tmp[['pp3']])

# ARDf
tmp <- getDataPpBPIC('arw-RL-mag-SAT-BV02', 'exp3')
BPICARD <- tmp$BPIC
qRTsARD <- getqRTsByCue(tmp[['data3']], tmp[['pp3']])

# Combine -----------------------------------------------------------------
allqRTs <- list(qRTsDDM, qRTsARD)
allBPICs <- cbind(BPICDDMa[,2],  BPICDDMm[,2], BPICDDMam[,2], BPICARD[,2])
sBPICs <- apply(allBPICs, 2, sum)
sBPICs - sBPICs[4]


# Plot --------------------------------------------------------------------
layoutM <- matrix(1:25, nrow=5, byrow=TRUE)
layoutM[c(1, 5),1:2] <- 1
layoutM[c(1, 5),4:5] <- 9
layoutM[2:4,1:2] <- 2:7 #matrix(c(2:13), nrow=3, byrow=TRUE)
layoutM[2:4,4:5] <- 10:15 #matrix(c(2:13), nrow=3, byrow=TRUE)
layoutM[,3] <- 8
layoutM

if(savePlot) pdf('./figures/exp3-SAT.pdf', width=7, height=7/4*3)
layout(layoutM, heights = c(0.01, .8, 1, 1, 0.01), widths=c(1,1,.1,1,1))
par(oma=c(3,4,2,0), mar=c(0, 0, 1, 0.5) + 0.1, #mfcol=c(3,4), 
    mgp=c(2.75,.75,0), las=1, bty='l')
i <- 0
data.cex=1.5
corrRTylim <- errRTylim <- c(.35,1.1)
for(qRTs in allqRTs) {
  plot.new()
  if(i == 0) { mtext('RL-DDM', side=3, cex=.66*1.2, font=2, line=2); mtext(paste0('BPIC = ', round(apply(allBPICs, 2, sum)[1])), cex=.66, line=1)}
  if(i == 2) {plot.new(); mtext('RL-ARD', side=3, cex=.66*1.2, font=2, line=2); mtext(paste0('BPIC = ', round(apply(allBPICs, 2, sum)[4])), cex=.66, line=1)}
  for(cue in c('SPD', 'ACC')) {
    i <- i+1
    idxD <- qRTs$meanAccByCue[[1]]$cue==cue
    idxM <- qRTs$meanAccByCue[[2]]$cue==cue
    
    plotDataPPBins(data=qRTs$meanAccByCue[[1]][idxD,], pp=qRTs$meanAccByCue[[2]][idxM,],
                   xaxt='n', draw.legend = i==1, data.cex = data.cex,
                   dep.var='acc', ylab='', xlab = '', yaxt='n',
                   legend.pos='bottomright', ylim=c(0.5, 0.9), hline.by=0.1)
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
