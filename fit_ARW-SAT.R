################################################################################
rm(list=ls())
source ("dmc/dmc.R")
source('utils.R')
load_model ("RW", "arw-RL.R")

# some names for saving
modelName = 'arw-RL-SAT-B'
samplesDir <- 'samplesNew'

for(dataName in c('SAT')) {
  #dataName = 'calibration'
  fn = paste0('model-', modelName, '_data-', dataName)
  
  #### Model ----
  model <- model.dmc(p.map=list(A="1",t0="1",st0="1",s="1",
                                B0="cue",
                                SR="1", aV="1",
                                V0="1", wV="1"),
                     match.map=list(M=list(s1=1, s1=2)),
                     constants=c(st0=0, s=1,
                                 SR=-10, A=0),
                     factors=list(S=c("s1"), cue=c('ACC', 'SPD')), 
                     responses=c("r1","r2"),
                     type="norm")
  
  p.vector  <- c(t0=.2, aV=-1.6,
                 wV=1, B0.SPD=1, B0.ACC=1, V0=1.5)
  
  #### Data ----
  tmp <- loadData(dataName, subset=FALSE)
  data <- tmp[['data']]
  dat <- tmp[['dat']]
  
  data <- data.model.dmc(data, model)
  
  # append cvs by sub as well
  for(sub in unique(dat$sub)) {
    d <- prepareForFitting(dat[dat$sub==sub,])
    attr(data[[sub]], 'cvs') <- d$outcomes
    attr(data[[sub]], 'VVchoiceIdx') <- d$VVchoiceIdx
    attr(data[[sub]], 'startingValues') <- d$values
    attr(data[[sub]], 'trialsToIgnore') <- d$trialsToIgnore
  }
  
  # Check model, only need this when you are developing
  # likelihood.dmc(p.vector, data[[1]])
  # simulate.dmc(p.vector, model=model, n=nrow(data[[1]]), adapt=TRUE, cvs=attr(data[[1]], 'cvs'))
  
  #### Priors ----
  pp.prior <- getPriors(p.vector)
  p.prior <- pp.prior[[1]]
  
  #### Sample  -----------------------------------------------------------------
  doSample(data, p.prior, pp.prior, nmcBurn=150, nCores=15, restart=FALSE, fileName=file.path(samplesDir, fn))
}


### Plots & checks ----------
samples <- loadSamples(fn, samplesDir)

#samples <- h.samples.dmc(nmc=0,add=TRUE,samples=samples,remove=1:400)
plot.dmc(samples, hyper=TRUE, density=TRUE, layout=c(4,4))


# simulate posterior predictives
calculateByBin <- function(df) {
  df$acc <- as.integer(df$R)==2
  
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
pp = h.post.predict.dmc(samples = samples, adapt=TRUE,save.simulation = TRUE, cores=15)
ppNoSim <- h.pp.summary(pp, samples=samples)
library(moments)
#### Append stimulus set info to data & model --------
nBins <- 10
pp2 <- lapply(1:length(pp), addStimSetInfo, input=pp, orig_dat=dat)
data2 <- lapply(1:length(data), addStimSetInfo, input=data, orig_dat=dat)
pp3 <- lapply(pp2, calculateByBin)
data3 <- lapply(data2, calculateByBin)

##
meanRTsBySub <- list(lapply(data3, function(x) aggregate(RT~bin, x, mean)), lapply(pp3, function(x) aggregate(RT~bin*reps, x, mean)))
meanAccsBySub <- list(lapply(data3, function(x) aggregate(acc~bin, x, mean)), lapply(pp3, function(x) aggregate(acc~bin*reps, x, mean)))

meanRTsOverTime <- list(getDescriptives(data3, dep.var='RT', attr.name='RTsOverBins', id.var1='~bin*s', id.var2="~bin"),
                        getDescriptives(pp3, dep.var='RT', attr.name='RTsOverBins'))
meanAccOverTime <- list(getDescriptives(data3, dep.var='acc', attr.name='AccOverBins', id.var1='~bin*s', id.var2="~bin"),
                        getDescriptives(pp3, dep.var='acc', attr.name='AccOverBins'))

meanRTsByCue <- list(getDescriptives(data3, dep.var='RT', attr.name='RTsByCue', id.var1='~bin*cue*s', id.var2="~bin*cue"),
                     getDescriptives(pp3, dep.var='RT', attr.name='RTsByCue', id.var1='~reps*bin*cue*s', id.var2="~reps*bin*cue"))
meanAccByCue <- list(getDescriptives(data3, dep.var='acc', attr.name='AccByCue', id.var1='~bin*cue*s', id.var2="~bin*cue"),
                     getDescriptives(pp3, dep.var='acc', attr.name='AccByCue', id.var1='~reps*bin*cue*s', id.var2="~reps*bin*cue"))

# Start plotting
pdf(file=file.path('plots', paste0(fn, '.pdf')), width=12, height=4)
plot.dmc(samples, hyper=TRUE, density=TRUE, layout=c(2,4))

## Overall 3-panel plot ------------------------------------------------------------
# nf <- layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE), respect = TRUE)
par(mar=c(5, 4, 2, 2) + 0.1, mfcol=c(1,4))
plot.pp.dmc(ppNoSim, style='cdf', x.min.max = c(0, 3), no.layout = TRUE, do.main=TRUE,
            fits.lcol='blue', data.lwd.mult = 2, lwd=c(2,2), model.legend=FALSE)
plotDataPPBins(data=meanRTsByCue[[1]][meanRTsByCue[[1]]$cue=='SPD',], pp=meanRTsByCue[[2]][meanRTsByCue[[2]]$cue=='SPD',], 
               ylim=c(0.45, 0.75), colorD='black', colorM='red'); title('Mean RT over time')
plotDataPPBins(data=meanRTsByCue[[1]][meanRTsByCue[[1]]$cue=='ACC',], pp=meanRTsByCue[[2]][meanRTsByCue[[2]]$cue=='ACC',], 
               plot.new=FALSE, draw.legend=FALSE, colorD='black', colorM='blue')

plotDataPPBins(data=meanAccByCue[[1]][meanAccByCue[[1]]$cue=='SPD',], pp=meanAccByCue[[2]][meanAccByCue[[2]]$cue=='SPD',], 
               dep.var='acc', colorD='black', colorM='red'); title('Accuracy over time')
plotDataPPBins(data=meanAccByCue[[1]][meanAccByCue[[1]]$cue=='ACC',], pp=meanAccByCue[[2]][meanAccByCue[[2]]$cue=='ACC',], 
               dep.var='acc', plot.new=FALSE, draw.legend=FALSE, colorD='black', colorM='blue')

# ## 3-panel plot by ease ------------------------------------------------------------
# for(ease in unique(meanRTsByEase[[1]]$ease)) {
#   pp2Subset <- lapply(1:length(pp2), function(x) pp2[[x]][pp2[[x]]$ease==ease,])
#   dataSubset <- do.call(rbind, lapply(1:length(data2), function(x) data2[[x]][data2[[x]]$ease==ease,]))
#   ppNoSimSubset <- h.pp.summary(pp2Subset, samples=samples, dat=dataSubset, ignore.subjects=TRUE)
#   plot.pp.dmc(ppNoSimSubset, style='cdf', no.layout = TRUE, do.main=FALSE,
#               fits.lcol='blue', data.lwd.mult = 2, lwd=c(2,2)); title(paste0('Ease ', ease))
#   
#   idxD = meanRTsByEase[[1]]$ease == ease
#   idxM = meanRTsByEase[[2]]$ease == ease
#   plotDataPPBins(data=meanRTsByEase[[1]][idxD,], pp=meanRTsByEase[[2]][idxM,], draw.legend=FALSE); title('Mean RT over time')
#   plotDataPPBins(data=meanAccByEase[[1]][idxD,], pp=meanAccByEase[[2]][idxM,], dep.var='acc', draw.legend=FALSE); title('Mean RT over time')
# }

# Overall CDFs by bin ------------------------------------------------------------
layout(1)
par(mfrow=c(2,3))
for(i in 1:nBins) {
  pp2Subset <- lapply(1:length(pp2), function(x) pp2[[x]][pp2[[x]]$bin==i,])
  dataSubset <- do.call(rbind, lapply(1:length(data2), function(x) data2[[x]][data2[[x]]$bin==i,]))
  ppNoSimSubset <- h.pp.summary(pp2Subset, samples=samples, dat=dataSubset, ignore.subjects=TRUE)
  plot.pp.dmc(ppNoSimSubset, style='cdf', no.layout = TRUE, do.main=FALSE, aname='',
              fits.lcol='blue', data.lwd.mult = 2, lwd=c(2,2)); title(paste0('Bin ', i))
}

# Overall change in RT/Accuracy by ease ------------------------------------------------------------
par(mfcol=c(2,4))
for(ease in unique(meanRTsByEase[[1]]$ease)) {
  idxD = meanRTsByEase[[1]]$ease == ease
  idxM = meanRTsByEase[[2]]$ease == ease
  plotDataPPBins(data=meanRTsByEase[[1]][idxD,], pp=meanRTsByEase[[2]][idxM,], draw.legend=FALSE, ylim=c(0.6, 0.8))
  title(ease)
  plotDataPPBins(data=meanAccByEase[[1]][idxD,], pp=meanAccByEase[[2]][idxM,], dep.var='acc', draw.legend=FALSE, ylim=c(0.4, 1))
}

### Lastly, plot by subject
# nf <- layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE), respect = TRUE)
par(mar=c(5, 4, 2, 2) + 0.1, mfrow=c(1,4))
for(i in 1:length(ppNoSim)) {
  plot.pp.dmc(ppNoSim[[i]], style='cdf', x.min.max = c(0, 3), no.layout = TRUE, do.main=FALSE,
              fits.lcol='blue', data.lwd.mult = 2, lwd=c(2,2)); title(paste0('Participant ', i)) #Overall def. CDFs')
  plotDataPPBins(data=meanRTsBySub[[1]][[i]], pp=meanRTsBySub[[2]][[i]]); title('Mean RT over time')
  plotDataPPBins(data=meanAccsBySub[[1]][[i]], pp=meanAccsBySub[[2]][[i]], dep.var='acc'); title('Accuracy over time')
}

dev.off()
