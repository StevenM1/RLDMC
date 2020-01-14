################################################################################
rm(list=ls())
source ("dmc/dmc.R")
source('utils.R')
load_model ("ddm", "ddm-RL.R")

# some names for saving
modelName = 'ddm-RL'
samplesDir <- 'samplesNew'

for(dataName in c('calibration', 'annie-chris', 'calibration-subset')) {
  #dataName = 'calibration'
  fn = paste0('model-', modelName, '_data-', dataName)
  
  #### Model ----
  model <- model.dmc(p.map=list(aV="1", SR='1', m='1',
                                t0='1', st0='1', d='1',
                                z='1', sz='1', a='1',
                                sv='1'),
                     match.map=list(M=list(s1=1, s1=2)),
                     constants=c(SR=-10, sv=0, sz=0, z=0.5, st0=0, d=0),
                     factors=list(S=c("s1")), 
                     responses=c("r1","r2"),
                     type="norm")
  
  p.vector  <- c(aV=-1.6, m=1, t0=0.25, a=1)
  
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
  #source('doSample.R')
  doSample(data, p.prior, pp.prior, nmcBurn=150, nCores=15, restart=FALSE, fileName=file.path(samplesDir, fn))
}

### Plots & checks ----------
samples <- loadSamples(fn, samplesDir)


# load('calibration_ddm_autoconverge.RData')
# # load('calibration_ddm_old.Rdata'); samples <- burn1
# samples <- h.samples.dmc(nmc=0,add=TRUE,samples=samples,remove=1:2000)
# plot.dmc(samples, hyper=TRUE, density=TRUE, layout=c(4,4))
# # simulate posterior predictives
# 
# simulate posterior predictives
pp = h.post.predict.dmc(samples = samples, adapt=TRUE,save.simulation = TRUE, cores=15)
ppNoSim <- h.pp.summary(pp, samples=samples)
library(moments)
#### Append stimulus set info to data & model --------
nBins <- 10
dat$ease <- as.factor(abs(dat$high_stim_prob-dat$low_stim_prob))
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

meanRTsByEase <- list(getDescriptives(data3, dep.var='RT', attr.name='RTsByEase', id.var1='~bin*ease*s', id.var2="~bin*ease"),
                      getDescriptives(pp3, dep.var='RT', attr.name='RTsByEase', id.var1='~reps*bin*ease*s', id.var2="~reps*bin*ease"))
meanAccByEase <- list(getDescriptives(data3, dep.var='acc', attr.name='AccByEase', id.var1='~bin*ease*s', id.var2="~bin*ease"),
                      getDescriptives(pp3, dep.var='acc', attr.name='AccByEase', id.var1='~reps*bin*ease*s', id.var2="~reps*bin*ease"))

# Start plotting
pdf(file=file.path('plots', paste0(fn, '.pdf')), width=8, height=4)
plot.dmc(samples, hyper=TRUE, density=TRUE, layout=c(2,4))

## Overall 3-panel plot ------------------------------------------------------------
# nf <- layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE), respect = TRUE)
par(mar=c(5, 4, 2, 2) + 0.1, mfrow=c(1,3))
plot.pp.dmc(ppNoSim, style='cdf', x.min.max = c(0, 3), no.layout = TRUE, do.main=FALSE,
            fits.lcol='blue', data.lwd.mult = 2, lwd=c(2,2)); title('Overall def. CDFs')
plotDataPPBins(data=meanRTsOverTime[[1]], pp=meanRTsOverTime[[2]]); title('Mean RT over time')
plotDataPPBins(data=meanAccOverTime[[1]], pp=meanAccOverTime[[2]], dep.var='acc'); title('Accuracy over time')

## 3-panel plot by ease ------------------------------------------------------------
for(ease in unique(meanRTsByEase[[1]]$ease)) {
  pp2Subset <- lapply(1:length(pp2), function(x) pp2[[x]][pp2[[x]]$ease==ease,])
  dataSubset <- do.call(rbind, lapply(1:length(data2), function(x) data2[[x]][data2[[x]]$ease==ease,]))
  ppNoSimSubset <- h.pp.summary(pp2Subset, samples=samples, dat=dataSubset, ignore.subjects=TRUE)
  plot.pp.dmc(ppNoSimSubset, style='cdf', no.layout = TRUE, do.main=FALSE,
              fits.lcol='blue', data.lwd.mult = 2, lwd=c(2,2)); title(paste0('Ease ', ease))
  
  idxD = meanRTsByEase[[1]]$ease == ease
  idxM = meanRTsByEase[[2]]$ease == ease
  plotDataPPBins(data=meanRTsByEase[[1]][idxD,], pp=meanRTsByEase[[2]][idxM,], draw.legend=FALSE); title('Mean RT over time')
  plotDataPPBins(data=meanAccByEase[[1]][idxD,], pp=meanAccByEase[[2]][idxM,], dep.var='acc', draw.legend=FALSE); title('Mean RT over time')
}

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
par(mar=c(5, 4, 2, 2) + 0.1, mfrow=c(1,3))
for(i in 1:length(ppNoSim)) {
  plot.pp.dmc(ppNoSim[[i]], style='cdf', x.min.max = c(0, 3), no.layout = TRUE, do.main=FALSE,
              fits.lcol='blue', data.lwd.mult = 2, lwd=c(2,2)); title(paste0('Participant ', i)) #Overall def. CDFs')
  plotDataPPBins(data=meanRTsBySub[[1]][[i]], pp=meanRTsBySub[[2]][[i]]); title('Mean RT over time')
  plotDataPPBins(data=meanAccsBySub[[1]][[i]], pp=meanAccsBySub[[2]][[i]], dep.var='acc'); title('Accuracy over time')
}

dev.off()



# ## Plot expected values over time ------------
# allEV1OverTimeM <- do.call(rbind, (lapply(1:length(pp3), function(x) {tmp <- attr(pp3[[x]], 'SR.r1OverBins'); tmp$s <- x; tmp})))
# allEV2OverTimeM <- do.call(rbind, (lapply(1:length(pp3), function(x) {tmp <- attr(pp3[[x]], 'SR.r2OverBins'); tmp$s <- x; tmp})))
# meanEV1OverTimeM <- aggregate(SR.r1~reps*bin*ease, aggregate(SR.r1~reps*bin*ease*s, allEV1OverTimeM, mean), mean)
# meanEV2OverTimeM <- aggregate(SR.r2~reps*bin*ease, aggregate(SR.r2~reps*bin*ease*s, allEV2OverTimeM, mean), mean)
# 
# allVOverTimeM <- do.call(rbind, (lapply(1:length(pp3), function(x) {tmp <- attr(pp3[[x]], 'vOverBins'); tmp$s <- x; tmp})))
# meanVOverTimeM <- aggregate(v~reps*bin*ease, aggregate(v~reps*bin*ease*s, allVOverTimeM, mean), mean)
# 
# allPEOverTimeM <- do.call(rbind, (lapply(1:length(pp3), function(x) {tmp <- attr(pp3[[x]], 'PEsOverBins'); tmp$s <- x; tmp})))
# meanPEsOverTimeM <- aggregate(PEs~reps*bin*ease, aggregate(PEs~reps*bin*ease*s, allPEOverTimeM, mean), mean)
# 
# # allRisk1OverTimeM <- do.call(rbind, (lapply(1:length(pp3), function(x) {tmp <- attr(pp3[[x]], 'risk.r1OverBins'); tmp$s <- x; tmp})))
# # risk1OverTimeM <- aggregate(risk.r1~reps*bin*ease, aggregate(risk.r1~reps*bin*ease*s, allRisk1OverTimeM, mean), mean)
# # allRisk2OverTimeM <- do.call(rbind, (lapply(1:length(pp3), function(x) {tmp <- attr(pp3[[x]], 'risk.r2OverBins'); tmp$s <- x; tmp})))
# # risk2OverTimeM <- aggregate(risk.r2~reps*bin*ease, aggregate(risk.r2~reps*bin*ease*s, allRisk2OverTimeM, mean), mean)
# 
# 
# allb1OverTimeM <- do.call(rbind, (lapply(1:length(pp3), function(x) {tmp <- attr(pp3[[x]], 'b.r1OverBins'); tmp$s <- x; tmp})))
# allb2OverTimeM <- do.call(rbind, (lapply(1:length(pp3), function(x) {tmp <- attr(pp3[[x]], 'b.r2OverBins'); tmp$s <- x; tmp})))
# meanb1OverTimeM <- aggregate(b.r1~reps*bin*ease, aggregate(b.r1~reps*bin*ease*s, allb1OverTimeM, mean), mean)
# meanb2OverTimeM <- aggregate(b.r2~reps*bin*ease, aggregate(b.r2~reps*bin*ease*s, allb2OverTimeM, mean), mean)
# 
# 
# 
# par(mfrow=c(1,1))
# plot(0,0, type='n', xlim=range(meanEV1OverTimeM$bin)+c(-.5, .5), ylim=c(0.3, 1), xlab='Bin', ylab='mean EV1')
# for(ease in unique(meanEV1OverTimeM$ease)) {
#   toplot <- meanEV1OverTimeM[meanEV1OverTimeM$ease==ease,]
#   points(toplot$bin, 
#          toplot$SR.r1, 
#          col=ifelse(ease=="0.6", 1, ifelse(ease=="0.4", 2, ifelse(ease=="0.3", 4, 5))), pch=20, cex=.25)
#   toplot <- meanEV2OverTimeM[meanEV2OverTimeM$ease==ease,]
#   points(toplot$bin, 
#          toplot$SR.r2, 
#          col=ifelse(ease=="0.6", 1, ifelse(ease=="0.4", 2, ifelse(ease=="0.3", 4, 5))), pch=20, cex=.25)
# }
# 
# tmp <- aggregate(SR.r1~bin, meanEV1OverTimeM, mean)
# tmp2 <- aggregate(SR.r2~bin, meanEV2OverTimeM, mean)
# lines(tmp$bin, tmp$SR.r1)
# lines(tmp$bin, tmp2$SR.r2)
# tmp <- cbind(tmp, tmp2$SR.r2)
# colnames(tmp) <- c('bin', 'r1', 'r2')
# 
# tmp$r1+tmp$r2
# 
# plot(0,0, type='n', xlim=range(risk1OverTimeM$bin)+c(-.5, .5), ylim=c(0., 1), xlab='Bin', ylab='Risk')
# for(ease in unique(risk1OverTimeM$ease)) {
#   toplot <- risk1OverTimeM[risk1OverTimeM$ease==ease,]
#   points(toplot$bin,
#          toplot$risk.r1,
#          col=ifelse(ease=="0.6", 1, ifelse(ease=="0.4", 2, ifelse(ease=="0.3", 4, 5))), pch=20, cex=.25)
#   toplot <- risk2OverTimeM[risk2OverTimeM$ease==ease,]
#   points(toplot$bin, 
#          toplot$risk.r2, 
#          col=ifelse(ease=="0.6", 1, ifelse(ease=="0.4", 2, ifelse(ease=="0.3", 4, 5))), pch=20, cex=.25)
# }
# 
# 
# plot(0,0, type='n', xlim=range(meanV1OverTimeM$bin)+c(-.5, .5), ylim=c(1.5, 2.5), xlab='Bin', ylab='mean v')
# for(ease in unique(meanb1OverTimeM$ease)) {
#   toplot <- meanb1OverTimeM[meanb1OverTimeM$ease==ease,]
#   points(toplot$bin, 
#          toplot$b.r1, 
#          col=ifelse(ease=="0.6", 1, ifelse(ease=="0.4", 2, ifelse(ease=="0.3", 4, 5))), pch=20, cex=.25)
#   toplot <- meanb2OverTimeM[meanb2OverTimeM$ease==ease,]
#   points(toplot$bin, 
#          toplot$b.r2, 
#          col=ifelse(ease=="0.6", 1, ifelse(ease=="0.4", 2, ifelse(ease=="0.3", 4, 5))), pch=20, cex=.25)
# }
# 
# plot(0,0, type='n', xlim=range(meanPEsOverTimeM$bin)+c(-.5, .5), ylim=c(0, .7), xlab='Bin', ylab='b')
# for(ease in unique(meanPEsOverTimeM$ease)) {
#   toplot <- meanPEsOverTimeM[meanPEsOverTimeM$ease==ease,]
#   points(toplot$bin, 
#          toplot$PEs, 
#          col=ifelse(ease=="0.6", 1, ifelse(ease=="0.4", 2, ifelse(ease=="0.3", 4, 5))), pch=20, cex=.25)
# }
# 
# 
# 
# ###
# for(pp in unique(dat$pp)) {
#   for(stimulus_set in unique(dat$stimulus_set)) {
#     idx <- dat$pp == pp & dat$stimulus_set == stimulus_set
#     dat[idx, 'trialNthisStimSet'] <- 1:sum(idx)
#   }
# }
# 
# 
# library(lme4)
# library(lmerTest)
# lm <- lmer(RT~ ease*trialNthisStimSet + ease|pp, dat)
# 
# summary(lm)
# 
# 
# 
# 
# 
# 
# par(mfrow=c(1,1))
# i=1
# plot(0,0, xlim=c(1,10), ylim=c(0.5, 0.85), xlab='bin', ylab='RT', type='n')
# for(ease in unique(meanRTsByEase[[1]]$ease)) {
#   idxD = meanRTsByEase[[1]]$ease == ease
#   toPlot <- meanRTsByEase[[1]][idxD,]
#   
#   offset = 0.8 - toPlot$RT[1]
#   toPlot$RT <- toPlot$RT+offset
#   lines(toPlot$bin, toPlot$RT, col=i)
#   i=i+1
#   #  plotDataPPBins(data=meanRTsByEase[[1]][idxD,], pp=meanRTsByEase[[2]][idxM,], draw.legend=FALSE)
#   #  plotDataPPBins(data=meanAccByEase[[1]][idxD,], pp=meanAccByEase[[2]][idxM,], dep.var='acc', draw.legend=FALSE)
#   #  title(ease)
#   #  plotDataPPBins(data=meanSkewByEase[[1]][idxD,], pp=meanSkewByEase[[2]][idxM,], dep.var='RT', draw.legend=FALSE)
# }
