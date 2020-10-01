setwd("/home/nsteven/RLDMC2")
source('./dmc/models/RW/dists.R')
source('./parameterRecoveries/utils.R')  # simulation function here

# Can we recover? ---------------------------------------------------------
source ("dmc/dmc.R")
source('utils.R')
load_model ("RW", "arw-RL-risk-mag-Niek.R")

# some names for saving
modelName = 'arw-RL-risk-mag-Niek'
samplesDir <- 'parameterRecoveries/samples'

dataName <- 'parameterRecovery-exp2-Risk'
fn = paste0('model-', modelName, '_data-', dataName)

#### Model set-up ----
model <- model.dmc(p.map=list(A="1",t0="1",st0="1",s="1",
                              B0="1",
                              SR="1", aV="1",
                              RR="1", aR="1",
                              V0="1", wV="1",
                              wS="1"),
                   match.map=list(M=list(s1=1, s1=2)),
                   constants=c(st0=0, s=1, RR=-10,
                               SR=-10, A=0),
                   factors=list(S=c("s1")), 
                   responses=c("r1","r2"),
                   type="norm")

p.vector  <- p.vector  <- c(t0=.2, aV=-1.6, wS=1, aR=-1.6,
                            wV=1, B0=1, V0=1.5)
oldwd <- getwd()
setwd("/home/nsteven/RLDMC/")
samples <- loadSamples("samples/model-arw-RL-CentralTendencyRisk-mag-Niek_data-exp2")
setwd(oldwd)
samplesSummary <- summary.dmc(samples)
medianPars <- data.frame(do.call(rbind, lapply(samplesSummary, function(x) x$quantiles[,3])))

#### Data ----
# Simulate data with the median estimated parameters and the exact same design
setProbability <- list(c(0.2, 0.8), c(0.2, 0.8), c(0.2, 0.8), c(0.2, 0.8),
                       c(0.3, 0.7), c(0.3, 0.7), c(0.3, 0.7), c(0.3, 0.7))
allData <- NULL
trueParameters <- list()
for(s in 1:nrow(medianPars)) {
  B <- medianPars$B0[s]
  v0 <- medianPars$V0[s]
  wV <- medianPars$wV[s]
  wS <- medianPars$wS[s]
  alpha <- runif(1, 0.001, 0.07)
  t0 <- medianPars$t0[s]
  aR <- runif(1, 0.01, 0.25)
  data <- simulate.full.reversal(nTrialsPerSet = 64, nSets=8, reversal_trialN = 32,
                                 B=B, v0=v0, wV=wV, wS=wS, t0=t0, alpha=alpha, aR = aR,
                                 setProbability = setProbability)
  trueParameters[[s]] <- c(B=B, v0=v0, wV=wV, wS=wS, alpha=alpha, t0=t0, aR = aR)
  data$s = s
  print(s)  # progress
  allData <- rbind(allData, data)
}
allData$s <- factor(allData$s)
save(allData, trueParameters, file = './parameterRecoveries/data/parameterRecovery-Risk-exp2.RData')

load('./parameterRecoveries/data/parameterRecovery-Risk-exp2.RData')
### ugly work-around, sorry about this
allData$rt <- allData$RT
allData$choiceIsHighP <- ifelse(allData$R==2, 1, 0)

cvs <- list()
choiceIdx <- list()
for(sub in unique(allData$s)) {
  d <- prepareForFitting(allData[allData$s==sub,])
  cvs[[sub]] <- d$outcomes
  choiceIdx[[sub]] <- d$VVchoiceIdx
}
data <- allData  # convert to DMC style data frame
data$S <- factor('s1')
data$R <- factor(data$R, levels=c(1,2), labels=c('r1', 'r2'))
data <- data.model.dmc(data[,c('s', 'S', 'R', 'RT', 'stimulus_set')], model)

# append cvs by sub as well
for(sub in unique(allData$s)) {
  d <- prepareForFitting(allData[allData$s==sub,])
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
doSample(data, p.prior, pp.prior, nmcBurn=250, nCores=30, restart=TRUE, fileName=file.path(samplesDir, fn))
# samples <- h.samples.dmc(nmc=1000, p.prior=p.prior, data=data, pp.prior=pp.prior)
# samples <- h.RUN.dmc(hsamples = samples, cores=20, cut.converge=1.03, thorough = TRUE, verbose=TRUE, saveFn=file.path(samplesDir, fn))
# save(samples, save=paste0(saveFn, '.RData'))


# # ## check fit
load('./parameterRecoveries/samples/model-arw-RL-risk-mag-Niek_data-parameterRecovery-exp2-Risk.RData')
load('./parameterRecoveries/data/parameterRecovery-exp2.RData')
samples <- hsamples; rm(hsamples)
plot.dmc(samples)

#samples <- h.samples.dmc(samples=samples, nmc=0, add=TRUE, remove=1:700)
posteriorSummary <- summary.dmc(samples)

trueParams <- data.frame(do.call(rbind, trueParameters))
posteriorMedian <- data.frame(do.call(rbind, lapply(posteriorSummary, function(x) x$quantiles[,3])))
posteriorMinCI <- data.frame(do.call(rbind, lapply(posteriorSummary, function(x) x$quantiles[,1])))
posteriorMaxCI <- data.frame(do.call(rbind, lapply(posteriorSummary, function(x) x$quantiles[,5])))

posteriorMinCI$B <- posteriorMinCI$B0
posteriorMinCI$alpha <- pnorm(posteriorMinCI$aV)
posteriorMinCI$v0 <- posteriorMinCI$V0
posteriorMinCI <- posteriorMinCI[,colnames(trueParams)]

posteriorMaxCI$B <- posteriorMaxCI$B0
posteriorMaxCI$alpha <- pnorm(posteriorMaxCI$aV)
posteriorMaxCI$v0 <- posteriorMaxCI$V0
posteriorMaxCI <- posteriorMaxCI[,colnames(trueParams)]

posteriorMedian$B <- posteriorMedian$B0
posteriorMedian$alpha <- pnorm(posteriorMedian$aV)
posteriorMedian$v0 <- posteriorMedian$V0
posteriorMedian <- posteriorMedian[,colnames(trueParams)]

rmse <- function(x, y) {
  sqrt(mean((x-y)^2))
}

pdf(file='figures/exp2_parrec_Niek.pdf', width=6, height=4)
par(oma=c(0,0,0,0), mar=c(3, 4, 2, 0.75) + 0.1, mfrow=c(2,3), mgp=c(2.75,.75,0), las=1)
for(parName in colnames(posteriorMedian)) {
  x <- trueParams[,parName]
  y <- posteriorMedian[,parName]
  cimin <- posteriorMinCI[,parName]
  cimax <- posteriorMaxCI[,parName]
  col=1
  
  # Rename parameters for plot. Note that a = threshold
  if(parName == 'B') { main <- expression(italic('a')); legend.pos='bottomright'}
  if(parName == 'v0') { main <- expression(italic('V'[0])); legend.pos='bottomright'}
  if(parName == 'wV') { main <- expression(italic('w'['D'])); legend.pos='bottomright'}
  if(parName == 'wS') { main <- expression(italic('w'['S'])); legend.pos='bottomright'}
  if(parName == 't0') { main <- expression(italic('t'[0])); legend.pos='bottomright'}
  if(parName == 'alpha') { main <- expression(italic(alpha)); legend.pos='bottomright'}
  
  plot(x, y, xlab='', ylab='Median posterior', main=main, type='n')
  mtext('Data-generating value', side=1, line=2, cex=0.66)
  abline(a=0, b=1)
  points(x, y, col=col)
  #  arrows(x, cimin, x, cimax, length=0.05, angle=90, code=3, col=col)
  #  legend(legend.pos, paste0('r = ', round(cor(x, y), 2), '\nRMSE = ', round(rmse(x,y), 2)), bty='n')
  tmp <- legend(legend.pos, c(" ", " "), bty='n', xjust=1, 
                text.width = strwidth("RMSE = 0.03"))
  text(tmp$rect$left + tmp$rect$w, tmp$text$y,
       c(paste0('r = ', round(cor(x, y), 2)), 
         paste0('RMSE = ', round(rmse(x,y), 2))), pos = 2)
}
dev.off()



# fit?
pp = h.post.predict.dmc(samples = samples, adapt=TRUE, save.simulation = TRUE)
ppNoSim <- h.pp.summary(pp, samples=samples, cores=30)
library(moments)
#### Append stimulus set info to data & model --------
nBins <- 10
dat <- allData
if(dataName == 'exp4') dat$ease <- dat$meanP  #   cheat a little bit and use 'meanP' as 'ease' for plotting
dat$ease <- 0.4
dat$stimulus_set <- rep(rep(1:4, each=64), 47)
dat$sub <- dat$s
pp2 <- lapply(1:length(pp), addStimSetInfo, input=pp, orig_dat=dat)
data2 <- lapply(1:length(data), addStimSetInfo, input=data, orig_dat=dat)
if(!sfIsRunning()) sfInit(parallel=TRUE, cpus =30); sfLibrary(moments)
pp3 <- sfLapply(pp2, calculateByBin)
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
#pdf(file=file.path('plots', paste0(fn, '.pdf')), width=8, height=4)
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
  plotDataPPBins(data=meanRTsByEase[[1]][idxD,], pp=meanRTsByEase[[2]][idxM,], draw.legend=FALSE, ylim=c(0.55, 0.90)); title('Mean RT over time')
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
  plotDataPPBins(data=meanRTsByEase[[1]][idxD,], pp=meanRTsByEase[[2]][idxM,], draw.legend=FALSE, ylim=c(0.5, 0.75))
  title(ease)
  plotDataPPBins(data=meanAccByEase[[1]][idxD,], pp=meanAccByEase[[2]][idxM,], dep.var='acc', draw.legend=FALSE, ylim=c(0.1, .9))
}

### Lastly, plot by subject
# nf <- layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE), respect = TRUE)
par(mar=c(5, 4, 2, 2) + 0.1, mfrow=c(1,3))
for(i in 1:length(ppNoSim)) {
  plot.pp.dmc(ppNoSim[[i]], style='cdf', x.min.max = c(0, 3), no.layout = TRUE, do.main=FALSE,
              fits.lcol='blue', data.lwd.mult = 2, lwd=c(2,2)); title(paste0('Participant ', i)) #Overall def. CDFs')
  plotDataPPBins(data=meanRTsBySub[[1]][[i]], pp=meanRTsBySub[[2]][[i]], ylim=c(0.35, 1)); title('Mean RT over time')
  plotDataPPBins(data=meanAccsBySub[[1]][[i]], pp=meanAccsBySub[[2]][[i]], dep.var='acc'); title('Accuracy over time')
}

#dev.off()


dat$trialNreversal <- rep(-32:32, 6)
dat$bin <- dat$trialNreversal
pp2 <- lapply(1:length(pp), addStimSetInfo, input=pp, orig_dat=dat, addColumns='bin')
data2 <- lapply(1:length(data), addStimSetInfo, input=data, orig_dat=dat, addColumns='bin')
if(!sfIsRunning()) sfInit(parallel=TRUE, cpus =30); sfLibrary(moments)
pp3 <- sfLapply(pp2, calculateByBin)
data3 <- lapply(data2, calculateByBin)

##
meanRTsBySub <- list(lapply(data3, function(x) aggregate(RT~bin, x, mean)), lapply(pp3, function(x) aggregate(RT~bin*reps, x, mean)))
meanAccsBySub <- list(lapply(data3, function(x) aggregate(acc~bin, x, mean)), lapply(pp3, function(x) aggregate(acc~bin*reps, x, mean)))

meanRTsOverTrials <- list(getDescriptives(data3, dep.var='RT', attr.name='RTsOverBins', id.var1='~bin*s', id.var2="~bin"),
                          getDescriptives(pp3, dep.var='RT', attr.name='RTsOverBins'))
meanAccOverTrials <- list(getDescriptives(data3, dep.var='acc', attr.name='AccOverBins', id.var1='~bin*s', id.var2="~bin"),
                          getDescriptives(pp3, dep.var='acc', attr.name='AccOverBins'))

par(oma=c(3,1,1,0), mar=c(3, 4, 2, 0.5) + 0.1, mfrow=c(2,3), mgp=c(2.75,.75,0), las=1)
layout(matrix(c(1,1,2,2,3,3,4,4,4,5,5,5), byrow=TRUE, ncol=6))
i <- 0
plot.pp.dmc(ppNoSim, style='cdf', x.min.max = c(0.2, 2), ylim=c(0, 0.65),
            no.layout = TRUE, do.main=FALSE,
            fits.lcol='blue', data.lwd.mult = 2,
            lwds=c(2,2), ltys=c(2,1), pos=NA, xlab = '',
            xaxt='s', axvlines=c(0.5, 1, 1.5, 2),
            axhlines=seq(0, 1, .2))
title('Defective CDF')
legend('bottomright', c('Data', 'Model', 'Error'), lwd=c(4,4,2), lty=c(1,1,2), bty='n', col=c('black', 'blue', 'black'))
#legend('topright', c('Error'), lwd=c(2), lty=c(2), bty='n', col=c('black'))
mtext('RT (s)', side=1, cex=.66, line=2)


plotDataPPBins(data=meanRTsOverTime[[1]], pp=meanRTsOverTime[[2]],  ylab='RT (s)',
               xaxt='s',
               xlab='', ylim=c(0.55, 0.75), hline.by=0.05)
abline(v=0, lty=2)
title('Mean RT over time')
mtext('Trial bin', side=1, cex=.66, line=2)

plotDataPPBins(data=meanAccOverTime[[1]], pp=meanAccOverTime[[2]],
               xaxt='s',
               dep.var='acc',
               ylab=expression('Proportion choice A'),
               xlab = '',
               legend.pos='topright', ylim=c(0.25, 0.85), hline.by=0.1)
mtext('Trial bin', side=1, cex=.66, line=2)
title('Choices over time')


##
plotDataPPBins(data=meanRTsOverTrials[[1]], pp=meanRTsOverTrials[[2]],  ylab='RT (s)',
               xaxt='s', xlim=c(-60, 40), plot.model.points = FALSE, data.lwd=1.5, data.cex=1.5,
               xlab='', ylim=c(0.55, 0.75), hline.by=0.05, axvlines=seq(-50, 50, 10))
title('Mean RT after reversal')
abline(v=0, lty=2)
mtext('Trial (relative to reversal point)', side=1, cex=.66, line=2)

plotDataPPBins(data=meanAccOverTrials[[1]], pp=meanAccOverTrials[[2]],
               xaxt='s', xlim=c(-60, 40), plot.model.points=FALSE,
               dep.var='acc', ylab=expression('Proportion choice A'),
               xlab = '', data.lwd=1.5, data.cex=1.5,
               legend.pos='topright', ylim=c(0.25, 0.85), hline.by=0.05, axvlines=seq(-50, 50, 10))
title('Choices after reversal')
axis(1, at=seq(-50, 50, 10), labels=rep(NA, 5))
abline(v=0, lty=2)
mtext('Trial (relative to reversal point)', side=1, cex=.66, line=2)

