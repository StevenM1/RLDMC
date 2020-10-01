source('./dmc/models/RW/dists.R')
source('./parameterRecoveries/utils.R')  # simulation function here

# Can we recover? ---------------------------------------------------------
source ("dmc/dmc.R")
source('utils.R')
load_model ("RW", "arw-RL-mag-mods.R")

# some names for saving
modelName = 'arw-RL-mag'
samplesDir <- 'parameterRecoveries/samples'

dataName <- 'parameterRecovery-exp3'
fn = paste0('model-', modelName, '_data-', dataName)

#### Model set-up ----
#### Model ----
model <- model.dmc(p.map=list(A="1",st0="1",s="1",
                              t0="1",
                              B0="1",
                              SR="1", aV="1",
                              V0="1", wV="1",
                              driftMod="1",
                              V0Mod="cue",
                              B0Mod="cue",
                              aVMod="1",
                              t0Mod="1",
                              wS="1"),
                   match.map=list(M=list(s1=1, s1=2)),
                   constants=c(st0=0, s=1, 
                               driftMod=0,
                               B0Mod.ACC=0,
                               V0Mod.ACC=0,
                               t0Mod=0,
                               aVMod=0,
                               SR=-10, A=0),
                   factors=list(S=c("s1"), cue=c('ACC', 'SPD')), 
                   responses=c("r1","r2"),
                   type="norm")

p.vector  <- c(t0=.2,# t0.ACC=.2,
               aV=-1.6,# aV.ACC=-1, 
               wS=1, # wS.ACC=1, 
               wV=1, # wV.ACC=1, 
               B0=1, # B0.ACC=1, 
               V0=1.5, #V0.ACC=1.5,
               B0Mod.SPD=0,
               V0Mod.SPD=0)

tmp <- loadData('exp3')
dat <- tmp$dat
data <- tmp$data

samples <- loadSamples(fn='model-arw-RL-mag-SAT-BV02_data-exp3', samplesDir='samples')
samplesSummary <- summary.dmc(samples)
medianPars <- data.frame(do.call(rbind, lapply(samplesSummary, function(x) x$quantiles[,3])))

#### Data ----
# Simulate participants
setProbability <- list(c(0.2, 0.8), c(0.2, 0.8), c(0.2, 0.8),
                       c(0.3, 0.7), c(0.3, 0.7), c(0.3, 0.7),
                       c(0.4, 0.6), c(0.4, 0.6), c(0.4, 0.6))
allData <- NULL
trueParameters <- list()
for(s in 1:nrow(medianPars)) {
  B0 = medianPars$B0[s]
  B0Mod.SPD = medianPars$B0Mod.SPD[s]
  V0Mod.SPD = medianPars$V0Mod.SPD[s]
  V0 = medianPars$V0[s]
  wV = medianPars$wV[s]
  wS = medianPars$wS[s]
  alpha = pnorm(medianPars$aV[s])
  t0 = medianPars$t0[s]
  data <- simulate.SAT(nTrialsPerSet = 36, nSets=9, 
                       driftMod.ACC=0, driftMod.SPD=0,
                       B0=B0, B0Mod.ACC=0, B0Mod.SPD=B0Mod.SPD,
                       V0=V0, V0Mod.ACC=0, V0Mod.SPD=V0Mod.SPD,
                       wV=wV, wS=wS, 
                       t0=t0, alpha=alpha, 
#                        B0=B0, B0Mod.SPD=B0Mod.SPD, V0=V0, V0Mod.SPD=V0Mod.SPD, wV=wV, wS=wS, t0=t0, alpha=alpha,
                        setProbability = setProbability)
  trueParameters[[s]] <- c(B0=B0, B0Mod.SPD=B0Mod.SPD, V0=V0, V0Mod.SPD=V0Mod.SPD, wV=wV, wS=wS, alpha=alpha, t0=t0)
  data$s = s
  print(s)  # progress
  allData <- rbind(allData, data)
}
plot(allData$RT, col=allData$s)
allData$s <- factor(allData$s)
allData$cue <- factor(as.character(allData$cue), levels=c('ACC', 'SPD'))
save(allData, trueParameters, file = './parameterRecoveries/data/parameterRecovery-exp3.RData')

load('./parameterRecoveries/data/parameterRecovery-exp3.RData')
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
data <- data.model.dmc(data[,c('s', 'S', 'R', 'RT', 'cue', 'stimulus_set')], model)

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
doSample(data, p.prior, pp.prior, nmcBurn=1000, nCores=30, restart=TRUE, fileName=file.path(samplesDir, fn))


# # ## check fit
load('./parameterRecoveries/samples/model-arw-RL-mag_data-parameterRecovery-exp3.RData')
load('./parameterRecoveries/data/parameterRecovery-exp3.RData')
samples <- hsamples
plot.dmc(samples)

posteriorSummary <- summary.dmc(samples)

trueParams <- data.frame(do.call(rbind, trueParameters))
posteriorMedian <- data.frame(do.call(rbind, lapply(posteriorSummary, function(x) x$quantiles[,3])))
posteriorMinCI <- data.frame(do.call(rbind, lapply(posteriorSummary, function(x) x$quantiles[,1])))
posteriorMaxCI <- data.frame(do.call(rbind, lapply(posteriorSummary, function(x) x$quantiles[,5])))

posteriorMinCI$B0Mod.SPD <- posteriorMinCI$B0Mod.SPD
posteriorMinCI$B0 <- posteriorMinCI$B0
posteriorMinCI$V0Mod.SPD <- posteriorMinCI$V0Mod.SPD
posteriorMinCI$V0 <- posteriorMinCI$V0
posteriorMinCI$alpha <- pnorm(posteriorMinCI$aV)
#posteriorMinCI$v0 <- posteriorMinCI$V0
posteriorMinCI <- posteriorMinCI[,colnames(trueParams)]
posteriorMinCI$V0.ACC <- posteriorMinCI$V0
posteriorMinCI$V0.SPD <- posteriorMinCI$V0 * (1+posteriorMinCI$V0Mod.SPD)
posteriorMinCI$B0.ACC <- posteriorMinCI$B0
posteriorMinCI$B0.SPD <- posteriorMinCI$B0 * (1+posteriorMinCI$B0Mod.SPD)


posteriorMaxCI$B0Mod.SPD <- posteriorMaxCI$B0Mod.SPD
posteriorMaxCI$B0 <- posteriorMaxCI$B0
posteriorMaxCI$V0Mod.SPD <- posteriorMaxCI$V0Mod.SPD
posteriorMaxCI$V0 <- posteriorMaxCI$V0
posteriorMaxCI$alpha <- pnorm(posteriorMaxCI$aV)
#posteriorMaxCI$v0 <- posteriorMaxCI$V0
posteriorMaxCI <- posteriorMaxCI[,colnames(trueParams)]
posteriorMaxCI$V0.ACC <- posteriorMaxCI$V0
posteriorMaxCI$V0.SPD <- posteriorMaxCI$V0 * (1+posteriorMaxCI$V0Mod.SPD)
posteriorMaxCI$B0.ACC <- posteriorMaxCI$B0
posteriorMaxCI$B0.SPD <- posteriorMaxCI$B0 * (1+posteriorMaxCI$B0Mod.SPD)

posteriorMedian$B0Mod.SPD <- posteriorMedian$B0Mod.SPD
posteriorMedian$B0 <- posteriorMedian$B0
posteriorMedian$V0Mod.SPD <- posteriorMedian$V0Mod.SPD
posteriorMedian$V0 <- posteriorMedian$V0
posteriorMedian$alpha <- pnorm(posteriorMedian$aV)
posteriorMedian <- posteriorMedian[,colnames(trueParams)]
posteriorMedian$V0.ACC <- posteriorMedian$V0
posteriorMedian$V0.SPD <- posteriorMedian$V0 * (1+posteriorMedian$V0Mod.SPD)
posteriorMedian$B0.ACC <- posteriorMedian$B0
posteriorMedian$B0.SPD <- posteriorMedian$B0 * (1+posteriorMedian$B0Mod.SPD)
#posteriorMedian$v0 <- posteriorMedian$V0

trueParams$V0.ACC <- trueParams$V0
trueParams$V0.SPD <- trueParams$V0 * (1+trueParams$V0Mod.SPD)
trueParams$B0.ACC <- trueParams$B0
trueParams$B0.SPD <- trueParams$B0 * (1+trueParams$B0Mod.SPD)

rmse <- function(x, y) {
  sqrt(mean((x-y)^2))
}

pdf(file='figures/exp3_parrec.pdf', width=8, height=4)
par(oma=c(0,0,0,0), mar=c(3, 4, 2, 0.75) + 0.1, mfrow=c(2,4), mgp=c(2.75,.75,0), las=1)
plot.arrows=TRUE
for(parName in c('B0.ACC', 'B0.SPD', 'V0.ACC', 'V0.SPD', 'wV', 'wS', 't0', 'alpha')) {
  x <- trueParams[,parName]
  y <- posteriorMedian[,parName]
  cimin <- posteriorMinCI[,parName]
  cimax <- posteriorMaxCI[,parName]
  col=1
  
  # Rename parameters for plot. Note that a = threshold
  if(parName == 'B0.SPD') { main <- expression(italic('a'['speed'])); legend.pos='bottomright'}
  if(parName == 'B0.ACC') { main <- expression(italic('a'['accuracy'])); legend.pos='bottomright'}
  if(parName == 'V0.SPD') { main <- expression(italic('V'['0,speed'])); legend.pos='bottomright'}
  if(parName == 'V0.ACC') { main <- expression(italic('V'['0,accuracy'])); legend.pos='bottomright'}
  if(parName == 'B0Mod.SPD') { main <- expression(italic('a'['0,spd'])); legend.pos='bottomright'}
  if(parName == 'V0') { main <- expression(italic('V'[0])); legend.pos='bottomright'}
  if(parName == 'V0Mod.SPD') { main <- expression(italic('V'['0,spd'])); legend.pos='bottomright'}
  if(parName == 'wV') { main <- expression(italic('w'['D'])); legend.pos='bottomright'}
  if(parName == 'wS') { main <- expression(italic('w'['S'])); legend.pos='bottomright'}
  if(parName == 't0') { main <- expression(italic('t'[0])); legend.pos='bottomright'}
  if(parName == 'alpha') { main <- expression(italic(alpha)); legend.pos='bottomright'}
  
  if(plot.arrows) {
    ylim <- c(min(cimin), max(cimax))
  } else {
    ylim <- c(min(y), max(y))+c(-.1, .1)
  }
  plot(x, y, xlab='', ylab='Median posterior', main=main, type='n', ylim=ylim)
  mtext('Data-generating value', side=1, line=2, cex=0.66)
  abline(a=0, b=1)
  points(x, y, col=col)
  if(plot.arrows) arrows(x, cimin, x, cimax, length=0.05, angle=90, code=3, col=col)
  #  legend(legend.pos, paste0('r = ', round(cor(x, y), 2), '\nRMSE = ', round(rmse(x,y), 2)), bty='n')
  tmp <- legend(legend.pos, c(" ", " "), bty='n', xjust=1, 
                text.width = strwidth("RMSE = 0.03"))
  text(tmp$rect$left + tmp$rect$w, tmp$text$y,
       c(paste0('r = ', round(cor(x, y), 2)), 
         paste0('RMSE = ', round(rmse(x,y), 2))), pos = 2)
}
dev.off()

#coverage?
apply((trueParams < posteriorMaxCI) & (trueParams > posteriorMinCI), 2, mean)

h.gelman.diag.dmc(samples)

plot.dmc(samples[[80]])


for(i in 1:length(p.prior)) plot.prior(i, p.prior)


sub = 80
# what are the likelihoods?
p.vector <- as.numeric(trueParams[sub,])
names(p.vector) <- c('B0', 'V0', 'wV', 'wS', 'aV', 't0')
p.vector[5] <- qnorm(p.vector[5])
sum(log(likelihood.dmc(p.vector, samples[[sub]]$data)))

p.vector2 <- as.numeric(posteriorMedian[sub,])
names(p.vector2) <- c('B0', 'V0', 'wV', 'wS', 'aV', 't0')
p.vector2[5] <- qnorm(p.vector2[5])
sum(log(likelihood.dmc(p.vector2, samples[[sub]]$data)))

log_likelihoods <- data.frame(true=NA, medianPosterior=NA)[c()]
for(row in 1:nrow(posteriorMedian)) {
  p.vector <- as.numeric(trueParams[2,])
  names(p.vector) <- c('B0', 'V0', 'wV', 'wS', 'aV', 't0')
  p.vector[5] <- qnorm(p.vector[5])
  log_likelihoods[row,'true'] = sum(log(likelihood.dmc(p.vector, samples[[row]]$data)))
  
  p.vector2 <- as.numeric(posteriorMedian[2,])
  names(p.vector2) <- c('B0', 'V0', 'wV', 'wS', 'aV', 't0')
  p.vector2[5] <- qnorm(p.vector2[5])
  sum(log(likelihood.dmc(p.vector2, samples[[2]]$data)))
  log_likelihoods[row,'medianPosterior'] = sum(log(likelihood.dmc(p.vector2, samples[[row]]$data)))
}

par(mfrow=c(1,1))
x=log_likelihoods[,'true']
y=log_likelihoods[,'medianPosterior']
plot(x, y, type='n', xlab='True LL', ylab='median posterior LL')
abline(a=1, b=1)
points(x,y)

which(x>y)

log_likelihoods[,'true'] < -2000



par(mfrow=c(3,2))
for(parName in colnames(posteriorMedian)) {
  x <- trueParams[,parName]
  y <- posteriorMedian[,parName]
  
  #  col = ifelse(log_likelihoods[,'true'] < log_likelihoods[,'medianPosterior'], 2, 1)
  #col = ifelse(log_likelihoods[,'true'] < -2000, 2, 1)
  #  col=1
  col = ifelse(posteriorMedian$B > 1.3*trueParams$B, 2, 1)
  
  plot(x, y, xlab='True', ylab='Median posterior', main=parName, type='n')
  abline(a=0, b=1)
  points(x, y, col=col)
  legend('topleft', paste0('cor: ', round(cor(x, y), 2), '\ncor2: ', round(cor(x[col==1], y[col==1]), 2)), bty='n')
}


par(mfrow=c(1,1))
x=log_likelihoods[,'true']
y=log_likelihoods[,'medianPosterior']
plot(x, y, type='n', xlab='True LL', ylab='median posterior LL')
abline(a=1, b=1)
points(x,y, col=ifelse(posteriorMedian$B > 1.3*trueParams$B, 2, 1))



subIdx <- which(posteriorMedian$B > 1.3*trueParams$B)
plot.dmc(samples[[subIdx[1]]])

# anything specific about the data of these subjects?
dataSubset <- allData[allData$s %in% subIdx,]
aggregate(RT~s, dataSubset, mean)
aggregate(R==2~s, dataSubset, mean)
aggregate(RT>1.5~s, dataSubset, mean)
aggregate(RT~s, dataSubset, IQR)
aggregate(RT~s, dataSubset, range)

dataSubset2 <- allData[!allData$s %in% subIdx,]
aggregate(RT~s, dataSubset2, mean)
aggregate(R==2~s, dataSubset2, mean)
aggregate(RT>1.5~s, dataSubset2, mean)
plot(aggregate(RT~s, dataSubset2, IQR)[,2])
aggregate(RT~s, dataSubset2, range)

plot(aggregate(RT~s, dataSubset2, IQR)[,2])
points(aggregate(RT~s, dataSubset, IQR)[,2], col='red')

plot(aggregate(RT~s, dataSubset2, mean)[,2])
points(aggregate(RT~s, dataSubset, mean)[,2], col='red')

plot(aggregate(R==2~s, dataSubset2, mean)[,2])
points(aggregate(R==2~s, dataSubset, mean)[,2], col='red')

plot(density(aggregate(RT~s, dataSubset2, skewness)[,2]))
lines(density(aggregate(RT~s, dataSubset, skewness)[,2]), col='red')


pp = h.post.predict.dmc(samples = samples, adapt=TRUE,save.simulation = TRUE, cores=15)
ppNoSim <- h.pp.summary(pp, samples=samples)

# library(moments)
# #### Append stimulus set info to data & model --------
nBins <- 10
dat <- do.call(rbind, lapply(1:length(samples), function(x) {data <- samples[[x]]$data; data$sub <- x; data}))
data <- lapply(1:length(samples), function(x) samples[[x]]$data)
dat$ease <- NA   #as.factor(abs(dat$high_stim_prob-dat$low_stim_prob))
dat$ease[dat$stimulus_set %in% c(1,2,3)] <- 0.6   #as.factor(abs(dat$high_stim_prob-dat$low_stim_prob))
dat$ease[dat$stimulus_set %in% c(4,5,6)] <- 0.4   #as.factor(abs(dat$high_stim_prob-dat$low_stim_prob))
dat$ease[dat$stimulus_set %in% c(7,8,9)] <- 0.2   #as.factor(abs(dat$high_stim_prob-dat$low_stim_prob))

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

# # Start plotting
# #pdf(file=file.path('plots', paste0(fn, '.pdf')), width=8, height=4)
# plot.dmc(samples, hyper=TRUE, density=TRUE, layout=c(2,4))
# 
## Overall 3-panel plot ------------------------------------------------------------
# nf <- layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE), respect = TRUE)
par(mar=c(5, 4, 2, 2) + 0.1, mfrow=c(1,3))
plot.pp.dmc(ppNoSim, style='cdf', x.min.max = c(0, 3), no.layout = TRUE, do.main=FALSE,
            fits.lcol='blue', data.lwd.mult = 2, lwd=c(2,2)); title('Overall def. CDFs')
plotDataPPBins(data=meanRTsOverTime[[1]], pp=meanRTsOverTime[[2]]); title('Mean RT over time')
plotDataPPBins(data=meanAccOverTime[[1]], pp=meanAccOverTime[[2]], dep.var='acc'); title('Accuracy over time')



# 
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
# par(mar=c(5, 4, 2, 2) + 0.1, mfrow=c(1,3))
# for(i in 1:length(ppNoSim)) {
#   plot.pp.dmc(ppNoSim[[i]], style='cdf', x.min.max = c(0, 3), no.layout = TRUE, do.main=FALSE,
#               fits.lcol='blue', data.lwd.mult = 2, lwd=c(2,2)); title(paste0('Participant ', i)) #Overall def. CDFs')
#   plotDataPPBins(data=meanRTsBySub[[1]][[i]], pp=meanRTsBySub[[2]][[i]]); title('Mean RT over time')
#   plotDataPPBins(data=meanAccsBySub[[1]][[i]], pp=meanAccsBySub[[2]][[i]], dep.var='acc'); title('Accuracy over time')
# }

par(mar=c(5, 4, 2, 2) + 0.1, mfrow=c(1,3))
for(i in subIdx) {
  plot.pp.dmc(ppNoSim[[i]], style='cdf', x.min.max = c(0, 3), no.layout = TRUE, do.main=FALSE,
              fits.lcol='blue', data.lwd.mult = 2, lwd=c(2,2)); title(paste0('Participant ', i)) #Overall def. CDFs')
  plotDataPPBins(data=meanRTsBySub[[1]][[i]], pp=meanRTsBySub[[2]][[i]]); title('Mean RT over time')
  plotDataPPBins(data=meanAccsBySub[[1]][[i]], pp=meanAccsBySub[[2]][[i]], dep.var='acc'); title('Accuracy over time')
}


pairs.dmc(samples[[subIdx[1]]], start=2250)
pairs.dmc(samples[[subIdx[2]]], start=2250)
pairs.dmc(samples[[subIdx[3]]], start=2250)
pairs.dmc(samples[[subIdx[4]]], start=2250)
pairs.dmc(samples[[subIdx[5]]], start=2250)
pairs.dmc(samples[[subIdx[6]]], start=2250)

for(sub in subIdx) {
  plot.dmc(samples[[sub]])
  mtext(sub, outer=TRUE, side=3, line=-3)
}

plot.dmc(samples[[subIdx[3]]]); mtext(subIdx[3], side=3, outer=TRUE, line=-3)
plot.dmc(samples[[subIdx[4]]]); mtext(subIdx[4], side=3, outer=TRUE, line=-3)
plot.dmc(samples[[160]]); mtext(160, side=3, outer=TRUE, line=-3)

plot.dmc(samples, hyper=TRUE, density=TRUE)

plot.dmc(samples[[180]]); mtext(180, side=3, outer=TRUE, line=-3)

ii <- 69
plot.dmc(samples[[ii]]); mtext(ii, side=3, outer=TRUE, line=-3)
