#rm(list=ls())
source('./dmc/models/RW/dists.R')
source('./parameterRecoveries/utils.R')  # simulation function here

# Can we recover? ---------------------------------------------------------
source ("dmc/dmc.R")
source('utils.R')
load_model ("RW", "arw-RL-mag.R")

# some names for saving
modelName = 'arw-RL-mag'
samplesDir <- 'parameterRecoveries/samples'

dataName <- 'parameterRecovery-exp-n360'
fn = paste0('model-', modelName, '_data-', dataName)

#### Model set-up ----
model <- model.dmc(p.map=list(A="1",t0="1",st0="1",s="1",
                              B0="1", wS='1',
                              SR="1", aV="1",
                              V0="1", wV="1"),
                   match.map=list(M=list(s1=1, s1=2)),
                   constants=c(st0=0, s=1,
                               SR=-10, A=0),
                   factors=list(S=c("s1")), 
                   responses=c("r1","r2"),
                   type="norm")

p.vector  <- c(t0=.2, aV=-1.6, wS=1,
               wV=1, B0=1, V0=1.5)

#samples <- loadSamples(fn='model-arw-RL-mag_data-exp1', samplesDir='samples')
#samplesSummary <- summary.dmc(samples)
#medianPars <- data.frame(do.call(rbind, lapply(samplesSummary, function(x) x$quantiles[,3])))

#### Data ----
# Simulate data with the median estimated parameters and the exact same design
nSubs <- 500
# Let's simulate a design with in total 9 stimuli sets (this would typically be three experimental blocks) with varying difficulty (reward contingencies)
setProbability <- list(c(0.2, 0.8), c(0.2, 0.8), c(0.2, 0.8),
                       c(0.3, 0.7), c(0.3, 0.7), c(0.3, 0.7),
                       c(0.4, 0.6), c(0.4, 0.6), c(0.4, 0.6))
nTrialsPerSet <- 40  # each stimulus set is shown 40 times, so we have 9*40 = 360 trials per participant

allData <- NULL
trueParameters <- list()
for(s in 1:nSubs) {
  while(TRUE) {  # this is generally bad coding practice
    B0 <- rnorm(1, mean=2, sd=1)
    V0 <- rnorm(1, mean=2, sd=1)
    wV <- rnorm(1, mean=3, sd=1)
    wS <- rnorm(1, mean=0.5, sd=1)
    aV <- pnorm(rnorm(n=1, mean=-2.3, sd=1))  # learning rate, sampled on a logit scale and then transformed to normal scale
    t0 <- max(rnorm(1, .3, .1), 0.025)
    if(any(c(B0, V0, wV)<0)) next  # these parameters cant be lower than 0, generate again
    
    data <- simulate.full(nTrialsPerSet = nTrialsPerSet, nSets=length(setProbability),   # function defined in parameterRecoveryUtils.R
                          B=B0, v0=V0, wV=wV, wS=wS, t0=t0, alpha=aV,
                          setProbability = setProbability)
    
    if(mean(data$RT) > 1) next  # too slow
    if(mean(data$RT) < 0.3) next  # too fast
    if(mean(data$R==2) < 0.6) next  # no learning / accuracy not high enough
    if(mean(data$R==2) > 0.9) next  # accuracy too high, too few observations of 'error' distribution
    
    ## in real data, error RTs are slower (in one of my RL datasets, ~1.081 (SD 0.098) times slower) than correct RTs 
    ## The asymmetry between error and correct RTs is important for recovery parameters when using lower trial numbers
    ## For datasets with (overall) slow RTs and symmetric RT distributions for correct and incorrect answers,
    ## B0 and t0 start to trade off strongly (I believe this is also the issue that Dora has with using a racing Wald as a stop-signal
    ## accumulator?
    ## Luckily, the error RTs are in real data typically slower than correct RTs. This asymmetry provides an additional constraint
    ## to disentangle B0 from t0 effects.
    ## anyway, maybe try to first generate and recover *without* this asymmetry constraint?
    # if(mean(data$RT[data$R==1])  <= 1.1*mean(data$RT[data$R==2])) next
    break
  }
  trueParameters[[s]] <- c(B0=B0, V0=V0, wV=wV, wS=wS, aV=aV, t0=t0)
  data$s = s
  print(s)  # progress
  allData <- rbind(allData, data)
}
allData$s <- factor(allData$s)
save(allData, trueParameters, file = './parameterRecoveries/data/parameterRecovery-exp-sim-n360.RData')

load('./parameterRecoveries/data/parameterRecovery-exp-sim-n360.RData')
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

#### Sample  -----------------------------------------------------------------
doSample(data, pp.prior[[1]], pp.prior, nmcBurn=250, nmc=1000, nCores=30, restart=FALSE, fileName=file.path(samplesDir, fn))


# 
# # # ## check fit
load('./parameterRecoveries/samples/model-arw-RL-mag_data-parameterRecovery-exp-n360.RData')
load('./parameterRecoveries/data/parameterRecovery-exp-sim-n360.RData')
samples <- hsamples; rm(hsamples)
# plot.dmc(samples)
# 
posteriorSummary <- summary.dmc(samples)

trueParams <- data.frame(do.call(rbind, trueParameters))
posteriorMedian <- data.frame(do.call(rbind, lapply(posteriorSummary, function(x) x$quantiles[,3])))
posteriorMinCI <- data.frame(do.call(rbind, lapply(posteriorSummary, function(x) x$quantiles[,1])))
posteriorMaxCI <- data.frame(do.call(rbind, lapply(posteriorSummary, function(x) x$quantiles[,5])))

posteriorMinCI$B <- posteriorMinCI$B0
posteriorMinCI$aV <- pnorm(posteriorMinCI$aV)
posteriorMinCI$v0 <- posteriorMinCI$V0
posteriorMinCI <- posteriorMinCI[,colnames(trueParams)]

posteriorMaxCI$B <- posteriorMaxCI$B0
posteriorMaxCI$aV <- pnorm(posteriorMaxCI$aV)
posteriorMaxCI$V0 <- posteriorMaxCI$V0
posteriorMaxCI <- posteriorMaxCI[,colnames(trueParams)]

posteriorMedian$B <- posteriorMedian$B0
posteriorMedian$aV <- pnorm(posteriorMedian$aV)
posteriorMedian$v0 <- posteriorMedian$V0
posteriorMedian <- posteriorMedian[,colnames(trueParams)]

bad <- (posteriorMedian$V0 > 1.2*trueParams$V0) & (posteriorMedian$t0 < 0.75*trueParams$t0)

rmse <- function(x, y) {
  sqrt(mean((x-y)^2))
}

# pdf(file='figures/exp1-parrec.pdf', width=6, height=4)
par(oma=c(0,0,0,0), mar=c(3, 4, 2, 0.75) + 0.1, mfrow=c(2,3), mgp=c(2.75,.75,0), las=1)
plot.arrows=FALSE
for(parName in colnames(posteriorMedian)) {
  x <- trueParams[,parName]
  y <- posteriorMedian[,parName]
  cimin <- posteriorMinCI[,parName]
  cimax <- posteriorMaxCI[,parName]
  col = ifelse(bad, 2, 1)

  # Rename parameters for plot. Note that a = threshold
  if(parName == 'B' | parName == 'B0') { main <- expression(italic('a')); legend.pos='bottomright'}
  if(parName == 'v0' | parName == 'V0') { main <- expression(italic('V'[0])); legend.pos='bottomright'}
  if(parName == 'wV') { main <- expression(italic('w'['D'])); legend.pos='bottomright'}
  if(parName == 'wS') { main <- expression(italic('w'['S'])); legend.pos='bottomright'}
  if(parName == 't0') { main <- expression(italic('t'[0])); legend.pos='bottomright'}
  if(parName == 'aV') { main <- expression(italic(alpha)); legend.pos='bottomright'}

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
  tmp <- legend(legend.pos, c(" ", " "), bty='n', xjust=1,
                text.width = strwidth("RMSE = 0.03"))
  text(tmp$rect$left + tmp$rect$w, tmp$text$y,
       c(paste0('r = ', round(cor(x, y), 2)),
         paste0('RMSE = ', round(rmse(x,y), 2))), pos = 2)
}

head(allData)
badIdx <- which(bad)

allData$isBad <- allData$s %in% badIdx

allData$recovery <- ifelse(allData$isBad, 'Bad recovery', 'Good recovery')
allData$accuracy <- ifelse(allData$choiceIsHighP, 'Correct', 'Error')

par(mfrow=c(2,2), las=1, bty='l')
plot(0,0,type='n', xlim=c(.5,3.5), ylim=c(0.3,1), main='Mean RT', ylab='RT', xaxt='n')
mRTs <- aggregate(RT~recovery*s, allData, mean)
boxplot(RT~recovery, aggregate(RT~recovery*s, allData, mean), main='Mean RT', col=c(2, NA), add=TRUE)
boxplot(aggregate(rt~pp, dat, mean)$rt, add=TRUE, at=3, col='cornflowerblue')
axis(1, 3, c('Real data'))

plot(0,0,type='n', xlim=c(.5,3.5), ylim=c(0.3,1), main='Accuracy', ylab='Accuracy', xaxt='n')
meanAccs <- aggregate(choiceIsHighP~recovery*s, allData, mean)
boxplot(choiceIsHighP~recovery, meanAccs, main='Accuracy', ylab='Accuracy', col=c(2, NA), add=TRUE)
boxplot(aggregate(choiceIsHighP~pp, dat, mean)$choiceIsHighP, add=TRUE, at=3, col='cornflowerblue')
axis(1, 3, c('Real data'))

plot(0,0,type='n', xlim=c(.5,6.5), ylim=c(0.3,1.3), main='Mean RT by choice', ylab='RT', xaxt='n')
rtByChoice <- aggregate(RT~recovery*accuracy*s, allData, mean)
boxplot(RT~recovery*accuracy, rtByChoice, main='mean RT by choice', xaxt='n', col=c(2, NA, 2, NA), add=TRUE)
axis(1, 1:4, c('Bad Rec.\nCorrect', 'Good Rec.\nCorrect', 'Bad Rec.\nError', 'Good Rec\nError'), mgp=c(4,2,0))
boxplot(rt~choiceIsHighP, aggregate(rt~pp*choiceIsHighP, dat, mean), add=TRUE, at=5:6, col='cornflowerblue', xaxt='n')
axis(1, 5:6, c('Real data\nError', 'Real data\nCorrect'), mgp=c(4,2,0))

plot(0,0,type='n', xlim=c(.5,6.5), ylim=c(0,.8), main='SD RT by choice', ylab='SD RT', xaxt='n')
sdRts <- aggregate(RT~recovery*accuracy*s, allData, sd)
boxplot(RT~recovery*accuracy, sdRts, main='SD RT by choice', ylim=c(0,.8), xaxt='n', ylab='SD RT', col=c(2,NA), add=TRUE)
axis(1, 1:4, c('Bad Rec.\nCorrect', 'Good Rec.\nCorrect', 'Bad Rec.\nError', 'Good Rec\nError'), mgp=c(4,2,0))

boxplot(rt~choiceIsHighP, aggregate(rt~pp*choiceIsHighP, dat, sd), add=TRUE, at=5:6, col='cornflowerblue', xaxt='n')
axis(1, 5:6, c('Real data\nError', 'Real data\nCorrect'), mgp=c(4,2,0))

library(moments)
aggregate(RT~isBad, aggregate(RT~isBad*s, allData, skewness), mean)
mean(aggregate(rt~pp, dat, skewness)$rt)

aggregate(RT~isBad, aggregate(RT~isBad*s, allData, mean), mean)
aggregate(choiceIsHighP~isBad, aggregate(choiceIsHighP~isBad*s, allData, mean), mean)
aggregate(RT~isBad*choiceIsHighP, aggregate(RT~isBad*s*choiceIsHighP, allData, mean), mean)
aggregate(RT~isBad*choiceIsHighP, aggregate(RT~isBad*s*choiceIsHighP, allData, sd), mean)

aggregate(RT~isBad*choiceIsHighP, aggregate(RT~isBad*s*choiceIsHighP, allData, mean), mean)


aggregate(RT~isBad*choiceIsHighP, aggregate(RT~isBad*s*choiceIsHighP, allData, skewness), mean)

plot(0,0,xlim=c(0,7), ylim=c(0,4), main='Skew by choice', ylab='Skewness', xaxt='n', xlab='')
boxplot(RT~recovery*choiceIsHighP, aggregate(RT~recovery*s*choiceIsHighP, allData, skewness), add=TRUE, xaxt='n', col=c(2, NA))
boxplot(rt~choiceIsHighP, aggregate(rt~pp*choiceIsHighP, dat, skewness), add=TRUE, at=5:6, xaxt='n', col='cornflowerblue')
axis(1, 1:4, c('Bad Rec.\nCorrect', 'Good Rec.\nCorrect', 'Bad Rec.\nError', 'Good Rec\nError'), mgp=c(4,2,0))
axis(1, 5:6, c('Real data\nError', 'Real data\nCorrect'), mgp=c(4,2,0))



#boxplot(RT~isBad*choiceIsHighP, aggregate(RT~isBad*s*choiceIsHighP, allData, skewness), add=TRUE)
aggregate(rt~choiceIsHighP, aggregate(rt~pp*choiceIsHighP, dat, skewness), mean)



tmp <- posteriorMedian
tmp$bad <- bad

tmp$ratio <- tmp$wV/tmp$V0
tmp$s <- 1:nrow(tmp)
boxplot(ratio~bad, aggregate(ratio~bad*s, tmp, mean))
aggregate(ratio~bad, tmp, mean)

load('./samples/model-arw-RL-mag_data-exp1.RData')
summReal <- summary.dmc(hsamples)
posteriorMedianReal <- data.frame(do.call(rbind, lapply(summReal, function(x) x$quantiles[,3])))
mean(posteriorMedianReal$wV/posteriorMedianReal$V0)
posteriorMedianReal$ratio <- posteriorMedianReal$wV/posteriorMedianReal$V0

par(mfrow=c(1,1))
plot(0,0,type='n', xlim=c(0,4), ylim=c(0,5), xaxt='n', yaxt='n', main='Ratio wD/V0', ylab='Ratio')
boxplot(ratio~bad, aggregate(ratio~bad*s, tmp, mean), add=TRUE, xaxt='n')
boxplot(posteriorMedianReal$ratio, add=TRUE, at=3)
axis(1, 1:3, c('Good rec.', 'Bad rec.', 'Real data'))


# dev.off()
# 
# apply((trueParams < posteriorMaxCI) & (trueParams > posteriorMinCI), 2, mean)
