rm(list=ls())
source('./dmc/models/RW/dists.R')
source('./dmc/models/RW/WAdists.R')  # Win-All dists
source('./parameterRecoveries/utils.R')  # simulation function here

 
# Simulate ----------------------------------------------------------------
simulate.WA <- function(nTrialsPerSet=100, nSets=9,
                        setProbability=c(0.70, 0.15, .15),
                        V0=1, wV=2, wS=0.3, t0=0.3, alpha=0.05, B=2,
                        startValues=c(0.0, 0.0, 0.0)) {
  # simulates an entire experiment for a single subject
  # some 'experimental settings' can be adjusted:
  # - the number of trials per stimulus set;
  # - the number of stimulus sets;
  # - the reward contingencies per stimulus set. Can be a list of length nSets, which allows for specifying varying reward contingencies

  # model parameters are:
  # - v0: "urgency" effect on drift rate
  # - wV: weight of the difference in values on drift rate
  # - wS: weight of the sum in values on drift rate
  # (such that v_i = v0 + wV*(v_i-v_j) + wS*(v_i+v_j) )
  # - t0: non-decision time
  # - alpha: learning rate
  # - B: threshold
  # finally: startValues are the values at trial=0

  df <- data.frame(S=1, R=NA, RT=NA, stimulus_set=NA, reward=NA, predictionError=NA, value1=NA, value2=NA, value3=NA)[c()]

  # loop over stimulus sets
  for(set in 1:nSets) {
    if(is.list(setProbability)) {
      rewardProbability <- setProbability[[set]]
    } else {
      rewardProbability <- setProbability
    }
    values <- startValues

    # loop over trials
    for(trialN in 1:nTrialsPerSet) {
      thisTrial <- WArnd(x=c('V0'=V0, 'wV'=wV, 'wS'=wS, 't0'=t0, 'B'=B, 'A'=0), 3, Trials=1, Inputs = values)
      # determine drift rate
      # v1 <- v0 + wV*(values[1]-values[2]) + wS*(values[1]+values[2])  # "incorrect"
      # v2 <- v0 + wV*(values[2]-values[1]) + wS*(values[1]+values[2])  # "correct"

      # simulate single trial based on current values
#      thisTrial <- rWaldRaceSM(n=1, A=0, B=B, t0=t0, v=c(v1, v2), s=1)
      choice <- thisTrial$R #ifelse(thisTrial$R==1, 2, 1)   # shortcut

      # sample reward
      reward <- rbinom(1, 1, rewardProbability[choice])

      # calculate prediction error
      predictionError <- reward - values[choice]

      # update value of chosen option with Rescorla-Wagner rule
      values[choice] <- values[choice] + alpha*predictionError

      # off-load to dataframe
      df <- rbind(df, cbind('S'=1, thisTrial, stimulus_set=set, reward=reward,
                            predictionError=predictionError, value1=values[1], value2=values[2], value3=values[3],
                            trialNthisSet=trialN))
    }
  }
  return(df)
}
#
# # Illustration of simulation, some plots ----------------------------------
# # nTrialsPerSet = 50
# # df = simulate.WA(nTrialsPerSet = nTrialsPerSet, nSets = 9,
# #                    B=2, v0=2, wV=2)
# #
# # par(mfrow=c(2,1))
# # plot(0, 0, xlim=c(0, nTrialsPerSet), ylim=c(0,1), type='n', xlab='Trial', ylab='Value')
# # for(stimSet in unique(df$stimulus_set)) {
# #   lines(df$trialNthisSet[df$stimulus_set==stimSet], df$value1[df$stimulus_set==stimSet], col=stimSet)
# #   lines(df$trialNthisSet[df$stimulus_set==stimSet], df$value2[df$stimulus_set==stimSet], col=stimSet, lty=2)
# #   lines(df$trialNthisSet[df$stimulus_set==stimSet], df$value3[df$stimulus_set==stimSet], col=stimSet, lty=3)
# #   points(df$trialNthisSet[df$stimulus_set==stimSet], df$R[df$stimulus_set==stimSet]-1, col=stimSet, pch=4)
# # }
# #
# # plot(0, 0, xlim=c(0, nTrialsPerSet), ylim=c(-1,1), type='n', xlab='Trial', ylab='Prediction error')
# # for(stimSet in unique(df$stimulus_set)) {
# #   lines(df$trialNthisSet[df$stimulus_set==stimSet], df$predictionError[df$stimulus_set==stimSet], col=stimSet)
# # }
#
# # Can we recover? ---------------------------------------------------------
source ("dmc/dmc.R")
source('utils.R')
load_model ("RW", "arw-RL-mag.R")

# some names for saving
# modelName = 'arw-RL-mag'
# samplesDir <- 'parameterRecoveries/samples'
# dataName <- 'parameterRecovery-expMultiAlternative2'
# fn = paste0('model-', modelName, '_data-', dataName)

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

samples <- loadSamples(fn='model-arw-RL-WA_data-expMA', samplesDir='samples')
samplesSummary <- summary.dmc(samples)
medianPars <- data.frame(do.call(rbind, lapply(samplesSummary, function(x) x$quantiles[,3])))

#### Data ----
# Simulate data with the median estimated parameters and the exact same design
setProbability <- list(c(0.2, .2, 0.6),
                       c(0.15, 0.15, 0.7),
                       c(0.25, 0.25, 0.8),
                       c(0.3, .3, 0.7),
                       c(0.2, .2, 0.6),
                       c(0.15, 0.15, 0.7),
                       c(0.25, 0.25, 0.8),
                       c(0.3, .3, 0.7),
                       c(0.2, .2, 0.6),
                       c(0.15, 0.15, 0.7),
                       c(0.25, 0.25, 0.8),
                       c(0.3, .3, 0.7))
allData <- NULL
trueParameters <- list()
for(s in 1:nrow(medianPars)) {
  B = medianPars$B0[s]
  V0 = medianPars$V0[s]
  wV = medianPars$wV[s]
  wS = medianPars$wS[s]
  alpha = pnorm(medianPars$aV[s])
  t0 = medianPars$t0[s]
  data <- simulate.WA(nTrialsPerSet = 36, nSets=12,
                      B=B, V0=V0, wV=wV, wS=wS, t0=t0, alpha=alpha,
                      setProbability = setProbability)
  trueParameters[[s]] <- c(B=B, V0=V0, wV=wV, wS=wS, alpha=alpha, t0=t0)
  data$s = s
  print(s)  # progress
  allData <- rbind(allData, data)
}
allData$s <- factor(allData$s)
save(allData, trueParameters, file = './parameterRecoveries/data/parameterRecovery-expMultiAlternative3.RData')



# Fit ---------------------------------------------------------------------
source ("dmc/dmc.R")
source('utils.R')
load_model ("RW", "arw-RL-WA.R")
model <- model.dmc(p.map=list(A="1",t0="1",st0="1",s="1",
                              B0="1", wS='1',
                              SR="1", aV="1",
                              V0="1", wV="1"),
                   match.map=list(M=list(s1=1, s1=2, s1=3)),
                   constants=c(st0=0, s=1,
                               SR=-10, A=0),
                   factors=list(S=c("s1")), 
                   responses=c("r1","r2","r3"),
                   type="norm")
load('./parameterRecoveries/data/parameterRecovery-expMultiAlternative3.RData')
### ugly work-around, sorry about this
allData$rt <- allData$RT
allData$choiceIsHighP <- ifelse(allData$R==3, 1, 0)

cvs <- list()
choiceIdx <- list()
for(sub in unique(allData$s)) {
  d <- prepareForFitting(allData[allData$s==sub,])
  cvs[[sub]] <- d$outcomes
  choiceIdx[[sub]] <- d$VVchoiceIdx
}
data <- allData  # convert to DMC style data frame
data$S <- factor('s1')
data$R <- factor(data$R, levels=c(1,2,3), labels=c('r1', 'r2', 'r3'))
data <- data.model.dmc(data[,c('s', 'S', 'R', 'RT', 'stimulus_set')], model)

allData$choice <- allData$R
# append cvs by sub as well
for(sub in unique(allData$s)) {
  d <- prepareForFitting(allData[allData$s==sub,])
  attr(data[[sub]], 'cvs') <- d$outcomes
  attr(attr(data[[sub]], 'cvs'), 'VVchoiceIdx') <- d$VVchoiceIdx
  attr(data[[sub]], 'VVchoiceIdx') <- d$VVchoiceIdx
  attr(data[[sub]], 'startingValues') <- d$values
  attr(data[[sub]], 'trialsToIgnore') <- d$trialsToIgnore
}

# Check model, only need this when you are developing
# p.vector <- as.numeric(trueParameters[[1]])
# names(p.vector) <- names(trueParameters[[1]])
# p.vector <- c(p.vector,
#               'aV'=qnorm(as.numeric(p.vector[which(names(p.vector)=='alpha')])),
#               'V0'=as.numeric(p.vector[which(names(p.vector)=='v0')]),
#               'B0'=as.numeric(p.vector[which(names(p.vector)=='B')]))
# likelihood.dmc((p.vector), data[[1]])
# simulate.dmc(p.vector, model=model, n=nrow(data[[1]]), adapt=TRUE, cvs=attr(data[[1]], 'cvs'))

#### Priors ----
pp.prior <- getPriors(attr(model, 'p.vector'))
p.prior <- pp.prior[[1]]

modelName = 'arw-RL-WA'
samplesDir <- 'parameterRecoveries/samples'
dataName <- 'parameterRecovery-expMultiAlternative3'
fn = paste0('model-', modelName, '_data-', dataName)
#### Sample  -----------------------------------------------------------------
doSample(data, p.prior, pp.prior, nmcBurn=250, nCores=15, restart=FALSE, fileName=file.path(samplesDir, fn))


# samples <- h.samples.dmc(nmc=1000, p.prior=p.prior, data=data, pp.prior=pp.prior)
# samples <- h.RUN.dmc(hsamples = samples, cores=20, cut.converge=1.03, thorough = TRUE, verbose=TRUE, saveFn=file.path(samplesDir, fn))
# save(samples, save=paste0(saveFn, '.RData'))


# analyze -----------------------------------------------------------------
source ("dmc/dmc.R")
source('utils.R')
load('./parameterRecoveries/samples/model-arw-RL-WA_data-parameterRecovery-expMultiAlternative3.RData')
load('./parameterRecoveries/data/parameterRecovery-expMultiAlternative3.RData')
samples <- hsamples; rm(hsamples)
plot.dmc(samples)

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

pdf(file='figures/MA-parrec-samplepars.pdf', width=6, height=4)
par(oma=c(0,0,0,0), mar=c(3, 4, 2, 0.75) + 0.1, mfrow=c(2,3), mgp=c(2.75,.75,0), las=1)
plot.arrows=FALSE
for(parName in colnames(posteriorMedian)) {
  x <- trueParams[,parName]
  y <- posteriorMedian[,parName]
  cimin <- posteriorMinCI[,parName]
  cimax <- posteriorMaxCI[,parName]
  col=1
  
  # Rename parameters for plot. Note that a = threshold
  if(parName == 'B') { main <- expression(italic('a')); legend.pos='bottomright'}
  if(parName == 'V0') { main <- expression(italic('V'[0])); legend.pos='bottomright'}
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
  tmp <- legend(legend.pos, c(" ", " "), bty='n', xjust=1,
                text.width = strwidth("RMSE = 0.03"))
  text(tmp$rect$left + tmp$rect$w, tmp$text$y,
       c(paste0('r = ', round(cor(x, y), 2)),
         paste0('RMSE = ', round(rmse(x,y), 2))), pos = 2)
}
dev.off()
# 
# apply((trueParams < posteriorMaxCI) & (trueParams > posteriorMinCI), 2, mean)

