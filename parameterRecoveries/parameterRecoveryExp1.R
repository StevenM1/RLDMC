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

dataName <- 'parameterRecovery-exp1'
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

samples <- loadSamples(fn='model-arw-RL-mag_data-exp1', samplesDir='samples')
samplesSummary <- summary.dmc(samples)
medianPars <- data.frame(do.call(rbind, lapply(samplesSummary, function(x) x$quantiles[,3])))

#### Data ----
# Simulate data with the median estimated parameters and the exact same design
setProbability <- list(c(0.2, 0.8),
                       c(0.3, 0.7),
                       c(0.35, 0.75),
                       c(0.4, 0.6))
allData <- NULL
trueParameters <- list()
for(s in 1:nrow(medianPars)) {
  B = medianPars$B0[s]
  v0 = medianPars$V0[s]
  wV = medianPars$wV[s]
  wS = medianPars$wS[s]
  alpha = pnorm(medianPars$aV[s])
  t0 = medianPars$t0[s]
  data <- simulate.full(nTrialsPerSet = 52, nSets=4, 
                        B=B, v0=v0, wV=wV, wS=wS, t0=t0, alpha=alpha,
                        setProbability = setProbability)
  trueParameters[[s]] <- c(B=B, v0=v0, wV=wV, wS=wS, alpha=alpha, t0=t0)
  data$s = s
  print(s)  # progress
  allData <- rbind(allData, data)
}
allData$s <- factor(allData$s)
save(allData, trueParameters, file = './parameterRecoveries/data/parameterRecovery-exp1.RData')

load('./parameterRecoveries/data/parameterRecovery-exp1.RData')
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
doSample(data, p.prior, pp.prior, nmcBurn=1000, nCores=30, restart=FALSE, fileName=file.path(samplesDir, fn))
# samples <- h.samples.dmc(nmc=1000, p.prior=p.prior, data=data, pp.prior=pp.prior)
# samples <- h.RUN.dmc(hsamples = samples, cores=20, cut.converge=1.03, thorough = TRUE, verbose=TRUE, saveFn=file.path(samplesDir, fn))
# save(samples, save=paste0(saveFn, '.RData'))


# # ## check fit
load('./parameterRecoveries/samples/model-arw-RL-mag_data-parameterRecovery-exp1.RData')
load('./parameterRecoveries/data/parameterRecovery-exp1.RData')
samples <- hsamples; rm(hsamples)
#plot.dmc(samples)

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

pdf(file='figures/exp1-parrec.pdf', width=6, height=4)
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
  if(parName == 'v0') { main <- expression(italic('V'[0])); legend.pos='bottomright'}
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

apply((trueParams < posteriorMaxCI) & (trueParams > posteriorMinCI), 2, mean)
