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

loadSamples <- function(fn, samplesDir=NULL) {
  if(!is.null(samplesDir)) fn <- file.path(samplesDir, fn)
  fnExt <- paste0(fn, '.RData')
  load(fnExt)
  return(allData)
}

#### Sample  -----------------------------------------------------------------
doSample(data, p.prior, pp.prior, nmcBurn=250, nCores=30, restart=TRUE, fileName=file.path(samplesDir, fn))
# samples <- h.samples.dmc(nmc=1000, p.prior=p.prior, data=data, pp.prior=pp.prior)
# samples <- h.RUN.dmc(hsamples = samples, cores=20, cut.converge=1.03, thorough = TRUE, verbose=TRUE, saveFn=file.path(samplesDir, fn))
# save(samples, save=paste0(saveFn, '.RData'))

