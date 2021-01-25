#rm(list=ls())
source('./dmc/models/RW/dists.R')
source('./modelRecoveries/utils.R')  # simulation function here

# Can we recover? ---------------------------------------------------------
source ("dmc/dmc.R")
source('utils.R')
load_model ("RW", "arw-RL-mag.R")

# some names for saving
samplesDir = './samples' # location of fit to *EMPIRICAL* data, not location to save fits to simulated data
modelName = 'arw-RL-mag'

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

samples <- loadSamples(fn='model-arw-RL-mag_data-exp1', samplesDir=samplesDir)
samplesSummary <- summary.dmc(samples)
medianPars <- data.frame(do.call(rbind, lapply(samplesSummary, function(x) x$quantiles[,3])))

#### Data ----
# Simulate data with the median estimated parameters and the exact same design
setProbability <- list(c(0.2, 0.8),
                       c(0.3, 0.7),
                       c(0.35, 0.75),
                       c(0.4, 0.6))

for(dataset_n in 1:100) {
  dataName <- paste0('modelRecovery_gen-RL-ARD_ds-', dataset_n)
  fn = paste0('model-', modelName, '_data-', dataName)

  allData <- NULL
  trueParameters <- list()
  for(s in 1:nrow(medianPars)) {
    B = medianPars$B0[s]
    v0 = medianPars$V0[s]
    wV = medianPars$wV[s]
    wS = medianPars$wS[s]
    alpha = pnorm(medianPars$aV[s])
    t0 = medianPars$t0[s]
    data <- simulate.rlard(nTrialsPerSet = 52, nSets=4, 
                          B=B, v0=v0, wV=wV, wS=wS, t0=t0, alpha=alpha,
                          setProbability = setProbability)
    trueParameters[[s]] <- c(B=B, v0=v0, wV=wV, wS=wS, alpha=alpha, t0=t0)
    data$s = s
    print(s)  # progress
    allData <- rbind(allData, data)
  }
  allData$s <- factor(allData$s)
  save(allData, trueParameters, file = file.path('./modelRecoveries/data/RL-ARD_simulated', paste0(dataName, '.RData')))
}
