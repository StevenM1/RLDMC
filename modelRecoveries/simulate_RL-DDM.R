rm(list=ls())
source('./dmc/models/RW/dists.R')
source('./modelRecoveries/utils.R')  # simulation function here

# Can we recover? ---------------------------------------------------------
source ("dmc/dmc.R")
source('utils.R')
load_model ("ddm", "ddm-RL.R")

# some names for saving
samplesDir = './samples' # location of fit to *EMPIRICAL* data, not location to save fits to simulated data
modelName = 'ddm-RL'

#### Model set-up ----
load_model ("ddm", "ddm-RL.R")
model <- model.dmc(p.map=list(aV="1", SR='1', m='1',
                              t0='1', st0='1', d='1',
                              z='1', sz='1', a='1',
                              sv='1'),
                   match.map=list(M=list(s1=1, s1=2)),
                   constants=c(SR=-10, z=0.5, d=0, sv=0, sz=0, st0=0),
                   factors=list(S=c("s1")), 
                   responses=c("r1","r2"),
                   type="norm")


samples <- loadSamples(fn='model-ddm-RL_data-exp1', samplesDir=samplesDir)
samplesSummary <- summary.dmc(samples)
medianPars <- data.frame(do.call(rbind, lapply(samplesSummary, function(x) x$quantiles[,3])))

#### Data ----
# Simulate data with the median estimated parameters and the exact same design
setProbability <- list(c(0.2, 0.8),
                       c(0.3, 0.7),
                       c(0.35, 0.75),
                       c(0.4, 0.6))

for(dataset_n in 1:100) {
  dataName <- paste0('modelRecovery_gen-RL-DDM_ds-', dataset_n)
  fn = paste0('model-', modelName, '_data-', dataName)
  
  allData <- NULL
  trueParameters <- list()
  for(s in 1:nrow(medianPars)) {
    a = medianPars$a[s]
    m = medianPars$m[s]
    alpha = pnorm(medianPars$aV[s])
    t0 = medianPars$t0[s]
    data <- simulate.rlddm(nTrialsPerSet = 52, nSets=4, 
                           a=a, m=m, t0=t0, alpha=alpha,
                           setProbability = setProbability)
    trueParameters[[s]] <- c(a=a, m=m, alpha=alpha, t0=t0)
    data$s = s
    print(s)  # progress
    allData <- rbind(allData, data)
  }
  allData$s <- factor(allData$s)
  save(allData, trueParameters, file = file.path('./modelRecoveries/data/RL-DDM_simulated', paste0(dataName, '.RData')))
}
