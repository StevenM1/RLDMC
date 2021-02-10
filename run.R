# Parse command line arguments (if provided) ------------------------------
args = commandArgs(trailingOnly=TRUE)
if(length(args)>0) {
  modelName <- args[1]
  dataName <- args[2]
  diagnosticsOnly <- as.logical(args[3])
} else {
  # manual run
  rm(list=ls())
  modelName <- 'arw-RL-WA' #'ddm-RL-nonlinear-svszst0'
  dataName <- 'expMA'     # NB: exp 3 = SAT; exp 2 = reversal learning
  diagnosticsOnly <- FALSE
}

# Source ------------------------------------------------------------------
source("dmc/dmc.R")
source('utils.R')
source('models.R')
samplesDir <- 'samples'

# get model specification
model <- setupModel(modelName)
fn <- paste0('model-', modelName, '_data-', dataName)

# Data ----
tmp <- loadData(dataName, removeBlock=NULL)
data <- data.model.dmc(tmp[['data']], model)
dat <- tmp[['dat']]

# append cvs by sub as well
for(sub in unique(dat$sub)) {
  d <- prepareForFitting(dat[dat$sub==sub,])
  attr(data[[sub]], 'cvs') <- d$outcomes
  attr(attr(data[[sub]], 'cvs'), 'VVchoiceIdx') <- d$VVchoiceIdx    # Trick for random.dmc
  attr(data[[sub]], 'VVchoiceIdx') <- d$VVchoiceIdx
  attr(data[[sub]], 'startingValues') <- d$values
  attr(data[[sub]], 'trialsToIgnore') <- d$trialsToIgnore
}

# Check model, only need this when you are developing
# likelihood.dmc(p.vector, data[[1]])
# simulate.dmc(p.vector, model=model, n=nrow(data[[1]]), adapt=TRUE, cvs=attr(data[[1]], 'cvs'))

# Priors ----
pp.prior <- getPriors(attr(model, 'p.vector'))

if(dataName == 'exp2' & modelName == 'ddm-RL-st02') {
  ## this model initially failed to converge; some participants got sz estimates of 1 (ie, sz covers the entire distance between both thresholds)
  # hence, we adjusted the upper bound of sz, after which the model did converge
  pp.prior[[1]]$sz$upper <- 0.5
}

# Sample  -----------------------------------------------------------------
if(!diagnosticsOnly) doSample(data, pp.prior[[1]], pp.prior, nmcBurn=250, nmc=1000, nCores=30, restart=FALSE, fileName=file.path(samplesDir, fn))

# Diagnostics -------------------------------------------------------------
runDiagnostics(fn, dat, plotFits=!grepl('softmax', fn))
