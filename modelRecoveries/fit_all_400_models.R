# Parse command line arguments (if provided) ------------------------------
args = commandArgs(trailingOnly=TRUE)
if(length(args)>0) {
  generatingModels <- c(args[1])
} else {
  # manual run
  rm(list=ls())
  generatingModels <- c('arw-RL-mag', 'ddm-RL') #'ddm-RL-nonlinear-svszst0'
}
source('./dmc/models/RW/dists.R')

# Can we recover? ---------------------------------------------------------
source ("dmc/dmc.R")
source('utils.R')
load_model ("RW", "arw-RL-mag.R")

# some names for saving
modelName = 'arw-RL-mag'
samplesDir <- 'modelRecoveries/samples'

# Load data
for(dataset_n in 1:100) {
  for(generatingModel in generatingModels) {  # loop over data-generating models
    for(modelName in c('arw-RL-mag', 'ddm-RL')) {   # loop over fitting models
      dataName <- paste0('modelRecovery_gen-', generatingModel, '_ds-', dataset_n)
      fn = paste0('model-', modelName, '_data-', dataName)
      print(paste0('Dataset: ', dataset_n, 
                   '\nGenerating model: ', generatingModel,
                   '\nFit model: ', modelName))
      if(dataset_n == 2 & generatingModel == 'ddm-RL') next  # skip this one, already gave it 100 tries
#      if(file.exists(file.path(samplesDir, paste0(fn, '.RData')))) next

      ## load model
      if(modelName == 'arw-RL-mag') {
        source('./dmc/models/RW/dists.R')
        source ("dmc/dmc.R")
        source('utils.R')
        load_model ("RW", "arw-RL-mag.R")
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
      } else {
        source('./dmc/models/ddm/dists.R')
        source ("dmc/dmc.R")
        source('utils.R')
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
        
      }
      
      # Load data
      load(file.path('./modelRecoveries/data', paste0(generatingModel, '_simulated'), paste0(dataName, '.RData')))
    
      ### ugly work-around, sorry about this
      if(!'rt' %in% colnames(allData)) {
        allData$rt <- allData$RT
      } else {
        allData$RT <- allData$rt
      }
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
      pp.prior <- getPriors(attr(model, 'p.vector'))
      p.prior <- pp.prior[[1]]
      
      #### Sample  -----------------------------------------------------------------
      doSample(data, p.prior, pp.prior, nmcBurn=250, nmc=1500, nCores=12, restart=FALSE, 
               cut.converge = 1.05, fileName=file.path(samplesDir, fn), max.try=20)
    }
  }
}

# samples <- h.samples.dmc(nmc=1000, p.prior=p.prior, data=data, pp.prior=pp.prior)
# samples <- h.RUN.dmc(hsamples = samples, cores=20, cut.converge=1.03, thorough = TRUE, verbose=TRUE, saveFn=file.path(samplesDir, fn))
# save(samples, save=paste0(saveFn, '.RData'))

# 
# # # ## check fit
# load('./parameterRecoveries/samples/model-arw-RL-mag_data-parameterRecovery-exp1.RData')
# load('./parameterRecoveries/data/parameterRecovery-exp1.RData')
# samples <- hsamples; rm(hsamples)
# #plot.dmc(samples)
# 
# posteriorSummary <- summary.dmc(samples)
# 
# trueParams <- data.frame(do.call(rbind, trueParameters))
# posteriorMedian <- data.frame(do.call(rbind, lapply(posteriorSummary, function(x) x$quantiles[,3])))
# posteriorMinCI <- data.frame(do.call(rbind, lapply(posteriorSummary, function(x) x$quantiles[,1])))
# posteriorMaxCI <- data.frame(do.call(rbind, lapply(posteriorSummary, function(x) x$quantiles[,5])))
# 
# posteriorMinCI$B <- posteriorMinCI$B0
# posteriorMinCI$alpha <- pnorm(posteriorMinCI$aV)
# posteriorMinCI$v0 <- posteriorMinCI$V0
# posteriorMinCI <- posteriorMinCI[,colnames(trueParams)]
# 
# posteriorMaxCI$B <- posteriorMaxCI$B0
# posteriorMaxCI$alpha <- pnorm(posteriorMaxCI$aV)
# posteriorMaxCI$v0 <- posteriorMaxCI$V0
# posteriorMaxCI <- posteriorMaxCI[,colnames(trueParams)]
# 
# posteriorMedian$B <- posteriorMedian$B0
# posteriorMedian$alpha <- pnorm(posteriorMedian$aV)
# posteriorMedian$v0 <- posteriorMedian$V0
# posteriorMedian <- posteriorMedian[,colnames(trueParams)]
# 
# rmse <- function(x, y) {
#   sqrt(mean((x-y)^2))
# }
# 
# pdf(file='figures/exp1-parrec.pdf', width=6, height=4)
# par(oma=c(0,0,0,0), mar=c(3, 4, 2, 0.75) + 0.1, mfrow=c(2,3), mgp=c(2.75,.75,0), las=1)
# plot.arrows=FALSE
# for(parName in colnames(posteriorMedian)) {
#   x <- trueParams[,parName]
#   y <- posteriorMedian[,parName]
#   cimin <- posteriorMinCI[,parName]
#   cimax <- posteriorMaxCI[,parName]
#   col=1
#   
#   # Rename parameters for plot. Note that a = threshold
#   if(parName == 'B') { main <- expression(italic('a')); legend.pos='bottomright'}
#   if(parName == 'v0') { main <- expression(italic('V'[0])); legend.pos='bottomright'}
#   if(parName == 'wV') { main <- expression(italic('w'['D'])); legend.pos='bottomright'}
#   if(parName == 'wS') { main <- expression(italic('w'['S'])); legend.pos='bottomright'}
#   if(parName == 't0') { main <- expression(italic('t'[0])); legend.pos='bottomright'}
#   if(parName == 'alpha') { main <- expression(italic(alpha)); legend.pos='bottomright'}
#   
#   if(plot.arrows) {
#     ylim <- c(min(cimin), max(cimax))
#   } else {
#     ylim <- c(min(y), max(y))+c(-.1, .1)
#   }
#   plot(x, y, xlab='', ylab='Median posterior', main=main, type='n', ylim=ylim)
#   mtext('Data-generating value', side=1, line=2, cex=0.66)
#   abline(a=0, b=1)
#   points(x, y, col=col)
#   if(plot.arrows) arrows(x, cimin, x, cimax, length=0.05, angle=90, code=3, col=col)
#   tmp <- legend(legend.pos, c(" ", " "), bty='n', xjust=1, 
#                 text.width = strwidth("RMSE = 0.03"))
#   text(tmp$rect$left + tmp$rect$w, tmp$text$y,
#        c(paste0('r = ', round(cor(x, y), 2)), 
#          paste0('RMSE = ', round(rmse(x,y), 2))), pos = 2)
# }
# dev.off()
# 
# apply((trueParams < posteriorMaxCI) & (trueParams > posteriorMinCI), 2, mean)


