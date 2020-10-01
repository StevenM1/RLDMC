################################################################################
rm(list=ls())
source ("dmc/dmc.R")
source('utils.R')
load_model('ddm', 'ddm-RL.R')
load_model ("ddm", "ddm.R")
samplesDir <- 'samples'

# some names for saving
modelName = 'ddm-RL'
dataName <- 'exp1'
fn = paste0('model-', modelName, '_data-', dataName, 'aggregated')

#### Model ----
model <- model.dmc(p.map=list(v='1',
                              t0='1', st0='1', d='1',
                              z='1', sz='1', a='1',
                              sv='1'),
                   match.map=list(M=list(s1=1, s2=2)),
                   constants=c(z=0.5, d=0),
                   factors=list(S=c("s1", 's2')), 
                   responses=c("r1","r2"),
                   type="norm")

#### Data ----
tmp <- loadData(dataName, removeBlock=NULL)
dat <- tmp[['dat']]
dat$S <- factor(ifelse(dat$p_win_left > dat$p_win_right, 's1', 's2'), levels=c('s1', 's2'))
dat$R <- factor(ifelse(dat$choiceDirection=='left', 'r1', 'r2'), levels=c('r1', 'r2'))
dat$RT <- dat$rt
data <- data.model.dmc(dat[,c('S', 'R', 'RT')], model)

# Check model, only need this when you are developing
# likelihood.dmc(p.vector, data[[1]])
# simulate.dmc(p.vector, model=model, n=nrow(data[[1]]), adapt=TRUE, cvs=attr(data[[1]], 'cvs'))

# Priors ----
pp.prior <- getPriors(attr(model, 'p.vector'))

# Sample  -----------------------------------------------------------------
burn <- run.unstuck.dmc(samples=samples.dmc(250, p.prior=pp.prior[[1]], data=data), p.migrate = 0.05, cores=25)
samples <- RUN.dmc(samples = burn, nmc=500, cores=25)

pp <- post.predict.dmc(burn)
plot.pp.dmc(pp)
#doSample(data, pp.prior[[1]], pp.prior, nmcBurn=1000, nCores=15, restart=FALSE, fileName=file.path(samplesDir, fn))

# Diagnostics -------------------------------------------------------------
runDiagnostics(fn, dat)









###### old stuff below
# ### Plots & checks ----------
# samples <- loadSamples(fn, samplesDir)
# #plot.dmc(samples, hyper=TRUE)
# 
# bpics = h.IC.dmc(samples)
# h.gelman.diag.dmc(samples)
# 
# summ <- summary.dmc(samples)
# medians <- data.frame(do.call(rbind, lapply(summ, function(x) x$quantiles[,3])))
# medians$aV <- pnorm(medians$aV)
# round(apply(medians, 2, mean), 5)
# round(apply(medians, 2, sd), 5)
# 
# 
# 
# #samples <- h.samples.dmc(nmc=0,add=TRUE,samples=samples,remove=1:400)
# plot.dmc(samples, hyper=TRUE, density=TRUE, layout=c(4,4))
# # # simulate posterior predictives
# # 
# # simulate posterior predictives
# pp = h.post.predict.dmc(samples = samples, adapt=TRUE,save.simulation = TRUE, cores=30)
# ppNoSim <- h.pp.summary(pp, samples=samples)
# library(moments)
# library(snowfall)
# #### Append stimulus set info to data & model --------
# nBins <- 10
# #dat$ease <- as.factor(abs(dat$high_stim_prob-dat$low_stim_prob))
# pp2 <- lapply(1:length(pp), addStimSetInfo, input=pp, orig_dat=dat)
# data2 <- lapply(1:length(data), addStimSetInfo, input=data, orig_dat=dat)
# if(!sfIsRunning()) sfInit(parallel=TRUE, cpus =30); sfLibrary(moments)
# pp3 <- sfLapply(pp2, calculateByBin)
# data3 <- lapply(data2, calculateByBin)
#  
# ##
# meanRTsBySub <- list(lapply(data3, function(x) aggregate(RT~bin, x, mean)), lapply(pp3, function(x) aggregate(RT~bin*reps, x, mean)))
# meanAccsBySub <- list(lapply(data3, function(x) aggregate(acc~bin, x, mean)), lapply(pp3, function(x) aggregate(acc~bin*reps, x, mean)))
# 
# meanRTsOverTime <- list(getDescriptives(data3, dep.var='RT', attr.name='RTsOverBins', id.var1='~bin', id.var2=NULL), #*s', id.var2="~bin"),
#                         getDescriptives(pp3, dep.var='RT', attr.name='RTsOverBins'))
# meanAccOverTime <- list(getDescriptives(data3, dep.var='acc', attr.name='AccOverBins', id.var1='~bin', id.var2=NULL), #*s', id.var2="~bin"),
#                         getDescriptives(pp3, dep.var='acc', attr.name='AccOverBins'))
# 
# meanRTsByEase <- list(getDescriptives(data3, dep.var='RT', attr.name='RTsByEase', id.var1='~bin*ease', id.var2 = NULL), #*s', id.var2="~bin*ease"),
#                       getDescriptives(pp3, dep.var='RT', attr.name='RTsByEase', id.var1='~reps*bin*ease', id.var2 = NULL)) #*s', id.var2="~reps*bin*ease"))
# meanAccByEase <- list(getDescriptives(data3, dep.var='acc', attr.name='AccByEase', id.var1='~bin*ease', id.var2=NULL),#, id.var2="~bin*ease"),
#                       getDescriptives(pp3, dep.var='acc', attr.name='AccByEase', id.var1='~reps*bin*ease', id.var2=NULL)) #', id.var2="~reps*bin*ease"))
# 
# # Start plotting
# pdf(file=file.path('plots', paste0(fn, '.pdf')), width=8, height=4)
# plot.dmc(samples, hyper=TRUE, density=TRUE, layout=c(4,4))
# 
# ## Overall 3-panel plot ------------------------------------------------------------
# # nf <- layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE), respect = TRUE)
# par(mar=c(5, 4, 2, 2) + 0.1, mfrow=c(1,3))
# plot.pp.dmc(ppNoSim, style='cdf', x.min.max = c(0, 3), no.layout = TRUE, do.main=FALSE,
#             fits.lcol='blue', data.lwd.mult = 2, lwd=c(2,2)); title('Overall def. CDFs')
# plotDataPPBins(data=meanRTsOverTime[[1]], pp=meanRTsOverTime[[2]]); title('Mean RT over time')
# plotDataPPBins(data=meanAccOverTime[[1]], pp=meanAccOverTime[[2]], dep.var='acc'); title('Accuracy over time')
# 
# ## 3-panel plot by ease ------------------------------------------------------------
# for(ease in unique(meanRTsByEase[[1]]$ease)) {
#   pp2Subset <- lapply(1:length(pp2), function(x) pp2[[x]][pp2[[x]]$ease==ease,])
#   dataSubset <- do.call(rbind, lapply(1:length(data2), function(x) data2[[x]][data2[[x]]$ease==ease,]))
#   ppNoSimSubset <- h.pp.summary(pp2Subset, samples=samples, dat=dataSubset, ignore.subjects=TRUE)
#   plot.pp.dmc(ppNoSimSubset, style='cdf', no.layout = TRUE, do.main=FALSE,
#               fits.lcol='blue', data.lwd.mult = 2, lwd=c(2,2)); title(paste0('Ease ', ease))
#   
#   idxD = meanRTsByEase[[1]]$ease == ease
#   idxM = meanRTsByEase[[2]]$ease == ease
#   plotDataPPBins(data=meanRTsByEase[[1]][idxD,], pp=meanRTsByEase[[2]][idxM,], draw.legend=FALSE); title('Mean RT over time')
#   plotDataPPBins(data=meanAccByEase[[1]][idxD,], pp=meanAccByEase[[2]][idxM,], dep.var='acc', draw.legend=FALSE); title('Mean RT over time')
# }
# 
# # Overall CDFs by bin ------------------------------------------------------------
# layout(1)
# par(mfrow=c(2,3))
# for(i in 1:nBins) {
#   pp2Subset <- lapply(1:length(pp2), function(x) pp2[[x]][pp2[[x]]$bin==i,])
#   dataSubset <- do.call(rbind, lapply(1:length(data2), function(x) data2[[x]][data2[[x]]$bin==i,]))
#   ppNoSimSubset <- h.pp.summary(pp2Subset, samples=samples, dat=dataSubset, ignore.subjects=TRUE)
#   plot.pp.dmc(ppNoSimSubset, style='cdf', no.layout = TRUE, do.main=FALSE, aname='',
#               fits.lcol='blue', data.lwd.mult = 2, lwd=c(2,2)); title(paste0('Bin ', i))
# }
# 
# # Overall change in RT/Accuracy by ease ------------------------------------------------------------
# par(mfcol=c(2,4))
# for(ease in unique(meanRTsByEase[[1]]$ease)) {
#   idxD = meanRTsByEase[[1]]$ease == ease
#   idxM = meanRTsByEase[[2]]$ease == ease
#   plotDataPPBins(data=meanRTsByEase[[1]][idxD,], pp=meanRTsByEase[[2]][idxM,], draw.legend=FALSE, ylim=c(0.6, 0.8))
#   title(ease)
#   plotDataPPBins(data=meanAccByEase[[1]][idxD,], pp=meanAccByEase[[2]][idxM,], dep.var='acc', draw.legend=FALSE, ylim=c(0.4, 1))
# }
# 
# ### Lastly, plot by subject
# # nf <- layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE), respect = TRUE)
# par(mar=c(5, 4, 2, 2) + 0.1, mfrow=c(1,3))
# for(i in 1:length(ppNoSim)) {
#   plot.pp.dmc(ppNoSim[[i]], style='cdf', x.min.max = c(0, 3), no.layout = TRUE, do.main=FALSE,
#               fits.lcol='blue', data.lwd.mult = 2, lwd=c(2,2)); title(paste0('Participant ', i)) #Overall def. CDFs')
#   plotDataPPBins(data=meanRTsBySub[[1]][[i]], pp=meanRTsBySub[[2]][[i]]); title('Mean RT over time')
#   plotDataPPBins(data=meanAccsBySub[[1]][[i]], pp=meanAccsBySub[[2]][[i]], dep.var='acc'); title('Accuracy over time')
# }
# 
# dev.off()
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# q10RTsOverTime <- list(getDescriptives(data3, dep.var='RT.10.', attr.name='qRTsCorrect', id.var1='~bin', id.var2=NULL), #*s', id.var2="~bin"),
#                        getDescriptives(pp3, dep.var='RT.10.', attr.name='qRTsCorrect'))
# q50RTsOverTime <- list(getDescriptives(data3, dep.var='RT.50.', attr.name='qRTsCorrect', id.var1='~bin', id.var2=NULL), #*s', id.var2="~bin"),
#                        getDescriptives(pp3, dep.var='RT.50.', attr.name='qRTsCorrect'))
# q90RTsOverTime <- list(getDescriptives(data3, dep.var='RT.90.', attr.name='qRTsCorrect', id.var1='~bin', id.var2=NULL), #*s', id.var2="~bin"),
#                        getDescriptives(pp3, dep.var='RT.90.', attr.name='qRTsCorrect'))
# q10RTsOverTimeE <- list(getDescriptives(data3, dep.var='RT.10.', attr.name='qRTsError', id.var1='~bin', id.var2=NULL), #*s', id.var2="~bin"),
#                         getDescriptives(pp3, dep.var='RT.10.', attr.name='qRTsError'))
# q50RTsOverTimeE <- list(getDescriptives(data3, dep.var='RT.50.', attr.name='qRTsError', id.var1='~bin', id.var2=NULL), #*s', id.var2="~bin"),
#                         getDescriptives(pp3, dep.var='RT.50.', attr.name='qRTsError'))
# q90RTsOverTimeE <- list(getDescriptives(data3, dep.var='RT.90.', attr.name='qRTsError', id.var1='~bin', id.var2=NULL), #*s', id.var2="~bin"),
#                         getDescriptives(pp3, dep.var='RT.90.', attr.name='qRTsError'))
# 
# meanAccOverTime <- list(getDescriptives(data3, dep.var='acc', attr.name='AccOverBins', id.var1='~bin', id.var2=NULL), #*s', id.var2="~bin"),
#                         getDescriptives(pp3, dep.var='acc', attr.name='AccOverBins'))
# 
# q10RTsByEase <- list(getDescriptives(data3, dep.var='RT.10.', attr.name='qRTsCorrectByEase', id.var1='~bin*ease', id.var2=NULL), #*s', id.var2="~bin*ease"),
#                      getDescriptives(pp3, dep.var='RT.10.', attr.name='qRTsCorrectByEase', id.var1='~reps*bin*ease', id.var2=NULL)) #*s', id.var2="~reps*bin*ease"))
# q50RTsByEase <- list(getDescriptives(data3, dep.var='RT.50.', attr.name='qRTsCorrectByEase', id.var1='~bin*ease', id.var2=NULL), #*s', id.var2="~bin*ease"),
#                      getDescriptives(pp3, dep.var='RT.50.', attr.name='qRTsCorrectByEase', id.var1='~reps*bin*ease', id.var2=NULL)) #', id.var2="~reps*bin*ease"))
# q90RTsByEase <- list(getDescriptives(data3, dep.var='RT.90.', attr.name='qRTsCorrectByEase', id.var1='~bin*ease', id.var2=NULL), #', id.var2="~bin*ease"),
#                      getDescriptives(pp3, dep.var='RT.90.', attr.name='qRTsCorrectByEase', id.var1='~reps*bin*ease', id.var2=NULL)) #', id.var2="~reps*bin*ease"))
# q10RTsByEaseE <- list(getDescriptives(data3, dep.var='RT.10.', attr.name='qRTsErrorByEase', id.var1='~bin*ease', id.var2=NULL), #*s', id.var2="~bin*ease"),
#                       getDescriptives(pp3, dep.var='RT.10.', attr.name='qRTsErrorByEase', id.var1='~reps*bin*ease', id.var2=NULL))#, #)*s', id.var2="~reps*bin*ease"))
# q50RTsByEaseE <- list(getDescriptives(data3, dep.var='RT.50.', attr.name='qRTsErrorByEase', id.var1='~bin*ease', id.var2=NULL), #*s', id.var2="~bin*ease"),
#                       getDescriptives(pp3, dep.var='RT.50.', attr.name='qRTsErrorByEase', id.var1='~reps*bin*ease', id.var2=NULL))#, #*s', id.var2="~reps*bin*ease"))
# q90RTsByEaseE <- list(getDescriptives(data3, dep.var='RT.90.', attr.name='qRTsErrorByEase', id.var1='~bin*ease', id.var2=NULL), #*s', id.var2="~bin*ease"),
#                       getDescriptives(pp3, dep.var='RT.90.', attr.name='qRTsErrorByEase', id.var1='~reps*bin*ease', id.var2=NULL))#, #*s', id.var2="~reps*bin*ease"))
# meanAccByEase <- list(getDescriptives(data3, dep.var='acc', attr.name='AccByEase', id.var1='~bin*ease', id.var2=NULL), #*s', id.var2="~bin*ease"),
#                       getDescriptives(pp3, dep.var='acc', attr.name='AccByEase', id.var1='~reps*bin*ease', id.var2=NULL)) #*s', id.var2="~reps*bin*ease"))
# 
# 
# # # # Experiment 1: Difficulty effect plot ------------------------------------
# if(dataName == 'exp1') {
#   pdf(file=paste0('./figures/exp1_difficulty_', modelName, '-QQ.pdf'), width=6, height=6)
#   par(oma=c(3,1,1,0), mar=c(0, 4, 1, 0.5) + 0.1, mfrow=c(4,3), mgp=c(2.75,.75,0), las=1)
#   i <- 0
#   for(ease in unique(meanAccByEase[[1]]$ease)) {
#     i <- i+1
#     # pp2Subset <- lapply(1:length(pp2), function(x) pp2[[x]][pp2[[x]]$ease==ease,])
#     # dataSubset <- do.call(rbind, lapply(1:length(data2), function(x) data2[[x]][data2[[x]]$ease==ease,]))
#     # ppNoSimSubset <- h.pp.summary(pp2Subset, samples=samples, dat=dataSubset, ignore.subjects=TRUE)
#   
#     # plot.pp.dmc(ppNoSimSubset, style='cdf', x.min.max = c(0.3, 2), ylim=c(0, 0.95),
#     #             no.layout = TRUE, do.main=FALSE,
#     #             fits.lcol='blue', data.lwd.mult = 2,
#     #             lwds=c(2,2), ltys=c(2,1), pos=NA, xlab = '',
#     #             xaxt='n', axvlines=c(0.5, 1, 1.5, 2),
#     #             axhlines=seq(0, 1, .2))
#     idxD = meanAccByEase[[1]]$ease == ease
#     idxM = meanAccByEase[[2]]$ease == ease
#   
#     plotDataPPBins(data=meanAccByEase[[1]][idxD,], pp=meanAccByEase[[2]][idxM,],
#                    xaxt=ifelse(i==4, 's', 'n'),
#                    dep.var='acc', ylab='Accuracy', xlab = '',
#                    legend.pos='bottomright', ylim=c(0.5, 0.95), hline.by=0.1)
#     if(i == 4) {
#       mtext('Trial bin', side=1, cex=.66, line=2)
#     } else {
#       axis(1, at=seq(2, 10, 2), labels=rep(NA, 5))
#     }
#     if(i == 1) title('Accuracy over time')
#     mtext(ifelse(round(ease,2) == 0.6, "0.8/0.2", ifelse(round(ease,2)==0.4, "0.7/0.3", ifelse(round(ease,2)==0.3, "0.65/0.35", "0.6/0.4"))),
#           side=2, cex=0.66*1.2, las=0, line=4, font=2)
#   
#     ##
#     plotDataPPBins(data=q10RTsByEase[[1]][idxD,], pp=q10RTsByEase[[2]][idxM,], dep.var='RT.10.', ylim=c(.4, 1.2), xaxt='n', ylab='RT (s)')#; title('Correct RTs over time')
#     plotDataPPBins(data=q50RTsByEase[[1]][idxD,], pp=q50RTsByEase[[2]][idxM,], dep.var='RT.50.', plot.new = FALSE, draw.legend=FALSE)
#     plotDataPPBins(data=q90RTsByEase[[1]][idxD,], pp=q90RTsByEase[[2]][idxM,], dep.var='RT.90.', plot.new = FALSE, draw.legend=FALSE)
#     if(i == 1) title('Correct RTs over time')
#     if(i == 4) {
#       axis(1, at=seq(2, 10, 2))
#       mtext('Trial bin', side=1, cex=.66, line=2)
#     } else {
#       axis(1, at=seq(2, 10, 2), labels=rep(NA, 5))
#     }
#   
#   
#     plotDataPPBins(data=q10RTsByEaseE[[1]][idxD,], pp=q10RTsByEaseE[[2]][idxM,], dep.var='RT.10.', ylim=c(.4, 1.2), xaxt='n', ylab='RT (s)')#; title('Error RTs over time')
#     plotDataPPBins(data=q50RTsByEaseE[[1]][idxD,], pp=q50RTsByEaseE[[2]][idxM,], dep.var='RT.50.', plot.new = FALSE, draw.legend=FALSE)
#     plotDataPPBins(data=q90RTsByEaseE[[1]][idxD,], pp=q90RTsByEaseE[[2]][idxM,], dep.var='RT.90.', plot.new = FALSE, draw.legend=FALSE)
#     # mtext(ifelse(i==1, 'RL-DDM', ifelse(i==2, 'RL-RD', 'RL-ARD')), side=2, cex=0.66*1.2, las=0, line=4, font=2)
#     if(i == 1) title('Error RTs over time')
#     if(i == 4) {
#       axis(1, at=seq(2, 10, 2))
#       mtext('Trial bin', side=1, cex=.66, line=2)
#     } else {
#       axis(1, at=seq(2, 10, 2), labels=rep(NA, 5))
#     }
#   
#     # idxD = meanRTsByEase[[1]]$ease == ease
#     # idxM = meanRTsByEase[[2]]$ease == ease
#     # plotDataPPBins(data=meanRTsByEase[[1]][idxD,], pp=meanRTsByEase[[2]][idxM,],  ylab='RT (s)',
#     #                xaxt=ifelse(i==4, 's', 'n'),
#     #                xlab='', ylim=c(0.55, 0.875), hline.by=0.05)
#     # if(i == 1) title('Mean RT over time')
#     # if(i == 4) {
#     #   mtext('Trial bin', side=1, cex=.66, line=2)
#     # } else {
#     #   axis(1, at=seq(2, 10, 2), labels=rep(NA, 5))
#     # }
#   
#     # plotDataPPBins(data=meanAccByEase[[1]][idxD,], pp=meanAccByEase[[2]][idxM,],
#     #                xaxt=ifelse(i==4, 's', 'n'),
#     #                dep.var='acc', ylab='Accuracy', xlab = '',
#     #                legend.pos='bottomright', ylim=c(0.5, 0.95), hline.by=0.1)
#     # if(i == 4) {
#     #   mtext('Trial bin', side=1, cex=.66, line=2)
#     # } else {
#     #   axis(1, at=seq(2, 10, 2), labels=rep(NA, 5))
#     # }
#     # if(i == 1) title('Accuracy over time')
#   }
#   dev.off()
# }
# 
# 
# if(dataName == 'exp2') {
#   dat$bin <- dat$trialNreversal
#   pp2 <- lapply(1:length(pp), addStimSetInfo, input=pp, orig_dat=dat, addColumns='bin')
#   data2 <- lapply(1:length(data), addStimSetInfo, input=data, orig_dat=dat, addColumns='bin')
#   if(!sfIsRunning()) sfInit(parallel=TRUE, cpus = 30); sfLibrary(moments)
#   pp3 <- sfLapply(pp2, calculateByBin)
#   data3 <- lapply(data2, calculateByBin)
#   
#   
#   q10RTsOverTrials <- list(getDescriptives(data3, dep.var='RT.10.', attr.name='qRTs', id.var1='~bin', id.var2=NULL),
#                          getDescriptives(pp3, dep.var='RT.10.', attr.name='qRTs'))
#   q50RTsOverTrials <- list(getDescriptives(data3, dep.var='RT.50.', attr.name='qRTs', id.var1='~bin', id.var2=NULL),
#                          getDescriptives(pp3, dep.var='RT.50.', attr.name='qRTs'))
#   q90RTsOverTrials <- list(getDescriptives(data3, dep.var='RT.90.', attr.name='qRTs', id.var1='~bin', id.var2=NULL),
#                          getDescriptives(pp3, dep.var='RT.90.', attr.name='qRTs'))
#   
#   meanAccOverTrials <- list(getDescriptives(data3, dep.var='acc', attr.name='AccOverBins', id.var1='~bin', id.var2=NULL),
#                           getDescriptives(pp3, dep.var='acc', attr.name='AccOverBins'))
#   
#   
#   
#   #####
#   pdf(file=paste0('./figures/exp2_reversal_', modelName, '-QQ.pdf'), width=6, height=4.5)
#   par(oma=c(3,1,1,0), mar=c(3, 4, 2, 0.5) + 0.1, mfrow=c(2,3), mgp=c(2.75,.75,0), las=1)
#   layout(matrix(c(1,1,2,2,3,3,4,4,4,5,5,5), byrow=TRUE, ncol=6))
#   i <- 0
#   plotDataPPBins(data=meanAccOverTime[[1]], pp=meanAccOverTime[[2]], dep.var='acc', ylab='Choice proportion A'); title('Accuracy over time')
#   mtext('Trial bin', side=1, cex=.66, line=2)
#   
#   plotDataPPBins(data=q10RTsOverTime[[1]], pp=q10RTsOverTime[[2]], dep.var='RT.10.', ylim=c(.3, 1.2), ylab='RT (s)'); title('Choice A RTs')
#   plotDataPPBins(data=q50RTsOverTime[[1]], pp=q50RTsOverTime[[2]], dep.var='RT.50.', plot.new = FALSE, draw.legend=FALSE)
#   plotDataPPBins(data=q90RTsOverTime[[1]], pp=q90RTsOverTime[[2]], dep.var='RT.90.', plot.new = FALSE, draw.legend=FALSE)
#   mtext('Trial bin', side=1, cex=.66, line=2)
#   
#   plotDataPPBins(data=q10RTsOverTimeE[[1]], pp=q10RTsOverTimeE[[2]], dep.var='RT.10.', ylim=c(.3, 1.2), ylab='RT (s)'); title('Choice B RTs')
#   plotDataPPBins(data=q50RTsOverTimeE[[1]], pp=q50RTsOverTimeE[[2]], dep.var='RT.50.', plot.new = FALSE, draw.legend=FALSE)
#   plotDataPPBins(data=q90RTsOverTimeE[[1]], pp=q90RTsOverTimeE[[2]], dep.var='RT.90.', plot.new = FALSE, draw.legend=FALSE)
#   mtext('Trial bin', side=1, cex=.66, line=2)
#   
#   
#   ##
#   plotDataPPBins(data=q10RTsOverTrials[[1]], pp=q10RTsOverTrials[[2]], dep.var='RT.10.', ylim=c(.4, 1.1), ylab='RT (s)',
#                  xaxt='s', xlim=c(-60, 40), plot.model.points = FALSE, data.lwd=1.5, data.cex=1.5,
#                  xlab='', hline.by=0.05, axvlines=seq(-50, 50, 10))
#   
#   plotDataPPBins(data=q50RTsOverTrials[[1]], pp=q50RTsOverTrials[[2]], dep.var='RT.50.', plot.new = FALSE, draw.legend=FALSE, plot.model.points = FALSE, data.lwd=1.5, data.cex=1.5)
#   plotDataPPBins(data=q90RTsOverTrials[[1]], pp=q90RTsOverTrials[[2]], dep.var='RT.90.', plot.new = FALSE, draw.legend=FALSE, plot.model.points = FALSE, data.lwd=1.5, data.cex=1.5)
#   abline(v=0, lty=2)
#   title('Choice A RTs after reversal')
#   mtext('Trial (relative to reversal point)', side=1, cex=.66, line=2)
#   
#   
#   plotDataPPBins(data=meanAccOverTrials[[1]], pp=meanAccOverTrials[[2]],
#                  xaxt='s', xlim=c(-60, 40), plot.model.points=FALSE,
#                  dep.var='acc', ylab=expression('Proportion choice A'),
#                  xlab = '', data.lwd=1.5, data.cex=1.5,
#                  legend.pos='topright', ylim=c(0.25, 0.85), hline.by=0.05, axvlines=seq(-50, 50, 10))
#   title('Choices after reversal')
#   axis(1, at=seq(-50, 50, 10), labels=rep(NA, 5))
#   abline(v=0, lty=2)
#   mtext('Trial (relative to reversal point)', side=1, cex=.66, line=2)
#   dev.off()
# }
# 


# ## Plot expected values over time ------------
# allEV1OverTimeM <- do.call(rbind, (lapply(1:length(pp3), function(x) {tmp <- attr(pp3[[x]], 'SR.r1OverBins'); tmp$s <- x; tmp})))
# allEV2OverTimeM <- do.call(rbind, (lapply(1:length(pp3), function(x) {tmp <- attr(pp3[[x]], 'SR.r2OverBins'); tmp$s <- x; tmp})))
# meanEV1OverTimeM <- aggregate(SR.r1~reps*bin*ease, aggregate(SR.r1~reps*bin*ease*s, allEV1OverTimeM, mean), mean)
# meanEV2OverTimeM <- aggregate(SR.r2~reps*bin*ease, aggregate(SR.r2~reps*bin*ease*s, allEV2OverTimeM, mean), mean)
# 
# allVOverTimeM <- do.call(rbind, (lapply(1:length(pp3), function(x) {tmp <- attr(pp3[[x]], 'vOverBins'); tmp$s <- x; tmp})))
# meanVOverTimeM <- aggregate(v~reps*bin*ease, aggregate(v~reps*bin*ease*s, allVOverTimeM, mean), mean)
# 
# allPEOverTimeM <- do.call(rbind, (lapply(1:length(pp3), function(x) {tmp <- attr(pp3[[x]], 'PEsOverBins'); tmp$s <- x; tmp})))
# meanPEsOverTimeM <- aggregate(PEs~reps*bin*ease, aggregate(PEs~reps*bin*ease*s, allPEOverTimeM, mean), mean)
# 
# # allRisk1OverTimeM <- do.call(rbind, (lapply(1:length(pp3), function(x) {tmp <- attr(pp3[[x]], 'risk.r1OverBins'); tmp$s <- x; tmp})))
# # risk1OverTimeM <- aggregate(risk.r1~reps*bin*ease, aggregate(risk.r1~reps*bin*ease*s, allRisk1OverTimeM, mean), mean)
# # allRisk2OverTimeM <- do.call(rbind, (lapply(1:length(pp3), function(x) {tmp <- attr(pp3[[x]], 'risk.r2OverBins'); tmp$s <- x; tmp})))
# # risk2OverTimeM <- aggregate(risk.r2~reps*bin*ease, aggregate(risk.r2~reps*bin*ease*s, allRisk2OverTimeM, mean), mean)
# 
# 
# allb1OverTimeM <- do.call(rbind, (lapply(1:length(pp3), function(x) {tmp <- attr(pp3[[x]], 'b.r1OverBins'); tmp$s <- x; tmp})))
# allb2OverTimeM <- do.call(rbind, (lapply(1:length(pp3), function(x) {tmp <- attr(pp3[[x]], 'b.r2OverBins'); tmp$s <- x; tmp})))
# meanb1OverTimeM <- aggregate(b.r1~reps*bin*ease, aggregate(b.r1~reps*bin*ease*s, allb1OverTimeM, mean), mean)
# meanb2OverTimeM <- aggregate(b.r2~reps*bin*ease, aggregate(b.r2~reps*bin*ease*s, allb2OverTimeM, mean), mean)
# 
# 
# 
# par(mfrow=c(1,1))
# plot(0,0, type='n', xlim=range(meanEV1OverTimeM$bin)+c(-.5, .5), ylim=c(0.3, 1), xlab='Bin', ylab='mean EV1')
# for(ease in unique(meanEV1OverTimeM$ease)) {
#   toplot <- meanEV1OverTimeM[meanEV1OverTimeM$ease==ease,]
#   points(toplot$bin, 
#          toplot$SR.r1, 
#          col=ifelse(ease=="0.6", 1, ifelse(ease=="0.4", 2, ifelse(ease=="0.3", 4, 5))), pch=20, cex=.25)
#   toplot <- meanEV2OverTimeM[meanEV2OverTimeM$ease==ease,]
#   points(toplot$bin, 
#          toplot$SR.r2, 
#          col=ifelse(ease=="0.6", 1, ifelse(ease=="0.4", 2, ifelse(ease=="0.3", 4, 5))), pch=20, cex=.25)
# }
# 
# tmp <- aggregate(SR.r1~bin, meanEV1OverTimeM, mean)
# tmp2 <- aggregate(SR.r2~bin, meanEV2OverTimeM, mean)
# lines(tmp$bin, tmp$SR.r1)
# lines(tmp$bin, tmp2$SR.r2)
# tmp <- cbind(tmp, tmp2$SR.r2)
# colnames(tmp) <- c('bin', 'r1', 'r2')
# 
# tmp$r1+tmp$r2
# 
# plot(0,0, type='n', xlim=range(risk1OverTimeM$bin)+c(-.5, .5), ylim=c(0., 1), xlab='Bin', ylab='Risk')
# for(ease in unique(risk1OverTimeM$ease)) {
#   toplot <- risk1OverTimeM[risk1OverTimeM$ease==ease,]
#   points(toplot$bin,
#          toplot$risk.r1,
#          col=ifelse(ease=="0.6", 1, ifelse(ease=="0.4", 2, ifelse(ease=="0.3", 4, 5))), pch=20, cex=.25)
#   toplot <- risk2OverTimeM[risk2OverTimeM$ease==ease,]
#   points(toplot$bin, 
#          toplot$risk.r2, 
#          col=ifelse(ease=="0.6", 1, ifelse(ease=="0.4", 2, ifelse(ease=="0.3", 4, 5))), pch=20, cex=.25)
# }
# 
# 
# plot(0,0, type='n', xlim=range(meanV1OverTimeM$bin)+c(-.5, .5), ylim=c(1.5, 2.5), xlab='Bin', ylab='mean v')
# for(ease in unique(meanb1OverTimeM$ease)) {
#   toplot <- meanb1OverTimeM[meanb1OverTimeM$ease==ease,]
#   points(toplot$bin, 
#          toplot$b.r1, 
#          col=ifelse(ease=="0.6", 1, ifelse(ease=="0.4", 2, ifelse(ease=="0.3", 4, 5))), pch=20, cex=.25)
#   toplot <- meanb2OverTimeM[meanb2OverTimeM$ease==ease,]
#   points(toplot$bin, 
#          toplot$b.r2, 
#          col=ifelse(ease=="0.6", 1, ifelse(ease=="0.4", 2, ifelse(ease=="0.3", 4, 5))), pch=20, cex=.25)
# }
# 
# plot(0,0, type='n', xlim=range(meanPEsOverTimeM$bin)+c(-.5, .5), ylim=c(0, .7), xlab='Bin', ylab='b')
# for(ease in unique(meanPEsOverTimeM$ease)) {
#   toplot <- meanPEsOverTimeM[meanPEsOverTimeM$ease==ease,]
#   points(toplot$bin, 
#          toplot$PEs, 
#          col=ifelse(ease=="0.6", 1, ifelse(ease=="0.4", 2, ifelse(ease=="0.3", 4, 5))), pch=20, cex=.25)
# }
# 
# 
# 
# ###
# for(pp in unique(dat$pp)) {
#   for(stimulus_set in unique(dat$stimulus_set)) {
#     idx <- dat$pp == pp & dat$stimulus_set == stimulus_set
#     dat[idx, 'trialNthisStimSet'] <- 1:sum(idx)
#   }
# }
# 
# 
# library(lme4)
# library(lmerTest)
# lm <- lmer(RT~ ease*trialNthisStimSet + ease|pp, dat)
# 
# summary(lm)
# 
# 
# 
# 
# 
# 
# par(mfrow=c(1,1))
# i=1
# plot(0,0, xlim=c(1,10), ylim=c(0.5, 0.85), xlab='bin', ylab='RT', type='n')
# for(ease in unique(meanRTsByEase[[1]]$ease)) {
#   idxD = meanRTsByEase[[1]]$ease == ease
#   toPlot <- meanRTsByEase[[1]][idxD,]
#   
#   offset = 0.8 - toPlot$RT[1]
#   toPlot$RT <- toPlot$RT+offset
#   lines(toPlot$bin, toPlot$RT, col=i)
#   i=i+1
#   #  plotDataPPBins(data=meanRTsByEase[[1]][idxD,], pp=meanRTsByEase[[2]][idxM,], draw.legend=FALSE)
#   #  plotDataPPBins(data=meanAccByEase[[1]][idxD,], pp=meanAccByEase[[2]][idxM,], dep.var='acc', draw.legend=FALSE)
#   #  title(ease)
#   #  plotDataPPBins(data=meanSkewByEase[[1]][idxD,], pp=meanSkewByEase[[2]][idxM,], dep.var='RT', draw.legend=FALSE)
# }
