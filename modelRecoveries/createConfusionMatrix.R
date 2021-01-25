rm(list=ls())
source ("dmc/dmc.R")
source('utils.R')

allSamples <- list.files('./modelRecoveries/samples/')
allSamples

getBPICs <- function(generatingModel, ds) {
  allSamples <- list.files('./modelRecoveries/samples/', pattern=paste0('gen-', generatingModel, '_ds-', ds, '.RData'))
  
  resultsDf <- data.frame(generatingModel=numeric(0), dataset_n=numeric(0), fitModel=numeric(0), subjectN=numeric(0), BPIC=numeric(0), mindev=numeric(0), converged=numeric(0))
  #resultsDf <- data.frame(generatingModel=numeric(0), dataset_n=numeric(0)) #, BPIC=numeric(0), mindev=numeric(0), converged=numeric(0))
  
  for(fn in allSamples) {
    fitModel = gsub('model-', '', strsplit(fn, '_')[[1]][1])
    generatingModel= gsub('gen-', '', strsplit(fn, '_')[[1]][3])
    dataset_n = gsub('.RData', '', gsub('ds-', '', strsplit(fn, '_')[[1]][4]))
    
    load(file.path('./modelRecoveries/samples', fn))
    if(fitModel == 'arw-RL-mag') {
      load_model("RW", "arw-RL-mag.R")
    } else {
      load_model("ddm", "ddm-RL.R")
    }
  #  gm <- h.gelman.diag.dmc(hsamples)
  #  converged = all(gm < 1.05)
    converged=NA
    ic <- h.IC.dmc(hsamples)
    bpic <- apply(ic, 2, sum)[2]
    mindev <- apply(ic, 2, sum)[1]
    
    resultsDf <- rbind(resultsDf, data.frame(generatingModel=generatingModel, dataset_n=ds, fitModel=fitModel, 
                                             subjectN=c(1:nrow(ic), 'dataset'), BPIC=c(ic[,2], bpic), mindev=c(ic[,1],mindev), converged=converged))
  }
  return(resultsDf)
}

allDatasets <- expand.grid(ds=1:50, generatingModel=c('RL-ARD', 'RL-DDM'))
allDatasets
#debugonce(getBPICs)
#lapply(1:nrow(allDatasets), function(x) getBPICs(allDatasets[x,'generatingModel'], allDatasets[x,'ds']))

library(snowfall)
sfInit(parallel = TRUE, cpus=30)
sfExportAll()
x <- sfLapply(1:nrow(allDatasets), function(x) getBPICs(allDatasets[x,'generatingModel'], allDatasets[x,'ds']))
x2 <- do.call(rbind, x)
library(reshape2)
x2wide <- reshape(x2, direction='wide', 
                  v.names=c('BPIC', 'mindev'), timevar=c('fitModel'), 
                  idvar=c('generatingModel', 'dataset_n', 'subjectN'))

x2wide$winner <- NA
x2wide$arwRLmagWinsBPIC <- x2wide$`BPIC.arw-RL-mag` < x2wide$`BPIC.ddm-RL`
x2wide$arwRLmagWinsMinDev <- x2wide$`mindev.arw-RL-mag` < x2wide$`mindev.ddm-RL`
x2wide$BPICwinner <- ifelse(x2wide$`BPIC.arw-RL-mag` < x2wide$`BPIC.ddm-RL`, 'RL-ARD', 'RL-DDM')
x2wide$mindevwinner <- ifelse(x2wide$`mindev.arw-RL-mag` < x2wide$`mindev.ddm-RL`, 'RL-ARD', 'RL-DDM')

x2wide$BPICwinner <- factor(x2wide$BPICwinner, levels=c('RL-ARD', 'RL-DDM'))
x2wide$mindevwinner <- factor(x2wide$mindevwinner, levels=c('RL-ARD', 'RL-DDM'))

# Subjects only
isDataset <- x2wide$subjectN == 'dataset'
# By BPIC
table(x2wide$BPICwinner[!isDataset], x2wide$generatingModel[!isDataset])
# By  min deviance
table(x2wide$mindevwinner[!isDataset], x2wide$generatingModel[!isDataset])

# Datasets
table(x2wide$BPICwinner[isDataset], x2wide$generatingModel[isDataset])
# By  min deviance
table(x2wide$mindevwinner[isDataset], x2wide$generatingModel[isDataset])


# Matrices
#TClass <- factor(c(0, 0, 1, 1))
#PClass <- factor(c(0, 1, 0, 1))
#Y      <- c(2816, 248, 34, 235)
#df <- data.frame(TClass, PClass, Y)
library(caret)
fourfoldplot(cm$table)


library(gridExtra)
library(ggplot2)

#cm <- confusionMatrix(x2wide$mindevwinner[isDataset], x2wide$generatingModel[isDataset])
plots <- list()
for(i in 1:4) {
  if(i == 2) {
    cm <- confusionMatrix(x2wide$mindevwinner[isDataset], x2wide$generatingModel[isDataset])
    main = 'Minimum deviance'
    ylabel = '' #Entire dataset'
  } else if(i == 1) {
    cm <- confusionMatrix(x2wide$BPICwinner[isDataset], x2wide$generatingModel[isDataset])
    main = 'BPIC'
    ylabel = 'Datasets\nWinning'
  } else if(i == 4) {
    cm <- confusionMatrix(x2wide$mindevwinner[!isDataset], x2wide$generatingModel[!isDataset])
    main = ''
    ylabel =''
  } else if(i == 3) {
    cm <- confusionMatrix(x2wide$BPICwinner[!isDataset], x2wide$generatingModel[!isDataset])
    main =''
    ylabel = 'Subjects\nWinning' #Entire dataset'
  }
  df <- data.frame(cm$table)
  df$Winning <- df$Prediction
  df$Generating <- df$Reference
  p1 <- ggplot(data =  df, mapping = aes(x = Generating, y = ordered(Winning, levels=rev(levels(df$Winning))))) +
    geom_tile(aes(fill = Freq), colour = "white") +
    geom_text(aes(label = sprintf("%1.0f", Freq)), vjust = 1) +
    scale_fill_gradient(low = "white", high = "cornflowerblue") +
    theme_bw() + theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) + labs(y = ylabel, title = main)
  plots[[i]] <- p1
}
pdf(file='./figures/exp1-confusionmatrices.pdf', width=6, height=5)
grid.arrange(plots[[1]],plots[[2]],plots[[3]],plots[[4]], nrow=2)
dev.off()


df <- data.frame(cm$table)
df$Winning <- df$Prediction
df$Generating <- df$Reference
p1 <- ggplot(data =  df, mapping = aes(x = Generating, y = ordered(Winning, levels=rev(levels(df$Winning))))) +
      geom_tile(aes(fill = Freq), colour = "white") +
      geom_text(aes(label = sprintf("%1.0f", Freq)), vjust = 1) +
      scale_fill_gradient(low = "white", high = "blue") +
      theme_bw() + theme(legend.position = "none") + labs(y = 'Winning')

#
cm <- confusionMatrix(x2wide$BPICwinner[isDataset], x2wide$BPICwinner[isDataset])
df <- data.frame(cm$table)
df$Winning <- df$Prediction
df$Generating <- df$Reference
p2 <- ggplot(data =  df, mapping = aes(x = Generating, y = ordered(Winning, levels=rev(levels(df$Winning))))) +
  geom_tile(aes(fill = Freq), colour = "white") +
  geom_text(aes(label = sprintf("%1.0f", Freq)), vjust = 1) +
  scale_fill_gradient(low = "black", high = "white") +
  theme_bw() + theme(legend.position = "none") + labs(y = 'Winning')


#
cm <- confusionMatrix(x2wide$mindevwinner[!isDataset], x2wide$generatingModel[!isDataset])
df <- data.frame(cm$table)
df$Winning <- df$Prediction
df$Generating <- df$Reference
p3 <- ggplot(data =  df, mapping = aes(x = Generating, y = ordered(Winning, levels=rev(levels(df$Winning))))) +
  geom_tile(aes(fill = Freq), colour = "white") +
  geom_text(aes(label = sprintf("%1.0f", Freq)), vjust = 1) +
  scale_fill_gradient(low = "black", high = "white") +
  theme_bw() + theme(legend.position = "none") + labs(y = 'Winning')

#
cm <- confusionMatrix(x2wide$BPICwinner[!isDataset], x2wide$BPICwinner[!isDataset])
df <- data.frame(cm$table)
df$Winning <- df$Prediction
df$Generating <- df$Reference
p4 <- ggplot(data =  df, mapping = aes(x = Generating, y = ordered(Winning, levels=rev(levels(df$Winning))))) +
  geom_tile(aes(fill = Freq), colour = "white") +
  geom_text(aes(label = sprintf("%1.0f", Freq)), vjust = 1) +
  scale_fill_gradient(low = "black", high = "white") +
  theme_bw() + theme(legend.position = "none") + labs(y = 'Winning')

grid.arrange(p1, p2, p3, p4, nrow=2)









resultsDf$dataset_n <- as.numeric(as.character(resultsDf$dataset_n))
winners <- aggregate(BPIC~generatingModel*dataset_n, resultsDf, 
                     function(x) ifelse(length(x) == 1, 'None yet',
                                        ifelse(which.min(x)==1, 'RL-ARD', 'RL-DDM')))
#winners$dataset_n <- as.numeric(as.character(winners$dataset_n))
winners50 <- winners[winners$dataset_n < 51,]
table(winners50$BPIC, winners50$generatingModel)
# ToDo: check parameter recoveries for all instances where generatingModel == fitModel (to ensure sampling went correctly)

## double-check

test <- data.frame(pred=c(runif(50,0,75),runif(50,25,100)), group=c(rep("A",50), rep("B",50)) )

# Data-generting model = DDM, which model fits better?
# generatingModel = 'RL-ARD'
# dataset_n = 1
# 
# load_model("RW", "arw-RL-mag.R")
# load(paste0('./modelRecoveries/samples/model-arw-RL-mag_data-modelRecovery_gen-', generatingModel, '_ds-', dataset_n, '.RData'))
# h.IC.dmc(hsamples)
# 
# load_model("ddm", "ddm-RL.R")
# load(paste0('./modelRecoveries/samples/model-ddm-RL_data-modelRecovery_gen-', generatingModel, '_ds-', dataset_n, '.RData'))
# h.IC.dmc(hsamples)







# fit? --------------------------------------------------------------------

# Fit?
modelName <- 'ddm-RL' #'arw-RL-mag'
generatingModel = 'RL-DDM'
dataset_n <- 1
dataName <- paste0('modelRecovery_gen-', generatingModel, '_ds-', dataset_n)
fn = paste0('model-', modelName, '_data-', dataName)
load(file.path('./modelRecoveries/data', paste0(generatingModel, '_simulated'), paste0(dataName, '.RData')))
load(paste0('./modelRecoveries/samples/', fn, '.RData'))

getDataPpBPIC <- function(samples, modelName, do.plot=FALSE, BPIConly=FALSE) {
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
#  model <- setupModel(modelName)  # calls load_model(), which loads transform.dmc() and transform2.dmc()
#  dat <- loadData(dataName, removeBlock = NULL)[['dat']]
  fn <- paste0('model-', modelName, '_data-', dataName)

  # Load, generate posterior preds -------------------------------------------
#  samples <- loadSamples(fn, samplesDir)
  data <- lapply(samples, function(x) x$data)

  # imidate 'original data'
  dat <- do.call(rbind, lapply(1:length(data),
                                function(x) {data[[x]]$sub <- x;
                                             data[[x]]$ease <- data[[x]]$stimulus_set; data[[x]]}))
  if(do.plot) plot.dmc(samples, hyper=TRUE, density=TRUE, layout=c(4,4))
  if(!BPIConly) {
    pp = h.post.predict.dmc(samples = samples, adapt=TRUE,save.simulation = TRUE, cores=30)
    ppNoSim <- h.pp.summary(pp, samples=samples)

    #### Append stimulus set info to data & model --------
    nBins <- 10
    pp2 <- lapply(1:length(pp), addStimSetInfo, input=pp, orig_dat=dat)
    data2 <- lapply(1:length(data), addStimSetInfo, input=data, orig_dat=dat)
    if(!sfIsRunning()) sfInit(parallel=TRUE, cpus =30); sfLibrary(moments)
    pp3 <- sfLapply(pp2, calculateByBin)
    data3 <- lapply(data2, calculateByBin)

    bpics <- h.IC.dmc(samples)
    return(list('pp3'=pp3, 'data3'=data3, 'BPIC'=bpics))
  } else {
    return(list('BPIC'=h.IC.dmc(samples)))
  }
}

# Load BPICs & quantiles per bin ---------------------------------------------------------------
tmp <- getDataPpBPIC(hsamples, modelName)
BPIC <- tmp$BPIC
data3 <- tmp[['data3']]
pp3 <- tmp[['pp3']]

q10RTsByEase <- list(getDescriptives(data3, dep.var='RT.10.', attr.name='qRTsCorrectByEase', id.var1='~bin*ease', id.var2=NULL),
                     getDescriptives(pp3, dep.var='RT.10.', attr.name='qRTsCorrectByEase', id.var1='~reps*bin*ease', id.var2=NULL))
q50RTsByEase <- list(getDescriptives(data3, dep.var='RT.50.', attr.name='qRTsCorrectByEase', id.var1='~bin*ease', id.var2=NULL),
                     getDescriptives(pp3, dep.var='RT.50.', attr.name='qRTsCorrectByEase', id.var1='~reps*bin*ease', id.var2=NULL))
q90RTsByEase <- list(getDescriptives(data3, dep.var='RT.90.', attr.name='qRTsCorrectByEase', id.var1='~bin*ease', id.var2=NULL),
                     getDescriptives(pp3, dep.var='RT.90.', attr.name='qRTsCorrectByEase', id.var1='~reps*bin*ease', id.var2=NULL))
q10RTsByEaseE <- list(getDescriptives(data3, dep.var='RT.10.', attr.name='qRTsErrorByEase', id.var1='~bin*ease', id.var2=NULL),
                      getDescriptives(pp3, dep.var='RT.10.', attr.name='qRTsErrorByEase', id.var1='~reps*bin*ease', id.var2=NULL))
q50RTsByEaseE <- list(getDescriptives(data3, dep.var='RT.50.', attr.name='qRTsErrorByEase', id.var1='~bin*ease', id.var2=NULL),
                      getDescriptives(pp3, dep.var='RT.50.', attr.name='qRTsErrorByEase', id.var1='~reps*bin*ease', id.var2=NULL))
q90RTsByEaseE <- list(getDescriptives(data3, dep.var='RT.90.', attr.name='qRTsErrorByEase', id.var1='~bin*ease', id.var2=NULL),
                      getDescriptives(pp3, dep.var='RT.90.', attr.name='qRTsErrorByEase', id.var1='~reps*bin*ease', id.var2=NULL))
meanAccByEase <- list(getDescriptives(data3, dep.var='acc', attr.name='AccByEase', id.var1='~bin*ease', id.var2=NULL),
                      getDescriptives(pp3, dep.var='acc', attr.name='AccByEase', id.var1='~reps*bin*ease', id.var2=NULL))

par(oma=c(3,3,1,0), mar=c(0, 1, 1, 0) + 0.1, mfcol=c(3,4), mgp=c(2.75,.75,0), las=1, bty='l')
i <- 0
corrRTylim <- c(0.45, 1.1)
errRTylim <- c(0.45, 1.1)
data.cex=1.5
for(ease in unique(meanAccByEase[[1]]$ease)) {
  i <- i+1
  idxD = meanAccByEase[[1]]$ease == ease
  idxM = meanAccByEase[[2]]$ease == ease

  plotDataPPBins(data=meanAccByEase[[1]][idxD,], pp=meanAccByEase[[2]][idxM,],
                 xaxt='n', draw.legend = i==1, data.cex = data.cex,
                 dep.var='acc', ylab='', xlab = '', yaxt='n',
                 legend.pos='topleft', ylim=c(0.4, 0.95), hline.by=0.1)
  axis(1, at=seq(2, 10, 2), labels=rep(NA, 5), lwd=2)
  if(i == 1) {
    mtext('Accuracy', side=2, cex=.66, line=3, las=0, font=1)
    axis(2, at=seq(.5, .9, .1), lwd=1.5)
  } else {
    axis(2, at=seq(.5, .9, .1), labels=rep(NA, 5), lwd=1.5)
  }
  if(i == 1) title('0.6/0.4 (Hardest)')
  if(i == 2) title('0.65/0.35')
  if(i == 3) title('0.7/0.3')
  if(i == 4) title('0.8/0.2 (Easiest)')

  ##
  plotDataPPBins(data=q10RTsByEase[[1]][idxD,], pp=q10RTsByEase[[2]][idxM,], dep.var='RT.10.',
                 ylim=corrRTylim, xaxt='n', ylab='', yaxt='n', draw.legend = FALSE, data.cex = data.cex)
  plotDataPPBins(data=q50RTsByEase[[1]][idxD,], pp=q50RTsByEase[[2]][idxM,], dep.var='RT.50.', plot.new = FALSE, draw.legend=FALSE, data.cex = data.cex)
  plotDataPPBins(data=q90RTsByEase[[1]][idxD,], pp=q90RTsByEase[[2]][idxM,], dep.var='RT.90.', plot.new = FALSE, draw.legend=FALSE, data.cex = data.cex)
  axis(1, at=seq(2, 10, 2), labels=NA, lwd=1.5)
  if(i == 1) {
    mtext('Correct RTs (s)', side=2, cex=.66, line=3, las=0, font=1)
    axis(2, seq(.4, 1.2, .2), lwd=1.5)
  } else {
    axis(2, seq(.4, 1.2, .2), labels=NA, lwd=1.5)
  }

  ##
  plotDataPPBins(data=q10RTsByEaseE[[1]][idxD,], pp=q10RTsByEaseE[[2]][idxM,], dep.var='RT.10.', ylim=errRTylim, xaxt='n', ylab='', yaxt='n', draw.legend = FALSE, data.cex = data.cex)
  plotDataPPBins(data=q50RTsByEaseE[[1]][idxD,], pp=q50RTsByEaseE[[2]][idxM,], dep.var='RT.50.', plot.new = FALSE, draw.legend=FALSE, data.cex = data.cex)
  plotDataPPBins(data=q90RTsByEaseE[[1]][idxD,], pp=q90RTsByEaseE[[2]][idxM,], dep.var='RT.90.', plot.new = FALSE, draw.legend=FALSE, data.cex = data.cex)
  if(i == 1) {
    mtext('Error RTs (s)', side=2, cex=.66, line=3, las=0, font=1)
    axis(2, seq(.4, 1.2, .2), lwd=1.5)
  } else {
    axis(2, seq(.4, 1.2, .2), labels=NA, lwd=1.5)
  }
  axis(1, at=seq(2, 10, 2), lwd=1.5)
  mtext('Trial bin', side=1, cex=.66, line=2)
}

