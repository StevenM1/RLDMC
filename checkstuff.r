rm(list=ls())
source("dmc/dmc.R")
source('utils.R')
source('models.R')
source('PlotUtils.R')

modelName <- "arw-RL-timing_B_Tmod_NoV0"
samples <- loadSamples('samples/model-arw-RL-timing_B_Tmod_NoV0_data-exp3')

model <- setupModel(modelName)
pp = h.post.predict.dmc(samples = samples, adapt=TRUE,save.simulation = TRUE, cores=1)
summary.dmc(samples)
pairs.dmc(samples, thin = 50)

plot.dmc(samples, start = 1500, hyper = T)
plot.dmc(samples, pll.chain = T, start = 1500, hyper = T)
pdf('./figures/pp_PLL_comp_Niek.pdf', width=7, height=7/4*3)
for(i in 1:length(samples)){
  plot.dmc(samples[[i]], start = 0, pll.chain = T)
}
dev.off()

pdf('./figures/pp_theta_comp_Niek.pdf', width=7, height=7/4*3)
for(i in 1:length(samples)){
  plot.dmc(samples[[i]], start = 0)
}
dev.off()
ppNoSim <- h.pp.summary(pp, samples=samples)
tmp <- getDataPpBPIC(model = model, samples = samples, dataName = 'exp1') 
BIC <- getDataPpBPIC(model = model, samples = samples, dataName = 'exp1', BPIConly=TRUE)$BPIC
BPICDDM <- tmp$BPIC
qRTs <- getqRTs(tmp[['data3']], tmp[['pp3']])

allqRTs <- list(qRTs)



pdf('./figures/modelcomparison-exp1_Niek.pdf', width=7, height=7/4*3)
par(oma=c(3,4,1,0), mar=c(0, 0, 1, 0.5) + 0.1, mfcol=c(3,4), mgp=c(2.75,.75,0), las=1, bty='l')
i <- 0
data.cex=1.5
corrRTylim <- errRTylim <- c(.35, 1.1)
for(qRTs in allqRTs) {
  i <- i+1
  
  plotDataPPBins(data=qRTs$meanAccOverTime[[1]], pp=qRTs$meanAccOverTime[[2]],
                 xaxt='n', draw.legend = i==1, data.cex = data.cex,
                 dep.var='acc', ylab='', xlab = '', yaxt='n',
                 legend.pos='topleft', ylim=c(0.5, 0.85), hline.by=0.1)
  axis(1, at=seq(2, 10, 2), labels=rep(NA, 5), lwd=2)
  if(i == 1) {
    mtext('Accuracy', side=2, cex=.66, line=3, las=0, font=1)
    axis(2, at=seq(.5, .9, .1), lwd=1.5)
  } else {
    axis(2, at=seq(.5, .9, .1), labels=rep(NA, 5), lwd=1.5)
  }
  if(i == 1) title('RL-timing-ARW')
  if(i == 2) title('RL-RD')
  if(i == 3) title('RL-ARD')
  if(i == 4) title('RL-fARD')
  
  ##
  plotDataPPBins(data=qRTs$q10RTsOverTime[[1]], pp=qRTs$q10RTsOverTime[[2]], dep.var='RT.10.', 
                 ylim=corrRTylim, xaxt='n', ylab='', yaxt='n', draw.legend = FALSE, data.cex = data.cex)
  plotDataPPBins(data=qRTs$q50RTsOverTime[[1]], pp=qRTs$q50RTsOverTime[[2]], dep.var='RT.50.', plot.new = FALSE, draw.legend=FALSE, data.cex = data.cex)
  plotDataPPBins(data=qRTs$q90RTsOverTime[[1]], pp=qRTs$q90RTsOverTime[[2]], dep.var='RT.90.', plot.new = FALSE, draw.legend=FALSE, data.cex = data.cex)
  axis(1, at=seq(2, 10, 2), labels=NA, lwd=1.5)
  if(i == 1) {
    mtext('Correct RTs (s)', side=2, cex=.66, line=3, las=0, font=1)
    axis(2, seq(.4, 1.2, .2), lwd=1.5)
  } else {
    axis(2, seq(.4, 1.2, .2), labels=NA, lwd=1.5)
  }
  
  ##
  plotDataPPBins(data=qRTs$q10RTsOverTimeE[[1]], pp=qRTs$q10RTsOverTimeE[[2]], dep.var='RT.10.', ylim=errRTylim, xaxt='n', ylab='', yaxt='n', draw.legend = FALSE, data.cex = data.cex)
  plotDataPPBins(data=qRTs$q50RTsOverTimeE[[1]], pp=qRTs$q50RTsOverTimeE[[2]], dep.var='RT.50.', plot.new = FALSE, draw.legend=FALSE, data.cex = data.cex)
  plotDataPPBins(data=qRTs$q90RTsOverTimeE[[1]], pp=qRTs$q90RTsOverTimeE[[2]], dep.var='RT.90.', plot.new = FALSE, draw.legend=FALSE, data.cex = data.cex)
  if(i == 1) {
    mtext('Error RTs (s)', side=2, cex=.66, line=3, las=0, font=1)
    axis(2, seq(.4, 1.2, .2), lwd=1.5)
  } else {
    axis(2, seq(.4, 1.2, .2), labels=NA, lwd=1.5)
  }
  axis(1, at=seq(2, 10, 2), lwd=1.5)
  mtext('Trial bin', side=1, cex=.66, line=2)
}
dev.off()




