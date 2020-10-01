## load dataset
source('dmc/models/RW/arw-RL-mag.R')
source('utils.R')


pdf(file='figures/across-block2.pdf', width=6, height=5)
par(oma=c(3,1,1,0), mar=c(0, 4, 1, 0.5) + 0.1, mfrow=c(2,2), mgp=c(2.75,.75,0), las=1, bty='l')


# Experiment 3 ------------------------------------------------------------
dataName <- 'exp3'
tmp <- loadData(dataName)
data <- tmp[['data']]
dat <- tmp[['dat']]
for(pp in unique(dat$pp)) {
  for(block in unique(dat$block_nr[dat$pp==pp])) {
    idx <- dat$pp==pp&dat$block_nr==block
    dat[idx,]$trialBin <- cut(dat[idx,'trial_nr'], breaks=10, labels = FALSE)
    
  }
}

## RTs
rtByBlock <- aggregate(rt~block_nr*trialBin, aggregate(rt~block_nr*trialBin*pp, dat, mean), mean)
sertByBlock <- aggregate(rt~block_nr*trialBin, aggregate(rt~block_nr*trialBin*pp, dat, mean), function(x) sd(x)/sqrt(length(x)))

nBlocks <- max(dat$block_nr)
plot(0,0,type='n', xlim=c(1, 10), ylim=c(.45, .8), xlab='Trial bin', ylab='Mean RT (s)', xaxt='n')
title('Mean RT (s)', line = .25)
mtext('Experiment 2',
      side=2, cex=0.66*1.2, las=0, line=4, font=2)
axis(side=1, labels=FALSE)
abline(h=seq(0, 1, .05), col='grey', lty=3)
abline(v=seq(0, 10, 1), col='grey', lty=3)
for(block in unique(rtByBlock$block_nr) ) {
  idx = rtByBlock$block_nr == block
  x <- rtByBlock[idx, 'trialBin'] - block/10 + (nBlocks/20)
  lines(x, rtByBlock[idx, 'rt'], col=block, lwd=2)
  
  arrows(x, 
         y0=rtByBlock[idx, 'rt']-sertByBlock[idx, 'rt'],
         y1=rtByBlock[idx, 'rt']+sertByBlock[idx, 'rt'], 
         length = .05, angle = 90, code=3, col=block, lwd=2)
}
legend('topright', paste0('Block ', 1:nBlocks), col=1:3, lty=c(1,1,1), lwd=2, bty='n')


## accuracies
accByBlock <- aggregate(choiceIsHighP~block_nr*trialBin, aggregate(choiceIsHighP~block_nr*trialBin*pp, dat, mean), mean)
seaccByBlock <- aggregate(choiceIsHighP~block_nr*trialBin, aggregate(choiceIsHighP~block_nr*trialBin*pp, dat, mean), function(x) sd(x)/sqrt(length(x)))

nBlocks <- max(dat$block_nr)
plot(0,0,type='n', xlim=c(1, 10), ylim=c(.45, .9), xlab='Trial bin', ylab='Accuracy', xaxt='n')
title('Accuracy', line = .25)
axis(side=1, labels=FALSE)
abline(h=seq(0, 1, .05), col='grey', lty=3)
abline(v=seq(0, 10, 1), col='grey', lty=3)
for(block in unique(rtByBlock$block_nr) ) {
  idx = accByBlock$block_nr == block
  x <- accByBlock[idx, 'trialBin'] - block/10 + (nBlocks/20)
  lines(x, accByBlock[idx, 'choiceIsHighP'], col=block, lwd=2)
  
  arrows(x, 
         y0=accByBlock[idx, 'choiceIsHighP']-seaccByBlock[idx, 'choiceIsHighP'],
         y1=accByBlock[idx, 'choiceIsHighP']+seaccByBlock[idx, 'choiceIsHighP'], 
         length = .05, angle = 90, code=3, col=block, lwd=2)
}
legend('bottomright', paste0('Block ', 1:nBlocks), col=1:3, lty=c(1,1,1), lwd=2, bty='n')



# Experiment 2 ------------------------------------------------------------
dataName <- 'exp2'
tmp <- loadData(dataName)
data <- tmp[['data']]
dat <- tmp[['dat']]
dat$block_nr <- dat$block
dat$choiceIsHighP <- dat$choiceIsHighP_orig

## RTs
rtByBlock <- aggregate(rt~block_nr*trialBin, aggregate(rt~block_nr*trialBin*pp, dat, mean), mean)
sertByBlock <- aggregate(rt~block_nr*trialBin, aggregate(rt~block_nr*trialBin*pp, dat, mean), function(x) sd(x)/sqrt(length(x)))

nBlocks <- max(dat$block_nr)
plot(0,0,type='n', xlim=c(1, 10), ylim=c(.525, .8), xlab='Trial bin', ylab='Mean RT (s)', xaxt='s')
axis(side=1, labels=FALSE)
abline(h=seq(0, 1, .05), col='grey', lty=3)
abline(v=seq(0, 10, 1), col='grey', lty=3)
mtext('Experiment 3',
      side=2, cex=0.66*1.2, las=0, line=4, font=2)
for(block in unique(rtByBlock$block_nr) ) {
  idx = rtByBlock$block_nr == block
  x <- rtByBlock[idx, 'trialBin'] - block/10 + (nBlocks/20)
  lines(x, rtByBlock[idx, 'rt'], col=block, lwd=2)
  
  arrows(x, 
         y0=rtByBlock[idx, 'rt']-sertByBlock[idx, 'rt'],
         y1=rtByBlock[idx, 'rt']+sertByBlock[idx, 'rt'], 
         length = .05, angle = 90, code=3, col=block, lwd=2)
}
legend('topright', paste0('Block ', 1:nBlocks), col=1:nBlocks, lty=rep(1, nBlocks), lwd=2, bty='n')


## accuracies
accByBlock <- aggregate(choiceIsHighP~block_nr*trialBin, aggregate(choiceIsHighP~block_nr*trialBin*pp, dat, mean), mean)
seaccByBlock <- aggregate(choiceIsHighP~block_nr*trialBin, aggregate(choiceIsHighP~block_nr*trialBin*pp, dat, mean), function(x) sd(x)/sqrt(length(x)))

nBlocks <- max(dat$block_nr)
plot(0,0,type='n', xlim=c(1, 10), ylim=c(.25, .9), xlab='Trial bin', ylab='Accuracy', xaxt='s')
axis(side=1, labels=FALSE)
abline(h=seq(0, 1, .05), col='grey', lty=3)
abline(v=seq(0, 10, 1), col='grey', lty=3)
for(block in unique(rtByBlock$block_nr) ) {
  idx = accByBlock$block_nr == block
  x <- accByBlock[idx, 'trialBin'] - block/10 + (nBlocks/20)
  lines(x, accByBlock[idx, 'choiceIsHighP'], col=block, lwd=2)
  
  arrows(x, 
         y0=accByBlock[idx, 'choiceIsHighP']-seaccByBlock[idx, 'choiceIsHighP'],
         y1=accByBlock[idx, 'choiceIsHighP']+seaccByBlock[idx, 'choiceIsHighP'], 
         length = .05, angle = 90, code=3, col=block, lwd=2)
}
legend('bottomleft', paste0('Block ', 1:nBlocks), col=1:nBlocks, lty=rep(1, nBlocks), lwd=2, bty='n')


dev.off()




# 
# # Experiment 4 ------------------------------------------------------------
# dataName <- 'exp4'
# tmp <- loadData(dataName)
# data <- tmp[['data']]
# dat <- tmp[['dat']]
# dat$block_nr <- dat$block
# 
# ## RTs
# rtByBlock <- aggregate(rt~block_nr*trialBin, aggregate(rt~block_nr*trialBin*pp, dat, mean), mean)
# sertByBlock <- aggregate(rt~block_nr*trialBin, aggregate(rt~block_nr*trialBin*pp, dat, mean), function(x) sd(x)/sqrt(length(x)))
# 
# nBlocks <- max(dat$block_nr)
# 
# plot(0,0,type='n', xlim=c(1, max(dat$trialBin)), ylim=c(.5, .85), xlab='Trial bin', ylab='Mean RT (s)')
# mtext('Experiment 4',
#       side=2, cex=0.66*1.2, las=0, line=4, font=2)
# abline(h=seq(0, 1, .05), col='grey', lty=3)
# abline(v=seq(0, 10, 1), col='grey', lty=3)
# mtext('Trial bin', side=1, cex=.66, line=2)
# for(block in unique(rtByBlock$block_nr) ) {
#   idx = rtByBlock$block_nr == block
#   x <- rtByBlock[idx, 'trialBin'] - block/10 + (nBlocks/20)
#   lines(x, rtByBlock[idx, 'rt'], col=block, lwd=2)
#   
#   arrows(x, 
#          y0=rtByBlock[idx, 'rt']-sertByBlock[idx, 'rt'],
#          y1=rtByBlock[idx, 'rt']+sertByBlock[idx, 'rt'], 
#          length = .05, angle = 90, code=3, col=block, lwd=2)
# }
# legend('topright', paste0('Block ', 1:3), col=1:3, lty=c(1,1,1), lwd=2, bty='n')
# legend('bottomleft', paste0('Block ', 4:nBlocks), col=4:nBlocks, lty=c(1,1), lwd=2, bty='n')
# 
# 
# ## accuracies
# accByBlock <- aggregate(choiceIsHighP~block_nr*trialBin, aggregate(choiceIsHighP~block_nr*trialBin*pp, dat, mean), mean)
# seaccByBlock <- aggregate(choiceIsHighP~block_nr*trialBin, aggregate(choiceIsHighP~block_nr*trialBin*pp, dat, mean), function(x) sd(x)/sqrt(length(x)))
# 
# nBlocks <- max(dat$block_nr)
# plot(0,0,type='n', xlim=c(1, max(dat$trialBin)), ylim=c(.45, .8), xlab='Trial bin', ylab='Accuracy')
# abline(h=seq(0, 1, .05), col='grey', lty=3)
# abline(v=seq(0, 10, 1), col='grey', lty=3)
# mtext('Trial bin', side=1, cex=.66, line=2)
# for(block in unique(rtByBlock$block_nr) ) {
#   idx = accByBlock$block_nr == block
#   x <- rtByBlock[idx, 'trialBin'] - block/20 + (nBlocks/40)
#   lines(x, accByBlock[idx, 'choiceIsHighP'], col=block, lwd=2)
#   
#   arrows(x, 
#          y0=accByBlock[idx, 'choiceIsHighP']-seaccByBlock[idx, 'choiceIsHighP'],
#          y1=accByBlock[idx, 'choiceIsHighP']+seaccByBlock[idx, 'choiceIsHighP'], 
#          length = .05, angle = 90, code=3, col=block, lwd=2)
# }
# legend('bottomright', paste0('Block ', 1:nBlocks), col=1:nBlocks, lty=c(1,1,1), lwd=2, bty='n')
# dev.off()
# 
# 
# 




# Statistics? -------------------------------------------------------------
library(lme4)
library(lmerTest)
library(report)

# Exp 2
dataName <- 'exp2'
tmp <- loadData(dataName)
data <- tmp[['data']]
dat <- tmp[['dat']]
dat$block_nr <- dat$block
dat$lTrialBin <- log(dat$trialBin)

mod2rt = lmer(rt~block_nr*trialBin + (1|pp), data=dat)
summary(mod2rt)
report(mod2rt)

mod2acc = glmer(choiceIsHighP_orig~block_nr*lTrialBin + (1|pp), data=dat, family='binomial')
summary(mod2acc)
report(mod2acc)


# Exp 3
dataName <- 'exp3'
tmp <- loadData(dataName)
data <- tmp[['data']]
dat <- tmp[['dat']]
for(pp in unique(dat$pp)) {
  for(block in unique(dat$block_nr[dat$pp==pp])) {
    idx <- dat$pp==pp&dat$block_nr==block
    dat[idx,]$trialBin <- cut(dat[idx,'trial_nr'], breaks=10, labels = FALSE)
    
  }
}
dat$lTrialBin <- log(dat$trialBin)

mod3rt = lmer(rt~block_nr*trialBin + (1|pp), data=dat)
summary(mod3rt)
report(mod3rt)

mod3acc = glmer(choiceIsHighP~block_nr*lTrialBin + (1|pp), data=dat, family='binomial')
summary(mod3acc)
report(mod3acc)

# Exp 4
dataName <- 'exp4'
tmp <- loadData(dataName)
data <- tmp[['data']]
dat <- tmp[['dat']]
dat$block_nr <- dat$block
dat$lTrialBin <- dat$trialBin

mod4rt = lmer(rt~block_nr*trialBin + (1|pp), data=dat)
summary(mod4rt)
report(mod4rt)

mod4acc = glmer(choiceIsHighP~block_nr*lTrialBin + (1|pp), data=dat, family='binomial')
summary(mod4acc)







summary(mod2rt)

