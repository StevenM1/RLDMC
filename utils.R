doSample <- function(data, p.prior, pp.prior, fileName, restart=FALSE,
                     nCores=31, nmcBurn=250, nmc=1000, cut.converge=1.03, thorough=TRUE) {
  if(restart==FALSE & file.exists(paste0(fileName, '.RData'))) {
    samples <- loadSamples(fileName)
  } else {
    ## run burn-in first - this should wash out most of the prior already, reduce burn-time of h.RUN.dmc
    samples <- h.samples.dmc(nmc=nmcBurn, p.prior=p.prior, data=data, pp.prior=pp.prior)
    samples <- h.run.unstuck.dmc(samples=samples, nmc=nmcBurn, cores=nCores, max.try=100, p.migrate=0.05, h.p.migrate=0.05, end.no.migrate=TRUE)
    save(samples, file=paste0(fileName, '.RData'))
    ## create new samples object, don't add (ie remove burn-in), which will be passed to h.RUN.dmc
    samples <- h.samples.dmc(samples=samples, nmc=nmc, add=FALSE)
  }
  if(samples[[1]]$nmc < nmc) {
    nAdd <- nmc - samples[[1]]$nmc
    print(paste0('Adding ', nAdd, ' samples...'))
    samples <- h.samples.dmc(samples=samples, nmc=nAdd, add=TRUE)
  }
  samples <- h.RUN.dmc(hsamples = samples, cores=nCores, cut.converge=cut.converge, thorough = thorough,
                       verbose=TRUE, saveFn=fileName)
}

dtnorm <- function (x, mean = 0, sd = 1, lower = -Inf, upper = Inf, log = FALSE) 
{
  ret <- numeric(length(x))
  ret[x < lower | x > upper] <- if (log) 
    -Inf
  else 0
  ret[upper < lower] <- NaN
  ind <- x >= lower & x <= upper
  if (any(ind)) {
    denom <- pnorm(upper, mean, sd) - pnorm(lower, mean, 
                                            sd)
    xtmp <- dnorm(x, mean, sd, log)
    if (log) 
      xtmp <- xtmp - log(denom)
    else xtmp <- xtmp/denom
    ret[x >= lower & x <= upper] <- xtmp[ind]
  }
  ret
}

resetPar <- function() {
  dev.new()
  op <- par(no.readonly = TRUE)
  dev.off()
  op
}

loadData <- function(name, exclude=TRUE, addBins=FALSE, removeBlock=NULL, subset=FALSE, dataRoot='./data') {
  otherCols <- c()
  if(name == 'exp1' | name == 'exp2' | name == 'exp3' | name == 'exp3b' | name == 'exp4') {
    load(file.path(dataRoot, paste0('data_', name, '.RData')))
    if(exclude) {
      # remove excluded subjects
      dat <- droplevels(dat[!dat$excl,])
    }
    if(subset) {
      dat <- dat[!dat$subset,]
    }
    dat$sub <- as.factor(as.integer(dat$pp))
  } 
  if(name == 'exp2') {
    # use accuracy coding, reversing accuracies after reversals
    dat$choiceIsHighP_orig <- dat$choiceIsHighP
    dat$choiceIsHighP <-  dat$choiceIsHighPpreRev
  }
  if(name == 'exp3' | name == 'exp3b') otherCols <- c(otherCols, 'cue')
  
  # remove RTs < .15 and null responses
  dat <- dat[dat$rt>.15 & !is.na(dat$rt),]
  if(!is.null(removeBlock)) dat <- dat[dat$block != removeBlock,]
  rownames(dat) <- NULL
  
  # Prepare for fitting
  cvs <- list()
  choiceIdx <- list()
  for(sub in unique(dat$sub)) {
    d <- prepareForFitting(dat[dat$sub==sub,])
    cvs[[sub]] <- d$outcomes
    choiceIdx[[sub]] <- d$VVchoiceIdx
    if(addBins) {
      for(stimSet in unique(dat$stimulus_set)) {
        idx = dat$sub==sub&dat$stimulus_set==stimSet
        dat[idx, 'bin'] <- paste0('b', as.numeric(cut(1:sum(idx), 6)))
      }
    }
  }
  
  # make DMC-style
  data <- dat
  data$s <- data$sub
  data$R <- factor(ifelse(data$choiceIsHighP==0, 'r1', 'r2'))  # accuracy-coded: r1 = error, r2 = correct
  data$RT <- dat$rt
  data$S <- factor('s1')
  
  if(!addBins) {
    data <- data[,c('s', 'S', 'R', 'RT', otherCols)]
  } else {
    data$bin <- factor(data$bin, levels=paste0('b', 1:6))
    data <- data[,c('s', 'S', 'R', 'RT', 'bin', otherCols)]
  }
  
  # Add covariates and other useful things as attributes
  attr(data, 'startingValues') <- d$values
  attr(data, 'trialsToIgnore') <- d$trialsToIgnore
  attr(data, 'VVchoiceIdx') <- choiceIdx
  attr(data, 'cvs') <- cvs
  
  return(list(data=data, dat=dat))
}


getPriors <- function(p.vector) {
  dists <- rep('tnorm', length(p.vector))
  p1 <- upper <- rep(NA, length(p.vector))
  p2 <- rep(5, length(p.vector))
  lower <- rep(0, length(p.vector))
  names(p1) <- names(p2) <- names(upper) <- names(lower) <- names(p.vector)
  
  #s
  p1[grepl(pattern='s', names(p1))] <- 1
  p2[grepl(pattern='s', names(p1))] <- 3 
  
  #s_T
  p1[grepl(pattern='s_T', names(p1))] <- 1
  p2[grepl(pattern='s_T', names(p1))] <- 3 
  
  
  #s_V
  p1[grepl(pattern='v_T', names(p1))] <- 1
  p2[grepl(pattern='v_T', names(p1))] <- 3 
  
  # t0
  p1[grepl(pattern='t0', names(p1))] <- 0.3
  p2[grepl(pattern='t0', names(p1))] <- 0.5
  upper[grepl(pattern='t0', names(p1))] <- 1
  lower[grepl(pattern='t0', names(p1))] <- 0.025
  
  #t0t
  p1[grepl(pattern='t0T', names(p1))] <- 0.3
  p2[grepl(pattern='t0T', names(p1))] <- 0.5
  upper[grepl(pattern='t0T', names(p1))] <- 1
  lower[grepl(pattern='t0T', names(p1))] <- 0.025
  
  # stimulus or risk representations SR/RR
  p1[names(p1)=='SR' | names(p1) == 'RR'] = 0
  p2[names(p1)=='SR' | names(p1) == 'RR'] = 1
  lower[names(p1)=='SR' | names(p1) == 'RR'] = NA
  upper[names(p1)=='SR' | names(p1) == 'RR'] = NA
  
  # # learning rate aV or aR
  # p1[names(p1)=='aV' | names(p1)=='aR'] <- -1.6  # corresponds to 0.05
  # lower[names(p1)=='aV' | names(p1)=='aR'] <- NA  # corresponds to 0
  p1[grepl(pattern='aV', names(p1)) | grepl(pattern='aR', names(p1))] <- -1.6  # corresponds to 0.05
  lower[grepl(pattern='aV', names(p1)) | grepl(pattern='aR', names(p1))] <- NA  # corresponds to 0
  
  # urgency v0 - may vary across conditions, so use grepl
  p1[grepl(pattern='V0', names(p1))] <- 2
  lower[grepl(pattern='V0', names(p1))] <- NA
  
  # weight on v
  p1[grepl(pattern='wV', names(p1))] = 9
  
  # weight on B
  p1[grepl(pattern='wB', names(p1))] = 0
  p2[grepl(pattern='wB', names(p2))] = 0.1
  lower[grepl(pattern='wB', names(p1))] = NA
  
  # threshold - this may vary also across conditions, so use grepl
  p1[grepl(pattern='B0', names(p1))] = 3
  p1[grepl(pattern='B_T', names(p1))] = 3
  
  # weight of magnitude effect / sum of EVs
  p1[grepl(pattern='wS', names(p1))] = 0
  p2[grepl(pattern='wS', names(p2))] = 3
  lower[grepl(pattern='wS', names(lower))] = NA
  
  # Potential SAT effects: Mod(ulators) for V0, B0, t0, aV
  p1[grepl(pattern='Mod.', names(p1))] = 0
  p2[grepl(pattern='Mod.', names(p2))] = 1
  lower[grepl(pattern='Mod.', names(lower))] = -1  # values < -1 would indicate make the absolute parameters negative, and they can't be negative
  
  # Potential SAT effects: Mod(ulators) for V0, B0, t0, aV
  p1[grepl(pattern='mod.', names(p1))] = 0
  p2[grepl(pattern='mod.', names(p2))] = 1
  lower[grepl(pattern='mod.', names(lower))] = -1  # values < -1 would indicate make the absolute parameters negative, and they can't be negative
  
  
  # threshold
  p1[names(p1) == 'a' | names(p1) == 'a.SPD' | names(p1) == 'a.ACC'] = 2
  p2[names(p1) == 'a' | names(p1) == 'a.SPD' | names(p1) == 'a.ACC'] = 0.5
  
  # linear effect m
  p1[names(p1) == 'm' | names(p1) == 'm.SPD' | names(p1) == 'm.ACC'] = 2
  p2[names(p1) == 'm' | names(p1) == 'm.SPD' | names(p1) == 'm.ACC'] = 5
  
  # non-linear effect vmax
  p1[names(p1) == 'vmax'] = 2
  p2[names(p2) == 'vmax'] = 5
  
  # sv
  p1[names(p1) == 'sv'] = 0.1
  p2[names(p2) == 'sv'] = 0.1
  
  # sz
  p1[names(p1) == 'sz'] = 0.1
  p2[names(p2) == 'sz'] = 0.1
  
  # st0
  p1[names(p1) == 'st0'] = 0.1
  p2[names(p1) == 'st0'] = 0.1
#  lower[names(lower) == 'st0'] = 0
  
  # softmax beta
  p1[names(p1) == 'beta' | names(p1) == 'beta.SPD' | names(p1) == 'beta.ACC'] = 1
  p2[names(p1) == 'beta' | names(p1) == 'beta.SPD' | names(p1) == 'beta.ACC'] = 5
  p1[names(p1) == 'Beta' | names(p1) == 'Beta.SPD' | names(p1) == 'Beta.ACC'] = 1
  p2[names(p1) == 'Beta' | names(p1) == 'Beta.SPD' | names(p1) == 'Beta.ACC'] = 5
  
  # lba sd_v
  p1[grepl(pattern='sd_v', x=names(p1))] = 1
  p2[grepl(pattern='sd_v', x=names(p1))] = 1
  lower[grepl(pattern='sd_v', x=names(p1))] = 0

  # lba A
  p1[grepl(pattern='A\\.', x=names(p1)) | names(p1) == 'A'] = 1
  p2[grepl(pattern='A\\.', x=names(p1)) | names(p1) == 'A'] = 1
  lower[grepl(pattern='A\\.', x=names(p1)) | names(p1) == 'A'] = 0
  
  # combine everything
  p.prior <- prior.p.dmc(
    dists = dists,
    p1   = p1,
    p2   = p2,
    lower= lower,
    upper= upper
  )
  
  # par(mfcol=c(2,3)); for (i in names(p.prior)) plot.prior(i,p.prior)
  sigma.prior <- prior.p.dmc(
    dists = rep("gamma",length(p.vector)),
    p1   =rep(1, length(p.vector)),
    p2   =rep(1, length(p.vector)),
    lower=rep(NA,length(p.vector)),
    upper=rep(NA,length(p.vector))
  )
  names(sigma.prior) <- names(p.prior)
  # par(mfcol=c(2,3)); for (i in names(p.prior)) plot.prior(i,sigma.prior)
  pp.prior <- list(p.prior, sigma.prior)
  return(pp.prior)
}

addStimSetInfo <- function(x, input, orig_dat, addColumns=NULL, nBins=10) {
  df <- input[[x]]
  df$stimSet <- orig_dat[orig_dat$sub==x, 'stimulus_set']
  df$ease <- orig_dat[orig_dat$sub==x, 'ease']
  df$trialNstimset <- NA
  df$bin <- NA
  
  # make data compatible with pp format
  if(!'reps' %in% colnames(df)) {
    df$reps = 1
  }
  
  # add trialN by stim set, and bin, for first rep
  idx = df$reps == 1
  for(stimSet in unique(df[idx, 'stimSet'])) {
    idx2 <- df$stimSet == stimSet
    df[idx & idx2, 'trialNstimSet'] <- seq(1, sum(idx & idx2))
    df[idx & idx2, 'bin'] <- as.numeric(cut(seq(1, sum(idx & idx2)), nBins))
  }
  
  # add other columns from original data if wanted/needed
  if(!is.null(addColumns)) {
    for(colName in addColumns) {
      df[idx, colName] <- orig_dat[orig_dat$sub==x, colName]
    }
  }

  # copy&paste for all other reps
  if(max(df$reps) > 1) {
    df$trialNstimSet <- rep(df[idx, 'trialNstimSet'], max(df$reps))
    df$bin <- rep(df[idx, 'bin'], max(df$reps))
    
    if(!is.null(addColumns)) {
      for(colName in addColumns) {
        df[, colName] <- rep(df[idx, colName], max(df$reps))
      }
    }
    
  }
  
  df
}

calculateByBin <- function(df) {
  df$acc <- as.integer(df$R)==2
  
  attr(df, 'qRTs') <- do.call(data.frame, aggregate(RT~reps*bin, df, quantile, probs=seq(.1, .9, .4))) #cbind(quants[quants$acc==1,c('reps', 'bin')], quants[quants$acc==1,'RT'][,1])
  attr(df, 'qRTsCorrect') <- do.call(data.frame, aggregate(RT~reps*bin, df[df$acc==1,], quantile, probs=seq(.1, .9, .4))) #cbind(quants[quants$acc==1,c('reps', 'bin')], quants[quants$acc==1,'RT'][,1])
  attr(df, 'qRTsError') <- do.call(data.frame, aggregate(RT~reps*bin, df[df$acc==0,], quantile, probs=seq(.1, .9, .4))) #cbind(quants[quants$acc==1,c('reps', 'bin')], quants[quants$acc==1,'RT'][,1])
  attr(df, 'qRTsCorrectByEase') <- do.call(data.frame, aggregate(RT~reps*bin*ease, df[df$acc==1,], quantile, probs=seq(.1, .9, .4))) #cbind(quants[quants$acc==1,c('reps', 'bin')], quants[quants$acc==1,'RT'][,1])
  attr(df, 'qRTsErrorByEase') <- do.call(data.frame, aggregate(RT~reps*bin*ease, df[df$acc==0,], quantile, probs=seq(.1, .9, .4))) #cbind(quants[quants$acc==1,c('reps', 'bin')], quants[quants$acc==1,'RT'][,1])
  
  attr(df, 'RTsOverBins') <- aggregate(RT~reps*bin, df, mean)
  attr(df, 'AccOverBins') <- aggregate(acc~reps*bin, df, mean)
  attr(df, 'SkewOverBins') <- aggregate(RT~reps*bin, df, skewness)
  attr(df, 'RTsByChoiceByEase') <- aggregate(RT~reps*bin*ease*R, df, mean)
  attr(df, 'RTsByEase') <- aggregate(RT~reps*ease*bin, df, mean)
  attr(df, 'AccByEase') <- aggregate(acc~reps*ease*bin, df, mean)
  attr(df, 'SkewByEase') <- aggregate(RT~reps*ease*bin, df, skewness)
  
  
  if(!is.null(attr(df, 'adapt'))) {
    adapted <- attr(df, 'adapt')
    adapted$ease <- df$ease
    adapted$bin <- df$bin
    for(colName in colnames(adapted)) {
      form <- as.formula(paste0(colName, '~reps*ease*bin'))
      attr(df, paste0(colName, 'OverBins')) <- aggregate(form, adapted, mean)
    }
  }
  df
}

getDescriptives <- function(x, dep.var='RT', attr.name='RTsOverBins', id.var1='~reps*bin*s', id.var2='~reps*bin') {
  #  print('huh')
  allOverTime <- do.call(rbind, (lapply(1:length(x), function(y) {tmp <- attr(x[[y]], attr.name); tmp$s <- y; tmp})))
  #  print(head(allOverTime))
  
  form1 <- as.formula(paste0(dep.var, id.var1))
  res <- aggregate(form1, allOverTime, mean)
  
  if(!is.null(id.var2)) {
    form2 <- as.formula(paste0(dep.var, id.var2))
    res <- aggregate(form2, aggregate(form1, res, mean), mean)
  }
  return(res)
  #   if(!is.null(res2)) {
  #     form2 <- as.formula(paste0(dep.var, id.var2))
  #     res <- aggregate(form2, aggregate(form1, res, mean), mean)
  # }
  # 
  # if('reps' %in% colnames(allOverTime)) {
  #   # is model
  #   form1 <- as.formula(paste0(dep.var, id.var1))
  #   form2 <- as.formula(paste0(dep.var, id.var2))
  #   meanOverTime <- aggregate(form2, aggregate(form1, allOverTime, mean), mean)
  # }
  # else {
  #   form1 <- as.formula(paste0(dep.var, id.var1))
  #   form2 <- as.formula(paste0(dep.var, id.var2))
  #   meanOverTime <- aggregate(form2, aggregate(form1, allOverTime, mean), mean)
  # }
  # return(meanOverTime)
}

plotDataPPBins <- function(data, pp, dep.var='RT', xaxis='bin', colorM='blue', colorD=1, 
                           draw.legend=TRUE, legend.pos='topright', plot.new=TRUE, draw.polygon=TRUE,
                           xlim=NULL, ylim=NULL, xlab=NULL, ylab=NULL, data.cex=2, data.lwd=2,
                           hline.by=0.1, axhlines=NULL, plot.model.points=TRUE,
                           vline.by=2, axvlines=NULL, xaxt='s', yaxt='s') {
  if(is.null(xlim))  xlim = range(data[,xaxis])+c(-.5, .5)
  if(is.null(ylim)) {
    ylimD = range(data[,dep.var])*c(.9, 1.1)
    ylimM = range(pp[,dep.var])*c(.9, 1.1)
    ylim = c(min(ylimD[1], ylimM[1]), max(ylimD[2], ylimM[2]))
  }
  if(is.null(xlab)) xlab <- xaxis
  if(is.null(ylab)) ylab <- dep.var
  
  # empty canvas
  if(plot.new) {
    plot(0,0, type='n', xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, xaxt=xaxt, yaxt=yaxt)
    if(!is.null(axhlines)) {
      abline(h=axhlines, col='lightgrey', lty=1)
    } else {
      abline(h=seq(0, 5, hline.by), col='lightgrey', lty=1)
    }
    if(!is.null(axvlines)) {
      abline(v=axvlines, col='lightgrey', lty=1)
    } else {
      abline(v=seq(0, 20, vline.by), col='lightgrey', lty=1)
    }
  }
  
  if(draw.polygon) {
    lowerQ <- aggregate(as.formula(paste0(dep.var, '~bin')), pp, quantile, .025)
    upperQ <- aggregate(as.formula(paste0(dep.var, '~bin')), pp, quantile, .975)
    xs <- c(lowerQ[,xaxis], rev(lowerQ[,xaxis]))
    ys <- c(lowerQ[,dep.var], rev(upperQ[,dep.var]))
    polygon(xs, ys, col=rgb(col2rgb(colorM)[1]/255, col2rgb(colorM)[2]/255, col2rgb(colorM)[3]/255, alpha=.3), lty = NULL, border=NA)
  }
  
  # model
  if(plot.model.points) points(pp[,xaxis], pp[,dep.var], pch=20, col=colorM, cex=.01)
  
  # data
  points(data[,xaxis], data[,dep.var], col=colorD, pch=20, cex=data.cex)
  lines(data[,xaxis], data[,dep.var], col=colorD, lwd=data.lwd)
  if(!draw.polygon) {
    if(draw.legend) legend(legend.pos, c('Data', 'Model'), lty=c(1, NA), lwd=c(2, 2), col=c(colorD,colorM), pch=c(20, 20), bty='n')
  } else {
    if(draw.legend) legend(legend.pos, c('Data', 'Model'), lty=c(1, NA), 
                           fill=c(NA, rgb(col2rgb(colorM)[1]/255, col2rgb(colorM)[2]/255, col2rgb(colorM)[3]/255, alpha=.3)),
                           border=c(NA, 'black'),
                           lwd=c(2, NA), col=c(colorD,colorM), pch=c(20, NA), bty='n')
  }
}

plotFits <- function(samples, dat, nCores=30) {
  # Fit: overall & per participant
  pp = h.post.predict.dmc(samples = samples, adapt=TRUE, save.simulation = TRUE, cores=nCores)
  library(moments)
  nBins <- 10
  data <- lapply(samples, function(x) x$data)
  pp2 <- lapply(1:length(pp), addStimSetInfo, input=pp, orig_dat=dat)
  data2 <- lapply(1:length(data), addStimSetInfo, input=data, orig_dat=dat)
  if(!sfIsRunning()) sfInit(TRUE, nCores); sfLibrary(moments);
  pp3 <- sfLapply(pp2, calculateByBin)
  data3 <- lapply(data2, calculateByBin)
  
  # get quantiles
  q10RTs <- list(getDescriptives(data3, dep.var='RT.10.', attr.name='qRTsCorrect', id.var1='~bin', id.var2=NULL),
                 getDescriptives(pp3, dep.var='RT.10.', attr.name='qRTsCorrect', id.var1='~reps*bin', id.var2=NULL))
  q50RTs <- list(getDescriptives(data3, dep.var='RT.50.', attr.name='qRTsCorrect', id.var1='~bin', id.var2=NULL),
                 getDescriptives(pp3, dep.var='RT.50.', attr.name='qRTsCorrect', id.var1='~reps*bin', id.var2=NULL))
  q90RTs <- list(getDescriptives(data3, dep.var='RT.90.', attr.name='qRTsCorrect', id.var1='~bin', id.var2=NULL),
                 getDescriptives(pp3, dep.var='RT.90.', attr.name='qRTsCorrect', id.var1='~reps*bin', id.var2=NULL))
  q10RTsE <- list(getDescriptives(data3, dep.var='RT.10.', attr.name='qRTsError', id.var1='~bin', id.var2=NULL),
                  getDescriptives(pp3, dep.var='RT.10.', attr.name='qRTsError', id.var1='~reps*bin', id.var2=NULL))
  q50RTsE <- list(getDescriptives(data3, dep.var='RT.50.', attr.name='qRTsError', id.var1='~bin', id.var2=NULL),
                  getDescriptives(pp3, dep.var='RT.50.', attr.name='qRTsError', id.var1='~reps*bin', id.var2=NULL))
  q90RTsE <- list(getDescriptives(data3, dep.var='RT.90.', attr.name='qRTsError', id.var1='~bin', id.var2=NULL),
                  getDescriptives(pp3, dep.var='RT.90.', attr.name='qRTsError', id.var1='~reps*bin', id.var2=NULL))
  meanAcc <- list(getDescriptives(data3, dep.var='acc', attr.name='AccOverBins', id.var1='~bin', id.var2=NULL),
                  getDescriptives(pp3, dep.var='acc', attr.name='AccOverBins', id.var1='~reps*bin', id.var2=NULL))
  
  par(mfrow=c(1,3), mar=c(5,4,4,1)+.1, oma=c(0,0,1,0))
  plotDataPPBins(data=meanAcc[[1]], pp=meanAcc[[2]],
                 xaxt='s',
                 dep.var='acc', ylab='Accuracy', xlab = 'Trial bin',
                 legend.pos='bottomright', ylim=c(0.5, 0.95), hline.by=0.1)
  title('Accuracy')
  
  plotDataPPBins(data=q10RTs[[1]], pp=q10RTs[[2]], dep.var='RT.10.', ylim=c(.4, 1.2), xaxt='s', ylab='RT (s)', xlab='Trial bin')
  plotDataPPBins(data=q50RTs[[1]], pp=q50RTs[[2]], dep.var='RT.50.', plot.new = FALSE, draw.legend=FALSE)
  plotDataPPBins(data=q90RTs[[1]], pp=q90RTs[[2]], dep.var='RT.90.', plot.new = FALSE, draw.legend=FALSE)
  title('Correct RTs')
  
  plotDataPPBins(data=q10RTsE[[1]], pp=q10RTsE[[2]], dep.var='RT.10.', ylim=c(.4, 1.2), xaxt='s', ylab='RT (s)', xlab='Trial bin')
  plotDataPPBins(data=q50RTsE[[1]], pp=q50RTsE[[2]], dep.var='RT.50.', plot.new = FALSE, draw.legend=FALSE)
  plotDataPPBins(data=q90RTsE[[1]], pp=q90RTsE[[2]], dep.var='RT.90.', plot.new = FALSE, draw.legend=FALSE)
  title('Error RTs')
  title('Across subjects', outer=TRUE)
  
  
  ## by subject
  q10RTs <- list(getDescriptives(data3, dep.var='RT.10.', attr.name='qRTsCorrect', id.var1='~bin*s', id.var2=NULL),
                 getDescriptives(pp3, dep.var='RT.10.', attr.name='qRTsCorrect', id.var1='~reps*bin*s', id.var2=NULL))
  q50RTs <- list(getDescriptives(data3, dep.var='RT.50.', attr.name='qRTsCorrect', id.var1='~bin*s', id.var2=NULL),
                 getDescriptives(pp3, dep.var='RT.50.', attr.name='qRTsCorrect', id.var1='~reps*bin*s', id.var2=NULL))
  q90RTs <- list(getDescriptives(data3, dep.var='RT.90.', attr.name='qRTsCorrect', id.var1='~bin*s', id.var2=NULL),
                 getDescriptives(pp3, dep.var='RT.90.', attr.name='qRTsCorrect', id.var1='~reps*bin*s', id.var2=NULL))
  q10RTsE <- list(getDescriptives(data3, dep.var='RT.10.', attr.name='qRTsError', id.var1='~bin*s', id.var2=NULL),
                  getDescriptives(pp3, dep.var='RT.10.', attr.name='qRTsError', id.var1='~reps*bin*s', id.var2=NULL))
  q50RTsE <- list(getDescriptives(data3, dep.var='RT.50.', attr.name='qRTsError', id.var1='~bin*s', id.var2=NULL),
                  getDescriptives(pp3, dep.var='RT.50.', attr.name='qRTsError', id.var1='~reps*bin*s', id.var2=NULL))
  q90RTsE <- list(getDescriptives(data3, dep.var='RT.90.', attr.name='qRTsError', id.var1='~bin*s', id.var2=NULL),
                  getDescriptives(pp3, dep.var='RT.90.', attr.name='qRTsError', id.var1='~reps*bin*s', id.var2=NULL))
  meanAcc <- list(getDescriptives(data3, dep.var='acc', attr.name='AccOverBins', id.var1='~bin*s', id.var2=NULL),
                  getDescriptives(pp3, dep.var='acc', attr.name='AccOverBins', id.var1='~reps*bin*s', id.var2=NULL))
  
  par(mfrow=c(1,3), mar=c(5,4,4,1)+.1)
  for(subject in unique(q10RTs[[1]]$s)) {
    idxD <- q10RTs[[1]]$s == subject
    idxM <- q10RTs[[2]]$s == subject
    
    plotDataPPBins(data=meanAcc[[1]][idxD,], pp=meanAcc[[2]][idxM,],
                   xaxt='s',
                   dep.var='acc', ylab='Accuracy', xlab = 'Trial bin',
                   legend.pos='bottomright', ylim=c(0.5, 1), hline.by=0.1)
    title('Accuracy')
    
    plotDataPPBins(data=q10RTs[[1]][idxD,], pp=q10RTs[[2]][idxM,], dep.var='RT.10.', ylim=c(.3, 1.5), xaxt='s', ylab='RT (s)', xlab='Trial bin')
    plotDataPPBins(data=q50RTs[[1]][idxD,], pp=q50RTs[[2]][idxM,], dep.var='RT.50.', plot.new = FALSE, draw.legend=FALSE)
    plotDataPPBins(data=q90RTs[[1]][idxD,], pp=q90RTs[[2]][idxM,], dep.var='RT.90.', plot.new = FALSE, draw.legend=FALSE)
    title('Correct RTs')
    
    # some subs didnt make errors in some bin(s) - get idx again..
    idxD <- q10RTsE[[1]]$s == subject
    idxM <- q10RTsE[[2]]$s == subject
    plotDataPPBins(data=q10RTsE[[1]][idxD,], pp=q10RTsE[[2]][idxM,], dep.var='RT.10.', ylim=c(.3, 1.5), xaxt='s', ylab='RT (s)', xlab='Trial bin')
    plotDataPPBins(data=q50RTsE[[1]][idxD,], pp=q50RTsE[[2]][idxM,], dep.var='RT.50.', plot.new = FALSE, draw.legend=FALSE)
    plotDataPPBins(data=q90RTsE[[1]][idxD,], pp=q90RTsE[[2]][idxM,], dep.var='RT.90.', plot.new = FALSE, draw.legend=FALSE)
    title('Error RTs')
    title(paste0('Subject ', subject), outer=TRUE)
  }
}

runDiagnostics <- function(fn, dat, samples=NULL, nCores=30, verbose=TRUE, chainsFileType='jpeg',
                           plotChains=TRUE, plotPosteriorPriors=TRUE, calculateESS=TRUE, plotFits=TRUE, calculateParametersBPIC=TRUE) {
  if(is.null(samples)) samples <- loadSamples(fn, samplesDir='samples')
  
  saveDir <- file.path('diagnostics', fn)
  dir.create(saveDir, showWarnings = FALSE)
  if(chainsFileType == 'jpeg') dir.create(file.path(saveDir, 'chains'), showWarnings = FALSE)

  # Chains: Plot hyper, append Gelman's
  if(plotChains) {
    if(verbose) print('Plotting chains...')
    gms <- h.gelman.diag.dmc(samples)
    if(chainsFileType == 'pdf') {
      pdf(file=file.path(saveDir, 'chains.pdf'))
    } else {
      jpeg(filename=file.path(saveDir, 'chains', 'chains-%03d.jpeg'), width=7, height=7, units='in', quality=100, res=200)
    }
    par(mpg=c(1,1,0), mar=c(4,3,3,1)+.1, oma=c(0,0,1,0), mfrow=c(4,4))
    plot.dmc(samples, hyper=TRUE, layout=c(4,4))
    title(paste0('Hyper, Gelmans diag = ', round(gms['hyper'], 4)), outer=TRUE, line=0)
    
    # Chains: Plot per participant, append Gelman's
    for(subject in names(samples)) {
      plot.dmc(samples, subject=subject, layout=c(4,4)); 
      title(paste0('Subject ', subject, ' Gelmans diag = ', round(gms[which(names(gms)==as.character(subject))], 4)), outer=TRUE, line=0)
    }
    
    
    dev.off()
  }
  
  # Posterior vs Prior
  if(plotPosteriorPriors) {
    if(!grepl('softmax', fn)) {
      pdf(file=file.path(saveDir, 'posteriorPrior.pdf'))
      if(verbose) print('Plotting posterior vs prior...')
      par(mpg=c(1,1,0), mar=c(4,3,3,1)+.1, oma=c(0,0,1,0), mfrow=c(4,4))
      plot.dmc(samples, hyper=TRUE, p.prior = attr(samples, 'hyper')$pp.prior, layout=c(4,4))
      for(subject in names(samples)) {
        plot.dmc(samples, subject=subject, p.prior = samples[[subject]]$p.prior, layout=c(4,4))
        title(paste0('Subject ', subject, ' Gelmans diag = ', round(gms[which(names(gms)==as.character(subject))], 4)), outer=TRUE, line=0)
      } 
      dev.off()
    } else {
      warning('Skipping posterior vs prior plot, this crashes for softmax... (too few parameters?)')
    }
  }
  
  # Effective sample size: save to csv
  if(calculateESS) {
    if(verbose) print('Calculating effective samples sizes...')
    effS <- data.frame(do.call(rbind, effectiveSize.dmc(samples)))
    effSH <- effectiveSize.dmc(samples, hyper=TRUE)
    effS <- rbind(effS, effSH[grepl('.h1', names(effSH))])
    effS <- rbind(effS, effSH[grepl('.h2', names(effSH))])
    effS <- rbind(effS, apply(effS, 2, min))
    row.names(effS) <- c(row.names(effS)[1:(nrow(effS)-3)], 'h1', 'h2', 'minimum')
    write.csv(effS, file=file.path(saveDir, 'effectiveSizes.csv'))
  }
  
  # Fits
  if(plotFits) {
    if(verbose) print('Plotting posterior predictives...')
    pdf(file=file.path(saveDir, 'posteriorPredictives.pdf'))
    plotFits(samples, dat, nCores)
    dev.off()
  }
  
  # Parameters: save to csv, append BPIC per participant
  if(calculateParametersBPIC) {
    if(verbose) print('Getting parameters, BPIC...')
    bpics <- h.IC.dmc(samples)
    summ <- summary.dmc(samples)
    medians <- data.frame(do.call(rbind, lapply(summ, function(x) x$quantiles[,3])))
    medians$aV <- pnorm(medians$aV)
    medians <- rbind(medians, apply(medians, 2, mean))
    medians <- rbind(medians, apply(medians, 2, sd))
    
    medians$minimum.deviances <- c(bpics[,1], sum(bpics[,1]), NA)
    medians$BPIC <- c(bpics[,2], sum(bpics[,2]), NA)
    row.names(medians) <- c(row.names(medians)[1:(nrow(medians)-2)], 'mean/sum', 'sd')
    write.csv(medians, file=file.path(saveDir, 'medianParametersBPICs.csv'))
  }
}

pp.summary <- function(sim, samples, dat=NULL, n.post=100, probs=c(1:99)/100,random=TRUE,
                       bw="nrd0",report=10,save.simulation=FALSE,factors=NA, cvs=NULL,adapt=FALSE,
                       save.simulation.as.attribute=FALSE,ignore.R2=FALSE,censor=c(NA,NA),
                       gglist=TRUE, probs.gglist=c(0.1, 0.5, 0.9),CI.gglist=c(0.025, 0.975))
{
  model <- attributes(samples$data)$model
  facs <- names(attr(model,"factors"))
  if(is.null(cvs)) { 
    if(!is.null(attr(samples$data, 'cvs'))) {
      cvs <- attr(samples$data, 'cvs')
    } else {
      cvs <- samples$data[,attr(model,"cvs"),drop=FALSE]
    }
  }
  attr(cvs,"row.facs") <- apply(apply(
    samples$data[,facs,drop=FALSE],2,as.character),1,paste,collapse=".")
  if ( ignore.R2 & any(names(samples$data)=="R2") )
    samples$data <- samples$data[,names(samples$data)[names(samples$data)!="R2"]]
  if (!is.null(factors) ) {
    if (any(is.na(factors))) factors <- facs
    if (!all(factors %in% facs)) 
      stop(paste("Factors argument must contain one or more of:",paste(facs,collapse=",")))
  }
  resp <- names(attr(model,"responses"))
  ns <- table(samples$data[,facs],dnn=facs)
  n.par <- dim(samples$theta)[2]
  thetas <- matrix(aperm(samples$theta,c(3,1,2)),ncol=n.par)
  colnames(thetas) <- dimnames(samples$theta)[[2]]
  if (is.na(n.post)) use <- c(1:dim(thetas)[1]) else {
    if (random) use <- sample(c(1:dim(thetas)[1]),n.post,replace=F) else
      use <- round(seq(1,dim(thetas)[1],length.out=n.post))
  }
  n.post <- length(use)
  posts <- thetas[use,]
  n.rep <- sum(ns)
  # sim <- data.frame(matrix(nrow=n.post*n.rep,ncol=dim(samples$data)[2]))
  #  names(sim) <- names(samples$data)  
  # Tweaks for Stop Signal
  if ( !any(names(samples$data)=="SSD") ) {
    SSD <- rep(Inf,sum(ns)) 
    leave.out <- -c(1:dim(samples$data)[2])[names(samples$data) %in% c("RT",names(cvs),"R2")]
    TRIALS <- NULL
  } else {
    # Assumes last two are SSD and RT! FIX ME. EG WONT WORK IF THERE ARE CVS
    if ( is.null(facs) ) {
      SSD <- samples$data$SSD 
      if (any(names(samples$data)=="TRIALS")) 
        TRIALS <- samples$data$TRIALS else 
          TRIALS <- NULL
    } else {
      SSD <- unlist(tapply(samples$data$SSD,samples$data[,facs],identity)) 
      if ( any(names(samples$data)=="TRIALS") ) 
        TRIALS <- unlist(tapply(samples$data$TRIALS,samples$data[,facs],identity)) else 
          TRIALS <- NULL
    }
    leave.out <- -c(1:dim(samples$data)[2])[names(samples$data) %in% c("RT","SSD")]
    # leave.out <- -c((dim(samples$data)[2]-1):dim(samples$data)[2])
  }
  # cat("\n")
  # cat(paste("Simulating (\'.\'=",report,"): ",sep=""))
  # for (i in names(samples$data)[leave.out])
  #   sim[[i]] <- factor(rep(NA,n.post*n.rep),levels=levels(samples$data[[i]]))
  # if(adapt & save.simulation) {
  #   adapted <- data.frame()[1:nrow(sim),]
  # }
  # for (i in 1:n.post) {
  #   tmp <- simulate.dmc(posts[i,],model,n=ns,SSD=SSD,cvs=cvs,adapt=adapt,TRIALS=TRIALS)
  #   if (ignore.R2) tmp <- tmp[,names(tmp)[names(tmp)!="R2"]]
  #   sim[(1+(i-1)*n.rep):(i*n.rep),names(tmp)] <- tmp
  #   if(adapt & save.simulation) {
  #     adapt.out <- attr(tmp, 'adapt')
  #     adapted[(1+(i-1)*n.rep):(i*n.rep),names(adapt.out)] <- adapt.out
  #   }
  #   if ( (i %% report) == 0) cat(".")
  # }
  # cat("\n")
  # if ( any(names(sim)=="R2") ) { # MTR model
  #   levs <- outer(levels(samples$data$R),sort(unique(samples$data$R2)),"paste",sep="")
  #   if (attributes(model)$type=="normDK") levs[,2] <- rep("DK",dim(levs)[1])
  #   levs <- sort(unique(as.vector(levs)))  
  #   sim$R <- paste(as.character(sim$R),as.character(sim$R2),sep="") 
  #   if (attributes(model)$type=="normDK") sim$R[sim$R2=="2"] <- "DK"
  #   sim$R <- factor(sim$R,levels=levs)
  #   samples$data$R <- paste(as.character(samples$data$R),as.character(samples$data$R2),sep="") 
  #   if (attributes(model)$type=="normDK") 
  #     samples$data$R[samples$data$R2=="2"] <- "DK"
  #   samples$data$R <- factor(samples$data$R,levels=levs)
  # }
  # reps <- rep(1:n.post,each=dim(samples$data)[1])
  # 
  # if (!is.na(censor[1])) {
  #   fast <- sim[,"RT"] < censor[1] 
  #   cat(paste("Removing",sum(fast),"(",round(100*mean(fast),2))," %) fast RTs\n")
  # } else fast <- rep(FALSE,dim(sim)[1])
  # if (!is.na(censor[2])) {
  #   slow <- sim[,"RT"] > censor[2] 
  #   cat(paste("Removing",sum(slow),"(",round(100*mean(slow),2))," %) slow RTs\n")
  # } else slow <- rep(FALSE,dim(sim)[1])
  # ok <- !fast & !slow 
  # sim <- sim[ok,]
  # reps <- reps[ok]
  # 
  
  reps <- sim[,1]
  sim <- sim[,2:ncol(sim)]
  # if ( save.simulation ) {
  #   sim <- cbind(reps,sim)
  #   if(adapt) {
  #     adapted <- cbind(reps, adapted)
  #     attr(sim, 'adapt') <- adapted
  #   }
  #   attr(sim,"data") <- samples$data
  #   sim
  # } else {
  sim.dqp <- get.dqp(sim,facs=factors,probs,n.post,ns=ns,bw=bw,reps=reps)
  dat.dqp <- get.dqp(sim=samples$data,factors,probs,bw=bw)
  names(dat.dqp) <- paste("data",names(dat.dqp),sep=".")
  out <- c(sim.dqp,dat.dqp)
  dpqs <- vector(mode="list",length=length(n.post))
  for (i in 1:n.post) {
    simi <- sim[reps==i,]
    dpqs[[i]] <- get.dqp(simi,factors,probs,1)
  }
  attr(out,"dpqs") <- dpqs
  if (save.simulation.as.attribute) 
    attr(out,"sim") <- cbind(reps,sim)
  if (gglist) attr(out, "gglist") <- 
    get.fitgglist.dmc(sim=cbind(reps,sim),data=samples$data,factors=factors, noR=FALSE, 
                      quantiles.to.get= probs.gglist, CI = CI.gglist)
  out
  # }
}

h.pp.summary <- function(sim, samples, dat=NULL, ignore.subjects=FALSE,
                         n.post=100,probs=c(1:99)/100,bw="nrd0",
                         save.simulation=FALSE,factors=NA,save.subject.posts=FALSE,adapt=TRUE,
                         cores=1,ignore.R2=FALSE, 
                         gglist= FALSE, probs.gglist =c(0.1, 0.5, 0.9), CI.gglist =c(0.025, 0.975),
                         censor=c(NA,NA))
  # apply post.predict to each subject
{
  # out = sim
  if(is.null(dat)) {
    dat <- lapply(samples,function(x){x$data})
  }
  
  if(!ignore.subjects) {
    if(cores == 1) {
      out = lapply(1:length(sim), function(x) pp.summary(sim=sim[[x]], samples=samples[[x]], dat=dat[[x]],
                                                         n.post=n.post,probs=probs,bw=bw,
                                                         factors=factors,save.simulation=save.simulation,gglist=FALSE,
                                                         adapt=adapt,
                                                         save.simulation.as.attribute=FALSE,ignore.R2=ignore.R2,censor=censor))
    } else if(cores > 1) {
      if(!sfIsRunning()) {
        sfInit(parallel = TRUE, cpus=cores)
        wasRunning = TRUE
      } else {
        wasRunning = FALSE
      }
      sfExport('pp.summary', 'get.dqp')
      out = sfLapply(1:length(sim), function(x) pp.summary(sim=sim[[x]], samples=samples[[x]], dat=dat[[x]],
                                                         n.post=n.post,probs=probs,bw=bw,
                                                         factors=factors,save.simulation=save.simulation,gglist=FALSE,
                                                         adapt=adapt,
                                                         save.simulation.as.attribute=FALSE,ignore.R2=ignore.R2,censor=censor))
      if(!wasRunning) sfStop()
    }
    names(out) <- names(sim)
    sim <- do.call(rbind, sim)
    dat <- do.call(rbind, dat)
  } else {
    out <- sim <- do.call(rbind, sim)
  }
  
  # sim <- do.call(rbind,lapply(out,function(x){attr(x,"sim")}))
  if ( (any(names(samples[[1]]$data)=="R2")) && !ignore.R2 ) for (i in 1:length(samples)) {
    levs <- outer(levels(samples[[i]]$data$R),sort(unique(samples[[i]]$data$R2)),"paste",sep="")
    if (attributes(attributes(samples[[i]]$data)$model)$type=="normDK") 
      levs[,2] <- rep("DK",dim(levs)[1])
    levs <- sort(unique(as.vector(levs)))  
    
    samples[[i]]$data$R <- 
      paste(as.character(samples[[i]]$data$R),as.character(samples[[i]]$data$R2),sep="") 
    if (attributes(attributes(samples[[i]]$data)$model)$type=="normDK") 
      samples[[i]]$data$R[samples[[i]]$data$R2=="2"] <- "DK"
    samples[[i]]$data$R <- factor(samples[[i]]$data$R,levels=levs)
  }
  facs <- names(attr(attributes(samples[[1]]$data)$model,"factors"))
  if (!is.null(factors)) {
    if (any(is.na(factors))) factors <- facs
    if (!all(factors %in% facs)) 
      stop(paste("Factors argument must contain one or more of:",paste(facs,collapse=",")))
  }
  sim.dqp <- get.dqp(sim[,-1],factors,probs,n.post=1,bw=bw)
  dat.dqp <- get.dqp(sim=dat,factors,probs,bw=bw)
  names(dat.dqp) <- paste("data",names(dat.dqp),sep=".")
  av <- c(sim.dqp,dat.dqp)
  dpqs <- vector(mode="list",length=length(n.post))
  for (i in 1:n.post) {
    simi <- sim[sim$reps==i,-1]
    dpqs[[i]] <- get.dqp(simi,factors,probs,1)
  }
  attr(av,"dpqs") <- dpqs
  # Strip global sim attribute and dpqs for each participant
  out <- lapply(out,function(x){
    attr(x,"sim") <- NULL
    if (!save.subject.posts) attr(x,"dpqs") <- NULL
    x
  })
  if (gglist) attr(av, "gglist") <- get.fitgglist.dmc(sim,dat,factors=factors, noR=FALSE, 
                                                      quantiles.to.get = probs.gglist, CI= CI.gglist)
  attr(out,"av") <- av
  out
}



loadSamples <- function(fn, samplesDir=NULL) {
  if(!is.null(samplesDir)) fn <- file.path(samplesDir, fn)
  fnExt <- paste0(fn, '.RData')
  load(fnExt)
  return(hsamples)
  # fnAutoConverge <- gsub('.RData', '_autoconverge.RData', fnExt)
  # hasLoaded <- FALSE
  # if(file.exists(fnAutoConverge)) {
  #   if(file.info(fnAutoConverge)$mtime > file.info(fnExt)$mtime) {
  #     warning('Autoconverge is newest samples available. Sampling probably hasn\'t finished yet')
  #     load(fnAutoConverge)
  #     loadFinal <- TRUE
  #   }}
  # if(!hasLoaded) {
  #   load(fnExt)
  #   if('samples' %in% ls()) {
  #     print('Loaded final samples')
  #   } else {
  #     warning('samples not found, only burn found...')
  #     return(burn1)
  #   }
  # }
  # return(samples)
}