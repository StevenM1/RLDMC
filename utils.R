doSample <- function(data, p.prior, pp.prior, fileName, restart=FALSE,
                     nCores=31, nmcBurn=150, nmc=250, nIter=100, useRUN=False) {
  if(useRUN) {
    samples <- h.samples.dmc(nmc=nmcBurn, p.prior=p.prior, data=data, pp.prior=pp.prior)
    samples <- h.RUN.dmc(hsamples = samples, cores=nCores, verbose=TRUE, saveFn=fileName)
    save(samples, save=paste0(saveFn, '.RData'))
  } else {
    fileNameFull <- paste0(fileName, '.RData')
    if(file.exists(fileNameFull)) {
      load(fileNameFull)
    }
    
    if(restart | !exists('burn1')) {
      burn <- h.samples.dmc(nmc=nmcBurn, p.prior=p.prior, data=data, pp.prior=pp.prior)
      # burn1 <- h.run.dmc(burn, cores=nCores, p.migrate=0.05, h.p.migrate = 0.05, report=10)
      burn1 <- h.run.unstuck.dmc(samples=burn, nmc=nmcBurn, cores=nCores, p.migrate=0.05, h.p.migrate=0.05, end.no.migrate=TRUE)
      save(burn1, file=fileNameFull)
    }
    
    if(restart | !exists('samples')) {
      samples <- burn1
    }
    
    # Sample from posterior, starting from the end result of burn1
    samples <- h.run.converge.dmc(samples=h.samples.dmc(nmc=nmc, samples=samples, add=TRUE),
                                  thorough=TRUE, nmc = nmc, cores=nCores,
                                  save=paste0(fileName, "_autoconverge"))
    save(burn1, samples, file=fileNameFull)
  }  
}

loadData <- function(name, subset=FALSE, addBins=FALSE, addStimSets=FALSE) {
  otherCols <- c()
  if(name == 'calibration_subset' | name == 'calibration' | name == 'calibration-subset') {
    load('../RLDDM/data/data_exp-all.Rdata')
    #  dat <- droplevels(dat[dat$pp==sub,])
    dat <- dat[dat$task=='Calibration' & !is.na(dat$pp),]
    dat$ease <- abs(dat$high_stim_prob-dat$low_stim_prob)
    
    if(subset | name == 'calibration_subset' | name == 'calibration-subset') {
      ppsNoLearning <- aggregate(choiceIsHighP~pp*ease, dat, mean)
      ppsNoLearning <- ppsNoLearning[ppsNoLearning$choiceIsHighP>0.95 | ppsNoLearning$choiceIsHighP<0.05,]
      idx <- rep(FALSE, nrow(dat))
      for(i in 1:nrow(ppsNoLearning)) {
        idx <- idx | (dat$pp == ppsNoLearning[i, 'pp'] & dat$ease == ppsNoLearning[i, 'ease'])
        print(dat[(dat$pp == ppsNoLearning[i, 'pp'] & dat$ease == ppsNoLearning[i, 'ease']), 'reward'][1:7])
      }
      dat <- droplevels(dat[!idx,])
    }
    #dat <- droplevels(dat[!dat$pp %in% removePps,])
    dat$sub <- as.factor(as.integer(dat$pp))
    # Remove RTs > 1.5 and RTs < .125
    dat <- droplevels(dat[dat$rt>.125,])
    dat <- dat[!is.na(dat$pp),]
    dat$sub <- as.factor(as.integer(as.factor(dat$pp)))
  } else if(name == 'annie') {
    load('../RLDDM/data/data_annie_N74.Rdat')
    dat$sub <- as.factor(as.numeric(as.factor(dat$workerId)))
    dat$reward <- dat$outcome_num
    dat$stimulus_set <- as.numeric(dat$stimulus_set)
    dat$ease <- round(abs(dat$prob_win_correct_rev - (1-dat$prob_win_correct_rev)),1)
    dat$RT <- dat$rt <- dat$rt/1000
    dat <- dat[!is.na(dat$RT),]
    dat <- dat[!is.na(dat$choices),]
    dat$choiceIsHighP <- dat$choices-1
  } else if(name == 'annie-chris' | name == 'annie_chris') {
    dat <- read.csv('../RLDDM/data/data_annie_Chris.csv')
    dat$sub <- as.factor(as.numeric(as.factor(dat$SubjID)))
    
    
    dat$reward <- dat$rewards
    dat$stimulus_set <- as.numeric(dat$stimulus_set)
    dat$ease <- round(abs(dat$prob_win_correct_rev - (1-dat$prob_win_correct_rev)),2)
    dat$choiceIsHighP <- dat$choices-1
    
    dat$block <- ifelse(dat$stimulus_set<2, 1, ifelse(dat$stimulus_set<4, 2, ifelse(dat$stimulus_set<6,3,4)))
    dat$RT <- dat$rt <-  as.numeric(as.character(dat$rt))
    dat <- dat[!is.na(dat$RT),]
    dat <- dat[!is.na(dat$choices),]
    
    # always exclude subjects 1 & 38 - experiment crashed
    dat <- dat[dat$sub != 38,]
    
    # exclude block 4 for subject 1
    dat <- dat[!(dat$sub==1 & dat$stimulus_set>5),]
    dat$sub <- as.factor(as.numeric(as.factor(dat$SubjID)))  # re-number
    dat <- droplevels(dat)
    
  } else if(name == 'SAT' | name=='sat') {
    load('../RLDDM/data/data_exp-pilot3.Rdata')
    #  dat <- droplevels(dat[dat$pp==sub,])
    #dat <- dat[dat$task=='Calibration' & !is.na(dat$pp),]
    # dat$ease <- abs(dat$high_stim_prob-dat$low_stim_prob)
    #dat <- droplevels(dat[!dat$pp %in% removePps,])
    dat$sub <- as.factor(as.integer(dat$pp))
    # Remove RTs > 1.5 and RTs < .125
    dat <- droplevels(dat[dat$rt>.125,])
    dat <- dat[!is.na(dat$pp),]
    dat$sub <- as.factor(as.integer(as.factor(dat$pp)))
    otherCols <- c('cue')
  }
  
  # remove RTs < .15
  dat <- dat[dat$rt>.15,]
  
  rownames(dat) <- NULL
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
  data <- dat
  data$s <- data$sub
  data$R <- factor(ifelse(data$choiceIsHighP==0, 'r1', 'r2'))
  data$RT <- dat$rt
  if(addStimSets) {
    data$S <- factor(paste('s', data$stimulus_set, sep=''), levels=paste0('s', 1:4))
  } else {
    data$S <- factor('s1')
  }
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

loadSamples <- function(fn, samplesDir) {
  fnExt <- file.path(samplesDir, paste0(fn, '.RData'))
  fnAutoConverge <- gsub('.RData', '_autoconverge.RData', fnExt)
  hasLoaded <- FALSE
  if(file.exists(fnAutoConverge)) {
    if(file.info(fnAutoConverge)$mtime > file.info(fnExt)$mtime) {
      warning('Autoconverge is newest samples available. Sampling probably hasn\'t finished yet')
      load(fnAutoConverge)
      loadFinal <- TRUE
    }}
  if(!hasLoaded) {
    load(fnExt)
    if('samples' %in% ls()) {
      print('Loaded final samples')
    } else {
      warning('samples not found, only burn found...')
      return(burn1)
    }
  }
  return(samples)
}

getPriors <- function(p.vector) {
  dists <- rep('tnorm', length(p.vector))
  p1 <- upper <- rep(NA, length(p.vector))
  p2 <- rep(5, length(p.vector))
  lower <- rep(0, length(p.vector))
  names(p1) <- names(p2) <- names(upper) <- names(lower) <- names(p.vector)
  
  # t0
  p1[names(p1)=='t0'] <- 0.2
  p2[names(p1)=='t0'] <- 0.5
  upper[names(p1)=='t0'] <- 1
  lower[names(p1)=='t0'] <- 0.025
  
  # stimulus or risk representations SR/RR
  p1[names(p1)=='SR' | names(p1) == 'RR'] = 0
  p2[names(p1)=='SR' | names(p1) == 'RR'] = 1
  lower[names(p1)=='SR' | names(p1) == 'RR'] = NA
  upper[names(p1)=='SR' | names(p1) == 'RR'] = NA
  
  # learning rate aV or aR
  p1[names(p1)=='aV' | names(p1)=='aR'] <- -1.6  # corresponds to 0.05
  lower[names(p1)=='aV' | names(p1)=='aR'] <- NA  # corresponds to 0
  
  # urgency v0 - may vary across conditions, so use grepl
  p1[grepl(pattern='V0', names(p1))] <- 2
  lower[grepl(pattern='V0', names(p1))] <- NA
  
  # weight on v
  p1[names(p1) == 'wV'] = 9
  
  # threshold - this may vary also across conditions, so use grepl
  p1[grepl(pattern='B0', names(p1))] = 3
  
  # weight of magnitude effect / sum of EVs
  p1[names(p1) == 'wS'] = 0
  p2[names(p2) == 'wS'] = 3
  lower[names(p2) == 'wS'] = NA
  
  # DDM parameters
  # threshold
  p1[names(p1) == 'a'] = 2
  p2[names(p1) == 'a'] = 0.5
  
  # linear effect m
  p1[names(p1) == 'm'] = 2
  p2[names(p1) == 'm'] = 5
  
  # sv
  p1[names(p1) == 'sv'] = 0.1
  p2[names(p1) == 'sv'] = 0.1
  
  # st0
  p1[names(p1) == 'st0'] = 0.1
  p2[names(p1) == 'st0'] = 0.1
  
  
  # softmax beta
  p1[names(p1) == 'beta'] = 1
  p2[names(p1) == 'beta'] = 1
  
  # lba sd_v
  p1[grepl(pattern='sd_v', x=names(p1))] = 1
  p2[grepl(pattern='sd_v', x=names(p1))] = 1
  lower[grepl(pattern='sd_v', x=names(p1))] = 0

  # lba A
  p1[grepl(pattern='A', x=names(p1))] = 1
  p2[grepl(pattern='A', x=names(p1))] = 1
  lower[grepl(pattern='A', x=names(p1))] = 0
  
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

addStimSetInfo <- function(x, input, orig_dat, addColumns=NULL) {
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

plotDataPPBins <- function(data, pp, dep.var='RT', xaxis='bin', colorM='blue', colorD=1, draw.legend=TRUE, plot.new=TRUE, draw.polygon=TRUE,
                           xlim=NULL, ylim=NULL) {
  if(is.null(xlim))  xlim = range(data[,xaxis])+c(-.5, .5)
  if(is.null(ylim)) {
    ylimD = range(data[,dep.var])*c(.9, 1.1)
    ylimM = range(pp[,dep.var])*c(.9, 1.1)
    ylim = c(min(ylimD[1], ylimM[1]), max(ylimD[2], ylimM[2]))
  }
  
  # empty canvas
  if(plot.new) {
    plot(0,0, type='n', xlim=xlim, ylim=ylim, xlab=xaxis, ylab=dep.var)
    abline(h=seq(0, 5, .1), col='grey', lty=2)
  }
  
  if(draw.polygon) {
    lowerQ <- aggregate(as.formula(paste0(dep.var, '~bin')), pp, quantile, .025)
    upperQ <- aggregate(as.formula(paste0(dep.var, '~bin')), pp, quantile, .975)
    xs <- c(lowerQ[,xaxis], rev(lowerQ[,xaxis]))
    ys <- c(lowerQ[,dep.var], rev(upperQ[,dep.var]))
    polygon(xs, ys, col=rgb(col2rgb(colorM)[1]/255, col2rgb(colorM)[2]/255, col2rgb(colorM)[3]/255, alpha=.3), lty = NULL, border=NA)
  }
  
  # model
  points(pp[,xaxis], pp[,dep.var], pch=20, col=colorM, cex=.01)
  
  # data
  points(data[,xaxis], data[,dep.var], col=colorD, pch=20, cex=2)
  lines(data[,xaxis], data[,dep.var], col=colorD, lwd=2)
  if(!draw.polygon) {
    if(draw.legend) legend('topright', c('Data', 'Model'), lty=c(1, NA), lwd=c(2, 2), col=c(colorD,colorM), pch=c(20, 20), bty='n')
  } else {
    if(draw.legend) legend('topright', c('Data', 'Model'), lty=c(1, NA), 
                           fill=c(NA, rgb(col2rgb(colorM)[1]/255, col2rgb(colorM)[2]/255, col2rgb(colorM)[3]/255, alpha=.3)),
                           border=c(NA, 'black'),
                           lwd=c(2, NA), col=c(colorD,colorM), pch=c(20, NA), bty='n')
  }
}

plotDataPPTrialN <- function(data, pp, dep.var='RT', xaxis='trialN_reversal', colorM=1, colorD=1, draw.legend=TRUE, plot.new=TRUE, draw.polygon=TRUE,
                           xlim=NULL, ylim=NULL) {
  if(is.null(xlim))  xlim = range(data[,xaxis])+c(-.5, .5)
  if(is.null(ylim)) {
    ylimD = range(data[,dep.var])*c(.9, 1.1)
    ylimM = range(pp[,dep.var])*c(.9, 1.1)
    ylim = c(min(ylimD[1], ylimM[1]), max(ylimD[2], ylimM[2]))
  }
  
  # empty canvas
  if(plot.new) {
    plot(0,0, type='n', xlim=xlim, ylim=ylim, xlab=xaxis, ylab=dep.var)
    abline(h=seq(0, 5, .1), col='grey', lty=2)
  }
  
  if(draw.polygon) {
    lowerQ <- aggregate(as.formula(paste0(dep.var, '~', xaxis)), pp, quantile, .025)
    upperQ <- aggregate(as.formula(paste0(dep.var, '~', xaxis)), pp, quantile, .975)
    xs <- c(lowerQ[,xaxis], rev(lowerQ[,xaxis]))
    ys <- c(lowerQ[,dep.var], rev(upperQ[,dep.var]))
    polygon(xs, ys, col=rgb(col2rgb('blue')[1]/255, col2rgb('blue')[2]/255, col2rgb('blue')[3]/255, alpha=.3), lty = NULL, border=NA)
  }
  
  # model
  points(pp[,xaxis], pp[,dep.var], pch=20, col=colorM, cex=.01)
  
  # data
  points(data[,xaxis], data[,dep.var], col=colorD, pch=20, cex=2)
  lines(data[,xaxis], data[,dep.var], col=colorD, lwd=2)
  
  if(!draw.polygon) {
    if(draw.legend) legend('topright', c('Data', 'Model'), lty=c(1, NA), lwd=c(2, 2), col=c(colorD,colorM), pch=c(20, 20), bty='n')
  } else {
    if(draw.legend) legend('topright', c('Data', 'Model'), lty=c(1, NA), 
                           fill=c(NA, rgb(col2rgb('blue')[1]/255, col2rgb('blue')[2]/255, col2rgb('blue')[3]/255, alpha=.3)),
                           border=c(NA, 'black'),
                           lwd=c(2, NA), col=c(colorD,colorM), pch=c(20, NA), bty='n')
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
    out = lapply(1:length(sim), function(x) pp.summary(sim=sim[[x]], samples=samples[[x]], dat=dat[[x]],
                                                       n.post=n.post,probs=probs,bw=bw,
                                                       factors=factors,save.simulation=save.simulation,gglist=FALSE,
                                                       adapt=adapt,
                                                       save.simulation.as.attribute=FALSE,ignore.R2=ignore.R2,censor=censor))
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