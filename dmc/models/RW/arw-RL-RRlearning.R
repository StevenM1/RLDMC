library(dmcAdapt)

# Template setup for n-choice LBA, B=b-A parameterization
#   External parameters types: A, B, t0, mean_v, sd_v, st0 = 0 (optional)
#   Internal parameters types: A, b, t0, mean_v, sd_v, st0 = 0 (optional)

# This ncl (no cell likelihood) version uses p.list.dmc in the likelihood to 
# avoid looping over design cells.

transform.dmc <- function(par.df,do.trans=TRUE) 
{
  par.df$aR[par.df$aR == Inf] <- par.df$aV[par.df$aR == Inf]
  par.df$aV = t(pnorm(t(par.df$aV)))  # value learning rate
  par.df$SR = t(pnorm(t(par.df$SR)))  # stimulus representations
  par.df$aR = t(pnorm(t(par.df$aR)))  # reward rate learning rate
  
  if (do.trans) {
    return(
      list(A=t(par.df$A),
           s=t(par.df$s),
           t0=t(par.df$t0),
           st0=t(par.df$st0),
           B0=t(par.df$B0),
           SR=t(par.df$SR),
           RR=t(par.df$RR),
           aV=t(par.df$aV),
           aR=t(par.df$aR),
           V0=t(par.df$V0),
           wV=t(par.df$wV),
           wR=t(par.df$wR)))
  } else {
    return(
      list(A=par.df$A,s=par.df$s,t0=par.df$t0,st0=par.df$st0,
           B0=par.df$B0,
           SR=par.df$SR,
           RR=par.df$RR,
           aV=par.df$aV,
           aR=par.df$aR,
           V0=par.df$V0,
           wV=par.df$wV,
           wR=par.df$wR))
  }
}

transform2.dmc <- function(pars, cvsValue, choiceIdxValue, cvsReward, choiceIdxReward) {
  ### SM
  
  # Create start point vector
  startValuesVal <- rep(pars[,1,'SR'], each=1, times=ncol(cvsValue)/2)
  startValuesRR <- rep(pars[,1,'RR'], each=1, times=ncol(cvsReward)/2)
  startValues <- c(startValuesVal, startValuesRR)
  
  # learning rates matrix
  learningRatesVal <- matrix(rep(pars[1,,'aV'], each=ncol(cvsValue)), ncol=ncol(cvsValue), byrow=TRUE)
  learningRatesRR <- matrix(rep(pars[1,,'aR'], each=ncol(cvsReward)), ncol=ncol(cvsReward), byrow=TRUE)
  learningRates <- cbind(learningRatesVal, learningRatesRR)
  
  # call C
  updated <- adapt.c.dmc(startValues = startValues, 
                         learningRates = learningRates, 
                         feedback = cbind(cvsValue, cvsReward),
                         learningRule='SARSA')
  
  # add back data to 'pars' array
  updatedValues <- updated$adaptedValues[,1:ncol(cvsValue)]
  updatedRR <- updated$adaptedValues[,1:ncol(cvsReward)+ncol(cvsValue)]
  pars[,,'SR'] <- matrix(updatedValues[choiceIdxValue], ncol=2, byrow=FALSE)
  RR <- matrix(updatedRR[choiceIdxReward], ncol=1, byrow=FALSE)
  
  # Determine parameter mean_v
  mean_v <- cbind(pars[1,,"V0"] + pars[1,,"wV"]*(pars[2,,"SR"]-pars[1,,"SR"]) + pars[1,,'wR'] * RR,
                  pars[2,,"V0"] + pars[2,,"wV"]*(pars[1,,"SR"]-pars[2,,"SR"]) + pars[2,,'wR'] * RR)
  
  # Determine parameter b
  b <- cbind(pars[1,,"B0"], pars[2,,"B0"])
  
  return(list(pars=pars, mean_v=mean_v, b=b, RR=RR, updated=updated))
}


random.dmc <- function(p.list,model,save.adapt=TRUE)
{
  pars <- array(unlist(p.list,use.names=FALSE),
                dim=c(dim(p.list[[1]]),length(p.list)),
                dimnames=list(NULL,NULL,names(p.list)))
  cvs <- attr(p.list,"cvs")
  n <- dim(p.list[[1]])[2]  # number of trials
  
  # use row.facs to see if pars need reordering
  if(!is.null(attr(cvs, 'row.facs'))) {
    row.facs = attr(cvs, 'row.facs')
    dfColnames <- names(attr(model, 'factors'))
    if('.' %in% row.facs[1]) {
      cues = sub('s1.', '', row.facs)
      newOrder = match(cues, attr(p.list, 'facs')$cue)
    } else {
      newOrder = match(row.facs, attr(p.list, 'facs')$S)
    }
    pars <- pars[,newOrder,]
  }
  
  cvsValue <- cvs[,1:(2/3*ncol(cvs))]
  cvsReward <- cvs[,(2/3*ncol(cvs)+1):ncol(cvs)]
  # re-generate choice idx (couldn't be passed...)
  choiceIdx <- apply(cvsValue, 1, function(x) which(!is.na(x)))
  choiceIdx <- ifelse(choiceIdx%%2==0, choiceIdx-1, choiceIdx)
  choiceIdx <- rep(choiceIdx, each=2)
  choiceIdx[seq(2,length(choiceIdx),2)] = choiceIdx[seq(2,length(choiceIdx),2)]+1
  VVchoiceIdx <- cbind(rep(1:(length(choiceIdx)/2),each=2), choiceIdx)
  
  RRchoiceIdx <- apply(cvsReward, 1, function(x) which(!is.na(x)))
  RRchoiceIdx <- cbind(1:length(RRchoiceIdx), RRchoiceIdx)
  
  # update
  tmp <- transform2.dmc(pars, cvsValue, VVchoiceIdx, cvsReward, RRchoiceIdx)
  mean_v = tmp$mean_v
  b = tmp$b
  pars = tmp$pars
  updated = tmp$updated
  RR = tmp$RR
  
  # save
  if(save.adapt) {
    adapt <- data.frame()[1:n,]
    for(i in 1:dim(pars)[1]) { # accumulator N
      accumulatorName <- paste0('r', i)
      for(parName in c('SR')) {
        ii <- which(dimnames(pars)[[3]]==parName)
        adapt[,paste0(parName, '.', accumulatorName)] <- pars[i,,ii]
      }
      # add mean_v, b
      adapt[,paste0('mean_v.', accumulatorName)] <- mean_v[,i]
      adapt[,paste0('b.', accumulatorName)] <- b[,i]
    }
    adapt[,'RR'] <- RR
    adapt[,'PEs'] <- apply(updated$predictionErrors[,1:ncol(cvsValue)], 1, sum, na.rm=TRUE)
  }
  
  # generate random trials
  out <- rWaldRaceSM(n=n,
                     A=list(pars[1,,'A'], pars[2,,'A']),
                     v=list(mean_v[,1], mean_v[,2]),
                     B=list(b[,1], b[,2]),
                     t0=list(pars[1,,'t0'], pars[2,,'t0']),
                     st0=0,
                     s=list(pars[1,,'s'], pars[2,,'s']),
                     silent=TRUE, simPerTrial=FALSE)
  colnames(out) <- c('RT', 'R')
  
  if(!is.null(attr(cvs, 'row.facs'))) {
    if('.' %in% attr(cvs, 'row.facs')[1]) {
      out$cue <- factor(cues, levels=attr(model,"factors")$cue)
    } else {
      out$S <- factor(row.facs, levels=attr(model,"factors")$S)
    }
  }
  
  if(save.adapt) attr(out, 'adapt') <- adapt
  out
}


likelihood.dmc <- function(p.vector,data,min.like=1e-10, use.c=TRUE)   
{
  
  p.list <- p.list.dmc(p.vector,model=attributes(data)$model,n1order=FALSE,
                       cells=attributes(data)$cells,
                       cvs=data[,attr(attributes(data)$model,"cvs")],
                       n1.index=attr(data,"n1.index"),
                       do.trans=TRUE)  # last line by SM
  
  pars <- array(unlist(p.list,use.names=FALSE),
                dim=c(dim(p.list[[1]]),length(p.list)),
                dimnames=list(NULL,NULL,names(p.list)))
  
  cvs <- attr(data,"cvs")
  cvsValue <- cvs[,1:(2/3*ncol(cvs))]
  cvsReward <- cvs[,(2/3*ncol(cvs)+1):ncol(cvs)]
  VVchoiceIdx <- attr(data, "VVchoiceIdx")
  RRchoiceIdx <- attr(data, "RRchoiceIdx")
  
  # update & unpack
  tmp = transform2.dmc(pars, cvsValue, VVchoiceIdx, cvsReward, RRchoiceIdx)
  pars = tmp$pars
  mean_v = tmp$mean_v
  b = tmp$b
  
  # likelihood
  pmax(dWaldRace(rt=data$RT, response=data$R,
                 A=list(pars[1,,'A'], pars[2,,'A']),
                 v=list(mean_v[,1], mean_v[,2]),
                 B=list(b[,1], b[,2]),
                 t0=list(pars[1,,'t0'], pars[2,,'t0']),
                 st0=pars[1,1,'st0'],
                 s=list(pars[1,,'s'], pars[2,,'s']),
                 silent=TRUE), min.like, na.rm=TRUE)
}



# Functions required for fitting RL-DDM
prepareForFitting <- function(dat, n_bins=5) {
  # some checks first
  if(!'reward' %in% colnames(dat)) stop('Column `reward` is missing')
  if(!'stimulus_set' %in% colnames(dat)) stop('Column `stimulus_set` is missing')
  if(!'rt' %in% colnames(dat)) stop('Column `rt` is missing (are you using `RT`?)')
  if(!'choiceIsHighP' %in% colnames(dat)) stop('Column `choiceIsHighP` is missing')
  if('choice' %in% colnames(dat)) warning('Warning! Column `choice` will be overwritten')
  
  # Remove non-responses / too fast / slow responses  
  df <- dat[!is.na(dat$rt),] # remove nonresponse trials
  trialsToIgnore <- df$rt<=.15 | df$rt > 2  # tag trials that were very fast or slow.
  #  df <- df[df$rt>.1,] # remove extremely fast responses
  
  # Fix stimulus set to be numerical and have no missing numbers
  df$stimulus_set <- df$stimulus_set-min(df$stimulus_set)
  df$stimulus_set <- match(df$stimulus_set, unique(df$stimulus_set))-1
  
  # Get response "correctness" - i.e., whether a choice corresponds to the optimal choice
  # code choice in terms of upper bound (2) or lower bound (1), using the DDM convention in `rtdists`
  df$choice <- ifelse(df$choiceIsHighP, 2, 1) 
  
  # check for outcome column, for compatibility with older data. This will be removed later.
  if('outcome' %in% colnames(df)) {
    if(max(df$outcome) == 100) {
      df$reward <- df$outcome / 100
    }
  }
  #df$reward <- df$outcome / 100  # re-code rewards to 0/1, assuming an identity value function
  
  # define bins per stimulus; useful for later checking model fit
  df$trialN_this_stim <- NA
  for(lvl in unique(df$stimulus_set)) {
    df$trialN_this_stim[df$stimulus_set==lvl] <- seq(1, sum(df$stimulus_set==lvl))
  }
  df$bin <- as.numeric(cut(df$trialN_this_stim, n_bins))
  
  # Set-up outcome matrix
  outcomes <- matrix(NA, nrow=nrow(df), ncol=length(unique(df$stimulus_set))*2)
  for(row in 1:nrow(df)) {
    cond = df$stimulus_set[row]
    outcomes[row,(cond)*2+ifelse(df$choice[row]==1, 2, 1)] <- df$reward[row]
  }
  
  rewardsByStim <- matrix(NA, nrow=nrow(df), ncol=length(unique(df$stimulus_set)))
  for(row in 1:nrow(df)) {
    cond = df$stimulus_set[row]
    rewardsByStim[row,cond+1] <- df$reward[row]
  }
  
  RRchoiceIdx <- cbind(1:nrow(rewardsByStim), apply(rewardsByStim, 1, function(x) which(!is.na(x))))
  
  # Set-up values vector
  values <- rep(0.5, ncol(outcomes))
  
  # On the basis of which alternatives is chosen? 
  # Make a matrix of nTrials x 2; first column = trialN, second column = choice number
  VVchoiceIdx <- matrix(FALSE, nrow=nrow(df), ncol=ncol(outcomes))
  for(tr in 1:nrow(df)) {
    stimulus_set <- df[tr, 'stimulus_set']
    VVchoiceIdx[tr, ((stimulus_set)*2+1):((stimulus_set)*2+2)] <- TRUE
  }
  VVchoiceIdx <- which(t(VVchoiceIdx), arr.ind = TRUE)[,2:1]  # Gives for every trial each column that is chosen
  choice <- df$choice
  
  return(list(df=df, VVchoiceIdx=VVchoiceIdx, RRchoiceIdx=RRchoiceIdx, 
              outcomes=outcomes, rewardsByStim=rewardsByStim,
              values=values, trialsToIgnore=trialsToIgnore))
}

