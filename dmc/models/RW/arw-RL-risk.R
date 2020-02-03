library(dmcAdapt)

# Template setup for n-choice LBA, B=b-A parameterization
#   External parameters types: A, B, t0, mean_v, sd_v, st0 = 0 (optional)
#   Internal parameters types: A, b, t0, mean_v, sd_v, st0 = 0 (optional)

# This ncl (no cell likelihood) version uses p.list.dmc in the likelihood to 
# avoid looping over design cells.

transform.dmc <- function(par.df,do.trans=TRUE) 
{
  par.df$aR[par.df$aR==Inf] <- par.df$aV[par.df$aR==Inf]
  par.df$aV = t(pnorm(t(par.df$aV)))
  par.df$aR = t(pnorm(t(par.df$aR)))
  
  par.df$SR <- t(pnorm(t(par.df$SR)))
  par.df$RR <- t(pnorm(t(par.df$RR))) / 2
  
  # if start point is Inf, assume risk is "known" / "correct"
  par.df$RR[par.df$RR == Inf] <- par.df$SR[par.df$RR == Inf]*(1-par.df$SR[par.df$RR == Inf])
  par.df$RR[par.df$RR == 0] <- 1e-10 # never 0, since sqrt is taken in update.
  
  if (do.trans) {
    return(
      list(A=t(par.df$A),
           s=t(par.df$s),
           t0=t(par.df$t0),
           st0=t(par.df$st0),
           B0=t(par.df$B0),
           SR=t(par.df$SR),
           aV=t(par.df$aV),
           aR=t(par.df$aR),
           RR=t(par.df$RR),
           V0=t(par.df$V0),
           wV=t(par.df$wV)))
  } else {
    return(
      list(A=par.df$A,
           s=par.df$s,
           t0=par.df$t0,
           st0=par.df$st0,
           B0=par.df$B0,
           SR=par.df$SR,
           aV=par.df$aV,
           aR=par.df$aR,
           RR=par.df$RR,
           V0=par.df$V0,
           wV=par.df$wV))
  }
}


transform2.dmc <- function(pars, cvs, choiceIdx) {
  ### SM
  
  # Create start point vector
  startValues <- rep(pars[1,,'SR'], each=2, times=ncol(cvs))
  riskStartValues <- rep(pars[1,,'RR'], each=2, times=ncol(cvs))
  
  # learning rates matrix
  learningRates <- matrix(rep(pars[1,,'aV'], each=ncol(cvs)), ncol=ncol(cvs), byrow=TRUE)
  riskLearningRates <- matrix(rep(pars[1,,'aR'], each=ncol(cvs)), ncol=ncol(cvs), byrow=TRUE)
  
  # call C
  updated <- adapt.c.dmc(startValues = startValues, 
                         learningRates = learningRates, 
                         feedback = cvs,
                         learningRule='SARSARisk',
                         riskLearningRates = riskLearningRates,
                         riskStartValues = riskStartValues)
  
  # add back data to 'pars' array
  pars[,,'SR'] <- matrix(updated$adaptedValues[choiceIdx], ncol=2, byrow=FALSE)
  pars[,,'RR'] <- matrix(updated$adaptedRiskValues[choiceIdx], ncol=2, byrow=FALSE)
  
  # Determine parameter mean_v
  mean_v <- cbind(pars[1,,"V0"] + pars[1,,"wV"]*(pars[2,,"SR"]-pars[1,,'SR']),
                  pars[2,,"V0"] + pars[2,,"wV"]*(pars[1,,"SR"]-pars[2,,'SR']))
  
  # s <- cbind(pars[1,,"s"] * (1 + pars[1,,"wS"]*pars[2,,"RR"]),
  #            pars[2,,"s"] * (1 + pars[2,,"wS"]*pars[1,,"RR"]))
  
  # if(any(is.na(s))) {warning('NA found!'); print(pars[1,1,])}
  # if(any(s<=0)) warning('s < 0 found!')
  # Determine parameter b
  b <- cbind(pars[1,,"B0"],
             pars[2,,"B0"])
  
  return(list(pars=pars, mean_v=mean_v, b=b, updated=updated))
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
    newOrder = match(row.facs, apply(attr(p.list, 'facs'), 1, paste, collapse='.'))
    pars <- pars[,newOrder,]
  }

  # re-generate choice idx (couldn't be passed...)
  choiceIdx <- apply(cvs, 1, function(x) which(!is.na(x)))
  choiceIdx <- ifelse(choiceIdx%%2==0, choiceIdx-1, choiceIdx)
  choiceIdx <- rep(choiceIdx, each=2)
  choiceIdx[seq(2,length(choiceIdx),2)] = choiceIdx[seq(2,length(choiceIdx),2)]+1
  choiceIdx <- cbind(rep(1:(length(choiceIdx)/2),each=2), choiceIdx)
  
  # update
  tmp <- transform2.dmc(pars, cvs, choiceIdx)
  mean_v = tmp$mean_v
  b = tmp$b
  pars = tmp$pars
  updated = tmp$updated
  
  # save
  if(save.adapt) {
    adapt <- data.frame()[1:n,]
    for(i in 1:dim(pars)[1]) { # accumulator N
      accumulatorName <- paste0('r', i)
      for(parName in c('SR', 'RR')) {
        ii <- which(dimnames(pars)[[3]]==parName)
        adapt[,paste0(parName, '.', accumulatorName)] <- pars[i,,ii]
      }
      # add mean_v, b
      adapt[,paste0('mean_v.', accumulatorName)] <- mean_v[,i]
      adapt[,paste0('b.', accumulatorName)] <- b[,i]
    }
    adapt[,'PEs'] <- apply(updated$predictionErrors[,1:(ncol(cvs)-1)], 1, sum, na.rm=TRUE)
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
    row.facs = do.call(rbind, strsplit(attr(cvs, 'row.facs'), split='.', fixed=TRUE))
    facNames <- colnames(attr(p.list, 'facs'))
    for(i in 1:length(facNames)) {
      out[,facNames[i]] <- row.facs[,i]
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
  choiceIdx <- attr(data, "VVchoiceIdx")
  
  # update & unpack
  tmp = transform2.dmc(pars, cvs, choiceIdx)
  pars = tmp$pars
  mean_v = tmp$mean_v
  b = tmp$b
  
  if(any(is.na(mean_v))) return(rep(min.like, nrow(data)))
  
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
  return(list(df=df, VVchoiceIdx=VVchoiceIdx, outcomes=outcomes, values=values, trialsToIgnore=trialsToIgnore))
}

