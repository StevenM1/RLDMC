library(dmcAdapt)

# Template setup for n-choice LBA, B=b-A parameterization
#   External parameters types: A, B, t0, mean_v, sd_v, st0 = 0 (optional)
#   Internal parameters types: A, b, t0, mean_v, sd_v, st0 = 0 (optional)

# This ncl (no cell likelihood) version uses p.list.dmc in the likelihood to 
# avoid looping over design cells.

transform.dmc <- function(par.df,do.trans=TRUE) 
{
  par.df$aV = t(pnorm(t(par.df$aV)))
#  par.df$SR = t(pnorm(t(par.df$SR)))
  par.df$d <- par.df$d*par.df$t0  # proportional to t0, bounded -1 to 1
  
  
  
  # par.df[,c("a","v","t0","z","d","sz","sv","st0")]
  
  if (do.trans) {
    return(
      list(a=t(par.df$a),
           m=t(par.df$m),
           t0=t(par.df$t0),
           z=t(par.df$z),
           d=t(par.df$d),
           sz=t(par.df$sz),
           sv=t(par.df$sv),
           st0=t(par.df$st0),
           aV=t(par.df$aV),
           SR=t(par.df$SR),
           vmax=t(par.df$vmax)))
  } else {
    return(
      list(a=par.df$a,
           m=par.df$m,
           t0=par.df$t0,
           z=par.df$z,
           d=par.df$d,
           sz=par.df$sz,
           sv=par.df$sv,
           st0=par.df$st0,
           aV=par.df$aV,
           SR=par.df$SR,
           vmax=par.df$vmax))
  }
}

transform2.dmc <- function(pars, cvs, choiceIdx) {
  ### SM
  
  # Create start point vector
  startValues <- rep(pars[1,,'SR'], each=2, times=ncol(cvs))
  
  # learning rates matrix
  learningRates <- matrix(rep(pars[1,,'aV'], each=ncol(cvs)), ncol=ncol(cvs), byrow=TRUE)
  
  # call C
  updated <- adapt.c.dmc(startValues = startValues, 
                         learningRates = learningRates, 
                         feedback = cvs,
                         learningRule='SARSA')
  
  # add back data to 'pars' array
  pars[,,'SR'] <- matrix(updated$adaptedValues[choiceIdx], ncol=2, byrow=FALSE)
  
  # ## probabilities per trial
  # PP <- t(exp(pars[,,'SR']*pars[,,'beta']))
  # PP <- PP/apply(PP, 1, sum)
  value_difference = (pars[1,,'SR']-pars[2,,'SR'])
  v = (2*pars[1,,'vmax'])/(1+exp(-pars[1,,'m']*value_difference)) - pars[1,,'vmax']
  # v = pars[1,,'m'] * (pars[1,,'SR']-pars[2,,'SR'])
  
  return(list(pars=pars, v=v, updated=updated))
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
  
  # # re-generate choice idx (couldn't be passed...)
  # choiceIdx <- apply(cvs, 1, function(x) which(!is.na(x)))
  # choiceIdx <- ifelse(choiceIdx%%2==0, choiceIdx-1, choiceIdx)
  # choiceIdx <- rep(choiceIdx, each=2)
  # choiceIdx[seq(2,length(choiceIdx),2)] = choiceIdx[seq(2,length(choiceIdx),2)]+1
  # choiceIdx <- cbind(rep(1:(length(choiceIdx)/2),each=2), choiceIdx)
  # re-generate choice idx (couldn't be passed...)
  if(sum(apply(cvs, 1, function(x) sum(!is.na(x))) > 1)) {
    # multi-feedback outcomes make it easy to define VVchoiceIdx
    choiceIdx <- !is.na(cvs)
    choiceIdx <- which(t(choiceIdx), arr.ind = TRUE)[,2:1]  # Gives for every trial each column that is chosen
  } else {
    choiceIdx <- apply(cvs, 1, function(x) which(!is.na(x)))
    choiceIdx <- ifelse(choiceIdx%%2==0, choiceIdx-1, choiceIdx)
    choiceIdx <- rep(choiceIdx, each=2)
    choiceIdx[seq(2,length(choiceIdx),2)] = choiceIdx[seq(2,length(choiceIdx),2)]+1
    choiceIdx <- cbind(rep(1:(length(choiceIdx)/2),each=2), choiceIdx)
  }  
  # update
  tmp <- transform2.dmc(pars, cvs, choiceIdx)
  pars = tmp$pars
  updated = tmp$updated
  v = tmp$v
  
  # save
  if(save.adapt) {
    adapt <- data.frame()[1:n,]
    for(i in 1:dim(pars)[1]) { # accumulator N
      accumulatorName <- paste0('r', i)
      for(parName in c('SR')) {
        ii <- which(dimnames(pars)[[3]]==parName)
        adapt[,paste0(parName, '.', accumulatorName)] <- pars[i,,ii]
      }
      # # add mean_v, b
      # adapt[,paste0('mean_v.', accumulatorName)] <- mean_v[,i]
      # adapt[,paste0('b.', accumulatorName)] <- b[,i]
    }
    adapt[,'v'] <- v
    adapt[,'PEs'] <- apply(updated$predictionErrors[,1:(ncol(cvs)-1)], 1, sum, na.rm=TRUE)
  }
  
  
  # generate random trials
  # out <- data.frame(RT=rnorm(n, 1000, 1), R=NA)  # generate some random RTs to trick dmc
  # out$R <- rbinom(n, 1, prob=PP[,1])+1
  out <- rdiffusion(n=n,
                    a=pars[1,,'a'],
                    v=v,
                    t0=pars[1,,'t0'],
                    z=pars[1,,'z']*pars[1,,'a'],
                    sz=pars[1,,'sz']*pars[1,,'a'],
                    sv=pars[1,,'sv'],
                    d=pars[1,,'d'],
                    st0=pars[1,,'st0'], s=1, precision=2.5)
  colnames(out) <- c('RT', 'R')
  out$R <- as.numeric(out$R) #factor(ifelse(out$R=='upper', 'r2', 'r1'), levels=c('r1', 'r2'))
  # out$R[out$R=='upper'] <- 'r2'
  # out$R[out$R=='upper'] <- 'r1'
  
  
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
  choiceIdx <- attr(data, "VVchoiceIdx")
  
  # update & unpack
  tmp = transform2.dmc(pars, cvs, choiceIdx)
  pars = tmp$pars
  v = tmp$v
  
  bad <- function(p) 
    # Stops ddiffusion crashing if given bad values.
  {
    (p[,'a']<0)      | (p[,'z'] <1e-6) | (p[,'z'] >.999999) | (p[,'t0']<1e-6)  | 
      (p[,'sz']<0) | (p[,'st0']<0)    | (p[,'sv']<0) |
      (p[,'sz']>1.999999*min(c(p[,'z'],1-p[,'z'])))
  }
  
  #  likes <- rep(min.like, nrow(data))
  if(any(bad(pars[1,,]))) return(rep(min.like, nrow(data)))
  
  # likelihood
  pmax(ddiffusion(rt=data$RT, response=as.numeric(data$R),
                  a=pars[1,,'a'],
                  v=v,
                  t0=pars[1,,'t0'],
                  z=pars[1,,'z']*pars[1,,'a'],
                  sz=pars[1,,'sz']*pars[1,,'a'],
                  sv=pars[1,,'sv'],
                  d=pars[1,,'d'],
                  st0=pars[1,,'st0'], s=1, precision=2.5), min.like, na.rm=TRUE)
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

