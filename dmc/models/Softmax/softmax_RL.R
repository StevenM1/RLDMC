library(dmcAdapt)

# Template setup for n-choice LBA, B=b-A parameterization
#   External parameters types: A, B, t0, mean_v, sd_v, st0 = 0 (optional)
#   Internal parameters types: A, b, t0, mean_v, sd_v, st0 = 0 (optional)

# This ncl (no cell likelihood) version uses p.list.dmc in the likelihood to 
# avoid looping over design cells.

transform.dmc <- function(par.df,do.trans=TRUE) 
{
  par.df$aV = t(pnorm(t(par.df$aV)))
  par.df$SR = t(pnorm(t(par.df$SR)))
  
  if (do.trans) {
    return(
      list(beta=t(par.df$beta),
           aV=t(par.df$aV),
           SR=t(par.df$SR)))
  } else {
    return(
      list(beta=par.df$beta,
           aV=par.df$aV,
           SR=par.df$SR))
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
  
  ## probabilities per trial
  PP <- t(exp(pars[,,'SR']*pars[,,'beta']))
  PP <- PP/apply(PP, 1, sum)
  
  return(list(pars=pars, updated=updated, PP=PP))
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
  PP = tmp$PP
  
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
    adapt[,'PEs'] <- apply(updated$predictionErrors[,1:(ncol(cvs)-1)], 1, sum, na.rm=TRUE)
  }
  
 
  
  # generate random trials
  out <- data.frame(RT=rnorm(n, 1000, 1), R=NA)  # generate some random RTs to trick dmc
  out$R <- rbinom(n, 1, prob=PP[,1])+1
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
  choiceIdx <- attr(data, "VVchoiceIdx")
  
  # update & unpack
  tmp = transform2.dmc(pars, cvs, choiceIdx)
  pars = tmp$pars

  # PP <- exp(pars[,,'SR']*pars[,,'beta'])
  # PP <- PP/apply(PP, 2, sum)
  PP <- tmp$PP
  ## note that choice == 2 is "correct", which corresponds to the *first* column in the PP matrix
  ## Sorry for that confusion
  # LL <- sum(log(PP[choice==1,2])) + sum(log(PP[choice==2,1]))	
  
  LL <- rep(min.like, nrow(data))
  LL[as.numeric(data$R)==1] <- PP[as.numeric(data$R)==1, 2]
  LL[as.numeric(data$R)==2] <- PP[as.numeric(data$R)==2, 1]
  LL
  # likelihood
  # pmax(dWaldRace(rt=data$RT, response=data$R,
  #                A=list(pars[1,,'A'], pars[2,,'A']),
  #                v=list(mean_v[,1], mean_v[,2]),
  #                B=list(b[,1], b[,2]),
  #                t0=list(pars[1,,'t0'], pars[2,,'t0']),
  #                st0=pars[1,1,'st0'],
  #                s=list(pars[1,,'s'], pars[2,,'s']),
  #                silent=TRUE), min.like, na.rm=TRUE)
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

