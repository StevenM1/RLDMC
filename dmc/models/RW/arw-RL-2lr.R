library(dmcAdapt)

# Template setup for n-choice LBA, B=b-A parameterization
#   External parameters types: A, B, t0, mean_v, sd_v, st0 = 0 (optional)
#   Internal parameters types: A, b, t0, mean_v, sd_v, st0 = 0 (optional)

# This ncl (no cell likelihood) version uses p.list.dmc in the likelihood to 
# avoid looping over design cells.

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


# save.adapt=TRUE
random.dmc <- function(p.list,model,save.adapt=TRUE)
{
  pars <- array(unlist(p.list,use.names=FALSE),
                dim=c(dim(p.list[[1]]),length(p.list)),
                dimnames=list(NULL,NULL,names(p.list)))
  cvs <- attr(p.list,"cvs")
  
  # use row.facs to see if pars need reordering
  if(!is.null(attr(cvs, 'row.facs'))) {
    row.facs = attr(cvs, 'row.facs')
    if('.' %in% row.facs[1]) {
      cues = sub('s1.', '', row.facs)
      newOrder = match(cues, attr(p.list, 'facs')$cue)
      pars <- pars[,newOrder,]
    }
  }
  
  n <- dim(p.list[[1]])[2]  # number of trials
  #  facs <- attr(p.list,"facs") # factor columns + CR = numeric correct response
  
  choiceIdx <- apply(cvs, 1, function(x) which(!is.na(x)))
  choiceIdx <- ifelse(choiceIdx%%2==0, choiceIdx-1, choiceIdx)
  choiceIdx <- rep(choiceIdx, each=2)
  choiceIdx[seq(2,length(choiceIdx),2)] = choiceIdx[seq(2,length(choiceIdx),2)]+1
  choiceIdx <- cbind(rep(1:(length(choiceIdx)/2),each=2), choiceIdx)
  
  startValues <- rep(pars[1,1,'SR'], ncol(cvs))
  #  startValues <- attr(data, "startingValues")
  learningRatesPos <- matrix(pars[1,,'aVP'], nrow=nrow(cvs), ncol=ncol(cvs))
  learningRatesNeg <- matrix(pars[1,,'aVN'], nrow=nrow(cvs), ncol=ncol(cvs))
  
  # call C
  updated <- adapt.c.dmc(startValues = startValues, 
                         learningRates = learningRatesPos, 
                         feedback = cvs, learningRatesNeg=learningRatesNeg)
  
  # add back data to 'pars' array
  pars[,,'SR'] <- matrix(updated$adaptedValues[choiceIdx], ncol=2, byrow=FALSE)
  mean_v <- cbind(pars[1,,"V0"] + pars[1,,"wV"]*(pars[2,,"SR"]-pars[1,,"SR"]),
                  pars[2,,"V0"] + pars[2,,"wV"]*(pars[1,,"SR"]-pars[2,,"SR"]))
  
  # "perceived difficulty"
  b <- cbind(pars[1,,"B0"] + pars[1,,"wB"]*(pars[2,,"SR"]-pars[1,,"SR"])*pars[1,,"B0"],
             pars[2,,"B0"] + pars[2,,"wB"]*(pars[2,,"SR"]-pars[1,,"SR"])*pars[2,,"B0"])
  
  if(save.adapt) {
    adapt <- data.frame()[1:n,]
    for(i in 1:dim(pars)[1]) { # accumulator N
      accumulatorName <- paste0('r', i)
      for(ii in 1:dim(pars)[3]) { # parameter N}
        parName <- dimnames(pars)[3][[1]][ii]
        adapt[,paste0(parName, '.', accumulatorName)] <- pars[i,,ii]
      }
      # add mean_v, b
      adapt[,paste0('mean_v.', accumulatorName)] <- mean_v[,i]
      adapt[,paste0('b.', accumulatorName)] <- mean_v[,i]
    }
  }
  
  
  #   for(i in 1:dim(pars)[2]) {
  #     # for(ii in 1:dim(pars)[3]) {
  #     #   adapt[[i]][[dimnames(pars)[[3]][ii]]] <- c('r1'=as.numeric(pars[1,i,ii]), 'r2'=as.numeric(pars[2,i,ii]))
  #     # }
  #     adapt[[i]][[dimnames(pars)[[3]]]]
  #     adapt[[i]][['mean_v']] = c('r1'=mean_v[i,1], 'r2'=mean_v[i,2])
  #     adapt[[i]][['b']] = c('r1'=b[i,1], 'r2'=b[i,2])
  #   }
  # }
  # S <- as.numeric(facs$S)-1 # use stimulus indicator as binary real stim value
  #  S <- cvs$stim # real valued stimulus
  # Get correct response as feedback and possibly corrupt with prbability 1-pcFB
  # FB <- as.numeric(facs$S)-1
  # FBresponse <- p.list$pcFB[1,] < 0
  # if ( any(!FBresponse) ) {
  #   flip <- rbinom(n,1,p.list$pcFB[1,!FBresponse])==0
  #   FB[!FBresponse][flip] <- as.numeric(!FB[!FBresponse][flip])
  # }
  # # Update SP and SR
  # p <- get.p(p.list,1,exclude=NULL)
  # p <- adapt.dmc(p,Si1=S[1]) # for first iteration
  # out <- matrix(nrow= n, ncol=2,dimnames=list(NULL,c("RT","R")))
  # if (save.adapt) {
  #   adapt <- vector(mode="list",length=n)
  # adapt[[1]] <- p[c("SP","SR","b","mean_v","sd_v","t0","A")] 
  # }
  
  
  out <- rWaldRaceSM(n=n,
                     A=list(pars[1,,'A'], pars[2,,'A']),
                     v=list(mean_v[,1], mean_v[,2]),
                     B=list(b[,1], b[,2]),
                     t0=list(pars[1,,'t0'], pars[2,,'t0']),
                     st0=0,  #pars[1,1,'st0'],#, pars[2,,'st0']),
                     #s=list(pars[1,,'s'], pars[2,,'s']),
                     silent=TRUE, simPerTrial=FALSE)
  colnames(out) <- c('RT', 'R')
  if(!is.null(attr(cvs, 'row.facs'))) {
    if('.' %in% attr(cvs, 'row.facs')[1]) {
      out$cue <- factor(cues, levels=attr(model,"factors")$cue)
    }    
  }
  #  if (save.adapt) attr(out,"adapt") <- do.call(rbind,lapply(adapt,unlist))
  if(save.adapt) attr(out, 'adapt') <- adapt
  out
}


transform.dmc <- function(par.df,do.trans=TRUE) 
{
  
  if (do.trans)
    list(A=t(par.df$A),s=t(par.df$s),t0=t(par.df$t0),st0=t(par.df$st0),
         B0=t(par.df$B0),
         SR=t(par.df$SR),
         wB=t(par.df$wB),
         aVP=pnorm(t(par.df$aVP)),
         aVN=pnorm(t(par.df$aVN)),
         V0=t(par.df$V0),wV=t(par.df$wV)) else
           list(A=par.df$A,s=par.df$s,t0=par.df$t0,st0=par.df$st0,
                B0=par.df$B0,
                SR=par.df$SR,
                wB=par.df$wB,
                aVP=t(pnorm(t(par.df$aVP))),
                aVN=t(pnorm(t(par.df$aVN))),
                V0=par.df$V0,wV=par.df$wV)       
}


likelihood.dmc <- function(p.vector,data,min.like=1e-10, use.c=TRUE)   
{
  
  do.n1 <- function(x) matrix(x[attr(data,"n1.index")],ncol=dim(x)[2])
  
  p.list <- p.list.dmc(p.vector,model=attributes(data)$model,n1order=FALSE,
                       cells=attributes(data)$cells,
                       cvs=data[,attr(attributes(data)$model,"cvs")],
                       n1.index=attr(data,"n1.index")
  )
  
  pars <- array(unlist(p.list,use.names=FALSE),
                dim=c(dim(p.list[[1]]),length(p.list)),
                dimnames=list(NULL,NULL,names(p.list)))
  
  cvs <- attr(data,"cvs")
  
  if(use.c) {
    ### SM
    # all that is required for adapting is:
    # 1. Trial-by-trial feedback
    # 2. Starting values
    # 3. Learning rate (trial-by-trial)
    
    # 1. Create trial-by-trial feedback matrix
    # this matrix is of size [nTrials x nAdapt], where nAdapt is the number of representations / probabilities that are adapted.
    # e.g., in this example, nAdapt is 3: two stimulus presentations are learnt, and 1 stimulus probability.
    # Example:
    # > head(feedback)
    #           [,1]     [,2] [,3]
    # [1,] 1.060133       NA    0
    # [2,] 1.248056       NA    0
    # [3,]       NA 3.997854    1
    # [4,]       NA 3.428983    1
    # [5,]       NA 5.182715    1
    # [6,] 2.336167       NA    0
    
    # NAs depict trials on which the given representation doesn't need to be updated (i.e., no feedback was given)
    # nacc = dim(pars)[2]
    # feedback <- matrix(NA, nrow=nrow(cvs), ncol=nacc+1)
    # for(accumulator in 1:nacc) {
    #   feedback[cvs$FB==accumulator, accumulator] = cvs$stim[cvs$FB==accumulator]
    # }
    # feedback[,nacc+1] = cvs$FB-1 # since probabilities must sum to 1, only one of the SPs need to be updated
    # 
    # Create start point vector
    #    startValues <- c(pars[1,,c('SR')], pars[1,1,c('SP')])
    startValues <- attr(data, "startingValues")
    #  startValues <- attr(data, "startingValues")
    learningRatesPos <- matrix(pars[1,,'aVP'], nrow=nrow(cvs), ncol=ncol(cvs))
    learningRatesNeg <- matrix(pars[1,,'aVN'], nrow=nrow(cvs), ncol=ncol(cvs))
    
    # call C
    updated <- adapt.c.dmc(startValues = startValues, 
                           learningRates = learningRatesPos, 
                           feedback = cvs, learningRatesNeg=learningRatesNeg)
    
    
    # add back data to 'pars' array
    pars[,,'SR'] <- matrix(updated$adaptedValues[attr(data, 'VVchoiceIdx')],
                           ncol=2, byrow=TRUE)
  }
  
  # Update parameters
  mean_v <- cbind(pars[,1,"V0"] + pars[,1,"wV"]*(pars[,2,"SR"]-pars[,1,"SR"]),
                  pars[,2,"V0"] + pars[,2,"wV"]*(pars[,1,"SR"]-pars[,2,"SR"]))
  # "perceived difficulty" effect
  b <- cbind(pars[,1,"B0"] + pars[,1,"wB"]*(pars[,2,"SR"]-pars[,1,"SR"])*pars[,1,"B0"],
             pars[,2,"B0"] + pars[,2,"wB"]*(pars[,2,"SR"]-pars[,1,"SR"])*pars[,2,"B0"])
  
  pmax(dWaldRace(rt=data$RT, response=data$R,
                 A=list(pars[,1,'A'], pars[,2,'A']),
                 v=list(mean_v[,1], mean_v[,2]),
                 B=list(b[,1], b[,2]), # b[,1], b[,2]),
                 t0=list(pars[,1,'t0'], pars[,2,'t0']),
                 st0=pars[1,1,'st0'],#, pars[2,,'st0']),
                 s=list(pars[,1,'s'], pars[,2,'s']),
                 silent=TRUE), min.like)
}



