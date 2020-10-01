simulate.full <- function(nTrialsPerSet=100, nSets=9, 
                          setProbability=c(0.25, 0.75), 
                          v0=1, wV=2, wS=0.3, t0=0.3, alpha=0.05, B=2,
                          startValues=c(0.0, 0.0)) {
  # simulates an entire experiment for a single subject
  # some 'experimental settings' can be adjusted: 
  # - the number of trials per stimulus set;
  # - the number of stimulus sets;
  # - the reward contingencies per stimulus set. Can be a list of length nSets, which allows for specifying varying reward contingencies
  
  # model parameters are:
  # - v0: "urgency" effect on drift rate
  # - wV: weight of the difference in values on drift rate
  # - wS: weight of the sum in values on drift rate
  # (such that v_i = v0 + wV*(v_i-v_j) + wS*(v_i+v_j) )
  # - t0: non-decision time
  # - alpha: learning rate
  # - B: threshold
  # finally: startValues are the values at trial=0
  
  df <- data.frame(S=1, R=NA, RT=NA, stimulus_set=NA, reward=NA, predictionError=NA, value1=NA, value2=NA)[c()]
  
  # loop over stimulus sets
  for(set in 1:nSets) {
    if(is.list(setProbability)) {
      rewardProbability <- setProbability[[set]]
    } else {
      rewardProbability <- setProbability
    }
    values <- startValues
    
    # loop over trials
    for(trialN in 1:nTrialsPerSet) {
      # determine drift rate
      v1 <- v0 + wV*(values[1]-values[2]) + wS*(values[1]+values[2])  # "incorrect"
      v2 <- v0 + wV*(values[2]-values[1]) + wS*(values[1]+values[2])  # "correct"
      
      # simulate single trial based on current values
      thisTrial <- rWaldRaceSM(n=1, A=0, B=B, t0=t0, v=c(v1, v2), s=1)
      choice <- thisTrial$R #ifelse(thisTrial$R==1, 2, 1)   # shortcut
      
      # sample reward
      reward <- rbinom(1, 1, rewardProbability[choice])
      
      # calculate prediction error
      predictionError <- reward - values[choice]
      
      # update value of chosen option with Rescorla-Wagner rule
      values[choice] <- values[choice] + alpha*predictionError
      
      # off-load to dataframe
      df <- rbind(df, cbind('S'=1, thisTrial, stimulus_set=set, reward=reward, 
                            predictionError=predictionError, value1=values[1], value2=values[2],
                            trialNthisSet=trialN))
    }
  }
  return(df)
}

# Illustration of simulation, some plots ----------------------------------
# nTrialsPerSet = 50
# df = simulate.full(nTrialsPerSet = nTrialsPerSet, nSets = 9,
#                    B=2, v0=2, wV=2)
# par(mfrow=c(2,1))
# plot(0, 0, xlim=c(0, nTrialsPerSet), ylim=c(0,1), type='n', xlab='Trial', ylab='Value')
# for(stimSet in unique(df$stimulus_set)) {
#   lines(df$trialNthisSet[df$stimulus_set==stimSet], df$value1[df$stimulus_set==stimSet], col=stimSet)
#   lines(df$trialNthisSet[df$stimulus_set==stimSet], df$value2[df$stimulus_set==stimSet], col=stimSet, lty=2)
#   points(df$trialNthisSet[df$stimulus_set==stimSet], df$R[df$stimulus_set==stimSet]-1, col=stimSet, pch=4)
# }
# 
# plot(0, 0, xlim=c(0, nTrialsPerSet), ylim=c(-1,1), type='n', xlab='Trial', ylab='Prediction error')
# for(stimSet in unique(df$stimulus_set)) {
#   lines(df$trialNthisSet[df$stimulus_set==stimSet], df$predictionError[df$stimulus_set==stimSet], col=stimSet)
# }





## let's simulate the ARW model
simulate.full.reversal <- function(nTrialsPerSet=100, nSets=9, 
                          setProbability=c(0.25, 0.75),  reversal_trialN=32,
                          v0=1, wV=2, wS=0.3, t0=0.3, alpha=0.05, B=2,
                          startValues=c(0.0, 0.0)) {
  # simulates an entire experiment for a single subject
  # some 'experimental settings' can be adjusted: 
  # - the number of trials per stimulus set;
  # - the number of stimulus sets;
  # - the reward contingencies per stimulus set. Can be a list of length nSets, which allows for specifying varying reward contingencies
  
  # model parameters are:
  # - v0: "urgency" effect on drift rate
  # - wV: weight of the difference in values on drift rate
  # - wS: weight of the sum in values on drift rate
  # (such that v_i = v0 + wV*(v_i-v_j) + wS*(v_i+v_j) )
  # - t0: non-decision time
  # - alpha: learning rate
  # - B: threshold
  # finally: startValues are the values at trial=0
  
  df <- data.frame(S=1, R=NA, RT=NA, stimulus_set=NA, reward=NA, predictionError=NA, value1=NA, value2=NA)[c()]
  
  # loop over stimulus sets
  for(set in 1:nSets) {
    if(is.list(setProbability)) {
      rewardProbability <- setProbability[[set]]
    } else {
      rewardProbability <- setProbability
    }
    values <- startValues
    
    # loop over trials
    for(trialN in 1:nTrialsPerSet) {
      if(trialN == reversal_trialN) rewardProbability <- rev(rewardProbability)  # reverse
      
      # determine drift rate
      v1 <- v0 + wV*(values[1]-values[2]) + wS*(values[1]+values[2])  # "incorrect"
      v2 <- v0 + wV*(values[2]-values[1]) + wS*(values[1]+values[2])  # "correct"
      
      # simulate single trial based on current values
      thisTrial <- rWaldRaceSM(n=1, A=0, B=B, t0=t0, v=c(v1, v2), s=1)
      choice <- thisTrial$R #ifelse(thisTrial$R==1, 2, 1)   # shortcut
      
      # sample reward
      reward <- rbinom(1, 1, rewardProbability[choice])
      
      # calculate prediction error
      predictionError <- reward - values[choice]
      
      # update value of chosen option with Rescorla-Wagner rule
      values[choice] <- values[choice] + alpha*predictionError
      
      # off-load to dataframe
      df <- rbind(df, cbind('S'=1, thisTrial, stimulus_set=set, reward=reward, 
                            predictionError=predictionError, value1=values[1], value2=values[2],
                            trialNthisSet=trialN))
    }
  }
  return(df)
}

## SAT
simulate.SAT <- function(nTrialsPerSet=100, nSets=9, 
                         setProbability=c(0.25, 0.75), 
                         driftMod.ACC=0, driftMod.SPD=0,
                         B0=1.5, B0Mod.ACC=0, B0Mod.SPD=0,
                         V0=1, V0Mod.ACC=0, V0Mod.SPD=0,
                         wV=2, wS=0.3, 
                         t0=0.3, alpha=0.05, 
                         startValues=c(0.0, 0.0)) {
  # simulates an entire experiment for a single subject
  # some 'experimental settings' can be adjusted: 
  # - the number of trials per stimulus set;
  # - the number of stimulus sets;
  # - the reward contingencies per stimulus set. Can be a list of length nSets, which allows for specifying varying reward contingencies
  
  # model parameters are:
  # - v0: "urgency" effect on drift rate
  # - wV: weight of the difference in values on drift rate
  # - wS: weight of the sum in values on drift rate
  # (such that v_i = v0 + wV*(v_i-v_j) + wS*(v_i+v_j) )
  # - t0: non-decision time
  # - alpha: learning rate
  # - B: threshold
  # finally: startValues are the values at trial=0
  
  df <- data.frame(S=1, R=NA, RT=NA, stimulus_set=NA, reward=NA, predictionError=NA, value1=NA, value2=NA)[c()]
  
  # loop over stimulus sets
  for(set in 1:nSets) {
    if(is.list(setProbability)) {
      rewardProbability <- setProbability[[set]]
    } else {
      rewardProbability <- setProbability
    }
    values <- startValues
    
    # loop over trials
    for(trialN in 1:nTrialsPerSet) {
      trialType <- ifelse(runif(1)<.5, 'SPD', 'ACC')
      if(trialType == 'SPD') {
        B0Mod <- B0Mod.SPD
        V0Mod <- V0Mod.SPD
        driftMod <- driftMod.SPD
      } else {
        B0Mod <- B0Mod.ACC
        V0Mod <- V0Mod.ACC
        driftMod <- driftMod.ACC
      }
#      B <- B0 * (1+B0Mod)  #ifelse(trialType == 'SPD', B.SPD, B.ACC)
#      V0 <- V0 * (1+V0Mod)
      
      # determine drift rate
      v1 <- (1+V0Mod)*V0 + wV*(values[1]-values[2]) + wS*(values[1]+values[2])  # "incorrect"
      v2 <- (1+V0Mod)*V0 + wV*(values[2]-values[1]) + wS*(values[1]+values[2])  # "correct"
      
      # simulate single trial based on current values
      thisTrial <- rWaldRaceSM(n=1, A=0, B=B0*(1+B0Mod), t0=t0, v=c(v1, v2), s=1)
      choice <- thisTrial$R #ifelse(thisTrial$R==1, 2, 1)   # shortcut
      
      # sample reward
      reward <- rbinom(1, 1, rewardProbability[choice])
      
      # calculate prediction error
      predictionError <- reward - values[choice]
      
      # update value of chosen option with Rescorla-Wagner rule
      values[choice] <- values[choice] + alpha*predictionError
      
      # off-load to dataframe
      df <- rbind(df, cbind('S'=1, thisTrial, stimulus_set=set, reward=reward, cue=trialType,
                            predictionError=predictionError, value1=values[1], value2=values[2],
                            trialNthisSet=trialN))
    }
  }
  return(df)
}
