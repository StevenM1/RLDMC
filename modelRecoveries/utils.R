simulate.rlard <- function(nTrialsPerSet=100, nSets=9, 
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

library(rtdists)

simulate.rlddm <- function(nTrialsPerSet=100, nSets=9, 
                           setProbability=c(0.25, 0.75), 
                           m=2, t0=0.3, alpha=0.05, a=2,
                           startValues=c(0.0, 0.0)) {
  # simulates an entire experiment for a single subject
  # some 'experimental settings' can be adjusted: 
  # - the number of trials per stimulus set;
  # - the number of stimulus sets;
  # - the reward contingencies per stimulus set. Can be a list of length nSets, which allows for specifying varying reward contingencies
  
  # model parameters are:
  # - wV: weight of the difference in values on drift rate
  # (such that v_i =  wV*(v_i-v_j) )
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
      v <- m*(values[2]-values[1]) # "correct"
      
      # simulate single trial based on current values
      thisTrial <- rdiffusion(n=1, v=v, a=a, t0=t0, z=0.5*a, s=1)
      thisTrial$R <- as.numeric(thisTrial$response)
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

# 
# nTrialsPerSet = 50
# df = simulate.rlddm(nTrialsPerSet = nTrialsPerSet, nSets = 9,
#                    a=2, m=3)
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

