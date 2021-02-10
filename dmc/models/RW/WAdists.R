## Functions from Van Ravenzwaaij et al. (2019) https://osf.io/2s6ax/files/

# produces contrast inputs
ToM = function (input)
{
  vMat = outer (input, input, "-")
  v = c(vMat[upper.tri(vMat)], vMat[lower.tri(vMat)])
  vNames = outer (names (input), names (input), "paste", sep = "-")
  names (v) = c(vNames[upper.tri(vNames)], vNames[lower.tri(vNames)])
  v
}

# produces sum inputs
ToP = function (input)
{
  vMat = outer (input, input, "+")
  v = c(vMat[upper.tri(vMat)], vMat[lower.tri(vMat)])
  vNames = outer (names (input), names (input), "paste", sep = "-")
  names (v) = c(vNames[upper.tri(vNames)], vNames[lower.tri(vNames)])
  v
}

# produces ratio inputs x/(x+y)
ToR = function (input, a0)
{
  vMat = outer (input, input, function (x, y) {x/(a0+x+y)})
  v = c(vMat[upper.tri(vMat)], vMat[lower.tri(vMat)])
  vNames = outer (names (input), names (input), "paste", sep = "-")
  names (v) = c(vNames[upper.tri(vNames)], vNames[lower.tri(vNames)])
  v
}


# determines the winning accumulator
get.winner = function (times, nms)
{
  times = sort (times)
  losers = substr (names (times), 3, 3)
  avail = rep (TRUE, length (nms))
  names (avail) = nms
  for (i in 1:length (times))
  {
    avail[losers[i]] = FALSE  
    if (sum (avail)==1)
    {
      return (list (winner = which (avail), rt = times[i]))
    }
  }
  stop ("No winner??")
}

# determines the second last set of accumulators to finish
pick.second.biggest = function (x)
{
  -sort (-x, partial = 1:2)[2]
}

# determines drift rates for winning accumulators
WinX = function (Xname, v)
{
  unlist (lapply (strsplit (names (v), "-"), function (x){x[[1]]})) == Xname
}

# determines drift rates for losing accumulators
LoseX = function (Xname, v)
{
  unlist (lapply (strsplit (names (v), "-"), function (x){x[[2]]})) == Xname
}


############################### PDF FOR WA

WApdf = function (x, Inputs, Resp, t0, A, B, v0, ws, wd, sl=0)
  ## CODE FROM Van Ravenzwaaij et al (2019) https://osf.io/2s6ax/files/
  # x = response time
  # Inputs = Q-value vector
  # R = response
{
  n = length (Inputs)
  Data = vector (mode = "list", length = n)
  Input = names (Inputs) = names (Data) = 1:n
  v = v0 + ws*ToP (Inputs) + wd*ToM (Inputs)
  s = 1 + sl * (n-3) * v
  
  # if (Fun == "LNR")
  # {
  #   pl = list (ter = rep (Ter, length (v)), v = -1 * v, s = rep (s, length (v)))
  #   PDF = LNRpdf
  #   CDF = LNRcdf
  # }
  # if (Fun == "LBA")
  # {
  #   pl = list (ter = rep (Ter, length (v)), v = v, s = rep (s, length (v)), A = rep (A, length (v)), b = rep (b, length (v)))
  #   PDF = LBApdf
  #   CDF = LBAcdf
  # }
  # if (Fun == "tLBA")
  # {
  #   pl = list (ter = rep (Ter, length (v)), v = v, s = rep (s, length (v)), A = rep (A, length (v)), b = rep (b, length (v)))
  #   PDF = tLBApdf
  #   CDF = tLBAcdf
  # }
  pl = list(t0=rep(t0, length(v)), v=v, B=rep(B, length(v)), s=rep(s, length(v)))
  
  for (i in Input)
  {
    Winner = WinX (i, pl$v) 
    aWinner = names (v)[Winner]			# all accs where Resp wins
    iWinner = c(1:length (Winner))[Winner]	# ind of all accs where Resp wins
    Data[[i]] = list (aWinner = aWinner, iWinner = iWinner)
  }
  
  DataR = Data[[Resp]]
  
  # store pdf of each accumulator in response i set
  pdf = matrix (NA, length (x), (n-1))
  for (j in 1:(n-1))
  {
    pdf[,j] = dWald(x-t0, v=v[DataR$iWinner[j]], B=B, A=A)  #PDF (x, pl = pl, i = DataR$iWinner[j])
  }
  colnames (pdf) = DataR$aWinner
  
  # store cdf for all accumulators
  cdf = matrix (NA, length (x), length (v))
  for (j in 1:length (v))
  { 
    cdf[,j] = pWald(x-t0, v=v[j], B=B, A=A) #CDF (x, pl = pl, i = j)
  }
  colnames (cdf) = names (v)
  
  # store Pr (at least one member of other response sets IS NOT done) = (1-Pr(member 1 done)*...*1-Pr(member r done))
  NotInput = Input[Input!=Resp]
  NotDone = matrix (NA, length (x), n-1)
  for (r in 1:(n-1))
  {
    NotDone[,r] = 1 - apply (matrix (cdf[,Data[[NotInput[r]]]$aWinner], length (x), n-1), 1, prod)
  }   
  colnames (NotDone) = NotInput
  
  # probability
  Prob = rep (0, length (x))
  for (j in DataR$aWinner)
  {
    NotLast = DataR$aWinner[DataR$aWinner!=j]
    Prob = Prob + 
      pdf[,j] *		# accumulator j finishes at time t as last of set i, so that response i wins
      apply (matrix ((length(NotLast)>0) * cdf[,NotLast], length (x), length (NotLast)), 1, prod) *		# all other winning accumulators done at t 
      apply (matrix (NotDone, length (x), n-1), 1, prod)		# at least one accumulator out of other response sets is not done at t
  }
  as.vector (Prob)
}


WAlikevec <- function(RT, R, Inputs, t0, A, B, v0, ws, wd, s=1)
  ## CODE based on Van Ravenzwaaij et al (2019) https://osf.io/2s6ax/files/
  # but vectorised (a bit, not exactly perfect - it's roughly 10x as fast as DvR's implementation, but *CAN ONLY DO 3 RESPONSES*)
  # RT = vector of RTs
  # Inputs = Q-value matrix (nTrials x nOptions)
  # R = responses
{
  nOptions <- ncol(Inputs)
  nTrials <- nrow(Inputs)
  Data = vector (mode = "list", length = nOptions)
  Input = names (Inputs) = names (Data) = 1:nOptions
  advantages <- t(apply(Inputs, 1, ToM))
  sums <- t(apply(Inputs, 1, ToP))
  v = v0 + ws*sums + wd*advantages
  #  s = 1 + sl * (n-3) * v
  nAccumulators <- ncol(v)
  
  pl = list(t0=matrix(t0, nrow=nTrials, ncol=nAccumulators), 
            v=v, 
            A=matrix(A, nrow=nTrials, ncol=nAccumulators),
            B=matrix(B, nrow=nTrials, ncol=nAccumulators), 
            s=matrix(s, nrow=nTrials, ncol=nAccumulators))
  
  # Create index matrix
  is.winner <- matrix(FALSE, nrow=nTrials, ncol=nAccumulators)
  respMat <- (matrix(R, nrow=nTrials, ncol=2, byrow=FALSE)-1)*2+1
  respMat <- respMat + matrix(seq(0,ncol(respMat)-1), nrow=nTrials, ncol=ncol(respMat), byrow=TRUE)  
  is.winner[cbind(1:nrow(respMat), respMat[,1])] <- TRUE
  is.winner[cbind(1:nrow(respMat), respMat[,2])] <- TRUE
  
  vWinner <- matrix(t(v)[t(is.winner)], nrow=nrow(v), byrow=TRUE)
  
  ## Get pdf of accumulators corresponding to WINNING responses
  k <- (pl$B+pl$A/2)[,1:2]
  l <- vWinner
  nu <- k[,1:2]/l
  nu.inf <- is.infinite(nu) | is.na(nu)
  lambda <- (k/pl$s[,1:2])^2
  #  t <- 0.5
  times <- matrix(RT-t0, nrow=nTrials, ncol=2)
  pdf <- matrix(0, nrow=nTrials, ncol=2)
  pdf[!nu.inf] <- pmax(t(dinvGauss(times[!nu.inf], 
                                   nu=nu[!nu.inf], 
                                   lambda=lambda[!nu.inf])), 0, na.rm=TRUE)
  
  # Get CDFs of *ALL* accumulators
  times <- matrix(RT-t0, nrow=nTrials, ncol=nAccumulators)
  cdf <- matrix(pWald(t=times, v=v, A=pl$A, B=pl$B, s=pl$s), nrow=nrow(v))
  
  # store Pr (at least one member of other response sets IS NOT done) = (1-Pr(member 1 done)*...*1-Pr(member r done))
  # Product of probabilities of all accumulators corresponding to each response
  PrNotDone <- 1-do.call(cbind, lapply(1:3, function(x) apply(cdf[,((x-1)*2+1):((x-1)*2+2)], 1, prod)))
  
  # Likelihood:
  # probability = probability * probability of other for response having finished * probability of, for the other responses, at least 1 accumulator not finished
  # probability_{1} = pdf_{1-2} * cdf_{1-3} * (1-(CDF_{2-3}*CDF_{2-1})) * (1-(CDF_{3-2}*CDF_{3-1}))) + 
  #                   pdf_{1-3} * cdf_{1-2} * (1-(CDF_{2-3}*CDF_{2-1})) * (1-(CDF_{3-2}*CDF_{3-1}))) =
  # probability_{1} = (pdf_{1-2} * cdf_{1-3} + 
  #                    pdf_{1-3} * cdf_{1-2}) * (1-(CDF_{2-3}*CDF_{2-1})) * (1-(CDF_{3-2}*CDF_{3-1}))) =
  
  # find CDFs of winners
  cdfWinner <- matrix(t(cdf)[t(is.winner)], nrow=nrow(cdf), byrow=TRUE)
  PrFinal = pdf[,1]*cdfWinner[,2] + pdf[,2]*cdfWinner[,1]
  
  # multiply by probability of at least one other accumulator per other response not having won
  indx = !is.winner[,c(1,3,5)]
  PrFinal = PrFinal * apply(matrix(t(PrNotDone)[t(indx)], ncol=2, byrow=TRUE),1,prod)
  
  PrFinal
}


############################### RANDOM DATA FOR WA

WArnd = function (x, N, Trials, Inputs, sl = 0, Fun = "tLBA")
  # x = parameter vector
  # N = number of inputs
  # Trials = nTrials
  # Inputs = values
  ## CODE FROM Van Ravenzwaaij et al (2019) https://osf.io/2s6ax/files/
{
  nms = 1:N
  names (Inputs) = nms
  tmp = outer (X = nms, Y = nms, FUN = paste, sep = "-")
  pair.nms = c(tmp[upper.tri(tmp)], tmp[lower.tri(tmp)]) # All pairs.
  
  watch.sets = array (dim = c(N-1, N, 1), dimnames = list (paste ("watch", 1:(N-1), sep = ""), nms, c("win")))
  for (k in 1:N) {watch.sets[,k,"win"] = grep (paste (nms[k], "-", sep = ""), pair.nms)}
  
  #v = x["v0"] + x["wS"] * ToP (Inputs) + x["wD"] * ToM (Inputs)
  v = x["V0"] + x["wS"] * ToP (Inputs) + x["wV"] * ToM (Inputs)
  params = array (dim = c(length(pair.nms), 5), dimnames = list (pair.nms, c("v", "B", "A", "s", "t0")))
  params[,"v"] = v
  params[,"B"] = rep (x["B"], length(pair.nms))
  params[,"A"] = rep (x["A"], length(pair.nms))
  params[,"s"] = rep(x['s'], length(pair.nms)) #1 + sl * (N-3) * v
  params[,"t0"] = rep (x["t0"], length(pair.nms))
  times = array (dim = c(length(pair.nms), Trials), dimnames = list (pair.nms, NULL))
  for (k in pair.nms) {times[k,] = rWald(n=Trials, B=params[k,'B'], v=params[k,'v'], A=params[k,'A'])+params[k,'t0'] }  #raccumulator (n = Trials, x = params[k,], dist = Fun)}
  
  # Times for each relevant watch.set of accumulators to finish.
  set.times = array (dim = c(N, Trials), dimnames = list (nms, 1:Trials))
  if(Trials > 1) {
    for (k in 1:N) if (N==2) set.times[k,] = times[watch.sets[,k,1],] else set.times[k,] = apply (times[watch.sets[,k,1],], 2, max)
  } else {
    for (k in 1:N) if (N==2) set.times[k,] = times[watch.sets[,k,1],] else set.times[k,] = max(times[watch.sets[,k,1],])
  }
  
  # Now find the final thing: decision, and RT, for the two models.
  data.frame (RT = apply (set.times, 2, min), R = apply (set.times, 2, which.min))
}
