library(dmcAdapt)

# Template setup for n-choice LBA, B=b-A parameterization
#   External parameters types: A, B, t0, mean_v, sd_v, st0 = 0 (optional)
#   Internal parameters types: A, b, t0, mean_v, sd_v, st0 = 0 (optional)

# This ncl (no cell likelihood) version uses p.list.dmc in the likelihood to 
# avoid looping over design cells.

# pi1=get.p(p.list,i+1);Si=S[i]; i1=S[i+1];a=p[c("SP","SR")];FBi=facs[i,"CR"]-1
adapt.dmc <- function(pi1,Si1=NULL,a=NULL,Si=NULL,FBi=NULL) 
  # pi1 is a list of parameters for current trial
  # a is a list of quantities to be adapted (if absent no update)
  # S/FB: stimulus/response feedback (0/1) on current/previous trial
{

  # Threhsold learning function
  # SP is stimulus probability, R is the response (0 or 1)
  dSP <- function(SP,FB,alpha=.1) SP + alpha*(FB-SP) 

  # Stimulus learning function
  # SR is the representation, R is reponse
  dSR <- function(SR,S,alpha=.1) SR + alpha*(S-SR) 

  if ( !is.null(a) ) { # Update
    # Learn bias
    if (pi1$aB[1] == 0) pi1$SP <- a$SP else  
      pi1$SP <- c(dSP(SP=a$SP[1],FB=FBi, alpha=pi1$aB[1]),
                dSP(SP=a$SP[2],FB=!FBi,alpha=pi1$aB[1]))
    # Learn stimulus representation
    if (pi1$aV[1] != 0) {
      if (FBi==0)
        pi1$SR <-         c(dSR(SR=a$SR[1],S=Si,alpha=pi1$aV[1]),a$SR[2]) else
        pi1$SR <- c(a$SR[1],dSR(SR=a$SR[2],S=Si,alpha=pi1$aV[1]))
    } else pi1$SR <- a$SR
  } 
  # wB=2 => b <- A + 2*B0*(1-SP), full proportional effect of prob.
  # When SP=.5 b = A + B0
  if (pi1$aB[1] == 0) pi1$b <- pi1$A + pi1$B0 else
    pi1$b <- pi1$A + pi1$B0*(1 + pi1$wB*(pi1$SP[2]-.5)*c(1,-1))
  # ALBA rates, difference only version 
  if (pi1$aV[1] == 0) pi1$mean_v <- pi1$V0 else
    pi1$mean_v <- c(pi1$V0[1] + pi1$wV[1]*(pi1$SR[2]-Si1), 
                    pi1$V0[2] + pi1$wV[2]*(Si1-pi1$SR[1])) 
  pi1
}


adapt.r.dmc <- function(pi1,Si1=NULL,a=NULL,Si=NULL,FBi=NULL) 
  # pi1 is a list of parameters for current trial
  # a is a list of quantities to be adapted 
  # S/FB: stimulus/response feedback (0/1) on current/previous trial
{

  # Threhsold learning function
  # SP is stimulus probability, R is the response (0 or 1)
  dSP <- function(SP,FB,alpha=.1) SP + alpha*(FB-SP) 

  # Stimulus learning function
  # SR is the representation, R is reponse
  dSR <- function(SR,S,alpha=.1) SR + alpha*(S-SR) 

  # Learn bias
  if (pi1[1,"aB"] != 0) pi1[,"SP"] <- 
                c(dSP(SP=a[1,"SP"],FB=FBi, alpha=pi1[1,"aB"]),
                  dSP(SP=a[2,"SP"],FB=!FBi,alpha=pi1[1,"aB"]))
  # Learn stimulus representation
  if (pi1[1,"aV"] != 0) {
    if (FBi==0)
      pi1[,"SR"] <- c(dSR(SR=a[1,"SR"],S=Si,alpha=pi1[1,"aV"]),a[2,"SR"]) else
      pi1[,"SR"] <- c(a[1,"SR"],dSR(SR=a[2,"SR"],S=Si,alpha=pi1[1,"aV"]))
  }
  pi1
}


# save.adapt=TRUE
random.dmc <- function(p.list,model,save.adapt=FALSE)
{
  
  get.p <- function(p.list,i,exclude=c("SP","SR"))
    lapply(p.list[!(names(p.list) %in% exclude)],function(x){x[,i]})
  
  n <- dim(p.list[[1]])[2]  # number of trials
  cvs <- attr(p.list,"cvs")
  facs <- attr(p.list,"facs") # factor columns + CR = numeric correct response
  # S <- as.numeric(facs$S)-1 # use stimulus indicator as binary real stim value
  S <- cvs$stim # real valued stimulus
  # Get correct response as feedback and possibly corrupt with prbability 1-pcFB
  FB <- as.numeric(facs$S)-1
  FBresponse <- p.list$pcFB[1,] < 0
  if ( any(!FBresponse) ) {
    flip <- rbinom(n,1,p.list$pcFB[1,!FBresponse])==0
    FB[!FBresponse][flip] <- as.numeric(!FB[!FBresponse][flip])
  }
  # Update SP and SR
  p <- get.p(p.list,1,exclude=NULL)
  p <- adapt.dmc(p,Si1=S[1]) # for first iteration
  out <- matrix(nrow= n, ncol=2,dimnames=list(NULL,c("RT","R")))
  if (save.adapt) {
    adapt <- vector(mode="list",length=n)
    adapt[[1]] <- p[c("SP","SR","b","mean_v","sd_v","t0","A")] 
  }
  for (i in 1:n) {

    out[i,] <- as.numeric(rlba.norm(1,
      A=p$A,b=p$b,t0=p$t0,st0=p$st0[1],mean_v=p$mean_v,sd_v=p$sd_v,
      posdrift = attr(model,"posdrift")))

# if (out[i,"R"]==1) s=1:2 else s=2:1
# A=matrix(p$A[s],nrow=1);b=matrix(p$b[s],nrow=1)
# sd_v=matrix(p$sd_v[s],nrow=1);mean_v=matrix(p$mean_v[s],nrow=1)
# like[i] <- n1PDFfixedt0.norm(dt=out[i,"RT"]-p$t0[1],A=A,b=b,sd_v=sd_v,mean_v=mean_v,
#   posdrift=attr(model,"posdrift"))

    if (i != n) { # update for next iteration
      # Guide adaptation with last response
      if ( FBresponse[i] & (runif(1) <= -p.list$pcFB[1,i]) ) FB[i] <- out[i,"R"]-1
      p <- adapt.dmc(pi1=get.p(p.list,i+1),Si1=S[i+1],
                     Si=S[i],a=p[c("SP","SR")],FBi=FB[i])
      if (save.adapt) adapt[[i+1]] <- p[c("SP","SR","b","mean_v","sd_v","t0","A")]
    }
  }
  attr(out,"cvs") <- cbind(cvs,FB=FB+1)
  if (save.adapt) attr(out,"adapt") <- do.call(rbind,lapply(adapt,unlist))
  out
}


transform.dmc <- function(par.df,do.trans=TRUE) 
{
  
  if (do.trans)
    list(A=t(par.df$A),sd_v=t(par.df$sd_v),t0=t(par.df$t0),st0=t(par.df$st0),
         SP=t(par.df$SP),aB=t(par.df$aB),B0=t(par.df$B0),wB=t(par.df$wB),
         SR=t(par.df$SR),aV=t(par.df$aV),V0=t(par.df$V0),wV=t(par.df$wV),
         pcFB=t(par.df$pcFB)) else
    list(A=par.df$A,sd_v=par.df$sd_v,t0=par.df$t0,st0=par.df$st0,
         SP=par.df$SP,aB=par.df$aB,B0=par.df$B0,wB=par.df$wB,
         SR=par.df$SR,aV=par.df$aV,V0=par.df$V0,wV=par.df$wV,
         pcFB=par.df$pcFB)       
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
    nacc = dim(pars)[2]
    feedback <- matrix(NA, nrow=nrow(cvs), ncol=nacc+1)
    for(accumulator in 1:nacc) {
      feedback[cvs$FB==accumulator, accumulator] = cvs$stim[cvs$FB==accumulator]
    }
    feedback[,nacc+1] = cvs$FB-1 # since probabilities must sum to 1, only one of the SPs need to be updated
  
    # Create start point vector
    startValues <- c(pars[1,,c('SR')], pars[1,1,c('SP')])
    
    # Create learning rates vector
    learningRates <- cbind(pars[,,'aV'], pars[,1,'aB'])
  
    # call C
    updated <- adapt.c.dmc(startValues = startValues, 
                           learningRates = learningRates, 
                           feedback = feedback)
    
    # add back data to 'pars' array
    pars[,,'SR'] <- updated$adaptedValues[,c(1,2)]
    pars[,,'SP'] <- cbind(updated$adaptedValues[,3], 1-updated$adaptedValues[,3])
    ## end SM
  } else {
    # Update representations (SP and SR)
    for (i in 1:(dim(data)[1]-1)) {
      p <- adapt.r.dmc(pi1=pars[i+1,,c("SP","SR","aB","aV")],
                       Si1=cvs$stim[i+1],Si=cvs$stim[i],a=pars[i,,c("SP","SR")],FBi=cvs[i,"FB"]-1)
      pars[i+1,,c("SP","SR")] <- as.matrix(data.frame(p)[,c("SP","SR")])
    }
  }
  
  # Update parameters
  b <- pars[,,"A"] + pars[,,"B0"]*(1 + pars[,,"wB"]*(pars[,2:1,"SP"]-.5))
  mean_v <- cbind(pars[,1,"V0"] + pars[,1,"wV"]*(pars[,2,"SR"]-cvs$stim),
                  pars[,2,"V0"] + pars[,2,"wV"]*(cvs$stim-pars[,1,"SR"]))
  
  
  # all(pars[,,"SR"]==adapt[,c("SR.r1","SR.r2")])
  # all(round(pars[,,"SP"],4)==round(adapt[,c("SP.r1","SP.r2")],4))
  # all(pars[,,"A"]==adapt[,c("A.r1","A.r2")])
  # all(pars[,,"sd_v"]==adapt[,c("sd_v.r1","sd_v.r2")])
  # all(pars[,,"t0"]==adapt[,c("t0.r1","t0.r2")])
  # all(mean_v==adapt[,c("mean_v.r1","mean_v.r2")])
  # all(round(b,4)==round(adapt[,c("b.r1","b.r2")],4))

  # all(pars[,,"mean_v"]==adapt[,c("mean_v.r1","mean_v.r2")])
  # all(pars[,,"b"]==adapt[,c("b.r1","b.r2")])
  
  pmax(n1PDFfixedt0.norm(dt=data$RT-pars[,1,"t0"],
                         A=do.n1(pars[,,"A"]),
                         b=do.n1(b),
                         mean_v=do.n1(mean_v),
                         sd_v=do.n1(pars[,,"sd_v"]),
                         posdrift=attr(attr(data,"model"),"posdrift")),min.like)
}



