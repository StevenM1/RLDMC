library(SuppDists)

check_n1_arguments <- function(arg, nn, n_v, dots = FALSE) {
  mc <- match.call()
  varname <- sub("dots$", "", deparse(mc[["arg"]]), fixed = TRUE)
  if (!is.list(arg)) {
    if ((!is.vector(arg, "numeric")) || (length(arg) < 1))
      stop(paste(varname, "needs to be a numeric vector of length >= 1!"))
    if (dots) {
      arg <- as.list(arg)
      arg <- lapply(arg, rep, length.out=nn)
    } else arg <- rep(arg, length.out=nn)
  } else {
    if (!dots && (length(arg) != n_v))
      stop(paste("if", varname, "is a list, its length needs to correspond to the number of accumulators."))
    for (i in seq_along(arg)) {
      if ((!is.vector(arg[[i]], "numeric")) || (length(arg[[i]]) < 1))
        stop(paste0(varname, "[[", i, "]] needs to be a numeric vector of length >= 1!"))
      arg[[i]] <- rep(arg[[i]], length.out=nn)
    }
  }
  return(unname(arg))
}

check_i_arguments <- function(arg, nn, n_v, dots = FALSE) {
  mc <- match.call()
  varname <- sub("dots$", "", deparse(mc[["arg"]]), fixed = TRUE)
  if (!is.list(arg)) {
    if ((!is.vector(arg, "numeric")) || (length(arg) < 1))
      stop(paste(varname, "needs to be a numeric vector of length >= 1!"))
    if (dots) {
      arg <- as.list(arg)
      arg <- lapply(arg, rep, length.out=nn)
    } else arg <- lapply(seq_len(n_v), function(x) rep(arg, length.out=nn))
  } else {
    if (!dots && (length(arg) != n_v))
      stop(paste("if", varname, "is a list, its length needs to correspond to the number of accumulators."))
    for (i in seq_along(arg)) {
      if ((!is.vector(arg[[i]], "numeric")) || (length(arg[[i]]) < 1))
        stop(paste0(varname, "[[", i, "]] needs to be a numeric vector of length >= 1!"))
      arg[[i]] <- rep(arg[[i]], length.out=nn)
    }
  }
  #if (length(arg) != n_v) stop(paste("size of", varname, "does not correspond to number of accumulators."))
  return(arg)
}


rWald <- function(n,B,v,A,s=1)
  # random function for single acumulator
{
  
  rwaldt <- function(n,k,l,s=1,tiny=1e-6) {
    # random sample of n from a Wald (or Inverse Gaussian)
    # k = criterion, l = rate, assumes sigma=1 Browninan motion
    # about same speed as statmod rinvgauss
    
    rlevy <- function(n=1, m=0, c=1) {
      if (any(c<0)) stop("c must be positive")
      c/qnorm(1-runif(n)/2)^2+m
    }
    
    flag <- l>tiny
    x <- rep(NA,times=n)
    
    x[!flag] <- rlevy(sum(!flag),0,k[!flag]^2)
    mu <- k/l
    lambda <- (k/s)^2
    
    y <- rnorm(sum(flag))^2
    mu.0 <- mu[flag]
    lambda.0 <- lambda[flag]
    
    x.0 <- mu.0 + mu.0^2*y/(2*lambda.0) -
      sqrt(4*mu.0*lambda.0*y + mu.0^2*y^2)*mu.0/(2*lambda.0)
    
    z <- runif(length(x.0))
    test <- mu.0/(mu.0+x.0)
    x.0[z>test] <- mu.0[z>test]^2/x.0[z>test]
    x[flag] <- x.0
    x[x<0] <- max(x)
    x
  }
  
  # Act as if negative v never terminates, cluge to do single accumulator
  # case by passing negative v
  if (length(v)!=n) v <- rep(v,length.out=n)
  if (length(B)!=n) B <- rep(B,length.out=n)
  if (length(A)!=n) A <- rep(A,length.out=n)
  if (length(s)!=n) s <- rep(s,length.out=n)
  
  # Kluge to return -Inf for negative rates, so can implment one accumulator case
  out <- numeric(n)
  ok <- v>0
  nok <- sum(ok)
  bs <- B[ok]+runif(nok,0,A[ok])
  out[ok] <- rwaldt(nok,k=bs,l=v[ok],s=s[ok])
  out[!ok] <- Inf
  out
}

# dWaldSuppDists <- function(t, k, l, alpha=1) {
#   ## transform k (criterion) / l (rate) parametrisation to nu/lambda parametrisation
#   nu = k / l
#   lambda = (k^2)/(alpha^2)
#   dinvGauss(t, nu=nu, lambda=lambda)
# }
# 

dWald <- function(t,v,B,A,s=1,useSuppDists=T)
  # density for single accumulator
{
  digt <- function(t,k=1,l=1,a=.1,s=1,tiny=1e-10) {
    # pdf of inverse gaussian at t with k +/- a/2 uniform variability
    # returns digt.0 if a<1e-10
    
    digt.0 <- function(t,k=1,l=1, s=1) {
      # pdf of inverse gaussian at t with no k variability
      # much faster than statmod's dinvgauss funciton
      
      lambda <- (k/s)^2
      l0 <- l==0
      e <- numeric(length(t))
      if ( any(!l0) ) {
        mu <- k[!l0]/l[!l0]
        e[!l0] <- -(lambda[!l0]/(2*t[!l0])) * (t[!l0]^2/mu^2 - 2*t[!l0]/mu  + 1)
      }
      if ( any(l0) )  e[l0] <- -.5*lambda[l0]/t[l0]
      x <- exp(e + .5*log(lambda) - .5*log(2*t^3*pi))
      x[t<=0] <- 0
      x
    }
    
    options(warn=-1)
    if(length(k)!=length(t)) k <- rep(k,length.out=length(t))
    if(length(l)!=length(t)) l <- rep(l,length.out=length(t))
    if(length(a)!=length(t)) a <- rep(a,length.out=length(t))
    if(length(s)!=length(t)) s <- rep(s,length.out=length(t))
    tpos <- t <=0
    atiny <- a <= tiny & !tpos
    a[atiny] <- 0
    
    ltiny <- (l<=tiny) & !atiny & !tpos
    notltiny <- (l>tiny) & !atiny & !tpos
    l[l<=tiny] <- 0
    
    x <- numeric(length(t))
    # No threshold variability
    if ( any(atiny) )
      if(useSuppDists) {
        nu <- k[atiny]/l[atiny]
        lambda <- (k[atiny]/s[atiny])^2
        nu.inf <- is.infinite(nu) | is.na(nu)
        z <- x[atiny]
        z[nu.inf] <- 0
        z[!nu.inf] <- dinvGauss(t[atiny][!nu.inf], 
                                nu=nu[!nu.inf], 
                                lambda=lambda[!nu.inf])
        x <- z
      } else {
        x[atiny] <- digt.0(t=t[atiny],k=k[atiny],l=l[atiny],s=s[atiny])
      }
    
    # Threshold variability. CANNOT DO TRIAL-BY-TRIAL s!!
    if ( any(!atiny) ) {
      
      if ( any(notltiny) ) { # rate non-zero
        
        sqr.t <- sqrt(t[notltiny])
        
        term.1a <- -(a[notltiny]-k[notltiny]+t[notltiny]*l[notltiny])^2/(2*t[notltiny])
        term.1b <- -(a[notltiny]+k[notltiny]-t[notltiny]*l[notltiny])^2/(2*t[notltiny])
        term.1 <- (exp(term.1a) - exp(term.1b))/sqrt(2*pi*t[notltiny])
        
        term.2a <- log(.5)+log(l[notltiny])
        term.2b <- 2*pnorm((-k[notltiny]+a[notltiny])/sqr.t+sqr.t*l[notltiny])-1
        term.2c <- 2*pnorm((k[notltiny]+a[notltiny])/sqr.t-sqr.t*l[notltiny])-1
        term.2d <- term.2b+term.2c
        term.2 <- exp(term.2a)*term.2d
        
        term.3 <- term.1+term.2
        term.4 <- log(term.3)-log(2)-log(a[notltiny])
        x[notltiny] <- exp(term.4)
      }
      
      if ( any(ltiny) ) {  # rate zero
        log.t <- log(t[ltiny])
        term.1 <- -.5*(log(2)+log(pi)+log.t)
        term.2 <- (k[ltiny]-a[ltiny])^2/(2*t[ltiny])
        term.3 <- (k[ltiny]+a[ltiny])^2/(2*t[ltiny])
        term.4 <- (exp(-term.2)-exp(-term.3))
        term.5 <- term.1+log(term.4) - log(2) - log(a[ltiny])
        x[ltiny] <- exp(term.5)
      }
      
    }
    
    x[x<0 | is.nan(x) ] <- 0
    x
  }
  out <- numeric(length(t))
  ok <- v>0
  B <- unlist(B)
  A <- unlist(A)
  out[ok] <- digt(t[ok],k=B[ok]+A[ok]/2,l=v[ok],a=A[ok]/2,s=s[ok])
  out[!ok] <- 0
  out
}






pWald <- function(t,v,B,A,s=1)
  # cumulative density for single accumulator
{
  pigt <- function(t,k=1,l=1,a=.1,s=1,tiny=1e-10) {
    # cdf of inverse gaussian at t with k +/- a/2 uniform variability
    # returns pigt.0 if a<=0
    
    pigt.0 <- function(t,k=1,l=1,s=1) {
      # cdf of inverse gaussian at t with no k variability
      # much faster than statmod's pinvgauss funciton
      
      mu <- k/l
      lambda <- (k/s)^2
      
      e <- exp(log(2*lambda) - log(mu))
      add <- sqrt(lambda/t) * (1 + t/mu)
      sub <- sqrt(lambda/t) * (1 - t/mu)
      
      p.1 <- 1 - pnorm(add)
      p.2 <- 1 - pnorm(sub)
      x <- exp(e + log(p.1)) + p.2
      
      x[t<0] <- 0
      x
    }
    
    options(warn=-1)
    if(length(k)!=length(t)) k <- rep(k,length.out=length(t))
    if(length(l)!=length(t)) l <- rep(l,length.out=length(t))
    if(length(a)!=length(t)) a <- rep(a,length.out=length(t))
    if(length(s)!=length(t)) s <- rep(s,length.out=length(t))
    
    tpos <- t<=0
    
    atiny <- a<=tiny & !tpos
    a[atiny] <- 0
    
    ltiny <- (l<=tiny) & !atiny & !tpos
    notltiny <- (l>tiny) & !atiny & !tpos
    l[l<=tiny] <- 0
    
    x <- numeric(length(t))
    
    # No threshold variability
    if ( any(atiny) )
      x[atiny] <- pigt.0(t[atiny],k[atiny],l[atiny],s=s[atiny])
    
    # Threshold variability
    if ( any(!atiny) ) {
      
      if ( any(notltiny) ) { # rate non-zero
        
        log.t <- log(t[notltiny])
        sqr.t <- sqrt(t[notltiny])
        
        term.1a <- .5*log.t-.5*log(2*pi)
        term.1b <- exp(-((k[notltiny]-a[notltiny]-t[notltiny]*l[notltiny])^2/t[notltiny])/2)
        term.1c <- exp(-((k[notltiny]+a[notltiny]-t[notltiny]*l[notltiny])^2/t[notltiny])/2)
        term.1 <- exp(term.1a)*(term.1b-term.1c)
        
        term.2a <- exp(2*l[notltiny]*(k[notltiny]-a[notltiny]) +
                         log(pnorm(-(k[notltiny]-a[notltiny]+t[notltiny]*l[notltiny])/sqr.t)))
        term.2b <- exp(2*l[notltiny]*(k[notltiny]+a[notltiny]) +
                         log(pnorm(-(k[notltiny]+a[notltiny]+t[notltiny]*l[notltiny])/sqr.t)))
        term.2 <- a[notltiny] + (term.2b-term.2a)/(2*l[notltiny])
        
        term.4a <- 2*pnorm((k[notltiny]+a[notltiny])/sqr.t-sqr.t*l[notltiny])-1
        term.4b <- 2*pnorm((k[notltiny]-a[notltiny])/sqr.t-sqr.t*l[notltiny])-1
        term.4c <- .5*(t[notltiny]*l[notltiny] - a[notltiny] - k[notltiny] + .5/l[notltiny])
        term.4d <- .5*(k[notltiny] - a[notltiny] - t[notltiny]*l[notltiny] - .5/l[notltiny])
        term.4 <- term.4c*term.4a + term.4d*term.4b
        
        x[notltiny] <- (term.4 + term.2 + term.1)/(2*a[notltiny])
      }
      
      if ( any(ltiny) ) {  # rate zero
        sqr.t <- sqrt(t[ltiny])
        log.t <- log(t[ltiny])
        term.5a <- 2*pnorm((k[ltiny]+a[ltiny])/sqr.t)-1
        term.5b <- 2*pnorm(-(k[ltiny]-a[ltiny])/sqr.t)-1
        term.5 <- (-(k[ltiny]+a[ltiny])*term.5a - (k[ltiny]-a[ltiny])*term.5b)/(2*a[ltiny])
        
        term.6a <- -.5*(k[ltiny]+a[ltiny])^2/t[ltiny] - .5*log(2) -.5*log(pi) + .5*log.t - log(a[ltiny])
        term.6b <- -.5*(k[ltiny]-a[ltiny])^2/t[ltiny] - .5*log(2) -.5*log(pi) + .5*log.t - log(a[ltiny])
        term.6 <- 1 + exp(term.6b) - exp(term.6a)
        
        x[ltiny] <- term.5 + term.6
      }
      
    }
    
    x[x<0 | is.nan(x) ] <- 0
    x
  }
  
  out <- numeric(length(t))
  ok <- v>0
  B <- unlist(B)
  A <- unlist(A)
  out[ok] <- pigt(t[ok],k=B[ok]+A[ok]/2,l=v[ok],a=A[ok]/2,s=s[ok])
  out[!ok] <- 0
  out
  
}

rWaldRace <- function(n,v,B,A,t0,s,gf=0,return.ttf=FALSE)
  # random function for Wald race.
{
  B[B<0] <- 0 # Protection for negatives
  A[A<0] <- 0
  n_v  <- ifelse(is.null(dim(v)), length(v), dim(v)[1])
  ttf <- matrix(t0 + rWald(n*n_v,B=B,v=v,A=A,s=s), nrow=n_v)
  if (return.ttf) return(ttf)
  resp <- apply(ttf, 2, which.min)
  out <- data.frame(RT = ttf[cbind(resp,1:n)], R = apply(ttf, 2, which.min))
  
  if (gf[1] > 0) {
    is.gf <- as.logical(rbinom(dim(out)[1],1,gf))
    out$RT[is.gf] <- NA
    out$R[is.gf] <- 1
  }
  
  out
}

rWaldRaceSM <- function(n, A, B, t0, v, s, st0 = 0, silent = FALSE, args.dist=NULL, simPerTrial=FALSE, timer = F)
{
  if(is.list(v)) {
    v <- matrix(unlist(v), ncol=length(v))
    if(nrow(v) != n & !silent) warning('Number of trials does not equal number of per-trial drift rates...')
    if(length(B) == 1) {
      B <- rep(B, n)
    }
    if(length(s) == 1) {
      s <- rep(s, n)
    }
    out <- data.frame(RT=NA, R=NA)
    if(simPerTrial) {
      for(i in 1:nrow(v)) {
        out[i,] <- rWaldRace(n=n, v=v[i,], B=B[i], A=A, t0=t0, s=s)
      }    
    } else {
      if(is.list(B)) B <- matrix(unlist(B), ncol=length(B))
      if(is.list(A)) A <- matrix(unlist(A), ncol=length(A))
      if(is.list(t0)) t0 <- matrix(unlist(t0), ncol=length(t0))
      if(is.list(s)) s <- matrix(unlist(s), ncol=length(s))
      B[B<0] <- 0 # Protection for negatives
      A[A<0] <- 0
      bs <- B + runif(length(B), 0, A)
      n_v  <- ifelse(is.null(dim(v)), length(v), dim(v)[2])
      n <- nrow(v)
      ttf <- matrix(NA, nrow=ncol(v), ncol=nrow(v))
      if (timer){n_v <- 1}
      for(i in 1:n_v) {
        ttf[i,] <- rinvGauss(n, nu=bs[,i]/v[,i], lambda=(bs[,i]/s[,i])^2)
      }
      ttf[is.na(ttf)] <- Inf
      if (any(ttf<0)) browser()
      if(timer){
        ttf <- ttf + t0[1,]
      }
      else{
        ttf <- ttf + t(t0)
      }
      
      if (!timer){
        resp <- apply(ttf, 2, which.min)
        out <- data.frame(RT = ttf[cbind(resp,1:n)], R = apply(ttf, 2, which.min))
      } else{
        out <- data.frame(RT = ttf[1,])
      }
      
      
    }
    
  } else {
    # all equal drifts, nice and fast
    out <- rWaldRace(n=n, v=v, B=B, A=A, t0=t0, s=s)
  }
  out
}

rWaldRace_Timing_NS <- function(n, t0T, t0E, v_E, v_T, s_E, s_T, B_T, B_E, A, st0 = 0, silent=TRUE, simPerTrial=FALSE){
  out <- rWaldRaceSM(n= n, t0 = t0E,A = A, B = B_E, v = v_E, s = s_E, st0 = 0, silent=TRUE, simPerTrial=FALSE)
  out$timer.response = F
  ttf_T <- rWaldRaceSM(n= n, t0 = t0T, A = A, B = B_T, v = v_T, s = s_T, st0 = 0, silent=TRUE, simPerTrial=FALSE, timer = T)
  timer.replace <- ttf_T$RT < out$RT
  n.acc <- ifelse(is.null(dim(v_E)), length(v_E), dim(v_E)[1])
  if(any(timer.replace)) {
    out$timer.response[timer.replace] <- TRUE
    out$RT[timer.replace] <- ttf_T$RT[timer.replace]
    # then determine which response was given
    # unbiased guessing probablity of 1/N
    p.guess <- rep(1/n.acc, n.acc)
    out$R[timer.replace] <- sample(1:n.acc, size=sum(timer.replace), replace=TRUE, prob=p.guess)
  }
  out
}


n1WaldTiming_NS <- function(dt,B,A,t0, B_T, v_T, s_T, A_T, t0T, p.guess = 0.5, gf=0, ...)
  # Generates defective PDF for responses on node=1, dt (decison time) is a vector of times
{
  ### Added by SM/NS
  dots <- list(...)
  n_v <- max(vapply(dots, length, 0))  # Number of responses
  #  if(!silent) message(paste("Results based on", n_v, "accumulators/drift rates."))
  if (n_v < 2) stop("There need to be at least two accumulators/drift rates.")
  #  distribution <- match.arg(distribution)
  #check_single_arg(t0 = t0)
  nn <- length(dt) #dt <- decision times
  B_E <- check_n1_arguments(B, nn=nn, n_v=n_v)
  v_E <- check_n1_arguments(dots[[1]], nn=nn, n_v=n_v)
  s_E <- check_n1_arguments(dots[[2]], nn=nn, n_v=n_v)
  A_E <- check_n1_arguments(A, nn=nn, n_v=n_v)
  t0E <- check_n1_arguments(t0, nn=nn, n_v=n_v)
  B_T <- check_n1_arguments(B_T, nn=nn, n_v=n_v)
  v_T <- check_n1_arguments(v_T, nn=nn, n_v=n_v)
  s_T <- check_n1_arguments(s_T, nn=nn, n_v=n_v)
  A_T <- check_n1_arguments(A_T, nn=nn, n_v=n_v)
  t0T <- check_n1_arguments(t0T, nn=nn, n_v=n_v)
  #  B[B<0] <- 0 # Protection for negatives
  #  A[A<0] <- 0
  # n_acc <- ifelse(is.null(dim(v)),length(v),dim(v)[1])
  # if (is.null(dim(dt))) dt <- matrix(rep(dt,each=n_acc),nrow=n_acc)
  if (is.null(dim(dt))) dt <- matrix(rep(dt, each=n_v), nrow=n_v)
  nonNa <- !is.na(dt[1,]) #Should all be non-NA right? (as long as response is self-timed)
  #  if (!is.matrix(v)) v <- matrix(rep(v,n.go),nrow=n_v)
  if (!is.matrix(v_E)) v_E <- do.call(rbind, v_E) #matrix(unlist(v), nrow=n_v)
  if (!is.matrix(A_E)) A_E <- do.call(rbind, A_E) # matrix(unlist(A), nrow=n_v)
  if (!is.matrix(B_E)) B_E <- do.call(rbind, B_E) #matrix(unlist(B), nrow=n_v)
  if (!is.matrix(s_E)) s_E <- do.call(rbind, s_E)
  if (!is.matrix(v_T)) v_T <- do.call(rbind, v_T) #matrix(unlist(v), nrow=n_v)
  if (!is.matrix(A_T)) A_T <- do.call(rbind, A_T) # matrix(unlist(A), nrow=n_v)
  if (!is.matrix(B_T)) B_T <- do.call(rbind, B_T) #matrix(unlist(B), nrow=n_v)
  if (!is.matrix(s_T)) s_T <- do.call(rbind, s_T)
  #  if (!is.matrix(B)) B <- matrix(rep(B,n.go),nrow=n_v)
  #  if (!is.matrix(A)) A <- matrix(rep(A,n.go),nrow=n_v)
  
  ################ Timing accumulator from Hawkins, Heathcote 2020 Racing Against the Clock, added by NS
  # E = Evidence Accumulator, T = Timing Accumulator
  ### timing process starts earlier than evidence process ###
  if(t0T[[1]] < t0E[[1]]) {
    dt[1, nonNa] <- dt[1, nonNa] - t0T[[1]]
    # assume evidence process starts at d=t0_E-t0_T after the timing process
    d <- t0E[[1]] - t0T[[1]] #Assume this is the new t0
    u <- dt[1, nonNa] < d
    if(any(u)) {  # evidence process hasn't started so must be a timer response with a guess
      #Density of Timer accumulator * guessing probability
      dt[1, u] <- dWald(dt[1,u], v=v_T[1,u], B=B_T[1,u], A=A_T[1,u], s=s_T[1,u]) * p.guess
    }
    if(any(!u)) { # could be a response from the timer or evidence process
      #Evidence Accumulator is density accumulator * cumulative distribution function other Evidence Accumulators
      pdf_E <- dWald(dt[1,!u]-d[!u], v=v_E[1,!u], B=B_E[1,!u], A=A_E[1,!u], s=s_E[1,!u]) 
      #Timing accumulator = density Timing * pguess * cumulative distribution function of E accumulators
      pdf_T <- dWald(dt[1,!u], v=v_T[1,!u], B=B_T[1,!u], A=A_T[1,!u], s=s_T[1,!u]) * p.guess
      
      for (i in 1:n_v){
        cdf_E <- (1-pWald(dt[i,!u]-d[!u], A=A_E[i,!u], v=v_E[i,!u], B=B_E[i,!u], s=s_E[i,!u]))
        #For E accums, cdf of current accum is not included
        if (i != 1) pdf_E <- pdf_E *  cdf_E
        pdf_T <- pdf_T * cdf_E
      }
      cdf_T <- 1-pWald(dt[1,!u], A=A_T[1,!u], v=v_T[1,!u], B=B_T[1,!u], s=s_T[1,!u])
      pdf_E <- pdf_E * cdf_T
      dt[1, !u] <- pdf_E + pdf_T
    }
    
    ### evidence process starts earlier than timing process ###
  } else if(t0E[[1]] < t0T[[1]]) {
    dt <- dt - t0E[[1]]
    # assume timing process starts at d=t0_T-t0_E after the evidence process
    d <- t0T[[1]] - t0E[[1]]
    u <- dt[1,] < d
    if(any(u)) {  # timing process hasn't started so must be an evidence response
      dt[1, u] <- dWald(dt[1,u], v=v_E[1,u], B=B_E[1,u], A=A_E[1,u], s=s_E[1,u]) 
      for (i in 2:n_v){
        dt[1, u] <- dt[1, u] * (1-pWald(dt[i,u], A=A_E[i,u], v=v_E[i,u], B=B_E[i,u], s=s_E[i,u]))
      }
    }
    if(any(!u)) { # could be a response from the timer or evidence process
      #Evidence Accumulator is density accumulator * cumulative distribution function other Evidence Accumulators
      pdf_E <- dWald(dt[1,!u], v=v_E[1,!u], B=B_E[1,!u], A=A_E[1,!u], s=s_E[1,!u]) 
      
      #Timing accumulator = density Timing * pguess * cumulative distribution function of E accumulators
      pdf_T <- dWald(dt[1,!u]-d[!u], v=v_T[1,!u], B=B_T[1,!u], A=A_T[1,!u], s=s_T[1,!u]) * p.guess
      
      for (i in 1:n_v){
        cdf_E <- (1-pWald(dt[i,!u], A=A_E[i,!u], v=v_E[i,!u], B=B_E[i,!u], s=s_E[i,!u]))
        #For E accums, cdf of current accum is not included
        if (i != 1) pdf_E <- pdf_E *  cdf_E
        pdf_T <- pdf_T * cdf_E
      }
      cdf_T <- 1-pWald(dt[1,!u]-d[!u], A=A_T[1,!u], v=v_T[1,!u], B=B_T[1,!u], s=s_T[1,!u])
      pdf_E <- pdf_E * cdf_T
      
      dt[1, !u] <- pdf_E + pdf_T
    }
    
    ### evidence process and timing process start at the same time ###
  } else if(t0E[[1]] == t0T[[1]]) {
    dt <- dt - t0E[[1]]
    #Evidence Accumulator is density accumulator * cumulative distribution function other Evidence Accumulators
    pdf_E <- dWald(dt[1,], v=v_E[1,], B=B_E[1,], A=A_E[1,], s=s_E[1,]) 
    #Timing accumulator = density Timing * pguess * cumulative distribution function of E accumulators
    pdf_T <- dWald(dt[1,], v=v_T[1,], B=B_T[1,], A=A_T[1,], s=s_T[1,]) * p.guess
    
    for (i in 1:n_v){
      cdf_E <- (1-pWald(dt[i,], A=A_E[i,], v=v_E[i,], B=B_E[i,], s=s_E[i,]))
      #For E accums, cdf of current accum is not included
      if (i != 1) pdf_E <- pdf_E *  cdf_E
      pdf_T <- pdf_T * cdf_E
    }
    cdf_T <- 1-pWald(dt[1,], A=A_T[1,], v=v_T[1,], B=B_T[1,], s=s_T[1,])
    pdf_E <- pdf_E * cdf_T
    dt[1, ] <- pdf_E+ pdf_T
  }
  
  dt[1,]
}

dWaldRace <- function(rt, response, A, B, t0, B_T, v_T, s_T, A_T, t0T, ..., 
                      st0 = 0, args.dist = list(), silent = FALSE)
{
  dots <- list(...)
  if (is.null(names(dots)))
    stop("... arguments need to be named.")
  if (is.data.frame(rt)) {
    response <- rt$response
    rt <- rt$rt
  }
  response <- as.numeric(response)
  nn <- length(rt)
  n_v <- max(vapply(dots, length, 0))
  if (!silent)
    message(paste("Results based on", n_v, "accumulators/drift rates."))
  if (!is.numeric(response) || max(response) > n_v)
    stop("response needs to be a numeric vector of integers up to number of accumulators.")
  if (any(response < 1))
    stop("the first response/accumulator must have value 1.")
  if (n_v < 2)
    stop("There need to be at least two accumulators/drift rates.")
  #  distribution <- match.arg(distribution)
  response <- rep(response, length.out = nn)
  A <- check_i_arguments(A, nn = nn, n_v = n_v)
  A_T <- check_i_arguments(A_T, nn = nn, n_v = n_v)
  B <- check_i_arguments(B, nn = nn, n_v = n_v)
  B_T <- check_i_arguments(B_T, nn = nn, n_v = n_v)
  t0 <- check_i_arguments(t0, nn = nn, n_v = n_v)
  t0T <- check_i_arguments(t0T, nn = nn, n_v = n_v)
  dots$v <- check_i_arguments(dots$v, nn = nn,
                              n_v = n_v, dots = TRUE)
  v_T <- check_i_arguments(v_T, nn = nn, n_v = n_v)
  dots$s <- check_i_arguments(dots$s, nn=nn, n_v=n_v, dots=TRUE)
  s_T <- check_i_arguments(s_T, nn = nn, n_v = n_v)
  for (i in seq_len(length(dots))) {
    if (length(dots[[i]]) < n_v) 
      dots[[i]] <- rep(dots[[i]], length.out = n_v)
  }
  out <- vector("numeric", nn)
  for (i in unique(response)) {
    sel <- response == i #make sure that only matrix entries belonging to this specific response are selected. 
    out[sel] <- do.call(n1WaldTiming_NS, args = c(dt = list(rt[sel]),
                                                  A = list(lapply(A, '[', i=sel)[c(i, seq_len(n_v)[-i])]),
                                                  A_T = list(lapply(A_T, '[', i=sel)[c(i, seq_len(n_v)[-i])]),
                                                  B = list(lapply(B, '[', i=sel)[c(i, seq_len(n_v)[-i])]),
                                                  B_T = list(lapply(B_T, '[', i=sel)[c(i, seq_len(n_v)[-i])]),
                                                  t0 = list(lapply(t0, '[', i=sel)[c(i, seq_len(n_v)[-i])]),
                                                  t0T = list(lapply(t0T, '[', i=sel)[c(i, seq_len(n_v)[-i])]),
                                                  s_T = list(lapply(s_T, '[', i=sel)[c(i, seq_len(n_v)[-i])]),
                                                  v_T = list(lapply(v_T, '[', i=sel)[c(i, seq_len(n_v)[-i])]),
                                                  gf=0,
                                                  lapply(dots, function(x) lapply(x, "[", i = sel)[c(i, 
                                                                                                     seq_len(n_v)[-i])])))
  }
  return(out)
}



