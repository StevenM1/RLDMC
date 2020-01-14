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

dWald <- function(t,v,B,A,s=1,useSuppDists=TRUE)
  # density for single accumulator
{

  digt <- function(t,k=1,l=1,a=.1,s=1,tiny=1e-10) {
    # pdf of inverse gaussian at t with k +/- a/2 uniform variability
    # returns digt.0 if a<1e-10

    digt.0 <- function(t,k=1,l=1) {
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

    tpos <- t<=0

    atiny <- a<=tiny & !tpos
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
        x[atiny][nu.inf] <- 0
        x[atiny][!nu.inf] <- dinvGauss(t[atiny][!nu.inf], 
                                       nu=nu[!nu.inf], 
                                       lambda=lambda[!nu.inf])
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

rWaldRaceSM <- function(n, A, B, t0, v, s, st0 = 0, silent = FALSE, args.dist=NULL, simPerTrial=FALSE)
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
      for(i in 1:n_v) {
        ttf[i,] <- rinvGauss(n, nu=bs[,i]/v[,i], lambda=(bs[,i]/s[,i])^2)
      }
      ttf <- ttf + t(t0)
      resp <- apply(ttf, 2, which.min)
      out <- data.frame(RT = ttf[cbind(resp,1:n)], R = apply(ttf, 2, which.min))
    }

  } else {
    # all equal drifts, nice and fast
    out <- rWaldRace(n=n, v=v, B=B, A=A, t0=t0, s=s)
  }
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

n1Wald <- function(dt,B,A,t0,gf=0, ...)
  # Generates defective PDF for responses on node=1, dt (decison time) is a vector of times
{
  ### Added by SM
  dots <- list(...)
  n_v <- max(vapply(dots, length, 0))  # Number of responses
  #  if(!silent) message(paste("Results based on", n_v, "accumulators/drift rates."))
  if (n_v < 2) stop("There need to be at least two accumulators/drift rates.")
  #  distribution <- match.arg(distribution)
  #check_single_arg(t0 = t0)
  nn <- length(dt)
  B <- check_n1_arguments(B, nn=nn, n_v=n_v)
  A <- check_n1_arguments(A, nn=nn, n_v=n_v)
  t0 <- check_n1_arguments(t0, nn=nn, n_v=n_v)
  v <- check_n1_arguments(dots[[1]], nn=nn, n_v=n_v)
  s <- check_n1_arguments(dots[[2]], nn=nn, n_v=n_v)
  
  #  B[B<0] <- 0 # Protection for negatives
  #  A[A<0] <- 0
  # n_acc <- ifelse(is.null(dim(v)),length(v),dim(v)[1])
  # if (is.null(dim(dt))) dt <- matrix(rep(dt,each=n_acc),nrow=n_acc)
  if (is.null(dim(dt))) dt <- matrix(rep(dt, each=n_v), nrow=n_v)
  dt <- dt-t0[[1]]  ## assume single t0

  is.go <- !is.na(dt[1,])
  n.go <- sum(is.go)

  #  if (!is.matrix(v)) v <- matrix(rep(v,n.go),nrow=n_v)
  if (!is.matrix(v)) v <- do.call(rbind, v) #matrix(unlist(v), nrow=n_v)
  if (!is.matrix(A)) A <- do.call(rbind, A) # matrix(unlist(A), nrow=n_v)
  if (!is.matrix(B)) B <- do.call(rbind, B) #matrix(unlist(B), nrow=n_v)
  if (!is.matrix(s)) s <- do.call(rbind, s)
  #  if (!is.matrix(B)) B <- matrix(rep(B,n.go),nrow=n_v)
  #  if (!is.matrix(A)) A <- matrix(rep(A,n.go),nrow=n_v)

  # Winner
  dt[1,is.go] <- (1-gf[1])*dWald(dt[1,is.go],A=A[1,],v=v[1,],B=B[1,],s=s[1,])
  if (n_v > 1) for (i in 2:n_v)
    dt[1,is.go] <- dt[1,is.go]*(1-pWald(dt[i,is.go],A=A[i,],v=v[i,],B=B[i,],s=s[i,]))

  dt[1,!is.go] <- gf[1]

  dt[1,]
}

dWaldRace <- function(rt, response, A, B, t0, ..., st0 = 0, args.dist = list(), silent = FALSE)
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
  B <- check_i_arguments(B, nn = nn, n_v = n_v)
  t0 <- check_i_arguments(t0, nn = nn, n_v = n_v)
  dots$v <- check_i_arguments(dots$v, nn = nn,
                               n_v = n_v, dots = TRUE)
  dots$s <- check_i_arguments(dots$s, nn=nn, n_v=n_v, dots=TRUE)

  for (i in seq_len(length(dots))) {
    if (length(dots[[i]]) < n_v) 
      dots[[i]] <- rep(dots[[i]], length.out = n_v)
  }
  out <- vector("numeric", nn)
  for (i in unique(response)) {
    sel <- response == i
    out[sel] <- do.call(n1Wald, args = c(dt = list(rt[sel]),
                                         A = list(lapply(A, '[', i=sel)[c(i, seq_len(n_v)[-i])]),
                                         B = list(lapply(B, '[', i=sel)[c(i, seq_len(n_v)[-i])]),
                                         t0 = list(lapply(t0, '[', i=sel)[c(i, seq_len(n_v)[-i])]),
                                         gf=0,
                                         lapply(dots, function(x) lapply(x, "[", i = sel)[c(i, 
                                                                                            seq_len(n_v)[-i])])))
  }
  return(out)
}