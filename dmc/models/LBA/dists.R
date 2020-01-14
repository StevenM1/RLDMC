## LBA
require(rtdists) 

# The following functions are from the rtdists package but with some tweaks 
# to vectorize t0.

## Simulate LBA trials ----

rlba.norm <- function (n,A,b,t0,mean_v,sd_v,st0=0,posdrift=TRUE,return.ttf=FALSE) 
{
  if (any(b < A)) 
    stop("b cannot be smaller than A!")
  n_v  <- ifelse(is.null(dim(mean_v)), length(mean_v), dim(mean_v)[1])
  if (posdrift) 
    drifts <- matrix(rtnorm(n = n * n_v, mean = mean_v, sd = sd_v, 
                            lower = 0), nrow = n_v)
  else drifts <- matrix(rnorm(n = n * n_v, mean = mean_v, sd = sd_v), 
                        nrow = n_v)
  make.r(drifts = drifts, b = b, A = A, n_v = n_v, t0 = t0, st0 = st0, n = n,
         return.ttf=return.ttf)
}


make.r <- function (drifts, b, A, n_v, t0, st0 = 0, n, return.ttf=FALSE) 
{
  drifts[drifts < 0] <- 0
  starts <- matrix(runif(min = 0, max = A, n = n * n_v), nrow = n_v)
  ttf <- t0 + (b - starts)/drifts
  if (return.ttf) return(ttf)
  resp <- apply(ttf, 2, which.min)
  rt <- ttf[cbind(resp,1:n)]
  if (st0[1]>0) rt <- rt + runif(min = 0, max = st0[1], n = n)
  bad <- !is.finite(rt)
  if (any(bad)) {
    warning(paste(sum(bad), "infinite RTs removed and less than", 
                  n, "rts returned"))
    resp <- resp[!bad]
    rt <- rt[!bad]
  }
  data.frame(rt = rt, response = resp)
}

### LBA likelihood ----

dlba.norm.core <- function (t,A,b,mean_v,sd_v,posdrift=TRUE,robust=FALSE,nn) 
# like dlba_norm_core but t0 dealt with outside (removed from t), and
# default added to nn
{
  
  pnormP <- function (x, mean = 0, sd = 1, lower.tail = TRUE) 
    ifelse(abs(x) < 7, pnorm(x, mean = mean, sd = sd, lower.tail = lower.tail), 
      ifelse(x < 0, 0, 1))
  
  dnormP <- function (x, mean = 0, sd = 1) 
    ifelse(abs(x) < 7, dnorm(x, mean = mean, sd = sd), 0)

    if (robust) {
        pnorm1 <- pnormP
        dnorm1 <- dnormP
    }
    else {
        pnorm1 <- pnorm
        dnorm1 <- dnorm
    }
    if (posdrift) 
        denom <- pmax(pnorm1(mean_v/sd_v), 1e-10)
    else denom <- rep(1, nn)
    if (any(A < 1e-10)) {
        A_small <- A < 1e-10
        out_A <- pmax(0, ((b[A_small]/t[A_small]^2) * 
                            dnorm1(b[A_small]/t[A_small], 
            mean_v[A_small], sd = sd_v[A_small]))/denom[A_small])
        zs <- t[!A_small] * sd_v[!A_small]
        zu <- t[!A_small] * mean_v[!A_small]
        chiminuszu <- b[!A_small] - zu
        chizu <- chiminuszu/zs
        chizumax <- (chiminuszu - A[!A_small])/zs
        out_o <- pmax(0, (mean_v[!A_small] * (pnorm1(chizu) - 
            pnorm1(chizumax)) + sd_v[!A_small] * (dnorm1(chizumax) - 
            dnorm1(chizu)))/(A[!A_small] * denom[!A_small]))
        out <- numeric(nn)
        out[!A_small] <- out_o
        out[A_small] <- out_A
        return(out)
    }
    else {
        zs <- t * sd_v
        zu <- t * mean_v
        chiminuszu <- b - zu
        chizu <- chiminuszu/zs
        chizumax <- (chiminuszu - A)/zs
        return(pmax(0, (mean_v * (pnorm1(chizu) - pnorm1(chizumax)) + 
            sd_v * (dnorm1(chizumax) - dnorm1(chizu)))/(A * denom)))
    }
}


plba.norm.core <- function (t,A,b,mean_v,sd_v,posdrift=TRUE,robust=FALSE,nn) 
# like plba_norm_core but t0 dealt with outside (removed from t), and
# default added to nn
{
  
    pnormP <- function (x, mean = 0, sd = 1, lower.tail = TRUE) 
      ifelse(abs(x) < 7, pnorm(x, mean = mean, sd = sd, lower.tail=lower.tail), 
        ifelse(x < 0, 0, 1))
  
    dnormP <- function (x, mean = 0, sd = 1) 
      ifelse(abs(x) < 7, dnorm(x, mean = mean, sd = sd), 0)

    if (robust) {
        pnorm1 <- pnormP
        dnorm1 <- dnormP
    }
    else {
        pnorm1 <- pnorm
        dnorm1 <- dnorm
    }
    if (posdrift) 
        denom <- pmax(pnorm1(mean_v/sd_v), 1e-10)
    else denom <- 1
    if (any(A < 1e-10)) {
        A_small <- A < 1e-10
        out_A <- pmin(1, pmax(0, (pnorm1(b[A_small]/t[A_small], 
            mean = mean_v[A_small], sd = sd_v[A_small], 
            lower.tail = FALSE))/denom[A_small]))
        zs <- t[!A_small] * sd_v[!A_small]
        zu <- t[!A_small] * mean_v[!A_small]
        chiminuszu <- b[!A_small] - zu
        xx <- chiminuszu - A[!A_small]
        chizu <- chiminuszu/zs
        chizumax <- xx/zs
        tmp1 <- zs * (dnorm1(chizumax) - dnorm1(chizu))
        tmp2 <- xx * pnorm1(chizumax) - chiminuszu * pnorm1(chizu)
        out_o <- pmin(pmax(0, (1 + (tmp1 + tmp2)/A[!A_small])/denom[!A_small]), 
            1)
        out <- numeric(length(mean_v))
        out[!A_small] <- out_o
        out[A_small] <- out_A
        return(out)
    }
    else {
        zs <- t * sd_v
        zu <- t * mean_v
        chiminuszu <- b - zu
        xx <- chiminuszu - A
        chizu <- chiminuszu/zs
        chizumax <- xx/zs
        tmp1 <- zs * (dnorm1(chizumax) - dnorm1(chizu))
        tmp2 <- xx * pnorm1(chizumax) - chiminuszu * pnorm1(chizu)
        return(pmin(pmax(0, (1 + (tmp1 + tmp2)/A)/denom), 1))
    }
}


dlba.norm <- function (t, A, b, mean_v, sd_v, posdrift = TRUE, robust = FALSE) 
{
    nn <- length(t)
    tpos <- (t>0) & (b >= A)
    out <- numeric(length(t))
    out[tpos] <- dlba.norm.core(t = t[tpos], A = A[tpos], b = b[tpos], 
      mean_v = mean_v[tpos], sd_v = sd_v[tpos], 
      posdrift = posdrift, robust = robust, nn = nn)
    out
}


plba.norm <- function (t, A, b, mean_v, sd_v, posdrift = TRUE, robust = FALSE) 
{
    nn <- length(t)
    tpos <- (t>0) & (b >= A)
    out <- numeric(length(t))
    out[tpos] <- plba.norm.core(t = t[tpos], A = A[tpos], b = b[tpos], 
      mean_v = mean_v[tpos], sd_v = sd_v[tpos], 
      posdrift = posdrift, robust = robust, nn = nn)
    out
}


n1PDFfixedt0.norm=function(dt,A,b,mean_v,sd_v, 
                           posdrift=TRUE,robust = FALSE) 
# Generates defective PDF for responses on node= 1
# dt (decison time) is a matrix with length(mean_v) rows, one row for
# each accumulator to allow for different start times
{
  
  n_acc <- dim(mean_v)[2]
  dt <- matrix(rep(dt,times=n_acc),ncol=n_acc)

  dt[,1] <- dlba.norm(dt[,1],A=A[,1],b=b[,1],mean_v=mean_v[,1],sd_v=sd_v[,1],
                      posdrift=posdrift,robust=robust)
  if (n_acc>1) for (i in 2:n_acc)
    dt[,1] <- dt[,1]*(1-plba.norm(dt[,i],
      A=A[,i],b=b[,i],mean_v=mean_v[,i],sd_v=sd_v[,i],
      posdrift=posdrift,robust=robust))
  dt[,1]
}

n1PDFfixedt0.norm.t=function(dt,A,b,mean_v,sd_v, 
                           posdrift=TRUE,robust = FALSE) 
# Same as n1PDFfixedt0.norm but with parameter matrices transposed
{
  
  n_acc <- dim(mean_v)[1]
  dt <- matrix(rep(dt,times=n_acc),ncol=n_acc)

  dt[,1] <- dlba.norm(dt[1,],A=A[1,],b=b[1,],mean_v=mean_v[1,],sd_v=sd_v[1,],
                      posdrift=posdrift,robust=robust)
  if (n_acc>1) for (i in 2:n_acc)
    dt[,1] <- dt[,1]*(1-plba.norm(dt[,i],
      A=A[i,],b=b[i,],mean_v=mean_v[i,],sd_v=sd_v[i,],
      posdrift=posdrift,robust=robust))
  dt[,1]
}

### LBA returning non-responses when both accumulators have rate < 0 ----

make.r.vGF <- function (drifts, b, A, n_v, t0, censor = Inf, st0 = 0, n, return.ttf = FALSE, gf=0) 
{
    nr <- apply(drifts,2,function(x){all(x<0)}) # both rates negative
    resp.nr <- apply(drifts[,nr,drop=FALSE],2,which.max)
    drifts[drifts < 0] <- 0
    starts <- matrix(runif(min = 0, max = A, n = n * n_v), nrow = n_v)
    ttf <- t0 + (b - starts)/drifts
    if (return.ttf) return(ttf)
    resp <- apply(ttf, 2, which.min)
    if (any(nr)) resp[nr] <- resp.nr
    rt <- ttf[cbind(resp, 1:n)]
    if (st0[1] > 0) 
        rt <- rt + runif(min = 0, max = st0[1], n = n)
    
    rt[!is.finite(rt) | (rt > censor) ] <- NA
    # resp[is.na(rt)] <- NA
    out <- data.frame(RT = rt, R = resp)
    
    if (gf > 0) {
      is.gf <- as.logical(rbinom(dim(out)[1],1,gf))
      out$RT[is.gf] <- NA
      out$R[is.gf] <- 1
    }
    
    out
}


rlba.norm.vGF <- function (n, A, b, t0, mean_v, sd_v, censor=Inf, return.ttf = FALSE, gf=0) 
{
    if (any(b < A)) 
        stop("b cannot be smaller than A!")
    n_v <- ifelse(is.null(dim(mean_v)), length(mean_v), dim(mean_v)[1])
    drifts <- matrix(rnorm(n = n * n_v, mean = mean_v, sd = sd_v), 
        nrow = n_v)
    make.r.vGF(drifts = drifts, b = b, A = A, n_v = n_v, t0 = t0, censor = censor, 
        st0 = 0, n = n, return.ttf = return.ttf, gf=gf)
}



n1PDFvGF <- function (rt, A, b, t0, mean_v, sd_v, censor=Inf, st0 = 0, silent = FALSE, gf=0) {
  
  
  n1PDFfixedt0=function(dt,A,b,mean_v,sd_v,posdrift=TRUE,robust = FALSE) 
  {
  
    n_acc <- length(mean_v)
    n <- length(dt)
    dt <- matrix(rep(dt,times=n_acc),ncol=n_acc)
    A <- matrix(rep(A,n),nrow=n_acc)
    b <- matrix(rep(b,n),nrow=n_acc)
    mean_v <- matrix(rep(mean_v,n),nrow=n_acc)
    sd_v <- matrix(rep(sd_v,n),nrow=n_acc)

    dt[,1] <- dlba.norm(dt[,1],A=A[1,],b=b[1,],mean_v=mean_v[1,],sd_v=sd_v[1,],
                      posdrift=posdrift,robust=robust)
    if (n_acc>1) for (i in 2:n_acc)
      dt[,1] <- dt[,1]*(1-plba.norm(dt[,i],
        A=A[i,],b=b[i,],mean_v=mean_v[i,],sd_v=sd_v[i,],
        posdrift=posdrift,robust=robust))
    dt[,1]
  }

  pN <- prod(pnorm(0,mean_v,sd_v))
  if ( is.infinite(censor) ) {
      pRespond <- try(integrate(n1PDFfixedt0,0,Inf,A=A,b=b,mean_v=mean_v,sd_v=sd_v,posdrift=FALSE)$value,silent=TRUE)
      pCensor <- 0
      bad <- pRespond=="try-error"
  } else {
      pRespond <- try(integrate(n1PDFfixedt0,0,censor-t0,A=A,b=b,mean_v=mean_v,sd_v=sd_v,posdrift=FALSE)$value,silent=TRUE)
      pCensor <- try(integrate(n1PDFfixedt0,censor-t0,Inf,A=A,b=b,mean_v=mean_v,sd_v=sd_v,posdrift=FALSE)$value,silent=TRUE)
      bad <- any(c(class(pRespond)=="try-error",class(pCensor)=="try-error"))
  }
  if (bad) return(rep(0,length(rt)))
  if (is.na(pN) || pN==1) p_omission <- pCensor else p_omission <- pCensor + pN*(pRespond+pCensor)/(1-pN)
  out <- numeric(length(rt))
  rt <- rt-t0
  nart <- is.na(rt)
  is.rt <- rt >= 0; is.rt[nart] <- FALSE 
  out[nart] <- rep(p_omission,sum(nart)) 
  out[is.rt] <- n1PDFfixedt0(dt=rt[is.rt],
          A=A,
          b=b,
          mean_v=mean_v,
          sd_v=sd_v,
          posdrift=FALSE)
  
  if (gf>0) {
    out[is.rt] <- (1-gf)*out[is.rt]
    out[nart] <- gf + (1-gf)*out[nart]
  }
  
  out
} 
  
# # Check go failure
# n=1e5
# v=c(2,1); sd_v = c(1,1); B=c(1,1); A=c(2,2);t0=1; gf=.2; censor=2
# b=A+B
# sim <- rlba.norm.vGF(n=n,A=A,b=b,mean_v=v,sd_v=sd_v,t0=t0,gf=gf,censor=censor)
# par(mfrow=c(1,3))
# dns <- plot.cell.density(sim,C=1,xlim=c(0,4),save.density=TRUE)
# # # DOESNT WORK
# # pc <- integrate(n1PDFvGF,lower=0,upper=Inf,A=A,b=b,mean_v=v,sd_v=sd_v,t0=t0,gf=gf)$value/(1-gf)
# # print(pc)
# dt=dns$correct$x
# d <- n1PDFvGF(dt,A=A,b=b,mean_v=v,sd_v=sd_v,gf=gf,t0=t0,censor=censor)
# # red=black?
# plot(dns$correct$x,dns$correct$y,type="l")
# lines(dns$correct$x,d,col="red")
# dt=dns$error$x
# d <- n1PDFvGF(dt,A=A[2:1],b=b[2:1],mean_v=v[2:1],sd_v=sd_v[2:1],gf=gf,t0=t0,censor=censor)
# plot(dns$error$x,dns$error$y,lty=2,type="l")
# lines(dns$error$x,d,col="red",lty=2)

