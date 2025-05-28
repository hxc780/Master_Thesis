# ------ PEAKS AND VARYING AMPLITUDE CORRECTION ------ #
require(truncnorm)

# estimate peak magnitude
estimatePeakSize = function(Y.arg, Valleys.arg, Peaks.arg){
  K = length(Peaks.arg)
  Y.size = length(Y.arg)
  MagnitudeVector <- rep(0,K)
  
  for(k in 1:K){ # for each top
    sign_vec = sign(Valleys.arg-Peaks.arg[k]) # below = -1 or above = +1
    below <- sum(sign_vec<0) # valley below
    above <- below + 1 # valley above
    
    if(below==0){# if no valley below
      Obs_below <- 1
      Obs_Above <- ceiling(Valleys.arg[above])
    }else if(below==length(Valleys.arg)){ # if all valleys below
      Obs_below <- ceiling(Valleys.arg[below])
      Obs_Above <- Y.size
    }else{ # one valley below, one valley above
      Obs_Above <- ceiling(Valleys.arg[above])
      Obs_below <- ceiling(Valleys.arg[below])
    }
    MagnitudeVector[k] <- median(sort(Y.arg[Obs_below:Obs_Above], decreasing=TRUE)[1:5])
  }
  
  return(list(idx=Peaks.arg, size=MagnitudeVector))
}


EstimateBx <-  function(Y.arg, xi.arg, A.arg, B0.arg, b.arg, delta.arg=1){
  
  # get valleys
  valleys.c <- sort(c(FindXstar(g(xi.arg, delta.arg), b.arg, 0),FindXstar(g(xi.arg, delta.arg), b.arg, pi)))
  
  # get peaks
  peaks.c <- sort(c(FindXstar(g(xi.arg, delta.arg),b.arg, pi/2),FindXstar(g(xi.arg, delta.arg), b.arg, 3*pi/2)))
  
  if(length(peaks.c) > 1){ # if we have at least 2 peaks ==> assume varying amplitude
    
    # estimate peak size
    yfit_star <- estimatePeakSize(Y.arg, valleys.c, peaks.c)$size
    
    n.obs <- length(xi.arg)
    # Ensure non-intersecting envelopes, i.e. that B(x) > -A/2 ==> r(x) > -A/2 - B0
    # Also ensure that slow primary wave dominate, i.e. that B(x) < A ==> r(x) < A - B0
    rk.c <- pmin(pmax(yfit_star - f(xi.arg, rep(B0.arg, n.obs), delta.arg,  A.arg, b.arg)[peaks.c], -A.arg/2-B0.arg), A.arg-B0.arg)
    
    Bx.c <- get_Bx(1:length(Y.arg), peaks.c, B0.arg + rk.c)
    
    
  }else{ # if we have only 1 or less peak ==> assume constant amplitude
    
    Bx.c <- rep(B0.arg, length(Y.arg))
    rk.c <- 0
    
  }
  
  return(list(rk.c=rk.c, Bx.c=Bx.c))

}

# compute varying amplitude correction B(x) = B0 + h(x)
get_Bx <- function(x.arg, xstar.id.arg,Bk.arg){
  
  N = length(x.arg)
  K = length(xstar.id.arg)
  Bx = rep(0,N)
  
  # if only 1 peak:
  if(K<=1){
    return(Bx)
  }
  
  # handle boundaries
  Bx[1:xstar.id.arg[1]] <- Bk.arg[1]
  Bx[xstar.id.arg[K-1]:N] <- Bk.arg[K]
  
  # interpolationg between peaks
  slope = (Bk.arg[-K] - Bk.arg[-1])/(xstar.id.arg[-K] - xstar.id.arg[-1]) # Slope of linear interpolation
  intercept = (Bk.arg[-K] * xstar.id.arg[-1] - Bk.arg[-1] * xstar.id.arg[-K])/(xstar.id.arg[-1] - xstar.id.arg[-K]) # Intercept of linear interpolation
  
  for(k in 1:(K-1)){
    Bx[xstar.id.arg[k]:xstar.id.arg[k+1]] <- intercept[k] + slope[k]*x.arg[xstar.id.arg[k]:xstar.id.arg[k+1]]
  }
  
  return(Bx)
  
}


# smooth signal and identify potential valleys
findvalleys <- function(x, y, w=1, ...) {
  require(zoo)
  n <- length(y)
  y.smooth <- loess(y ~ x, ...)$fitted
  y.min <- rollapply(zoo(y.smooth), 2*w+1, min, 
                     align="center")
  delta <- y.min - y.smooth[-c(1:w, n+1-1:w)]
  i.min <- which(delta >= 0) + w
  list(x=x[i.min], i=i.min, y.hat=y.smooth)
}

# Find all recurrent points, fulfiling
# g(x) + b = n*2*pi+offset
FindXstar <- function(g,b,offset=pi/2){
  NoMorePeaks <- FALSE
  Xstar <- numeric()
  n <- 0
  i <- 1
  while(!NoMorePeaks){
    
    NewPeakStart <- sum(g + b - offset - n*2*pi < 0)
    
    if(NewPeakStart==length(g)){
      NoMorePeaks <- TRUE
    }else if(NewPeakStart > 0){
      Xstar[i] <- NewPeakStart
      i <- i + 1
    }
    n <- n + 1
  }
  return(Xstar)
}

complete.log.lik <- function(Y.obs,  xi.arg, delta.arg = delta, Bx.arg, Bk.arg, sigma_eps.arg, sigma_B.arg, psi.arg, 
                             omega2dBeta.arg, a.arg, A0.arg, B0.arg, b.arg){
  #  --- w/ CIR-exact transition density --- #
  c = 2/(1-psi.arg)/omega2dBeta.arg
  lambda = 2*c*xi.arg[1:(length(xi.arg) - 1)]*psi.arg
  return(sum(dnorm(Y.obs, f(xi.arg, Bx.arg , delta.arg, A0.arg = A0.arg, b.arg = b.arg), sigma_eps.arg, log = TRUE)) +
           sum(dchisq( c * xi.arg[2:length(xi.arg)], df = 4 * a.arg  / omega2dBeta.arg, ncp = lambda, log = TRUE))+
           sum(dnorm(Bk.arg,B0.arg,sigma_B.arg,log=TRUE)))
}

# ------ END ------- #

# ----- MARTINGALE ESTIAMTING FUNCTIONS ----- #
MSF.a <- function(xi.arg, psi.arg){
  bias.cor <- 1
  n.xi = length(xi.arg)
  term <- mean(xi.arg[2:n.xi]) + psi.arg/((n.xi-1)*(1-psi.arg))*(xi.arg[n.xi]-xi.arg[1])
  return(bias.cor*term)
}

MSF.psi <- function(xi.arg){
  n.xi = length(xi.arg)
  num = (n.xi-1) * sum(xi.arg[2:n.xi]/xi.arg[1:(n.xi-1)])-sum(xi.arg[2:n.xi])*sum(xi.arg[1:(n.xi-1)]^(-1))
  den = (n.xi-1)^2 - sum(xi.arg[1:(n.xi-1)])*sum(xi.arg[1:(n.xi-1)]^(-1))
  return(num/den)
}

MSF.omega2 <- function(xi.arg, psi.arg, a.arg, delta.arg){
  n.xi = length(xi.arg)
  num = sum(xi.arg[1:(n.xi-1)]^(-1)*(xi.arg[2:n.xi]-xi.arg[1:(n.xi-1)]*psi.arg-a.arg*(1-psi.arg))^2)
  den = sum(xi.arg[1:(n.xi-1)]^(-1)*((1/2*a.arg-xi.arg[1:(n.xi-1)])*psi.arg^2-(a.arg-xi.arg[1:(n.xi-1)])*psi.arg +1/2*a.arg) )
  return(-1/delta.arg*log(psi.arg)*num/den)
}

# ------ END ------- #

# ----- SIMULATION AUXILIARY FUNCTIONS ----- #

rxi.exact <- function(a.arg=a,beta.arg=beta,omega.arg=omega,n.arg=N,delta.arg=delta){
  omega2dBeta.arg = omega.arg^2/beta.arg
  psi.arg = exp(-delta.arg*beta.arg)
  c = 2/(1-psi.arg)/omega2dBeta.arg
  xi.arg = rep(a.arg, n.arg)
  for(i in 2:length(xi.arg)){
    xi.arg[i]  <- 1/(2*c)*rchisq(1, df = 4 * a.arg  / omega2dBeta.arg, ncp = 2*c*xi.arg[i-1]*psi.arg)
  }
  return(xi.arg)
}

transfer.xi.noncentralchisq <- function(xi.arg, psi.arg = psi, a.arg=a, omega2dBeta.arg = omega2dBeta){
  c = 2/(1-psi.arg)/omega2dBeta.arg
  lambda = 2*c*xi.arg*psi.arg
  #if(any(lambda >= 1e5)) print("Unstable Noncentral Chisquare")
  return(1/(2*c)*rchisq(length(xi.arg), df = 4 * a.arg  / omega2dBeta.arg, ncp = lambda))
}

g <- function(xi.arg, delta.arg = delta) {
  return(cumsum(xi.arg)*delta.arg)
}

f <- function(xi.arg, Bx.arg=Bx, delta.arg = delta, A.arg = A, b.arg = b) {
  A.arg * sin(g(xi.arg, delta.arg) + b.arg) -
    Bx.arg * cos(2 * (g(xi.arg, delta.arg) + b.arg))
}

# ------ END ------- #


# ---- SIMULATE SIGNAL ---- #

SimSignal <- function(a,b,A,B0,beta,omega,sigma_B,sigma, N, delta){
  require(truncnorm)
  
  psi <- exp(-delta*beta)
  x = seq(delta,N*delta,delta) 
  
  xi <- rxi.exact(a.arg=a,beta.arg=beta,omega.arg=omega,n.arg=N,delta.arg=delta)
  gxi <- g(xi,delta)
  
  hightops <- FindXstar(gxi,b)
  lowtops <- FindXstar(gxi,b,3*pi/2)
  
  xkstar.idx  <- sort(c(hightops,lowtops))
  xkstar <- x[xkstar.idx]
  
  K <- length(xkstar)
  
  # ensure that B(x) < A, otherwise the secondary wave will dominate
  # also ensure that upper env = A + Bx > -Bx = lower env ==> Bx > -A/2
  Bk <- rtruncnorm(K,mean=B0, a=-A/2, b=A, sd=sigma_B) 
  Bx <- get_Bx(1:N, xkstar.idx, Bk)
  
  noise <- rnorm(N, 0, sigma)
  Y0 <- f(xi.arg=xi,Bx.arg=Bx,delta.arg=delta, A.arg=A,b.arg=b)
  Y <- f(xi.arg=xi,Bx.arg=Bx,delta.arg=delta, A.arg=A,b.arg=b) + noise
  
  return(list(x=x, Y=Y, Y0=Y0, xi=xi,gxi=gxi, Bx=Bx, Bk=Bk, ExpGLG=gxi[N]/(2*pi), xkstar=xkstar, xkstar.idx = xkstar.idx))
  
}
# ------ END ------- #
