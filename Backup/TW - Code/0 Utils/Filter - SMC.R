#########################################
#########  SMC helper functions #########
#########################################

###############################################
#### --- SMC ENSEMBLE (INITIALIZATION) --- ####
###############################################
smc.init.xi.ens <- function(n.time, n.part) {
  return(list(xi.part = matrix(nrow = n.part, ncol = n.time),
              w.part = matrix(nrow = n.part, ncol = n.time),
              xi.c = rep(0, n.time),
              gxi.c = rep(0,n.time),
              psi.c = 0,
              omega2dBeta.c = 0,
              a.c = 0))
}


smc.xi.ensemble <- function (Y.obs, n.part, total.n.part.min, sigma.arg = sigma, A.arg = A0, Bx.arg = Bx, a.arg = a, b.arg = b, delta.arg=delta, FixOmega=FALSE, track=FALSE) {
  
  n.time <- length(Y.obs)
  Rx.input <- rep(0,n.time)
  choptime <- 50 # every 50th iteration
  
  # initialize matrix of parameters: beta, psi, omega2dBeta giving rise to "n.ens" ensembles
  beta.cand <- -1/delta.arg*log(seq(0.35,0.95,0.1))
  #a.cand <- c(a.arg/2, a.arg, 2*a.arg)
  a.cand <- a.arg
  if(FixOmega){ # assume omega^2/beta  = a
    pars <- expand.grid(beta=beta.cand,a=a.cand)
    pars$omega2dBeta = pars$a
  }else{ # or free omega^2
    pars <- expand.grid(beta = beta.cand, a=a.cand, lambda=seq(0.2,0.8,0.2))
    omega2.func <- function(beta.arg, lambda.arg) lambda.arg*2*a.arg*beta.arg # let omega = lambda*sqrt(2ab) < sqrt(2ab) for number 0 < lambda < 1
    pars$omega2dBeta = omega2.func(pars$beta, pars$lambda)/pars$beta
  }
  pars$psi <- exp(-delta.arg*pars$beta)
  n.ens <- nrow(pars)
  pars_ext <- pars[rep(row.names(pars), each = n.part), ]
  row.names(pars_ext) <- NULL
  pars_ext$ens <- rep(1:n.ens, each=n.part)
  
  # avoid implausible (evaluated to -Inf) candidates
  lower_trim <- 30
  f.smc <- function(gxi.arg, Bx.arg=Bx, delta.arg = delta, A.arg = A0, b.arg = b) {
    A.arg * sin(gxi.arg + b.arg) -
      Bx.arg * cos(2 * (gxi.arg + b.arg))
  }
  
  # initialize particles
  total.n.part <- n.ens*n.part
  
  smc.obj <- smc.init.xi.ens(n.time, total.n.part)
  sample.size = n.part
  
  smc.obj$xi.part[, 1] <- rep(a.arg, total.n.part)
  smc.obj$w.part[, 1] <- rep(1 / n.part, total.n.part)
  
  # particles at each iteration (halving at each i%%choptime iteration)
  particles_at_it <- function(N, n.max=total.n.part, n.min=total.n.part.min) {
    vector <- numeric(N)  # Initialize an empty numeric vector of length N
    
    for (i in 1:N) {
      if (i %% choptime == 1 && i != 1 && n.max > n.min) {  # Check if iteration is a multiple of choptime and not the first iteration
        n.max <- round(n.max / 2)  # Halve n.max
      }
      if (n.max >= n.min) {
        vector[i] <- n.max  # Assign n.max to the current index if it's greater than or equal to n.min
      } else {
        vector[i] <- n.min  # Otherwise, assign n.min
      }
    }
    
    return(vector)
  }
  
  # Auxilliary function
  which_idx_min <- function(v, k) {
    sorted_indices <- order(v)  # Get the indices that would sort v in ascending order
    smallest_indices <- sorted_indices[1:k]  # Select the first k indices
    
    return(smallest_indices)
  }
  
  particles <- particles_at_it(n.time)
  particles_to_remove <- particles[1:(n.time-1)]-particles[2:n.time]
  
  for (i in 2:n.time) {
    
    if(track) print(i)
    
    xi.part <- transfer.xi.noncentralchisq(smc.obj$xi.part[,i-1], pars_ext$psi, pars_ext$a, pars_ext$omega2dBeta)
    
    # Save particle draws
    smc.obj$xi.part[, i] <- xi.part
    
    # Compute weight and resample particles with high likelihood paths
    smc.obj$w.part[, i] <- lower_trim + 1e-1 + pmax(dnorm(rep(Y.obs[i],n.part), f.smc(rowSums(smc.obj$xi.part[, 1:i])*delta.arg, Bx.arg[i], delta.arg, A.arg = A.arg, b.arg = b.arg), sigma.arg,log=TRUE),-lower_trim) 
    
    # Resample
    likely.particles <- sample(total.n.part, size=total.n.part, replace=T, prob=smc.obj$w.part[,i])
    smc.obj$xi.part[, 1:i] <- smc.obj$xi.part[likely.particles, 1:i]
    smc.obj$w.part[, 1:i] <- smc.obj$w.part[likely.particles, 1:i]
    pars_ext <- pars_ext[likely.particles,]
    rownames(pars_ext) <- NULL
    
    if(i < n.time && particles_to_remove[i]!=0){
      idx_to_remove <- which_idx_min(smc.obj$w.part[,i], particles_to_remove[i]) #  get small weight particle and remove
      smc.obj$xi.part <- smc.obj$xi.part[-idx_to_remove, ]
      smc.obj$w.part <- smc.obj$w.part[-idx_to_remove, ]
      pars_ext <- pars_ext[-idx_to_remove,]
      total.n.part <- total.n.part - length(idx_to_remove) # update total number of particles
    }
    
  }
  
  # # trajectory sampling
  #draw <- sample(total.n.part,1,prob=smc.obj$w.part[,n.time])
  draw <- which.max(smc.obj$w.part[, n.time])
  smc.obj$xi.c <- smc.obj$xi.part[draw, ]
  smc.obj$gxi.c <- g(smc.obj$xi.c,delta.arg)
  smc.obj$psi.c <- pars_ext[draw,]$psi
  smc.obj$omega2dBeta.c <- pars_ext[draw,]$omega2dBeta
  smc.obj$a.c <- pars_ext[draw,]$a
  smc.obj$pargrid <- unique(pars_ext)
  
  return(smc.obj)
}


#############################################
### ----     SMC - Estimation of Xi ---- ####
#############################################
smc.init.xi <- function(n.time, n.part) {
  return(list(xi.part = matrix(nrow = n.part, ncol = n.time),
              w.part = matrix(nrow = n.part, ncol = n.time),
              xi.c = rep(0, n.time),
              gxi.c = rep(0,n.time)))
}


smc.xi <- function (Y.obs, n.part,  xkstar.arg, Bx.arg, sigma.arg = sigma, psi.arg = psi, 
                    omega2dBeta.arg = omega2dBeta, A.arg = A0, a.arg = a, b.arg = b, delta.arg=delta, track=FALSE) {
  
  n.time <- length(Y.obs)
  
  if(length(xkstar.arg)==1){
    Bx.input <- rep(mean(Bx.arg), n.time)
  }else{
    Bx.input <- Bx.arg
  }
  
  smc.obj <- smc.init.xi(n.time, n.part)
  sample.size = n.part
  
  # exact likelihood:
  lower_trim <- 30
  weight.ll.exact <- function(xi.arg, idx){
    return(dnorm(Y.obs[idx],  f(xi.arg[1:idx], Bx.input[1:idx], delta.arg, A.arg = A.arg, b.arg = b.arg)[idx], sigma.arg, log = TRUE))
  }
  
  
  # initialize particles
  smc.obj$xi.part[, 1] <- rep(a.arg, n.part)
  smc.obj$w.part[, 1] <- rep(1 / n.part, n.part)
  
  for (i in 2:n.time) {
    
    if(track) print(i)
    
    xi.part <- transfer.xi.noncentralchisq(smc.obj$xi.part[,i-1], psi.arg, a.arg, omega2dBeta.arg)
    
    # Save particle draws
    smc.obj$xi.part[, i] <- xi.part
    
    # Compute weight and resample particles with high likelihood paths
    smc.obj$w.part[, i] <- lower_trim+0.1+pmax(apply(smc.obj$xi.part[, 1:i], 1, weight.ll.exact, i), -lower_trim) #-apply(smc.obj$xi.part[, (i-1):i],1, density.sum.xi.noncentralchisq, psi.arg, a.arg, omega2dBeta.arg)
    
    # Resample:
    likely.particles <- sample(n.part, size=n.part, replace=T, prob=smc.obj$w.part[,i])
    smc.obj$xi.part[, 1:i] <- smc.obj$xi.part[likely.particles, 1:i]
    smc.obj$w.part[, 1:i] <- smc.obj$w.part[likely.particles, 1:i]
    
    
  }
  
  # trajectory sampling
  draw <- sample(n.part,1,prob=smc.obj$w.part[,n.time])
  #draw <- which.max(smc.obj$w.part[, n.time])
  smc.obj$xi.c <- smc.obj$xi.part[draw, ]
  smc.obj$gxi.c <- g(smc.obj$xi.c,delta.arg)
  
  return(smc.obj)
}


#################################
#### --- SAEM - SMC  --- ########
#################################
warp.smc <- function(Y.arg,  delta.arg=delta, InputData, n.saem, n.mem, n.part ,track=TRUE, stopepsilon = 1e-1, fixedOmega=FALSE, fast=FALSE){
  require(minpack.lm)
  require(nnls)
  require(stats)
  
  n <- length(Y.arg) # number of observations
  x = 1:n 
  
  # small epsilon to avoid division by zero in rel differences for stopping criteria
  SmallEps <- 1e-10
  
  # SMC particles at each SAEM-iteration
  particles <- rep(n.part, n.saem) 
  particles[2:min(6,n.saem)] <- n.part*2 # the first at least five (note 1 = initialization) iterations use double particles
  
  # SAEM memory
  alpha <- rep(1, n.saem)
  steps <- ((1:(n.saem - n.mem + 1))^0.95)
  if(n.saem>n.mem) alpha[n.mem:n.saem] <- 1 / steps # add memory when n > n.mem
  
  # Initial parameters
  a0 <- InputData$pars$a0/delta.arg
  a0_min <- max(a0/2, 2*pi/(n*delta)) # period / 2, cannot be lower than 1 cycle
  a0_max <- 2*a0 # period *2
  b0 <- InputData$pars$b0
  A0 <- InputData$pars$A.0
  B0.0 <- InputData$pars$B0.0
  sigma0 <- InputData$pars$sigma0
  xkstar0 <- InputData$xkstar
  Bk0 <- InputData$rk + B0.0 # Bk = rk + B0
  Bx0 <- InputData$Bx
  
  # INITIALIZATION of trajectories
  a.est <- rep(a0, n.saem)
  b.est <- rep(b0, n.saem)
  A.est <- rep(A0, n.saem)
  B0.est <- rep(B0.0, n.saem)
  sigma.est <- rep(sigma0,n.saem)
  
  # Range of signal
  l <- max(Y.arg)-min(Y.arg)
  
  # Sufficient and approximated sufficient statistics
  s1.c <- 0 # sigma
  s2.c <- 0 # psi
  s3.c <- 0 # a
  s4.c <- 0 # omega^2
  s5.c <- 0 # sigma_B
  
  k <- 2 # first run (k = 1 is initialization)
  StoppingCriteriaMet <- FALSE # initiate stopping criteria bool
  r21 <- 1 # stopping criteria measure (calculated between two iterations)
  stopTracker <- 0 # number of consecutive r21 < stopthreshold (stop SAEM after 5)
  
  while(k <= n.saem && !StoppingCriteriaMet){
    
    # track SAEM ("print")
    if(track){
      cat("Current iteration:", k-1, " | Relative change:", r21, "\n")
    }
    
    if(k==2){
      # SMC ensemble: multiple SMC over random grid of (psi, omega^2)
      xi.obj <- smc.xi.ensemble(Y.arg, particles[[k]], 5000,
                                sigma.arg = sigma0,
                                A.arg = A0,
                                Bx.arg = Bx0,
                                a.arg = a0,
                                b.arg = b0,
                                delta.arg=delta.arg,
                                fixedOmega)
      
      # (New) initial estimate for: Xi, Psi (Beta) and Omega^2
      xi0 <- xi.obj$xi.c
      psi0 <- xi.obj$psi.c
      beta0 <- -1/delta.arg*log(psi0)
      omega2.0 <- xi.obj$omega2dBeta.c*beta0
      
      psi.est <- rep(psi0, n.saem) 
      beta.est <- rep(beta0, n.saem)
      omega2.est <- rep(omega2.0,n.saem)
      
      # (New) initial estimate for: Bk, sigma_B
      xkstar0 <- sort(c(FindXstar(g(xi0,delta.arg),b0), FindXstar(g(xi0,delta.arg),b0,3*pi/2)))
      BxData <- EstimateBx(Y.arg, xi0, A0, B0.0, b0, delta.arg)
      Bk0 <- B0.0 + BxData$rk.c
      Bx0 <- BxData$Bx.c
      sigma_B0 <- max(1e-1,sqrt(1/length(Bk0)*sum((B0.0-Bk0)^2))) # max in case |Bk0| = 1 ==> sigma_B = 0
      sigma_B.est <- rep(sigma_B0,n.saem)
      
      # INITIALIZE sufficient statistics (first run of SAEM)
      S1 <- sigma0
      S2 <- psi0
      S3 <- a0
      S4 <- omega2.0
      S5 <- sigma_B0 # std. dev of amplitude displacement random effects
      
    }else{
      # SMC_plus
      xi.obj <- smc.xi(Y.arg, particles[[k]],
                       xkstar.arg = xkstar0, # full model or simple model: xkstar0=0
                       Bx.arg=Bx0,  # full model
                       sigma.arg = sigma0,
                       psi.arg=psi0,
                       omega2dBeta.arg = omega2.0/beta0,
                       A.arg = A0,
                       a.arg= a0,
                       b.arg = b0,
                       delta.arg=delta.arg)
      xi0 <- xi.obj$xi.c
      xkstar0 <- sort(c(FindXstar(g(xi0,delta.arg),b0), FindXstar(g(xi0,delta.arg),b0,3*pi/2)))
      
      # Compute varying amplitude
      BxData <- EstimateBx(Y.arg, xi0, A0, B0.0, b0, delta.arg)
      Bx0 <- BxData$Bx.c
      Bk0 <- B0.0 + BxData$rk.c
    }
    
    # Sufficient statistics
    S1 <- mean((Y.arg - f(xi0, Bx0, delta.arg, A0, b0))^2)
    S2 <- MSF.psi(xi0) # Martingale Estimation
    S3 <- MSF.a(xi0, psi0) # Martingale Estimation
    S4 <- MSF.omega2(xi0, psi0, a0, delta.arg) # Martingale Estimation
    # Assume Bk ~ gaussian (but Bk ~ truncated gaussian)
    if(length(Bk0)>1){
      S5 <- sqrt(1/length(Bk0)*sum((Bk0 - mean(Bk0))^2))
    }else{
      S5 <- 0 # no information? set to zero
    }
    
    # SA 
    if(is.finite(S1)) s1.c <- s1.c + alpha[[k]] * (S1 - s1.c)
    if(is.finite(S2)) s2.c <- s2.c + alpha[[k]] * (S2 - s2.c)
    if(is.finite(S3)) s3.c <- s3.c + alpha[[k]] * (S3 - s3.c)
    if(is.finite(S4)) s4.c <- s4.c + alpha[[k]] * (S4 - s4.c)
    if(is.finite(S5)) s5.c <- s5.c + alpha[[k]] * (S5 - s5.c)
    
    # Update parameters
    # Only accept if inside parameter space
    sigma.0 <- sqrt(s1.c)
    psi0 <- if(s2.c<.98 && s2.c>0) s2.c else psi0 
    beta0 <- -1/delta.arg*log(psi0)
    a0 <- if(s3.c>a0_min && s3.c<a0_max) s3.c else a0 # a*N*delta/(2*pi) > 1 ==> a > 2*pi/(N*delta) && a*N*delta/2*pi < 12 ==> a < 24*pi/(N*delta)
    if(fixedOmega){ # omega^2 = a*beta (default: false)
      omega2.0 <- a0*beta0
    }else{
      omega2.0 <- if(s4.c < 2*a0*beta0) s4.c else omega2.0 # ensure omega^2 < 2*a*beta 
    }
    sigma_B0 <-  if(s5.c > 0) s5.c else sigma_B0
    
    # ### --- START: non-linear least squares ---###
    maxRx <- max(Bx0 - B0.0) # introduced to ensure that A > B(x) = B0 + R(x) ==> A - B0 > max(R(x))
    nls.obj <- nlsLM( # parametrization: B0.c = A.c - k - Beps
      Y.arg ~ f(xi0, A.c - maxRx - Beps + Bx0 - B0.0, delta.arg, A.c, b.c),
      start = list(A.c = A0, Beps = 0.01, b.c = b0 - floor(b0/(2*pi))*(2*pi)),
      lower = c(1/100*l, 0, 0),  # Lower bounds for A.c, Beps, and b.c
      upper = c(l, l-maxRx, 2*pi)   # Upper bounds for A.c, Beps, and b.c
    )
    A0 <- coef(nls.obj)[[1]]
    B0.0 <- A0 - maxRx  - coef(nls.obj)[[2]]
    b0 <- coef(nls.obj)[[3]]
    b0 <- b0-floor(b0/(2*pi))*(2*pi)
    # ###--- END:  non-linear least squares ---### 
    
    
    A.est[[k]] <- A0
    B0.est[[k]] <- B0.0
    b.est[[k]] <- b0
    a.est[[k]] <- a0
    psi.est[[k]] <- psi0
    beta.est[[k]] <- beta0
    sigma.est[[k]] <- sigma.0
    omega2.est[[k]] <- omega2.0
    sigma_B.est[[k]] <- sigma_B0
    
    
    ################################################################################################################
    # --- Stopping criteria: max of relative differences of estimates < threshold for 5 consecutive iterations --- #
    ################################################################################################################
    
    r21 <- min(
      sort( # or max. mean = avoid single unidentifiable parameters to slow down convergence
        c(abs(A.est[[k]]-A.est[[k-1]])/pmax(SmallEps,abs(A.est[[k]])), 
          abs(B0.est[[k]]-B0.est[[k-1]])/pmax(SmallEps,abs(B0.est[[k]])),
          abs(b.est[[k]]-b.est[[k-1]])/pmax(SmallEps,abs(b.est[[k]])),
          abs(a.est[[k]]-a.est[[k-1]])/pmax(SmallEps,abs(a.est[[k]])),
          abs(psi.est[[k]]-psi.est[[k-1]])/pmax(SmallEps,abs(psi.est[[k]])),
          abs(sigma.est[[k]]-sigma.est[[k-1]])/pmax(SmallEps,sigma.est[[k]]),
          abs(omega2.est[[k]]-omega2.est[[k-1]])/pmax(SmallEps,omega2.est[[k]]),
          abs(sigma_B.est[[k]]-sigma_B.est[[k-1]])/pmax(SmallEps,sigma_B.est[[k]])
        ),decreasing=TRUE)[2],    1)
    
    
    # Check if ratios is below threshold, otherwise reset StopTracker
    if(r21 < stopepsilon){
      stopTracker <- stopTracker + 1
    }else{
      stopTracker <- 0
    }
    
    # If 5 consecutive iterations has ratio < threshold, stoppingcriteria is met.
    if(k > 5 & stopTracker > 4) StoppingCriteriaMet <- TRUE
    
    #new run
    k <- k+1
    
  }
  print("Done.")
  print(paste("Ratio:",r21))
  
  # Before EXIT: Final estimate of Xi and B(x)
  print("Exit SMC")
  if(!fast){ # if fast: skip final estimation
    xi.obj <- smc.xi(Y.arg, ifelse(n.saem==2, 7000, 4000),
                     xkstar.arg=xkstar0, # full model
                     Bx.arg=Bx0,  # full model
                     sigma.arg = sigma0,
                     psi.arg=psi0,
                     omega2dBeta.arg = omega2.0/beta0,
                     A.arg = A0,
                     a.arg= a0,
                     b.arg = b0,
                     delta.arg=delta.arg)
    xi0 <- xi.obj$xi.c
    xkstar0 <- sort(c(FindXstar(g(xi0,delta.arg),b0), FindXstar(g(xi0,delta.arg),b0,3*pi/2)))
    
    # Compute varying amplitude
    BxData <- EstimateBx(Y.arg, xi0, A0, B0.0, b0, delta.arg)
    Bx0 <- BxData$Bx.c
    Bk0 <- B0.0 + BxData$rk.c
  }
  
  
  return(list(xi.obj = xi.obj, Bk.c = Bk0, Bx.c = Bx0, xkstar.c = xkstar0, n.final = k,
              a.c = a0, b.c = b0, A.c = A0, B0.c = B0.0, sigma_B.c = sigma_B0, 
              psi.c = psi0, beta.c = beta0, sigma.c = sigma0, omega2.c = omega2.0,
              a.est = a.est, b.est = b.est, A.est = A.est,
              B0.est = B0.est, sigma_B.est = sigma_B.est, psi.est = psi.est, beta.est = beta.est, 
              sigma.est = sigma.est, omega2.est = omega2.est))
}




