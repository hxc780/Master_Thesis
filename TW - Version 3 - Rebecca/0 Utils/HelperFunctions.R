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

# ---- VERIFY CONSTRAINTS ---- #
verify_constraints <- function(a.cand, rho.cand, omega2.cand, delta, amin, amax) {
  
  # rho must be between 0 and 1
  rho_ok <- (rho.cand > 0) && (rho.cand < 1)
  
  if(rho_ok){
    
    # Feller condition: 2 * a * beta > omega^2
    feller_ok <- (2 * a.cand * (-1/delta*log(rho.cand)) > omega2.cand)
    
    # a must be between amin and amax
    a_ok <- (a.cand > amin) && (a.cand < amax)
  }
  
  # Only accept if all conditions are met
  all_constraints_met <- feller_ok && a_ok && rho_ok
  
  return(all_constraints_met)
}

perturb_BetaAndOmega <- function(a.arg, beta.arg, prev_attempts, gamma=0.8, beta_lower = 0.01, beta_upper = 0.5, max_iter = 1000) {
  # cooling factor gamma for perturbation of beta and omega
  # *** strategy: for every iteration, that constraints are not met
  # increase the variance of the random walk. Whenever constraints are met
  # reset variance to lowest (rw) standard deviaton.
  require(truncnorm)
  
  # sd_fac goes to 1 when number of iterations of unmet constraints goes towards Inf
  sd_fac <- 1 - (1 - 1/20)*gamma^(prev_attempts + 1)  
  
  beta_new <- rtruncnorm(1, a=beta_lower, b=beta_upper, mean=beta.arg, sd = sd_fac*beta.arg)
  omega_frac_of_maxfeller <- runif(1, 0.2, 0.95) 
  omega_new <- omega_frac_of_maxfeller*sqrt(2*a.arg*beta_new)

  return(list(beta = beta_new, omega = omega_new))
}


# ----- Complete or exact log likelihood ----- #
tw.LogLik <- function(Y.obs, gx.arg, theta.arg, delta.arg, complete=FALSE){
  #  --- Exact CIR likelihood --- #
  
  beta.arg <- theta.arg[1]
  rho.arg <- exp(-delta.arg*beta.arg)
  a.arg <- theta.arg[2]
  omega.arg <- theta.arg[3]
  A.arg <- theta.arg[4]
  B.arg <- theta.arg[5]
  b.arg <- theta.arg[6]
  sigma.arg <- theta.arg[7]
  xi.arg <- diff(gx.arg)
  xi.arg <- c(xi.arg, xi.arg[length(xi.arg)])
  
  c <- (omega.arg^2 * (1 - rho.arg)) / (4*beta.arg)
  d <- 4 * beta.arg * a.arg / (omega.arg^2)
  lambda <- (4 * beta.arg * rho.arg * xi.arg[1:(length(xi.arg)-1)]) / (omega.arg^2 * (1 - rho.arg))
  
  ll <- sum(dnorm(Y.obs, f(gx.arg,  A.arg, b.arg), sigma.arg, log = TRUE))
  
  if(complete){
    ll + sum(dchisq( c * xi.arg[2:length(xi.arg)], df = d, ncp = lambda, log = TRUE))
  }
  
  return(ll)
}


# ----- SIMULATION AUXILIARY FUNCTIONS ----- #

g <- function(xi.arg, delta.arg = delta) {
  return(cumsum(xi.arg)*delta.arg)
}

f <- function(gx.arg, A.arg = A, b.arg = b) {
  A.arg * sin(gx.arg + b.arg)
}

# ------ END ------- #

# ------- INITIALIZATION TECHNIQUES  ------- #

# Estimate A and B, and sigma:
estimate_ABsigma <- function(Y.arg, span_start = 0.05){
  # THEORY: Let z = g(x) + b, i.e. f(z) = Asin(z) - cos(2z)
  # f'(z) = Acos(z)+2Bsin(2z) = cos(z)(A+4Bsinz) = 0 implies...
  # cosz = 0: solutions z = pi/ + k*pi, which yields f(z) = A+B (ymax) and f(z) = B-A
  # sinz = A/(4B) (only if A/(4B) <= 1) implies cos(2z) = 1- A^2/(8*B) resulting in
  # f(z) = -B - A^2/(8B) (ymin)
  # NOW: solve ymin = -B - A^2/(8B) and ymax = A+B for A and B. This gives:
  # B = 1/9*(ymax-4*ymin + c(-1,1)*sqrt((ymax-4*ymin)^2-9*ymax^2))
  # A = ymax - B
  
  n.obs <- length(Y.arg)
  X <- 1:n.obs
  span <- span_start  # initial spanning for loess smoother
  
  repeat{
    # Smooth the signal
    Ysm <- predict(loess(Y.arg ~ X, span = span), newdata = X)
    
    # Compute values
    ymax <- max(Ysm)
    ymin <- min(Ysm)
    d <- (ymax - 4 * ymin)
    D <- d^2 - 9 * ymax^2
    
    # Check condition
    if (D > 0 || span <= 0.01) {
      break  # Stop if D is positive or if span reaches 0.01
    }
    
    span <- span - 0.01  # Decrease span
  }
  
  # Compute A and B based on final values
  if (D > 0) {
    B.cand <- 1/9 *(d + c(-1, 1) * sqrt(D))
    A.cand <- ymax - B.cand
  } else {
    A.cand <- ymax
    B.cand <- 1/10 * ymax
  }
  
  # Compute also Sigma
  sigma.cand <- sqrt(mean((Ysm - Y.arg)^2))
  
  return(list(A.cand = A.cand, B.cand = B.cand, sigma.cand = sigma.cand))
}

# Estimate A and sigma:
estimate_Asigma <- function(Y.arg, span = 0.05){
  # THEORY: ymax = max(A*sin(z)) = A
  
  n.obs <- length(Y.arg)
  X <- 1:n.obs
  Ysm <- predict(loess(Y.arg ~ X, span = span), newdata = X)
  
  A.cand <- max(Ysm)
  sigma.cand <- sqrt(mean((Ysm - Y.arg)^2))
  
  return(list(A.cand = A.cand, sigma.cand = sigma.cand))
}


# Estimate a:
estimate_a <- function(Y.arg, delta.arg){
  require(signal)
  require(gsignal)
  
  # Compute analytic signal
  analytic_signal <- gsignal::hilbert(Y.arg)
  phase <- atan2(Im(analytic_signal), Re(analytic_signal))
  
  # Avoid discontinuities by "unwrapping" phase
  phase_unwrapped <- signal::unwrap(phase)
  
  # dphase/dt = freq.
  instantaneous_frequency <- c(0, diff(phase_unwrapped)) / delta.arg 
  
  return(mean(instantaneous_frequency))
}

# Estimate b: TODO: make this initial estimation
estimate_b <- function(Y.arg){
  return(0)
}


# ----- PIECEWISE LINEAR FUNCTION ----- #
pwlinear <- function(x, slopes, breakpoints) {
  return(
    ifelse(
      x < breakpoints[1],
      slopes[1] * x,
      ifelse(
        x < breakpoints[2],
        slopes[2] * (x - breakpoints[1]) + slopes[1] * breakpoints[1],
        slopes[3] * (x - breakpoints[2]) + slopes[2] * (breakpoints[2] - breakpoints[1]) + slopes[1] * breakpoints[1]
      )
    )
  )
}

# ----- MARTINGALE ESTIMATING FUNCTIONS ----- #
MSF.a <- function(S1.arg, S2.arg, S3.arg, rho.arg, n.xi = N){
  1/n.xi*S1.arg + rho.arg/(n.xi-n.xi*rho.arg)*(S2.arg - S3.arg)
}

MSF.rho <- function(S1.arg, S4.arg, S5.arg, S6.arg, n.xi=N){
  (n.xi*S4.arg - S1.arg*S5.arg)/(n.xi^2-S5.arg*S6.arg)
}

MSF.omega2 <- function(S1.arg, S4.arg, S5.arg, S6.arg, S7.arg, rho.arg, a.arg, n.xi=N){
  num <- (S7.arg + rho.arg^2*S6.arg + a.arg^2*(1-rho.arg)^2*S5.arg-2*rho.arg*S1.arg-2*a.arg*(1-rho.arg)*S4.arg+2*rho.arg*a.arg*(1-rho.arg))
  den <- (1/2*a.arg*rho.arg^2*S5.arg-rho.arg^2-a.arg*rho.arg*S5.arg+rho.arg+1/2*a.arg*S5.arg/rho.arg)
  return(num/den)
}

# --- Rejection Sampling of Beta / Omega under Fellers constraint --- #
Beta_Omega_Candidate_Generator <- function(a.arg, 
                                           beta_range = c(0.01, 0.5), 
                                           omega_range = c(sqrt(0.004), 0.1)) {
  repeat {
    beta_candidate <- runif(1, beta_range[1], beta_range[2])
    omega_candidate <- runif(1, omega_range[1], omega_range[2])
    
    # Compute the ratio beta/omega^2
    ratio <- beta_candidate / (omega_candidate^2)
    
    # Check Feller constraint and empirical constraint
    if (ratio > 1 / (2 * a.arg) && ratio < 62) {
      return(list(beta = beta_candidate, omega = omega_candidate))
    }
  }
}

# -- Rejection Sampling of Beta, Omega, a under Feller constraint using Parallel processing -- #
generate_feller_parameters <- function(n, 
                                       beta_range = c(0.01, 0.5), 
                                       omega_range = c(0.01, 0.1), 
                                       a_range = c(0.02, 0.1)) {
  
  g <- expand.grid(beta=seq(beta_range[1], beta_range[2], length.out=n^(1/3)),
                   omega=seq(omega_range[1], omega_range[2], length.out=n^(1/3)),
                   a=seq(a_range[1], a_range[2], length.out=n^(1/3)))
  
  feller_candidates <- which(2*g[,3]*g[,1]>2*g[,2]^2)
  g <- g[sample(feller_candidates, n, replace=TRUE),]
  row.names(g) <- NULL
  
  return(as.data.frame(g))
}

