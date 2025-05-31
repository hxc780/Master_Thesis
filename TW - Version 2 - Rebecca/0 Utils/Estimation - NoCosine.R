###############################################
#### --- MAIN Estimation Function  --- ########
###############################################
warp.alg <- function(Y.arg=NULL, params = NULL, delta = 1,  n.saem = 10, n.particles = 2000, n.mcmc = 50, n.memory = 5, n.siman = 5000, SimAn.MinPeaks = .5, SimAn.MaxPeaks = 12, SimAn.Pen = 1, Csyntax = FALSE, justPlot = FALSE){
  
  # input vector "params" is either NULL or contain:
  # params[1] = beta
  # params[2] = omega
  # params[3] = a
  # params[4] = A
  # params[5] = b
  # params[6] = sigma
  # params[7] = N
  
  # INITIALIZATION
  # Note: init.data is a list of 5 items:
  # item 1: pars      list (a0, b0, A.0, sigma0)
  # item 2: gx        vector g(x) integrated CIR process
  # item 3: xi        vector xi(x) CIR process
  # item 4: xkstar    vector peak position indices
  # item 5: Y.an.fit  vector \hat Y(x) fitted with sim.an parameters
  # note: rho0, beta0, omega2.0 are initialized directly in the estimation/warp.alg
  require(pomp)
  require(minpack.lm)
  require(truncnorm)
  require(foreach)
  
  ########################################################################
  #################          PART 1           ############################
  ########################################################################
  ###                                                                  ###
  ### SIMULATE SIGNAL AND/OR DEFINE STATE SPACE MODEL AS A POMP OBJECT ###
  ###                                                                  ###
  ########################################################################

  # --- Define Basic Components for POMP-object --- #
  # these includes: dmeasure: marginal likelihood
  # rmeasure: random generator for f(x)
  # rprocess: random step generator for stochastic process
  
  if(Csyntax){ # USE C-syntax (computational efficient)
    
    par_names <- c("bta", "a", "omega", "A", "b", "sigma", "X.0", "gX.0")
    
    # process simulator (rprocess)
    rproc_tw <- Csnippet("
  // Compute auxiliary constants
  double rho = exp(-bta*dt);
  double gamma_val = omega*omega*(1-rho)/bta;
  double mean_X  = X*rho+a*(1-rho);
  double var_X = gamma_val*(X*rho + a*(1-rho)/2.0);
  double epsilon = 1e-6;  // Small threshold to avoid exactly zero
  
  // Update state variables
  // Note: non-central chi squared approximated by normal density in C
  // Truncate to ensure positivity
  X = (var_X > 0.0) ? fmax(rnorm(mean_X, sqrt(var_X)), epsilon) : 0.0; ;
  
  gX = gX + X * dt;
")
    
    # observation simulator (rmeasure)
    rmeas_tw <- Csnippet("
  // Compute the signal f(x)
  double f_x = A * sin(gX + b);
  // Generate observation with Gaussian error
  Y = rnorm(f_x, sigma);
")
    
    # measurement density (dmeasure)
    dmeas_tw <- Csnippet("
  lik = dnorm(Y, A * sin(gX + b), sigma, give_log);
") 
  }else{ # Use R syntax (less computationally efficient)
    
    par_names <- c("beta", "a", "omega", "A", "b", "sigma", "X_0", "gX_0")
    
    rproc_tw <- function(X, gX, beta, a, omega, ..., delta.t) {
      
      c <- (omega^2 * (1 - exp(-beta*delta.t))) / (4*beta)
      d <- 4 * beta * a / (omega^2)
      lambda <- (4 * beta * exp(-beta*delta.t) * X) / (omega^2 * (1 - exp(-beta*delta.t)))
      
      # Sample from the non-central chi-squared distribution
      X <- c*rchisq(1, df = d, ncp = lambda)
      gX = gX + X*delta.t
      
      c(X = X, gX = gX)
    }
      
      # observation simulator (rmeasure)
      rmeas_tw <- function(X, gX, A, b, sigma, ...) {
        
        # Compute f(x)
        f_x <- A * sin(gX + b)
        
        # Add normal noise
        c(Y = rnorm(n = 1, mean = f_x, sd = sigma))
        
      }
      
      # measurement density (dmeasure)
      dmeas_tw <- function(A, b, sigma, X, gX, Y, ... , log = FALSE) {
        
        # Compute f(x)
        f_x <- A * sin(gX + b)
        
        dnorm(x=Y, mean = f_x, sd=sigma,log=log)
      }
  }
  
  if(is.null(Y.arg) & is.null(params)){ # draw random parameters, simulate, then estimate
    print("Generating parameters and simulating signal...")
    N <- sample(c(200, 400, 600, 800),1) # divide into 4 groups of different length
    Cmin <- 2 # min cycles
    Cmax <- pmax(6, N/100*1.5) # max cycles
    beta <- runif(1, 0.1, 3)
    rho <- exp(-beta*delta)
    sigma <- runif(1, 0.1, 0.5)
    A <- runif(1, 0.4, 0.6)
    a <- runif(1,min=2*pi*Cmin/N,max=2*pi*Cmax/N) # from approx 1~12 cycles in a signal
    b <- runif(1,min=0,max=2*pi)
    omega <- runif(1,0.3,0.9)*sqrt(2*a*beta)
    
    if(Csyntax){
      theta <- c(bta = beta, omega = omega, a = a, A = A, 
                b = b, sigma = sigma, X.0 = a, gX.0 = 0)
    }else{ # R syntax
      theta <- c(beta = beta, omega = omega, a = a, A = A, 
                b = b, sigma = sigma, X_0 = a, gX_0 = 0)
    }
    
    tw.pomp <- simulate(
      times = 1:N, t0 = 0,
      params=theta,
      
      rprocess = discrete_time(step.fun = rproc_tw, delta.t = delta),
      rmeasure = rmeas_tw,
      dmeasure = dmeas_tw,
      
      statenames = c("X", "gX"),
      obsnames = "Y",
      paramnames = par_names
    )
    
  }else if(is.null(Y.arg)){ # simulate signal from parameters,then estimate
    print("Simulating signal from input parameters...")
    
    if(Csyntax){
      theta <- c(bta = params[1], omega =  params[2], a =  params[3], A =  params[4], 
                b =  params[5], sigma =  params[6], X.0 = a, gX.0 = 0)
      N <- params[7]
    }else{ # R-syntax
      theta <- c(beta = params[1], omega =  params[2], a =  params[3], A =  params[4], 
                 b =  params[5], sigma =  params[6], X_0 = 0, gX_0 = 0)
      N <- params[7]
    }
    
    tw.pomp <- simulate(
      times = 1:N, t0 = 0,
      params = theta,
      
      rprocess = discrete_time(step.fun = rproc_tw, delta.t = delta),
      rmeasure = rmeas_tw,
      dmeasure = dmeas_tw,
      
      statenames = c("X", "gX"),
      obsnames = "Y",
      paramnames = par_names
    )
    
  }else if(is.null(params)){ # estimate parameters given signal
    
    print("Saving input signal in pomp-object...")
    tw.pomp <- pomp(
      data = data.frame(Y = Y.arg, time=1:length(Y.arg)),
      times="time",
      t0 = 0,
      
      rprocess = discrete_time(step.fun = rproc_tw, delta.t = delta),
      rmeasure = rmeas_tw,
      dmeasure = dmeas_tw,
      
      statenames = c("X", "gX"),
      obsnames = "Y",
      paramnames = par_names
    )
    
  }else{ # abort with error: nothing to estimate
    stop("Error: Nothing to estimate. Please check your input data or model parameters.")
  }
  
  Y0 = tw.pomp@data[1,]
  N <- length(Y0)
  
  if(justPlot){
    plot(tw.pomp)
    return(tw.pomp)
  }
  
  
  ########################################################################
  #################          PART 2           ############################
  ########################################################################
  ###                                                                  ###
  ###                    INITIALIZE PARAMETERS                         ###
  ###                                                                  ###
  ########################################################################
  
  # --- Initialization using Simulated Annealing --- #
  print("Initialization step.")
  
  # Simulated Annealing initialization:
  simAn.data <- SimAn.Init(tw.pomp@data[1,], iterations=n.siman, minCycles = SimAn.MinPeaks, maxCycles = SimAn.MaxPeaks, maxRatio = 5, penaltyWeight = SimAn.Pen, NoCosine=TRUE)
  a0 <- simAn.data$pars$a0
  A0 <- simAn.data$pars$A0
  b0 <- simAn.data$pars$b0
  sigma0 <- simAn.data$pars$sigma0
  
  ## Simple initialization:
  #a0 <- estimate_a(Y0, delta)
  #Asigma.0 <- estimate_Asigma(Y0)
  #A0 <- Asigma.0$A.cand
  #sigma0 <- Asigma.0$sigma.cand
  #b0 <-estimate_b(Y0)
  
  beta0 <- 0.1
  omega0 <- sqrt(a0*beta0) # omega^2/beta =  a < 2a (feller) and beta = 0.1
  
  # New guess at estimates:
  if(Csyntax){
    theta.guess <- c(bta = beta0, omega=omega0, a=a0, 
                   A = A0, 
                   b = b0, sigma = sigma0,
                   fellerfrac=0)
  }else{
    theta.guess <- c(beta = beta0,  omega=omega0, a=a0, 
                     A = A0,
                     b = b0, sigma = sigma0,
                     fellerfrac=0)
  }

  ########################################################################
  #################          PART 3           ############################
  ########################################################################
  ###                                                                  ###
  ###        SAEM-SMC: ESTIMATE PARAMETERS AND STOCHASTIC PROCESS      ###
  ###                                                                  ###
  ########################################################################
  
  # Initialize matrix of parameter estimates
  param.est <- matrix(ncol = length(theta.guess), 
                      nrow = n.saem + 1)
  param.est[1,] <- theta.guess
  colnames(param.est) <- c("beta",names(theta.guess)[-1])
  
  # Initialize rho estimate (auxilliary parameter)
  rho0 <- exp(-delta*beta0)
  
  # Sufficient and approximated sufficient statistics
  S1 <- s1 <- 0 #
  S2 <- s2 <- 0 # 
  S3 <- s3 <- 0 # 
  S4 <- s4 <- 0 # 
  S5 <- s5 <- 0 # 
  S6 <- s6 <- 0 # 
  S7 <- s7 <- 0 # 
  
  # squeezing/memory
  n.squeeze = min(75, round(n.saem/2))
  alpha <- rep(1, n.saem) # Kuhn and Lavielle (2004): alpha = (k-K)^-1 when k > K, else 1
  alpha[(n.squeeze+1):n.saem] <- ((n.squeeze+1):n.saem - n.squeeze)^(-0.8) # ændret fra 1 til 0.8 her

  # Range of signal
  l <- max(Y0)-min(Y0)
  
  # Min and max frequency based on cycles
  amin <- max(a0/2, 2*pi/length(Y0)) # period / 2, cannot be lower than 1 cycle
  amax <- 2*a0 # period *2
  
  # initialize iteration
  m <- 2
  m_eff <- 2 # effective time (counting also rejection steps)
  m_update <- 10 # update every ?? iteration
  
  print("Estimation step: SAEM.")
  
  #################################
  # --- SAEM Algorithm: Start --- #
  #################################
  while(m  < n.saem+1 && m_eff < 5*n.saem){
    
    if((m_eff-1)%%m_update == 0) print(paste(paste0(round(100*max(m/n.saem, m_eff/(5*n.saem)),2),"%"), "complete"))
    
    # bootstrap SMC filter
    pobj <- pfilter(tw.pomp,
                    params = c(theta.guess, X.0 = a0, gX.0 = 0),
                    Np = n.particles,
                    filter.traj=TRUE)
    xi0 <- pobj@filter.traj[1,,1:N] # \hat xi
    gx0 <- pobj@filter.traj[2,,1:N] # \hat g(x)
    
    # (2) Update sufficient statistics
    S1 <- sum(xi0[2:N]) # used for a, rho, omega estimator
    S2 <- xi0[N] # used for a-estimator
    S3 <- xi0[1] # used for a-estimator
    S4 <- sum(xi0[2:N]/xi0[1:(N-1)]) # used for omega, rho estimator
    S5 <- sum(xi0[1:(N-1)]^(-1)) # used for omega, rho estimator
    S6 <- sum(xi0[1:(N-1)]) # used for omega, rho estimator
    S7 <- sum(xi0[2:N]^2/xi0[1:(N-1)]) # used for omega estimator
    
    # (3) Ensure candidates fulfill constraints
    # Feller condition and a \in (amin, amax) as well as rho \in (0,1)
    rho.cand <- min(MSF.rho(S1, S4, S5, S6, N), .99)
    a.cand <- MSF.a(S1, S2, S3, rho.cand, N)
    omega2.cand <- MSF.omega2(S2, S4, S5, S6, S7, rho.cand, a.cand, N)
    all_constraints_met <- verify_constraints(a.cand, 
                                              rho.cand,
                                              omega2.cand,
                       delta,
                       amin,amax)
    
    if(all_constraints_met){ # only update "information" if constraints are met
      
      if(is.finite(S1)) s1 <- s1 + alpha[m-1] * (S1 - s1)
      if(is.finite(S2)) s2 <- s2 + alpha[m-1] * (S2 - s2)
      if(is.finite(S3)) s3 <- s3 + alpha[m-1] * (S3 - s3)
      if(is.finite(S4)) s4 <- s4 + alpha[m-1] * (S4 - s4)
      if(is.finite(S5)) s5 <- s5 + alpha[m-1] * (S5 - s5)
      if(is.finite(S6)) s6 <- s6 + alpha[m-1] * (S6 - s6)
      if(is.finite(S7)) s7 <- s7 + alpha[m-1] * (S7 - s7)
      
      # (4) Re-estimate parameters
      rho0 <- min(MSF.rho(s1, s4, s5, s6, N), .99)
      beta0 <- -1/delta*log(rho0)
      a0 <- MSF.a(s1, s2, s3, rho0, N)
      omega0 <- sqrt(MSF.omega2(s2, s4, s5, s6, s7, rho0, a0, N))
      
      
      # Non-linear least squares estimation of A and b
      nls.obj <- nlsLM( 
        Y0 ~ f(gx0, A.c, 0, b.c),
        start = list(A.c = A0, b.c = b0 - floor(b0/(2*pi))*(2*pi)),
        lower = c(1/100*l, 0),  # Lower bounds
        upper = c(l, 2*pi)   # Upper bounds 
      )
      
      A0 <- coef(nls.obj)[[1]]
      b0 <- coef(nls.obj)[[2]]
      b0 <- b0-floor(b0/(2*pi))*(2*pi)
      sigma0 <- sqrt(mean((Y0 - f(gx0, A0, 0, b0))^2))
      
      # New guess at estimates:
      if(Csyntax){
        theta.guess <- c(bta = beta0, omega = omega0, a = a0, 
                         A = A0,
                         b=b0, sigma=sigma0,
                         fellerfrac = omega0^2/(2*a0*beta0))
      }else{ # R syntax
        theta.guess <- c(beta = beta0,  omega = omega0, a = a0,
                         A = A0, 
                         b=b0, sigma=sigma0,
                         fellerfrac = omega0^2/(2*a0*beta0))
      }
      
      param.est[m,] <- theta.guess
      m <- m + 1
      m_eff <- m_eff + 1
    }else{
      ## Inform (in console) if constraints are met
      #print(paste("SAEM iteration:", m-1, "- constraints NOT met."))
      m_eff <- m_eff + 1
      # add a small perturbation to beta, omega
      feller_frac <- rtruncnorm(1, a=0, b=1, mean=omega0^2/(2*a0*beta0), sd=1/10)
      beta0 <- rtruncnorm(1, a=0.01, b=0.5, mean=beta0, sd=beta0/10)
      omega0 <- feller_frac*sqrt(2*a0*beta0) # feller_frac ensures feller condition: omega < sqrt(2*a*beta)
      
    }
      
  }
  
  names(theta.guess) <- c("beta", names(theta.guess)[-1])
  
  print("100% complete.")
  
  fit.data <- list(Y = Y0,
                   Y.hat = f(gx0, A0, 0, b0),
                   gx.hat = gx0,
                   xi.hat = xi0,
                   param.est = param.est,
                   theta.hat = theta.guess)
    
  return(fit.data)
  
}

