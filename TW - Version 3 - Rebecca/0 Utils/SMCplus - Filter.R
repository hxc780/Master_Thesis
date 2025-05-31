# Initialization SMC+: Multiple parameters 
# Lars Nørtoft Reiter, lrn@math.ku.dk, lars.reiter.nielsen@gmail.com

# ----- Base Components: Densities and Random Transition -----

# Transition simulator (CIR via noncentral chi-square)
rprocess <- function(xi_prev, a, beta, omega, dt = 1) {
  rho     <- exp(-beta * dt)
  omega2b <- omega^2 / beta
  c_const <- 2 / ((1 - rho) * omega2b)
  lambda  <- 2 * c_const * xi_prev * rho
  df      <- 4 * a / omega2b
  return((1 / (2 * c_const)) * rchisq(length(xi_prev), df = df, ncp = lambda))
}

# Helper: initialize SMC storage
create_smcdataobj <- function(N, n_part) {
  list(
    xi_part = matrix(NA, nrow = n_part, ncol = N),
    w_part  = matrix(NA, nrow = n_part, ncol = N),
    ess_track = numeric(N),
    xi.c  = numeric(N),
    gxi.c = numeric(N),
    b.c   = numeric(),
    beta.c  = numeric(),
    omega.c = numeric()
  )
}

# ----- SMC filter: beta, omega, b estimation -----
SMCplus <- function(y, params,
                       n_part     = 1500,
                       beta_range = c(0.01, 0.5),
                       omega_range= c(0.01, 0.1),
                       a_range = NULL,
                       dt          = 1,
                       Delta       = 60,
                       L           = 5,
                       xizero      = 0,
                       track       = FALSE) {
  
  # Unpack known parameters
  A     <- params[1]; sigma <- params[2]; 
  b     <- params[3];  N     <- length(y)
  
  if (is.null(a_range) || is.null(beta_range) || is.null(omega_range)) {
    stop("Must supply a_range, beta_range, and omega_range for SMCfilter initialization")
  }
  
  # Stage 1: Static parameter sampling
  feller   <- generate_feller_parameters(n_part,
                                         beta_range  = beta_range,
                                         omega_range = omega_range,
                                         a_range = a_range)
  a_vec <- feller$a
  beta_vec  <- feller$beta
  omega_vec <- feller$omega
  
  smc <- create_smcdataobj(N, n_part)
  smc$xi_part[,1] <- rep(xizero, n_part)
  smc$w_part[,1]  <- rep(1 / n_part, n_part)
  
  for (t in 2:N) {
    if (track) cat("t =", t, "\n")
    
    # Stage 2: Latent state filtering
    smc$xi_part[,t] <- rprocess(smc$xi_part[,t-1], a_vec, beta_vec, omega_vec, dt)
    gxi         <- rowSums(smc$xi_part[,1:t]) * dt
    mu_theta    <- f(gxi, A, b)
    
    # Log-weight trimming for numerical stability
    raw_logw <- dnorm(y[t], mean = mu_theta, sd = sigma, log = TRUE)
    raw_logw <- pmax(raw_logw, -Delta)  # Hard floor to avoid collapse
    
    w <- exp(raw_logw - max(raw_logw))
    w <- w / sum(w)
    smc$w_part[,t] <- w
    
    # Resample static parameters and particles
    if (t %% L == 0) {
      idx <- sample(seq_len(n_part), size = n_part, replace = TRUE, prob = w)
      smc$xi_part   <- smc$xi_part[idx,,drop=FALSE]
      smc$w_part    <- smc$w_part[idx,,drop=FALSE]
      a_vec         <- a_vec[idx]
      beta_vec      <- beta_vec[idx]
      omega_vec     <- omega_vec[idx]
    }
  }
  
  # Output: max-weight estimates and summaries
  idx_max <- which.max(smc$w_part[,N])
  smc$xi_est    <- smc$xi_part[idx_max,]
  smc$gxi_est   <- cumsum(smc$xi_est) * dt
  smc$a_est     <- a_vec[idx_max]
  smc$beta_est  <- beta_vec[idx_max]
  smc$omega_est <- omega_vec[idx_max]
  
  return(list(
    smc      = smc,
    logLik   = sum(log(colSums(smc$w_part))),
    summary  = list(
      a_est     = smc$a_est,
      beta_est  = smc$beta_est,
      omega_est = smc$omega_est
    )
  ))
}
