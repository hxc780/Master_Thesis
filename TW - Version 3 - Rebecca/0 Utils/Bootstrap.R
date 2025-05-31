# RESIDUAL BOOTSTRAPPING
resbs <- function(Y.obs, Y.fit) {
  # Calculate residuals
  residuals <- Y.obs - Y.fit
  
  # Bootstrap residuals by resampling with replacement
  boot_residuals <- sample(residuals, length(residuals), replace = TRUE)
  
  # Bootstrapped signal
  Y.boot <- Y.fit + boot_residuals
  
  return(Y.boot)
}

tw_bootstrap <- function(y.obs, y.fit, M=10, n.saem=80, n.particles=1200, minPeaks=5, maxPeaks=25){
  
  # Bootstrap initialization
  total_cycles <- list()
  fitting_objects <- list()
  
  m <- 1
  while(m <= M){
    print(paste("BS:",m))
    
    tryCatch({
      y.new <- resbs(y.obs, y.fit) # residual bootstrap replicate
      
      fit <- warp.alg(y.new,
                      n.saem=n.saem, 
                      n.particles = n.particles,
                      MinPeaks = minPeaks, # minimum number of cycles (= years)
                      MaxPeaks = maxPeaks, # maximum number of peaks (= years)
                      Csyntax = TRUE) # use C-syntax: faster, but normal approximated model
      
      total_cycles[[m]] <- fit$gx.hat[length(y.new)]/(2*pi)
      fitting_objects[[m]] <- fit
      m <- m + 1  # Increment only if successful
      
    }, error = function(e){
      message("Error encountered, retrying iteration ", m, ": ", e$message)
    })
  }
  
  bootstrap_data <- list(total_cycles = total_cycles, fitting_objects = fitting_objects)
  
  return(bootstrap_data)
  
}