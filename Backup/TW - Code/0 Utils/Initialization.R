# INITIALIZATION methods should always return a list element containing 7 items:
# item 1: pars      list (a0, b0, A.0, B0.0, sigma0)
# item 2: gx        vector g(x) integrated CIR process
# item 3: xi        vector xi(x) CIR process
# item 4: xkstar    vector peak position indices
# item 5: rk        vector r_k* peak k displacement
# item 6: Bx        vector B(x) = B0.0 + R(x) varying amplitude
# item 7: Y.an.fit  vector \hat Y(x) fitted with sim.an parameters
# note: psi0, beta0, omega2.0 are initialized directly in the estimation/warp.alg

#####################################################################################################
## INITIALIZATION METHOD 1: Fit model using Simulated Annealing and assume g(x) = piecewise linear ##
#####################################################################################################
SimAn.Init <- function(Y.arg, sim.its = 5000, min_peaks = .5, max_peaks = 10, max_ratio = 5, pen = 1, bp2 = TRUE){
  require(GenSA) # sim annealing package
  
  x_length <- length(Y.arg)
  x.arg <- seq(1, x_length)
  
  years_to_a <- function(years)
    2 * pi * years / x_length
  
  min_slope <- years_to_a(min_peaks / 2)
  max_slope <- years_to_a(max_peaks+5)
  
  l = max(Y.arg) - min(Y.arg)
  
  # penalization parameters
  max_y_dev = sqrt(sum(Y.arg^2)) # when fit = 0
  
  if(bp2){ # ALTERNATIVE OPTION: 2 change (frequency) points
    
    # parameters in order: b, level, HT, LT, slope 1 (alpha), slope 2 (beta), breakpoint 1 (nu)
    b_bounds <- c(0, 2*pi)
    dA_bounds <- c(1/100*l, l)
    B_bounds <- c(1/100*l, l)
    par_lower <- c(b_bounds[1], dA_bounds[1], B_bounds[1], rep(0,3), 0, 0)
    par_upper <- c(b_bounds[2], dA_bounds[2], B_bounds[2],  rep(max_slope,3), x_length, x_length)
    
    fx <- function(x_trans, v)
    {
      b <- v[1] # phase offset
      A <- v[3] + v[2]  # ensure that A > B (primary cycle dominates)
      B <- v[3] 
      
      A * sin(x_trans + b) - B * cos(2 * x_trans + 2 * b)
    }
    
    get_SS <- function(x, y)
    {
      function(v)
      {
        # restrict frequency jumps
        slope.ratios <- c(v[4]/v[5], v[5]/v[6])
        if (any(is.na(slope.ratios)) || any(slope.ratios<1/max_ratio) || any(slope.ratios> max_ratio))
          return(Inf)
        
        # restrict to minGLG and maxGLG
        gx <- linfunc2bp(x, v[4:6],v[7]+c(0,v[8]))
        agee <- gx[length(x)]/(2*pi)
        
        # penalize deviation from expected GLG = 1/2*(maxp+minp)
        penalty_term <- max_y_dev * ((min_peaks / agee - 1) * (agee < min_peaks) + (agee / max_peaks - 1) * (agee > max_peaks))^2
        
        # squared loss:
        sum((fx(gx, v) - y)^2) + pen*penalty_term
      }
    }
    
    SS <- get_SS(x.arg, Y.arg)
    
    gensa_opt <- GenSA::GenSA(fn = SS,
                              lower = par_lower,
                              upper = par_upper,
                              control = list(
                                smooth = TRUE,
                                simple.function = FALSE,
                                trace.mat = FALSE,
                                threshold.stop = 0.1,
                                maxit = sim.its
                              ))
    
    # save result
    gx.final = linfunc2bp(x.arg, gensa_opt$par[4:6],gensa_opt$par[7]+c(0,gensa_opt$par[8]))
    xi.final = diff(gx.final)
    xi.final[x_length] <- xi.final[x_length-1]
    xkstar.final <- sort(c(FindXstar(gx.final, gensa_opt$par[1], pi/2),FindXstar(gx.final, gensa_opt$par[1], 3*pi/2)))
    BxData <- EstimateBx(Y.arg, xi.final, gensa_opt$par[2], gensa_opt$par[3], gensa_opt$par[1])
    Bx.final <- BxData$Bx
    rk.final <- BxData$rk
    Y.an.fit <- f(xi.final, Bx.final, 1, gensa_opt$par[2], gensa_opt$par[1])
    
  }else{ # OPTION DEFAULT: 1 change (frequency) point
    
    # parameters in order: b, level, HT, LT, slope 1 (alpha), slope 2 (beta), breakpoint 1 (nu)
    b_bounds <- c(0, 2*pi)
    dA_bounds <- c(1/100*l, l)
    B_bounds <- c(1/100*l, l)
    par_lower <- c(b_bounds[1], dA_bounds[1], B_bounds[1], rep(0,2), 0)
    par_upper <- c(b_bounds[2], dA_bounds[2], B_bounds[2],  rep(max_slope,2), x_length)
    
    fx <- function(x_trans, v)
    {
      b <- v[1] # phase offset
      A <- v[3] + v[2]  # ensure that A > B (primary cycle dominates)
      B <- v[3] 
      
      A * sin(x_trans + b) - B * cos(2 * x_trans + 2 * b)
    }
    
    get_SS <- function(x, y)
    {
      function(v)
      {
        # restrict frequency jumps
        if (is.na(v[4]/v[5]) || v[4]/v[5]<1/max_ratio || v[4]/v[5]> max_ratio)
          return(Inf)
        
        # restrict to minGLG and maxGLG
        gx <- linfunc1bp(x, v[4:5],v[6])
        agee <- gx[length(x)]/(2*pi)
        
        # penalize deviation from expected GLG = 1/2*(maxp+minp)
        penalty_term <- max_y_dev * ((min_peaks / agee - 1) * (agee < min_peaks) + (agee / max_peaks - 1) * (agee > max_peaks))^2
        
        # squared loss:
        sum((fx(gx, v) - y)^2) + pen*penalty_term
        
      }
    }
    
    SS <- get_SS(x.arg, Y.arg)
    
    gensa_opt <- GenSA::GenSA(fn = SS,
                              lower = par_lower,
                              upper = par_upper,
                              control = list(
                                smooth = TRUE,
                                simple.function = FALSE,
                                trace.mat = FALSE,
                                threshold.stop = 0.1,
                                maxit = sim.its
                              ))
    
    
    gx.final = linfunc1bp(x.arg, gensa_opt$par[4:5],gensa_opt$par[6])
    xi.final = diff(gx.final)
    xi.final[x_length] <- xi.final[x_length-1]
    xkstar.final <- sort(c(FindXstar(gx.final, gensa_opt$par[1], pi/2),FindXstar(gx.final, gensa_opt$par[1], 3*pi/2)))
    BxData <- EstimateBx(Y.arg, xi.final, gensa_opt$par[2], gensa_opt$par[3], gensa_opt$par[1])
    Bx.final <- BxData$Bx
    rk.final <- BxData$rk
    Y.an.fit <- f(xi.final, Bx.final, 1, gensa_opt$par[2], gensa_opt$par[1])
    
  }
  
  # save and return initialization elements
  InitData <- list(pars = list(a0 = mean(diff(gx.final)),
                               b0 = gensa_opt$par[1],
                               A.0 = gensa_opt$par[3] + gensa_opt$par[2],
                               B0.0 = gensa_opt$par[3],
                               sigma0 = sqrt(mean((Y.arg - Y.an.fit)^2))),
                   gx = gx.final,
                   xi = xi.final,
                   xkstar = xkstar.final,
                   rk = rk.final,
                   Bx = Bx.final,
                   Yfit = Y.an.fit)
  
  return(InitData)
  
}






#####################################################################################################
################## INITIALIZATION METHOD 2: Fit model using AMPD - peak detection ###################
################## LINK TO PAPER:     https://www.mdpi.com/1999-4893/5/4/588      ###################
#####################################################################################################
AMPD.init <- function(signal, correction = TRUE, corr_neighbors = 0.01, beta=0.2){
  
  #####################################
  ###     PART 1: PEAK DETECTION    ###
  #####################################
  require(segmented)
  require(stats)
  require(minpack.lm)
  # Step 1: Detrend the signal (piecewise linear with 1 bp)
  x  <- seq_along(signal)
  y_sm <- stats::predict(loess(signal~x,span=0.15)) # rough smooth
  lm_fit <- lm(y_sm ~  x)
  seg_fit <- segmented(lm_fit, seg.Z = ~x, psi = median(x))  
  piecewise_trend <- predict(seg_fit)
  y_det <- signal - piecewise_trend
  N_org <- length(y_det)
  
  # pad y
  pad <- round(beta*N_org)
  y <- c(rep(y_det[1], pad), y_det, rep(y_det[N_org], pad)) # extended/padded y
  
  N <- length(y)
  L <- ceiling(N / 2) - 1
  
  # select uniform random number
  alpha <- 1
  
  # Step 2: Calculate Local Maxima Scalogram (LMS)
  LMS <- matrix(runif(N*L)+alpha, nrow = L, ncol = N)
  
  for (k in 1:L) {
    # Define window length
    for (i in (k + 2):(N - k + 1)){
      # k = 1, i = (k+2) = 3 then y2>y1 and y2>y3, then y1 is scale-max.
      # k = L = N/2 -1, i runs from  N/2+1 to N/2+2
      if (y[i - 1 ] > y[i - 1 - k] && y[i - 1] > y[i - 1 + k]){
        LMS[k,i] <- 0
        #print(0)
      }
    }
  }
  
  # Remove padding
  LMS <- LMS[,-c(1:pad, (N-pad+1):N)]
  
  # Step 3: Summarize Local Maxima
  gamma <- rowSums(LMS)
  lambda <- which.min(gamma)
  
  # Step 4: Reshape LMS matrix
  LMS_rescaled <- LMS[1:lambda, ]
  
  # Step 5: Detect Peaks
  column_sd <- apply(LMS_rescaled, 2, sd)
  peak_indices <- which(column_sd == 0)
  
  if(correction){
    peaks_cor <- peak_indices
    for (i in 1:length(peak_indices)) {
      p.idx <- peak_indices[i]
      lower <- max(p.idx-round(corr_neighbors*N_org),1)
      upper <- min(p.idx+round(corr_neighbors*N_org), N_org)
      max.neighbor <- which.max(signal[lower:upper])
      peaks_cor[i] <- max.neighbor + lower -1
    }
    peak_indices <- peaks_cor
  }
  
  
  ############################################
  ###     PART 2: Construct Init Element   ###
  ############################################
  NumberOfPeaks <- length(peak_indices)
  
  PeakIndicator <- numeric(length(peak_indices)) # vector = 1 if high peak, 0 otherwise.
  
  # CLUSTER PEAKS INTO LOW/HIGH PEAK CATEGORY 
  # AND ESTIMATE BASELINE FREQUENCY a
  if(NumberOfPeaks<2){ # only one spotted peak (assume max peak)

    HT <- signal[peak_indices]
    LT <- min(signal)
    
    fullperiod = 2*min(N-peak_indices, peak_indices) 
    f0 = 1/(fullperiod) # nominal frequency
    a0 = f0*2*pi # angular frequency = 2*pi*f0
    
    PeakIndicator[1] <- 1 # high peak
    
  }else if(NumberOfPeaks < 3){# only two spotted peaks

    minmax <- (signal[peak_indices]-min(signal))/(max(signal[peak_indices])-min(signal))
    
    PeakIndicator <- as.integer(minmax > 0.5) # high peaks ( = 1)
    HT.peaks <- peak_indices[minmax > 0.5]
    HT <- signal[HT.peaks] # if the smallest peak > 0.5 of the highest, then it is a high peak
    
    if(length(HT)>1){ # if both peaks are high
      
      LT <- min(signal) # then LT is the minimum value in the signal
      
      fullperiod = diff(HT.peaks) 
      f0 = 1/(fullperiod) # nominal frequency
      a0 = median(f0*2*pi) # angular frequency = 2*pi*f0
      
    }else{ 
      
      # only one high peak
      LT <- signal[peak_indices[minmax <= 0.5]] 
      
      fullperiod = 2*diff(peak_indices) 
      f0 = 1/(fullperiod) # nominal frequency
      a0 = f0*2*pi # angular frequency = 2*pi*f0
      
    }
  }else{ # > 2 spotted peaks

    minmax <- (signal[peak_indices]-min(signal))/(max(signal[peak_indices])-min(signal))
    threshold <- median(minmax) - sd(minmax) # assume median \approx true mean, and sd \approx true sd, and HT normal dist., then we want P(reject as HT | is HT) \approx 1/3
    PeakIndicator <- as.integer(minmax > threshold) # high peaks ( = 1)
    HT.peaks <- peak_indices[minmax > threshold]
    HT <- signal[HT.peaks] # peaks higher than quartile of minmax is high peak
    
    if(length(HT) == length(peak_indices)){ # if all peaks are high
      
      LT <- min(signal) # then LT is the minimum value in the signal
      
      fullperiod = diff(peak_indices) 
      f0 = 1/(fullperiod) # nominal frequency
      a0 = median(f0*2*pi) # angular frequency = 2*pi*f0
      
    }else if(length(HT) == 1){  # if only one peak is high

      high.peak <- sum(minmax > threshold)
      if(high.peak < length(peak_indices)){
        low.peak <- high.peak + 1
      }else{
        low.peak <- high.peak - 1 
      }
      fullperiod = 2*diff(abs(peak_indices[low.peak, high.peak])) 
      f0 = 1/(fullperiod) # nominal frequency
      a0 = median(f0*2*pi) # angular frequency = 2*pi*f0
      
    }else{ # if at least 2 highs and lows
      
      LT <- signal[peak_indices[minmax <= threshold]]
      
      fullperiod = diff(HT.peaks) 
      f0 = 1/(fullperiod) # nominal frequency
      a0 = median(f0*2*pi) # angular frequency = 2*pi*f0
      
    }
  }
  
  # Construct b
  if(NumberOfPeaks == 1){
    b0 <- runif(1, 0, pi/4) # b has to be < pi/2
  }else{
    Peak1IsHigh <- PeakIndicator[1] # = 1 (true) or 0 (false)
    Peak2IsHigh <- PeakIndicator[2] # = 1 (true) or 0 (false)
    Peak1Period <- peak_indices[1]
    Peak2Period <- peak_indices[2]
    
    if(sum(c(Peak1IsHigh,Peak2IsHigh))==1){ # if both are high or low
      b0 <- (1 - Peak1Period/Peak2Period)*2*pi
    }else{
      b0 <- (2 - Peak1Period/Peak2Period)*pi
    }
  }
  b0 <- b0-floor(b0/(2*pi))*(2*pi)
  
  # Etimate A and B0
  B0.0 = (mean(HT)+mean(LT))/2
  A.0 =  (mean(HT)-mean(LT))/2
  
  # noise estimate
  sigma0 <- sqrt(1/(N_org-1)*sum((signal - y_sm)^2))
  
  # Construct g(x), xi, Bx, rk, Y.fit and a, b
  gx.final <- constructgxi(peak_indices, PeakIndicator, N_org)
  xi.final <- diff(gx.final)
  xi.final[N_org]  <- xi.final[N_org - 1]
  a0 <- mean(xi.final)
  
  xkstar.final <- sort(c(FindXstar(gx.final, b0, pi/2),FindXstar(gx.final, b0, 3*pi/2)))
  AmpData <- EstimateBx(signal, xi.final, A.0, B0.0, b0)
  Bx.final <- AmpData$Bx
  rk.final <- AmpData$rk
  Y.fit <- f(xi.final, Bx.final, 1, A.0, b0)
  
  
  InitData <- list(pars = list(a0 = a0,
                               b0 = b0,
                               A.0 = A.0,
                               B0.0 = B0.0,
                               sigma0 = sigma0),
                   gx = gx.final,
                   xi = xi.final,
                   xkstar = xkstar.final,
                   rk = rk.final,
                   Bx = Bx.final,
                   Yfit = Y.fit)
  
  return(InitData)
  
}



### INITIALIZATION helper functions:


# likelihood loss
lloss <- function(Y.arg, Ysm.arg, fx.arg, penterm1, penterm1.sd=1){
  -sum(dnorm(Y.arg, fx.arg, sd=sqrt(mean((Ysm.arg-fx.arg)^2)),log=TRUE))
    -dnorm(penterm1, sd=penterm1.sd^(-1),log=TRUE)
}

# Define the linear piecewise function
linfunc1bp <- function(x, slopes, breakpoints) {
  return(
    ifelse(
      x < breakpoints[1],
      slopes[1] * x,
      slopes[2] * (x - breakpoints[1]) + slopes[1] * breakpoints[1])
  )
}

linfunc2bp <- function(x, slopes, breakpoints) {
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

##############################################################
######## --- Construct Growth Process from peaks --- #########
##############################################################
constructgxi <- function(peaks.arg, PeakIndicator.arg, signal.length){
  if(length(PeakIndicator.arg)==1){ # one high peak
    period <- diff(c(0,peaks.arg, signal.length))
    f_vec = 1/period # nominal frequency
    a_vec = f_vec*pi
    gxi_est <- cumsum(rep(a_vec, times=period))
  }else if(length(PeakIndicator.arg)==2){
    period <- diff(peaks.arg)
    f_vec = 1/period # nominal frequency
    if(PeakIndicator.arg[1] != PeakIndicator.arg[2]){
      a_vec = f_vec*pi
    }else{
      a_vec = f_vec*2*pi
    }
    gxi_est <- cumsum(rep(a_vec, times=signal.length))
  }else{
    period <- diff(peaks.arg)
    f_vec = 1/period # nominal frequency
    V <- as.integer(diff(PeakIndicator.arg) != 0) # = 1 if high-low peak in interval, 0 if low-low or high-high
    a_vec = f_vec*(pi*V + (1-V)*(2*pi))
    period_ext <- period
    period_ext[1] <- period_ext[1] + peaks.arg[1]
    period_ext[length(period_ext)] <- period_ext[length(period_ext)] + signal.length - peaks.arg[length(peaks.arg)]
    gxi_est <- cumsum(rep(a_vec, times=period_ext))
    
  }
  return(gxi_est)
}
