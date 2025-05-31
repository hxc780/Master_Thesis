AmpEstimation <- function(Y.arg, hspan = 0.2, denoise=FALSE) {
  # Computes the analytic signal using the Hilbert transform,  
  # extracting the instantaneous amplitude as the envelope of Y.  
  # The amplitude process is then smoothed using a moving average  
  # with a window size / fraction defined by hspan. Partial smoothing
  # is performed in the start and end.
  
  n <- length(Y.arg)
  amp <- numeric(n)
  if(!requireNamespace("gsignal", quietly = TRUE)){
    stop("Package 'gsignal' needed for the hilbert method. Please install it.")
  }
  if (!requireNamespace("zoo", quietly = TRUE)) {
    stop("Package 'zoo' needed for the rolling median filter. Please install it.")
  }
  
  if(denoise){
    Y_in <- zoo::rollapply(Y.arg, width = round(n/400), FUN = median, align = "center", fill = "extend")
  }else{
    Y_in <- Y.arg
  }
  
  analytic_signal <- gsignal::hilbert(Y_in)
  amp_raw <- abs(analytic_signal)
  
  # Smooth the amplitude envelope using a moving average.
  window_size <- max(3, round(hspan * n))
  # Use centered filter for smoothing
  smooth_filter <- rep(1/window_size, window_size)
  #amp_smoothed <- stats::filter(amp_raw, filter = smooth_filter, sides = 2)
  amp_smoothed <- zoo::rollapply(amp_raw, 
                                 width = window_size, 
                                 FUN = mean, 
                                 fill = NA, # "extend" instead of 2 lines below?
                                 align = "center", 
                                 partial = TRUE)
  # Replace NA values at boundaries with the raw amplitude values
  amp <- as.numeric(amp_smoothed)
  amp[is.na(amp)] <- amp_raw[is.na(amp)]
  
  return(amp)
}
