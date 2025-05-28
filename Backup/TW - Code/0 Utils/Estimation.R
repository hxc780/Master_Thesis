###############################################
#### --- MAIN Estimation Function  --- ########
###############################################
warp.alg <- function(Y.arg, delta.arg = 1, n.rep = 10, n.particles = 1000, n.mcmc = 50, n.memory = 5, init = "AMPD", SimAn.rep = 5000, SimAn.MinPeaks = .5, SimAn.MaxPeaks = 10, SimAn.Pen = 1, fast = FALSE){
  
  # INITIALIZATION
  # Note: init.data is a list of 7 items:
  # item 1: pars      list (a0, b0, A.0, B0.0, sigma0)
  # item 2: gx        vector g(x) integrated CIR process
  # item 3: xi        vector xi(x) CIR process
  # item 4: xkstar    vector peak position indices
  # item 5: rk        vector r_k* peak k displacement
  # item 6: Bx        vector B(x) = B0.0 + R(x) varying amplitude
  # item 7: Y.an.fit  vector \hat Y(x) fitted with sim.an parameters
  # note: psi0, beta0, omega2.0 are initialized directly in the estimation/warp.alg
  if(init == "SimAn") {
    print("Initialization: Simulated Annealing")
    init.data <- SimAn.Init(Y.arg, sim.its=SimAn.rep, min_peaks = SimAn.MinPeaks, max_peaks = SimAn.MaxPeaks, pen = SimAn.Pen)
  }else if (init == "AMPD") {
    print("Initialization: AMPD Peak Detection")
    init.data <- AMPD.init(Y.arg, correction = TRUE, beta=0.1, corr_neighbors = 0.1)
  }else{
    # Return an error if the input is invalid
    stop("Invalid value for 'init'. Please use either 'SimAn' or 'AMPD'.")
  }
  
  # ESTIMATION
  out.data <- warp.smc(Y.arg, delta.arg, init.data, n.rep, n.memory, n.particles, fast = fast)
  
  return(out.data)
  
}

